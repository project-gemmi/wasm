// Copyright 2019 Global Phasing Ltd.

#define POCKETFFT_CACHE_SIZE 0
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_VECTORS
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/fourier.hpp>  // for update_cif_block
#include <gemmi/math.hpp>     // for Variance
#include <emscripten/bind.h>

class MtzMap {
public:
  // char* or void* cannot be used
  // https://github.com/emscripten-core/emscripten/issues/9448
  MtzMap(int32_t data, size_t size)
    : data_((const char*) data) {
    try {
      mtz_.read_stream(gemmi::MemoryStream(data_, data_ + size), false);
    } catch (std::runtime_error& e) {
      (void) e;
    }
  }

  int32_t calculate_map(bool diff_map) {
    static const char* default_labels[] = {
      "FWT", "PHWT", "DELFWT", "PHDELWT",
      "2FOFCWT", "PH2FOFCWT", "FOFCWT", "PHFOFCWT",
    };
    const gemmi::Mtz::Column* f_col;
    const gemmi::Mtz::Column* phi_col;
    for (int i = (diff_map ? 2 : 0); i < 8; i += 4)
      if ((f_col = mtz_.column_with_label(default_labels[i], nullptr)) &&
          (phi_col = mtz_.column_with_label(default_labels[i+1], nullptr))) {
        break;
      }
    if (!f_col || !phi_col) {
      std::puts("Default map coefficient labels not found");
      return 0;
    }
    try {
      const float* raw_data = (const float*)(data_ + 80);
      gemmi::MtzExternalDataProxy proxy(mtz_, raw_data);
      auto size = gemmi::get_size_for_hkl(proxy, {{0, 0, 0}}, 3.);
      gemmi::Grid<std::complex<float>> coefs
        = gemmi::get_f_phi_on_grid<float>(proxy, f_col->idx, phi_col->idx,
                                          size, /*half_l=*/true,
                                          gemmi::HklOrient::LKH);
      gemmi::transform_f_phi_grid_to_map_(std::move(coefs), grid_);
    } catch (std::runtime_error& e) {
      std::puts(e.what());
      return 0;
    }
    gemmi::Variance grid_variance(grid_.data.begin(), grid_.data.end());
    rmsd_ = std::sqrt(grid_variance.for_population());
    return (int32_t) grid_.data.data();
  }

  int get_nx() const { return grid_.nu; }
  int get_ny() const { return grid_.nv; }
  int get_nz() const { return grid_.nw; }
  double get_rmsd() const { return rmsd_; }
  double get_cell_param(int n) const {
    const gemmi::UnitCell& cell = mtz_.cell;
    switch (n) {
      case 0: return cell.a;
      case 1: return cell.b;
      case 2: return cell.c;
      case 3: return cell.alpha;
      case 4: return cell.beta;
      case 5: return cell.gamma;
    }
    return 0.;
  }

private:
  gemmi::Mtz mtz_;
  const char* data_;
  gemmi::Grid<float> grid_;
  double rmsd_;
};

using emscripten::class_;
using emscripten::allow_raw_pointers;

EMSCRIPTEN_BINDINGS(GemmiMtz) {
  class_<MtzMap>("MtzMap")
    .constructor<int32_t, size_t>(allow_raw_pointers())
    .function("calculate_map", &MtzMap::calculate_map, allow_raw_pointers())
    .property("nx", &MtzMap::get_nx)
    .property("ny", &MtzMap::get_ny)
    .property("nz", &MtzMap::get_nz)
    .property("rmsd", &MtzMap::get_rmsd)
    .function("cell_param", &MtzMap::get_cell_param)
    ;
}
