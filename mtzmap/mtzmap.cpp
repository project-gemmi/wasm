// Copyright 2019 Global Phasing Ltd.

#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/fourier.hpp>  // for update_cif_block
#include <gemmi/math.hpp>     // for Variance
#include <emscripten/bind.h>

class MtzMap {
public:
  MtzMap(int32_t data_, size_t size) {
    // char* or void* cannot be used
    // https://github.com/emscripten-core/emscripten/issues/9448
    const char* data = (const char*) data_;
    try {
      mtz_.read_stream(gemmi::MemoryStream(data, data + size), true);
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
    if (!f_col || !phi_col)
      return 0; // "Default map coefficient labels not found."
    try {
      gemmi::Grid<std::complex<float>> coef_grid
        = gemmi::get_f_phi_on_grid<float>(gemmi::MtzDataProxy{mtz_},
                                          f_col->idx, phi_col->idx, true,
                                          {{0, 0, 0}}, 3.);
      grid_ = gemmi::transform_f_phi_grid_to_map(std::move(coef_grid));
    } catch (std::runtime_error& e) {
      (void) e;
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
