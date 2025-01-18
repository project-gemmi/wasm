#include <gemmi/it92.hpp>      // for IT92
#include <gemmi/c4322.hpp>     // for C4322
#include <gemmi/dencalc.hpp>   // for DensityCalculator
#include <gemmi/fourier.hpp>   // for transform_map_to_f_phi
#include <emscripten/emscripten.h>

struct DensityRoundTripResult {
  struct Point {
    float x, y;
    bool operator<(const Point& o) const { return x < o.x; }
  };
  double atom_radius;  // 0
  int grid_size[3];    // 8, 12, 16
  int used_points;     // 20
  int pre_len;         // 24
  float* density_pre;  // 28
  int post_len;        // 32
  float* density_post; // 36
  std::vector<Point> pre;
  std::vector<Point> post;
};
DensityRoundTripResult dr_result;

void remove_almost_identical(std::vector<DensityRoundTripResult::Point>& vec, double eps) {
  if (vec.empty())
    return;
  size_t prev = 0;
  for (size_t i = 1; i < vec.size(); ++i)
    if (vec[i].x - vec[prev].x >= eps)
      if (++prev != i)
        vec[prev] = vec[i];
  vec.resize(prev+1);
}

extern "C" {

void EMSCRIPTEN_KEEPALIVE run_srand() {
  std::srand(12345);
}

DensityRoundTripResult* EMSCRIPTEN_KEEPALIVE
density_roundtrip(int beam, int elem, double b_iso,
                  double a, double b, double c, double alpha, double beta, double gamma,
                  double dmin, double oversampling, double blur,
                  int smooth_cutoff,  // 0=no, 1=sinc, 2=cos taper
                  int add_more_points) {
  using TableX = gemmi::IT92<float>;
  using TableE = gemmi::C4322<float>;
  // DensityCalculator uses Table only in add_atom_density_to_grid which we don't call
  // so it doesn't matter if it's TableX or TableE
  gemmi::DensityCalculator<TableX, float> dencalc;
  dencalc.d_min = dmin;
  dencalc.rate = oversampling;
  dencalc.grid.unit_cell.set(a, b, c, alpha, beta, gamma);
  dencalc.initialize_grid();
  using CReal = float;
  // cf. DensityCalculator::do_add_atom_density_to_grid, isotropic case
  gemmi::Element el(elem);
  auto precal = beam == 2
      ? TableE::get(el, 0).precalculate_density_iso(b_iso)
      : TableX::get(el, 0).precalculate_density_iso(b_iso);
  if (!add_more_points) {
    // calculate pre for a line plot
    int n = 200;
    dr_result.pre.resize(n);
    for (size_t i = 0; i < dr_result.pre.size(); ++i) {
      double x = 5.0 * i / n;
      dr_result.pre[i].x = x;
      dr_result.pre[i].y = precal.calculate((CReal)x*x);
    }
    dr_result.post.clear();
  }

  // calculate post
  if (blur != 0.) {
    dencalc.blur = blur;
    precal = beam == 2
        ? TableE::get(el, 0).precalculate_density_iso(b_iso + blur)
        : TableX::get(el, 0).precalculate_density_iso(b_iso + blur);
  }
  CReal radius = dencalc.estimate_radius(precal, b_iso + blur);
  int counter = 0;
  auto random = []() { return 1.0 / (RAND_MAX + 1.0) * std::rand(); };
  gemmi::Fractional fpos(
      0.5 + random() / dencalc.grid.nu,
      0.5 + random() / dencalc.grid.nv,
      0.5 + random() / dencalc.grid.nw);
  dencalc.grid.template use_points_around<true>(fpos, radius, [&](float& point, double r2) {
      point += precal.calculate((CReal)r2);
      ++counter;
  }, /*fail_on_too_large_radius=*/false);
  gemmi::FPhiGrid<float> coefs = transform_map_to_f_phi(dencalc.grid, /*half_l=*/true);
  double max_inv_d2 = 1./(dmin*dmin);
  for (auto grid_point : coefs) {
    gemmi::Miller hkl = coefs.to_hkl(grid_point);
    double inv_d2 = coefs.unit_cell.calculate_1_d2(hkl);
    if (inv_d2 > max_inv_d2)
      *grid_point.value = 0;
    else {
      if (blur != 0.)
        *grid_point.value *= dencalc.reciprocal_space_multiplier(inv_d2);
      if (smooth_cutoff != 0) {
        double d2_ratio = inv_d2 / max_inv_d2;
        double factor = 1;
        if (smooth_cutoff == 1) {
          double k = gemmi::pi() * std::sqrt(d2_ratio);
          factor = sin(k) / k;
        } else if (smooth_cutoff == 2) {
          double k = std::sqrt(d2_ratio);
          double q1 = 0.5;
          if (k > q1)
            factor = 0.5 * (1 + cos(gemmi::pi() * (k - q1) / (1 - q1)));
        }
        *grid_point.value *= factor;
      }
    }
  }
  gemmi::transform_f_phi_grid_to_map_(std::move(coefs), dencalc.grid);
  dencalc.grid.template use_points_around<true>(fpos, radius, [&](float& point, double r2) {
      dr_result.post.push_back({std::sqrt((float)r2), point});
  }, /*fail_on_too_large_radius=*/false);
  std::sort(dr_result.post.begin(), dr_result.post.end());
  // remove almost identical values to not clutter the plot
  remove_almost_identical(dr_result.post, 0.002);

  // gather results
  dr_result.atom_radius = radius;
  dr_result.grid_size[0] = dencalc.grid.nu;
  dr_result.grid_size[1] = dencalc.grid.nv;
  dr_result.grid_size[2] = dencalc.grid.nw;
  dr_result.used_points = counter;
  dr_result.pre_len = (int) dr_result.pre.size();
  dr_result.density_pre = &dr_result.pre[0].x;
  dr_result.post_len = (int) dr_result.post.size();
  dr_result.density_post = &dr_result.post[0].x;

  return &dr_result;
}

} // extern "C"
