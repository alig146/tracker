/*
 * src/tracker/analysis/vertex.cc
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <tracker/analysis/vertex.hh>

#include <tracker/core/stat.hh>
#include <tracker/core/units.hh>

#include <tracker/util/algorithm.hh>
#include <tracker/util/error.hh>
#include <tracker/util/io.hh>
#include <tracker/util/math.hh>

#include "../helper/analysis.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////

//__Calculate Distance from Vertex to Track at Fixed Time_______________________________________
real _vertex_track_r3_distance(const real t,
                               const real x,
                               const real y,
                               const real z,
                               const track& track) {
  const auto track_point = track.at_t(t);
  return util::math::hypot(track_point.x - x, track_point.y - y, track_point.z - z);
}
//----------------------------------------------------------------------------------------------

//__Calculate Distance with Error from Vertex to Track at Fixed Time____________________________
const stat::type::uncertain_real _vertex_track_r3_distance_with_error(const real t,
                                                                      const real x,
                                                                      const real y,
                                                                      const real z,
                                                                      const track& track) {
  const auto track_point = track.at_t(t);
  const auto dx = track_point.x - x;
  const auto dy = track_point.y - y;
  const auto dz = track_point.z - z;
  const auto total_dt = t - track.t0_value();
  const auto distance = util::math::hypot(dx, dy, dz);
  const auto inverse_distance = 1.0L / distance;
  const auto dx_by_D = dx * inverse_distance;
  const auto dy_by_D = dy * inverse_distance;
  const auto dz_by_D = dz * inverse_distance;
  const real_array<6UL> gradient{
    -util::math::fused_product(track.vx_value(), dx_by_D,
                               track.vy_value(), dy_by_D,
                               track.vz_value(), dz_by_D),
    dx_by_D,
    dy_by_D,
    total_dt * dx_by_D,
    total_dt * dy_by_D,
    total_dt * dz_by_D};
  return stat::type::uncertain_real(
    distance,
    stat::error::propagate(gradient, track.covariance_matrix()));
}
//----------------------------------------------------------------------------------------------

//__Calculate Squared Residual of Vertex wrt Track______________________________________________
real _vertex_squared_residual(const stat::type::uncertain_real& distance) {
  return util::math::sum_squares(distance.value / distance.error);
}
//----------------------------------------------------------------------------------------------

//__Calculate Squared Residual of Vertex wrt Track______________________________________________
real _vertex_squared_residual(const real t,
                              const real x,
                              const real y,
                              const real z,
                              const track& track) {
  return _vertex_squared_residual(_vertex_track_r3_distance_with_error(t, x, y, z, track));
}
//----------------------------------------------------------------------------------------------

//__Fast Guess of Initial Track Parameters______________________________________________________
vertex::fit_parameters _guess_vertex(const track_vector& tracks) {
  const auto size = tracks.size();

  std::vector<full_hit> track_fronts;
  track_fronts.reserve(size);
  util::algorithm::back_insert_transform(tracks, track_fronts, [](const auto& track) {
    const auto front = track.full_front();
    const auto front_t = front.t;
    const auto point = track.at_t(front_t);
    const auto error = track.error_at_t(front_t);
    return full_hit{point.t, point.x, point.y, point.z,
             r4_point{front.width.t, error.x, error.y, error.z}};
  });

  real_vector t_errors, x_errors, y_errors, z_errors;
  t_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, t_errors,
    [](const auto& front) { return front.width.t; });
  x_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, x_errors,
    [](const auto& front) { return stat::error::uniform(front.width.x); });
  y_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, y_errors,
    [](const auto& front) { return stat::error::uniform(front.width.y); });
  z_errors.reserve(size);
  util::algorithm::back_insert_transform(track_fronts, z_errors,
    [](const auto& front) { return stat::error::uniform(front.width.z); });

  const auto average_point = std::accumulate(track_fronts.cbegin(), track_fronts.cend(), r4_point{},
    [](const auto sum, const auto& point) { return sum + reduce_to_r4(point); })
      / static_cast<real>(size);

  return {{average_point.t, stat::error::propagate_average(t_errors), 0, 0},
          {average_point.x, stat::error::propagate_average(x_errors), 0, 0},
          {average_point.y, stat::error::propagate_average(y_errors), 0, 0},
          {average_point.z, stat::error::propagate_average(z_errors), 0, 0}};
}
//----------------------------------------------------------------------------------------------

//__Gaussian Negative Log Likelihood Calculation________________________________________________
thread_local track_vector&& _nll_fit_tracks = {};
void _gaussian_nll(Int_t&, Double_t*, Double_t& out, Double_t* x, Int_t) {
  out = std::accumulate(_nll_fit_tracks.cbegin(), _nll_fit_tracks.cend(), 0.0L,
    [&](const auto sum, const auto& track) {
      const auto distance = _vertex_track_r3_distance_with_error(x[0], x[1], x[2], x[3], track);
      return sum + std::fma(0.5L, _vertex_squared_residual(distance), std::log(distance.error));
    });
}
//----------------------------------------------------------------------------------------------

//__MINUIT Gaussian Fitter______________________________________________________________________
bool _fit_tracks_minuit(const track_vector& tracks,
                        vertex::fit_parameters& parameters,
                        vertex::covariance_matrix_type& covariance_matrix) {
  using namespace helper::minuit;
  auto& t = parameters.t;
  auto& x = parameters.x;
  auto& y = parameters.y;
  auto& z = parameters.z;
  _nll_fit_tracks = tracks;

  TMinuit minuit;
  initialize(minuit, "T", t, "X", x, "Y", y, "Z", z);

  // FIXME: error handling for -> execute(minuit, _gaussian_nll);
  if (execute(minuit, _gaussian_nll) == error::diverged)
    return false;

  get_parameters(minuit, t, x, y, z);
  get_covariance<vertex::free_parameter_count>(minuit, covariance_matrix);
  return true;
}
//----------------------------------------------------------------------------------------------

} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(const track_vector& tracks) {
  reset(tracks);
}
//----------------------------------------------------------------------------------------------

//__Vertex Constructor__________________________________________________________________________
vertex::vertex(track_vector&& tracks) {
  reset(std::move(tracks));
}
//----------------------------------------------------------------------------------------------

//__Vertex Point________________________________________________________________________________
const r4_point vertex::point() const {
  return r4_point{t_value(), x_value(), y_value(), z_value()};
}
//----------------------------------------------------------------------------------------------

//__Error in Calculation of Vertex Point________________________________________________________
const r4_point vertex::point_error() const {
  return r4_point{t_error(), x_error(), y_error(), z_error()};
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter from Vertex_______________________________________________________________
const fit_parameter vertex::fit_of(const vertex::parameter p) const {
  switch (p) {
    case vertex::parameter::T: return _final.t;
    case vertex::parameter::X: return _final.x;
    case vertex::parameter::Y: return _final.y;
    case vertex::parameter::Z: return _final.z;
  }
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Value from Vertex_________________________________________________________
real vertex::value(const vertex::parameter p) const {
  return fit_of(p).value;
}
//----------------------------------------------------------------------------------------------

//__Get Fit Parameter Error from Vertex_________________________________________________________
real vertex::error(const vertex::parameter p) const {
  return fit_of(p).error;
}
//----------------------------------------------------------------------------------------------

//__Check if Track Fit Converged________________________________________________________________
bool vertex::fit_diverged() const noexcept {
  return _guess != _final && _final == vertex::fit_parameters{};
}
//----------------------------------------------------------------------------------------------

//__Get Distance from Each Track at Time of Vertex______________________________________________
real_vector vertex::distances() const {
  real_vector out;
  out.reserve(_tracks.size());
  util::algorithm::back_insert_transform(_tracks, out, [&](const auto& track) {
    return _vertex_track_r3_distance(t_value(), x_value(), y_value(), z_value(), track); });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Get Error in Distance of Each Track at Time of Vertex_______________________________________
real_vector vertex::distance_errors() const {
  real_vector out;
  out.reserve(_tracks.size());
  util::algorithm::back_insert_transform(_tracks, out, [&](const auto& track) {
    return _vertex_track_r3_distance_with_error(t_value(), x_value(), y_value(), z_value(), track).error; });
  return out;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared Test Statistic__________________________________________________________________
real vertex::chi_squared() const {
  return std::accumulate(_delta_chi2.cbegin(), _delta_chi2.cend(), 0.0L);
}
//----------------------------------------------------------------------------------------------

//__Vertex Degrees of Freedom___________________________________________________________________
std::size_t vertex::degrees_of_freedom() const {
  return 4UL;
}
//----------------------------------------------------------------------------------------------

//__Chi-Squared per Degree of Freedom___________________________________________________________
real vertex::chi_squared_per_dof() const {
  return chi_squared() / degrees_of_freedom();
}
//----------------------------------------------------------------------------------------------

//__Get Variance of a Vertex Parameter__________________________________________________________
real vertex::variance(const vertex::parameter p) const {
  return covariance(p, p);
}
//----------------------------------------------------------------------------------------------

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Get Shift Index of Vertex Parameters for Covariance Matrix__________________________________
constexpr std::size_t _shift_covariance_index(const vertex::parameter p) {
  switch (p) {
    case vertex::parameter::T: return 0;
    case vertex::parameter::X: return 1;
    case vertex::parameter::Y: return 2;
    case vertex::parameter::Z: return 3;
  }
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Get Covariance between Vertex Parameters____________________________________________________
real vertex::covariance(const vertex::parameter p,
                        const vertex::parameter q) const {
  return _covariance[4 * _shift_covariance_index(p) + _shift_covariance_index(q)];
}
//----------------------------------------------------------------------------------------------

//__Reset Vertex with Given Tracks______________________________________________________________
std::size_t vertex::reset(const track_vector& tracks) {
  _tracks = tracks;
  const auto new_size = _tracks.size();
  if (new_size > 1) {
    _guess = _guess_vertex(_tracks);
    _final = _guess;
    if (_fit_tracks_minuit(_tracks, _final, _covariance)) {
      util::algorithm::back_insert_transform(_tracks, _delta_chi2,
        [&](const auto& track) {
          return _vertex_squared_residual(t_value(), x_value(), y_value(), z_value(), track);
        });
      return size();
    }
  }
  _delta_chi2.resize(new_size, 0);
  _final = {};
  _covariance = {};
  return new_size;
}
//----------------------------------------------------------------------------------------------

//__Insert Track into Vertex and Refit__________________________________________________________
std::size_t vertex::insert(const track& track) {
  if (std::find(_tracks.cbegin(), _tracks.cend(), track) != _tracks.cend()) {
    _tracks.push_back(track);
    _tracks.shrink_to_fit();
    return reset(_tracks);
  }
  return size();
}
//----------------------------------------------------------------------------------------------

//__Insert Tracks into Vertex and Refit_________________________________________________________
std::size_t vertex::insert(const track_vector& tracks) {
  _tracks.reserve(size() + tracks.size());
  const auto begin = _tracks.cbegin();
  std::copy_if(tracks.cbegin(), tracks.cend(), std::back_inserter(_tracks),
    [&](const auto& t) { return std::find(begin, _tracks.cend(), t) != _tracks.cend(); });
  _tracks.shrink_to_fit();
  return reset(_tracks);
}
//----------------------------------------------------------------------------------------------

//__Remove Track from Vertex and Refit__________________________________________________________
std::size_t vertex::remove(const std::size_t index) {
  // TODO: improve efficiency
  const auto s = size();
  if (index >= s)
    return s;

  track_vector saved_tracks;
  saved_tracks.reserve(s - 1);
  const auto begin = _tracks.cbegin();
  const auto end = _tracks.cend();
  saved_tracks.insert(saved_tracks.cend(), begin, begin + index);
  saved_tracks.insert(saved_tracks.cend(), begin + index + 1, end);
  return reset(saved_tracks);
}
//----------------------------------------------------------------------------------------------

//__Remove Tracks from Vertex and Refit_________________________________________________________
std::size_t vertex::remove(const std::vector<std::size_t>& indices) {
  // TODO: improve efficiency
  const auto sorted = util::algorithm::copy_sort_range(indices);
  const auto s = size();
  track_vector saved_tracks;
  saved_tracks.reserve(s);
  for (std::size_t hit_index{}, removal_index{}; hit_index < s; ++hit_index) {
    if (sorted[removal_index] == hit_index) {
      ++removal_index;
    } else {
      saved_tracks.push_back(std::move(_tracks[hit_index]));
    }
  }
  saved_tracks.shrink_to_fit();
  return reset(saved_tracks);
}
//----------------------------------------------------------------------------------------------

//__Remove Track from Vertex Below Maximum Chi-Squared and Refit________________________________
std::size_t vertex::prune_on_chi_squared(const real max_chi_squared) {
  const auto s = size();
  std::vector<std::size_t> indices;
  indices.reserve(s);
  for (std::size_t i{}; i < s; ++i) {
    if (chi_squared_vector()[i] > max_chi_squared)
      indices.push_back(i);
  }
  return remove(indices);
}
//----------------------------------------------------------------------------------------------

//__Fill Plots with Vertexing Variables_________________________________________________________
void vertex::fill_plots(plot::histogram_collection& collection,
                        const vertex::plotting_keys& keys) const {
  if (collection.count(keys.t)) collection[keys.t].insert(t_value() / units::time);
  if (collection.count(keys.x)) collection[keys.x].insert(x_value() / units::length);
  if (collection.count(keys.y)) collection[keys.y].insert(y_value() / units::length);
  if (collection.count(keys.z)) collection[keys.z].insert(z_value() / units::length);
  if (collection.count(keys.t_error)) collection[keys.t_error].insert(t_error() / units::time);
  if (collection.count(keys.x_error)) collection[keys.x_error].insert(x_error() / units::length);
  if (collection.count(keys.y_error)) collection[keys.y_error].insert(y_error() / units::length);
  if (collection.count(keys.z_error)) collection[keys.z_error].insert(z_error() / units::length);

  if (collection.count(keys.distance)) {
    auto& distance_histogram = collection[keys.distance];
    for (const auto& distance : distances())
      distance_histogram.insert(distance / units::length);
  }
  if (collection.count(keys.distance_error)) {
    auto& distance_error_histogram = collection[keys.distance_error];
    for (const auto& distance_error : distance_errors())
      distance_error_histogram.insert(distance_error / units::length);
  }

  if (collection.count(keys.chi_squared_per_dof)) collection[keys.chi_squared_per_dof].insert(chi_squared_per_dof());
  if (collection.count(keys.size)) collection[keys.size].insert(size());
}
//----------------------------------------------------------------------------------------------

//__Draw Fit Vertex_____________________________________________________________________________
void vertex::draw(plot::canvas& canvas,
                  const real size,
                  const plot::color color,
                  const bool with_errors) const {
  // TODO: decide what to do with convergence
  if (fit_converged()) {
    canvas.add_point(point(), size, color);
    if (with_errors)
      canvas.add_box(point(), point_error().x, point_error().y, point_error().z, size, color);
  }
}
//----------------------------------------------------------------------------------------------

//__Draw Guess Vertex___________________________________________________________________________
void vertex::draw_guess(plot::canvas& canvas,
                        const real size,
                        const plot::color color,
                        const bool with_errors) const {
  // TODO: add with_errors argument
  canvas.add_point(r3_point{guess_fit().x.value, guess_fit().y.value, guess_fit().z.value}, size, color);
}
//----------------------------------------------------------------------------------------------


//__Vertex Output Stream Operator_______________________________________________________________
std::ostream& operator<<(std::ostream& os,
                         const vertex& vertex) {
  static const std::string bar(80, '-');
  os << bar << "\n";

  if (vertex.fit_diverged()) {
    os << "* Vertex Status: " << util::io::bold << "DIVERGED" << util::io::reset_font << "\n";
    const auto guess = vertex.guess_fit();
    os << "* Guess Parameters:\n"
       << "    T: " << guess.t.value << "  (+/- " << guess.t.error << ")\n"
       << "    X: " << guess.x.value << "  (+/- " << guess.x.error << ")\n"
       << "    Y: " << guess.y.value << "  (+/- " << guess.y.error << ")\n"
       << "    Z: " << guess.z.value << "  (+/- " << guess.z.error << ")\n";
  } else {
    os << "* Vertex Status: " << util::io::bold << "CONVERGED" << util::io::reset_font << "\n";
    os << "* Parameters:\n"
       << "    T: " << vertex.t_value() << "  (+/- " << vertex.t_error() << ")\n"
       << "    X: " << vertex.x_value() << "  (+/- " << vertex.x_error() << ")\n"
       << "    Y: " << vertex.y_value() << "  (+/- " << vertex.y_error() << ")\n"
       << "    Z: " << vertex.z_value() << "  (+/- " << vertex.z_error() << ")\n";

    os << "* Tracks: \n";
    const auto size = vertex.size();
    const auto& tracks = vertex.tracks();
    const auto& distances = vertex.distances();
    const auto& errors = vertex.distance_errors();
    for (std::size_t i{}; i < size; ++i) {
      const auto& track = tracks[i];
      os << "    " << distances[i] << "  (+/- " << errors[i]
         << ")\n      from (" << track.t0_value() << ", "
                              << track.x0_value() << ", "
                              << track.y0_value() << ", "
                              << track.z0_value() << ", "
                              << track.vx_value() << ", "
                              << track.vy_value() << ", "
                              << track.vz_value() << ")\n";
    }

    os.precision(7);
    os << "* Statistics: \n"
       << "    dof:      " << vertex.degrees_of_freedom()               << "\n"
       << "    chi2:     " << vertex.chi_squared() << " = ";
    util::io::print_range(vertex.chi_squared_vector(), " + ", "", os)   << "\n";
    os << "    chi2/dof: " << vertex.chi_squared_per_dof()              << "\n"
       << "    p-value:  " << stat::chi_squared_p_value(vertex)         << "\n"
       << "    cov mat:  | ";
    const auto matrix = vertex.covariance_matrix();
    for (size_t i = 0; i < 4; ++i) {
      if (i > 0) os << "              | ";
      for (size_t j = 0; j < 4; ++j) {
        const auto cell = matrix[4*i+j];
        if (i == j) {
          os << util::io::bold << util::io::underline
             << cell << util::io::reset_font << " ";
        } else {
          os << cell << " ";
        }
      }
      os << "|\n";
    }
  }

  return os << bar;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */