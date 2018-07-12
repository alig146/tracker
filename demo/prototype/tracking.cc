/*
 * demo/prototype/tracking.cc
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

#include <tracker/analysis/analysis.hh>
#include <tracker/analysis/monte_carlo.hh>
#include <tracker/analysis/track.hh>
#include <tracker/analysis/vertex.hh>

// TODO: remove
#include <tracker/analysis/tree.hh>

#include <tracker/geometry.hh>
#include <tracker/plot.hh>
#include <tracker/reader.hh>

#include <tracker/util/index_vector.hh>
#include <tracker/util/bit_vector.hh>
#include <tracker/util/io.hh>

#include "geometry.hh"
#include "logging.hh"

//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {

//__Find Primary Tracks for Prototype___________________________________________________________
const analysis::track_vector find_primary_tracks(const analysis::event& event,
                                                 const reader::tracking_options& options,
                                                 analysis::event& non_track_points) {
  analysis::event combined_rpc_hits;
  analysis::full_event original_rpc_hits;
  const auto optimized_event = combine_rpc_hits(event, combined_rpc_hits, original_rpc_hits);
  const auto layers          = analysis::partition(optimized_event, options.layer_axis, options.layer_depth);
  const auto seeds           = analysis::seed(options.seed_size, layers, options.line_width);
  const auto tracking_vector = reset_seeds(analysis::join_all(seeds), combined_rpc_hits, original_rpc_hits);
  const auto out = analysis::overlap_fit_seeds(tracking_vector, options.layer_axis, 1UL);

  // TODO: improve efficiency
  const auto size = event.size();
  util::bit_vector save_list(size);
  for (const auto& track : out) {
    for (const auto& point : track.event()) {
      const auto search = util::algorithm::range_binary_find_first(event,
                                                                   point,
                                                                   type::t_ordered<analysis::hit>{});
      if (search != event.cend())
        save_list.set(static_cast<std::size_t>(search - event.cbegin()));
    }
  }

  non_track_points.clear();
  non_track_points.reserve(size - save_list.count());
  save_list.unset_conditional_push_back(event, non_track_points);

  return out;
}
//----------------------------------------------------------------------------------------------

//__Find Secondary Tracks for Prototype_________________________________________________________
const analysis::track_vector find_secondary_tracks(const analysis::event& event,
                                                   const reader::tracking_options& options,
                                                   analysis::event& non_track_points) {
  analysis::event combined_rpc_hits;
  analysis::full_event original_rpc_hits;
  const auto optimized_event = combine_rpc_hits(event, combined_rpc_hits, original_rpc_hits);
  const auto layers          = analysis::partition(optimized_event, options.layer_axis, options.layer_depth);
  const auto seeds           = analysis::seed(2UL, layers, options.line_width);
  const auto tracking_vector = reset_seeds(seeds, combined_rpc_hits, original_rpc_hits);
  const auto out = analysis::overlap_fit_seeds(tracking_vector, options.layer_axis, 0UL);

  // TODO: improve efficiency
  const auto size = event.size();
  util::bit_vector save_list(size);
  for (const auto& track : out) {
    for (const auto& point : track.event()) {
      const auto search = util::algorithm::range_binary_find_first(event,
                                                                   point,
                                                                   type::t_ordered<analysis::hit>{});
      if (search != event.cend())
        save_list.set(static_cast<std::size_t>(search - event.cbegin()));
    }
  }

  non_track_points.clear();
  non_track_points.reserve(size - save_list.count());
  save_list.unset_conditional_push_back(event, non_track_points);

  return out;
}
//----------------------------------------------------------------------------------------------

//__Prototype Tracking Algorithm________________________________________________________________
int prototype_tracking(int argc,
                       char* argv[]) {
  const auto options = reader::parse_input(argc, argv);
  const auto detector_map = reader::import_detector_map(options.geometry_map_file);
  const auto time_resolution_map = reader::import_time_resolution_map(options.geometry_time_file);

  plot::init(options.draw_events);
  geometry::open(options.geometry_file, options.default_time_error, time_resolution_map);

  std::cout << "Begin Tracking in " << options.data_directory << ":\n\n";
  const auto statistics_path_prefix = options.statistics_directory + "/" + options.statistics_file_prefix;
  const plot::value_tag filetype_tag("FILETYPE", "MATHUSLA TRACKING STATFILE");
  const plot::value_tag project_tag("PROJECT", "Prototype");

  std::uint_fast64_t path_counter{};
  for (const auto& path : reader::root::search_directory(options.data_directory, options.data_file_extension)) {
    const auto path_counter_string = std::to_string(path_counter++);
    const auto statistics_save_path = statistics_path_prefix
                                    + path_counter_string
                                    + "." + options.statistics_file_extension;

    print_bar();
    std::cout << "Read Path: " << path << "\n";

    const auto event_bundle = reader::root::import_event_mc_bundle(path, options, detector_map);
    const auto imported_events = event_bundle.events;
    const auto import_size = imported_events.size();
    if (import_size == 0)
      continue;

    const auto mc_imported_events = event_bundle.true_events;

    auto histograms = generate_histograms();
    for (std::uint_fast64_t event_counter{}; event_counter < import_size; ++event_counter) {
      const auto& event = imported_events[event_counter];
      const auto event_size = event.size();
      const auto event_counter_string = std::to_string(event_counter);

      const auto compressed_event = analysis::compress(event, options.time_smearing);
      const auto compression_size = event_size / static_cast<type::real>(compressed_event.size());

      if (event.empty() || compression_size == event_size)
        continue;

      const auto event_density = modified_geometry_event_density(compressed_event);
      print_event_summary(event_counter, event_size, compression_size, event_density);

      if (event_density >= options.event_density_limit)
        continue;

      plot::canvas canvas("event" + event_counter_string, path + event_counter_string);
      if (options.draw_events) {
        draw_detector_centers(canvas);
        draw_mc_tracks(canvas, analysis::mc::convert(mc_imported_events[event_counter]));
        canvas.add_points(compressed_event, 0.8, plot::color::BLACK);
      }

      analysis::event non_primary_track_points, non_secondary_track_points;
      auto tracks = find_primary_tracks(compressed_event, options, non_primary_track_points);
      /*
      auto secondary_tracks = find_secondary_tracks(non_primary_track_points,
                                                    options,
                                                    non_secondary_track_points);
      tracks.reserve(tracks.size() + secondary_tracks.size());
      tracks.insert(tracks.cend(),
                    std::make_move_iterator(secondary_tracks.cbegin()),
                    std::make_move_iterator(secondary_tracks.cend()));
      */
      save_tracks(tracks, canvas, histograms, options);
      print_tracking_summary(event, tracks);

      save_vertex(analysis::vertex(tracks), canvas, histograms, options);

      canvas.draw();
    }
    plot::value_tag input_tag("DATAPATH", path);
    plot::value_tag event_tag("EVENTS", std::to_string(import_size));
    histograms.draw_all();
    plot::save_all(statistics_save_path, histograms, filetype_tag, project_tag, input_tag, event_tag);
  }

  print_bar();
  geometry::close();
  plot::end();
  return 0;
}
//----------------------------------------------------------------------------------------------

//__Silent Prototype Tracking Algorithm_________________________________________________________
int silent_prototype_tracking(int argc,
                              char* argv[]) {
  util::io::remove_buffer(std::cout, std::cerr, std::clog);
  return prototype_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------

} /* namespace MATHUSLA */

//__Main Function: Prototype Tracker____________________________________________________________
int main(int argc,
         char* argv[]) {
  return MATHUSLA::prototype_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------