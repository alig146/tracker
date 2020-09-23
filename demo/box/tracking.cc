/*
 * demo/box/tracking.cc
 *
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

#include <tracker/analysis.hh>

#include <tracker/geometry.hh>
#include <tracker/plot.hh>
#include <tracker/reader.hh>
#include <tracker/script.hh>

#include <tracker/analysis/type.hh>
#include <tracker/core/type.hh>
#include <tracker/core/units.hh>

#include "geometry.hh"
#include "io.hh"
#include "TTree.h"
#include <TFile.h>
#include "TSystem.h"
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class std::vector<double>+;
#pragma link C++ class std::vector<int>+;
#endif


//__Namespace Alias_____________________________________________________________________________
namespace analysis = MATHUSLA::TRACKER::analysis;
namespace mc       = analysis::mc;
namespace geometry = MATHUSLA::TRACKER::geometry;
namespace plot     = MATHUSLA::TRACKER::plot;
namespace reader   = MATHUSLA::TRACKER::reader;
namespace script   = MATHUSLA::TRACKER::script;
namespace reader   = MATHUSLA::TRACKER::reader;
namespace box      = MATHUSLA::box;
namespace util     = MATHUSLA::util;
//----------------------------------------------------------------------------------------------

namespace MATHUSLA {



//__Alter Event_________________________________________________________________________________
const analysis::full_event alter_event(const analysis::full_event& event,
                                       const script::tracking_options& options,
                                       const box::io::extension_parser& extension) {
  return box::geometry::restrict_layer_count(
           mc::use_efficiency(
             mc::add_noise<box::geometry>(
               event,
               options.simulated_noise_rate,
               options.event_time_window.begin,
               options.event_time_window.end,
               box::geometry::full(extension.layer_count)),
             options.simulated_efficiency),
           extension.layer_count);
}
//----------------------------------------------------------------------------------------------

//__Find Tracks for Box_________________________________________________________________________
const analysis::track_vector find_tracks(const analysis::full_event& event,
                                         const script::tracking_options& options,
                                         const type::real limit_chi_squared,
                                         const std::size_t overlap,
                                         plot::canvas& canvas,
                                         analysis::full_event& non_track_points) {
  namespace ash = analysis::seed_heuristic;
  const auto layers = analysis::partition(event, options.layer_axis, options.layer_depth);
  const auto seeds = analysis::seed(
    options.seed_size,
    layers,
    ash::all<ash::double_cone, ash::piecewise_speed>{
      {options.line_width},
      {0.9L * units::speed_of_light}});


  // for (const auto& seed : seeds) {
  // 	  //    for (std::size_t i{}; i < seed.size() - 1UL; ++i)
  // 	  //      canvas.add_line(type::reduce_to_r4(seed[i]), type::reduce_to_r4(seed[i + 1UL]), 1, plot::color::BLACK);
  //   for (const auto& point : seed)
  //     std::cout << type::reduce_to_r4(point) << " ";
  //   std::cout << "\n";
  // }


  const auto joined = analysis::join_all(seeds);
//std::cout << "inputttttt:" << joined << '\n';
  auto first_tracks = analysis::independent_fit_seeds(joined, options.layer_axis);
  for (auto& track : first_tracks)
    track.prune_on_chi_squared(limit_chi_squared);

  const auto out = analysis::overlap_fit_tracks(first_tracks, overlap);
  non_track_points = analysis::non_tracked_points(event, out, true);
  return out;
}
//----------------------------------------------------------------------------------------------

//__Track an Event Bundle_______________________________________________________________________
void track_event_bundle(const script::path_vector& paths,
                        const mc::full_event_vector_bundle& bundle,
                        const script::tracking_options& options,
                        const box::io::extension_parser& extension,
                        const script::path_type& save_path) {
  static const plot::value_tag filetype_tag("FILETYPE", "MATHUSLA TRACKING STATFILE");
  static const plot::value_tag project_tag("PROJECT", "Box");

  const auto energy_events = bundle.energy_events;
  const auto complete_events = bundle.complete_events;
  const auto imported_events = bundle.events;

  const auto import_size = imported_events.size();
  if (import_size == 0UL)
    return;

  const auto mc_imported_events = bundle.true_events;

  std::cout << "Event Count: " << import_size << "\n";

  // analysis::track::tree track_tree{"track_tree", "MATHUSLA Track Tree"};
  // analysis::vertex::tree vertex_tree{"vertex_tree", "MATHUSLA Vertex Tree"};
  // track_tree.set_file(save_path);
  // vertex_tree.set_file(save_path);
  // track_tree.add_friend(vertex_tree, "vertex");
  // vertex_tree.add_friend(track_tree, "track");

  std::vector<double> layer_center = {6002.5, 6105.5, 8002.5, 8105.5, 8502.5, 8605.5, 8708.5, 8811.5, 8914.5};
  std::vector<std::vector<double>> detector_x_limits = {{-4955, -4045}, {-3955, -3045}, {-2955, -2045}, {-1955, -1045}, {-955, -45}, {45, 955}, {1045, 1955}, {2045, 2955}, {3045, 3955}, {4045, 4955}};
  std::vector<std::vector<double>> detector_z_limits = {{6995, 7905}, {7995, 8905}, {8995, 9905}, {9995, 10905}, {10995, 11905}, {11995, 12905}, {12995, 13905}, {13995, 14905}, {14995, 15905}, {15995, 16905}};

  TTree integral_tree("integral_tree", "MATHUSLA Tree");

  //__Make Vertex Branches________________________________________________________________________
  gSystem->Load("libMylib.so");
  std::vector<double> vertex_numTracks;
  std::vector<double> vertex_chi2;
  std::vector<double> vertex_chi2_per_dof;
  std::vector<double> vertex_chi2_p_value;
  std::vector<double> vertex_t;
  std::vector<double> vertex_x;
  std::vector<double> vertex_y;
  std::vector<double> vertex_z;
  std::vector<double> vertex_t_error;
  std::vector<double> vertex_x_error;
  std::vector<double> vertex_y_error;
  std::vector<double> vertex_z_error;
  std::vector<std::vector<std::vector<std::vector<int>>>> vertex_indices;
  Double_t numvertices;

  auto v_numtracks  = integral_tree.Branch("NumVertices", &numvertices, "NumVertices/D");
  auto v_num_hits = integral_tree.Branch("Vertex_numTracks", "std::vector<double>", &vertex_numTracks, 32000, 99);
  auto v_t = integral_tree.Branch("Vertex_t", "std::vector<double>", &vertex_t, 32000, 99);
  auto v_x = integral_tree.Branch("Vertex_x", "std::vector<double>", &vertex_x, 32000, 99);
  auto v_y = integral_tree.Branch("Vertex_y", "std::vector<double>", &vertex_y, 32000, 99);
  auto v_z = integral_tree.Branch("Vertex_z", "std::vector<double>", &vertex_z, 32000, 99);
  auto v_t_error = integral_tree.Branch("Vertex_ErrorT", "std::vector<double>", &vertex_t_error, 32000, 99);
  auto v_x_error = integral_tree.Branch("Vertex_ErrorX", "std::vector<double>", &vertex_x_error, 32000, 99);
  auto v_y_error = integral_tree.Branch("Vertex_ErrorY", "std::vector<double>", &vertex_y_error, 32000, 99);
  auto v_z_error = integral_tree.Branch("Vertex_ErrorZ", "std::vector<double>", &vertex_z_error, 32000, 99);
  auto v_chi2 = integral_tree.Branch("Vertex_chi2", "std::vector<double>", &vertex_chi2, 32000, 99);
  auto v_chi2_per_dof = integral_tree.Branch("Vertex_chi2PerNdof", "std::vector<double>", &vertex_chi2_per_dof, 32000, 99);
  auto v_chi2_p_value = integral_tree.Branch("Vertex_chi2PValue", "std::vector<double>", &vertex_chi2_p_value, 32000, 99);
  auto v_indices = integral_tree.Branch("Vertex_indices", "std::vector<std::vector<std::vector<std::vector<int>>>>", &vertex_indices, 32000, 99);

  //__Make Track Branches_________________________________________________________________________
  std::vector<double> track_numHits;
  std::vector<double> track_chi2;
  std::vector<double> track_chi2_per_dof;
  std::vector<double> track_chi2_p_value;
  std::vector<double> track_beta;
  std::vector<double> track_beta_error;
  std::vector<double> track_angle;
  std::vector<double> track_angle_error;
  std::vector<double> unique_detector_count;
  std::vector<double> track_t;
  std::vector<double> track_x;
  std::vector<double> track_y;
  std::vector<double> track_z;
  std::vector<double> track_vx;
  std::vector<double> track_vy;
  std::vector<double> track_vz;
  std::vector<double> track_t_error;
  std::vector<double> track_x_error;
  std::vector<double> track_y_error;
  std::vector<double> track_z_error;
  std::vector<double> track_vx_error;
  std::vector<double> track_vy_error;
  std::vector<double> track_vz_error;
  std::vector<std::vector<double>> track_expected_hit_layer;
  std::vector<std::vector<double>> track_missing_hit_layer;
  std::vector<std::vector<std::vector<Int_t>>> track_indices;


  Double_t numtracks;

  auto t_numtracks  = integral_tree.Branch("NumTracks", &numtracks, "NumTracks/D");
  auto t_num_hits = integral_tree.Branch("Track_numHits", "std::vector<double>", &track_numHits, 32000, 99);
  auto t_t = integral_tree.Branch("Track_t0", "std::vector<double>", &track_t, 32000, 99);
  auto t_x = integral_tree.Branch("Track_x0", "std::vector<double>", &track_x, 32000, 99);
  auto t_y = integral_tree.Branch("Track_y0", "std::vector<double>", &track_y, 32000, 99);
  auto t_z = integral_tree.Branch("Track_z0", "std::vector<double>", &track_z, 32000, 99);
  auto t_vx = integral_tree.Branch("Track_velX", "std::vector<double>", &track_vx, 32000, 99);
  auto t_vy = integral_tree.Branch("Track_velY", "std::vector<double>", &track_vy, 32000, 99);
  auto t_vz = integral_tree.Branch("Track_velZ", "std::vector<double>", &track_vz, 32000, 99);
  auto t_t_error = integral_tree.Branch("Track_ErrorT0", "std::vector<double>", &track_t_error, 32000, 99);
  auto t_x_error = integral_tree.Branch("Track_ErrorX0", "std::vector<double>", &track_x_error, 32000, 99);
  auto t_y_error = integral_tree.Branch("Track_ErrorY0", "std::vector<double>", &track_y_error, 32000, 99);
  auto t_z_error = integral_tree.Branch("Track_ErrorZ0", "std::vector<double>", &track_z_error, 32000, 99);
  auto t_vx_error = integral_tree.Branch("Track_ErrorVx", "std::vector<double>", &track_vx_error, 32000, 99);
  auto t_vy_error = integral_tree.Branch("Track_ErrorVy", "std::vector<double>", &track_vy_error, 32000, 99);
  auto t_vz_error = integral_tree.Branch("Track_ErrorVz", "std::vector<double>", &track_vz_error, 32000, 99);
  auto t_chi2 = integral_tree.Branch("Track_chi2", "std::vector<double>", &track_chi2, 32000, 99);
  auto t_chi2_per_dof = integral_tree.Branch("Track_chi2PerNdof", "std::vector<double>", &track_chi2_per_dof, 32000, 99);
  auto t_chi2_p_value = integral_tree.Branch("Track_chi2PValue", "std::vector<double>", &track_chi2_p_value, 32000, 99);
  auto t_beta = integral_tree.Branch("Track_beta", "std::vector<double>", &track_beta, 32000, 99);
  auto t_beta_error = integral_tree.Branch("Track_ErrorBeta", "std::vector<double>", &track_beta_error, 32000, 99);
  auto t_angle = integral_tree.Branch("Track_angle", "std::vector<double>", &track_angle, 32000, 99);
  auto t_angle_error = integral_tree.Branch("Track_ErrorAngle", "std::vector<double>", &track_angle_error, 32000, 99);
  auto t_unique_detector_count = integral_tree.Branch("Track_detCount", "std::vector<double>", &unique_detector_count, 32000, 99);
  auto t_expected_hit_layer = integral_tree.Branch("Track_expectedHitLayer", "std::vector<std::vector<double>>", &track_expected_hit_layer, 32000, 99);
  auto t_missing_hit_layer = integral_tree.Branch("Track_missingHitLayer", "std::vector<std::vector<double>>", &track_missing_hit_layer, 32000, 99);
  auto t_indices = integral_tree.Branch("Track_indices", "std::vector<std::vector<std::vector<Int_t>>>", &track_indices, 32000, 99);

  //___Make Digi Branches_____________________________________________________________________
  std::vector<double> digi_hit_t;
  std::vector<double> digi_hit_x;
  std::vector<double> digi_hit_y;
  std::vector<double> digi_hit_z;
  std::vector<double> digi_hit_e;
  std::vector<double> digi_hit_px;
  std::vector<double> digi_hit_py;
  std::vector<double> digi_hit_pz;
  Int_t Digi_numHits;
  std::vector<std::vector<Int_t>> digi_hit_indices; ///

  auto branch_digi_num_hits  = integral_tree.Branch("Digi_numHits", &Digi_numHits, "Digi_numHits/I");
  auto branch_t  = integral_tree.Branch("Digi_time", "std::vector<double>", &digi_hit_t, 32000, 99);
  auto branch_x  = integral_tree.Branch("Digi_x", "std::vector<double>", &digi_hit_x, 32000, 99);
  auto branch_y  = integral_tree.Branch("Digi_y", "std::vector<double>", &digi_hit_y, 32000, 99);
  auto branch_z  = integral_tree.Branch("Digi_z", "std::vector<double>", &digi_hit_z, 32000, 99);
  auto branch_e  = integral_tree.Branch("Digi_energy", "std::vector<double>", &digi_hit_e, 32000, 99);
  auto branch_px = integral_tree.Branch("Digi_px", "std::vector<double>", &digi_hit_px, 32000, 99);
  auto branch_py = integral_tree.Branch("Digi_py", "std::vector<double>", &digi_hit_py, 32000, 99);
  auto branch_pz = integral_tree.Branch("Digi_pz", "std::vector<double>", &digi_hit_pz, 32000, 99);
  auto branch_indices = integral_tree.Branch("Digi_hitIndices", "std::vector<std::vector<Int_t>>", &digi_hit_indices, 32000, 99);

  //____________________________________________________________________________________________

 //___Make Sim Branches_________________________________________________________________________
 Double_t sim_numhits;
 std::vector<double> *sim_hit_e = nullptr; std::vector<double> **sim_hit_e_address = &sim_hit_e;
 std::vector<double> *sim_hit_t = nullptr;  std::vector<double> **sim_hit_t_address = &sim_hit_t;
 std::vector<double> *sim_hit_detId = nullptr;  std::vector<double> **sim_hit_detId_address = &sim_hit_detId;
 std::vector<double> *sim_hit_particlePdgId = nullptr;  std::vector<double> **sim_hit_particlePdgId_address = &sim_hit_particlePdgId;
 std::vector<double> *sim_hit_G4TrackId = nullptr;  std::vector<double> **sim_hit_G4TrackId_address = &sim_hit_G4TrackId;
 std::vector<double> *sim_hit_G4ParentTrackId = nullptr;  std::vector<double> **sim_hit_G4ParentTrackId_address = &sim_hit_G4ParentTrackId;
 std::vector<double> *sim_hit_x = nullptr;  std::vector<double> **sim_hit_x_address = &sim_hit_x;
 std::vector<double> *sim_hit_y = nullptr;  std::vector<double> **sim_hit_y_address = &sim_hit_y;
 std::vector<double> *sim_hit_z = nullptr;  std::vector<double> **sim_hit_z_address = &sim_hit_z;
 std::vector<double> *sim_hit_particleEnergy = nullptr;  std::vector<double> **sim_hit_particleEnergy_address = &sim_hit_particleEnergy;
 std::vector<double> *sim_hit_px = nullptr;  std::vector<double> **sim_hit_px_address = &sim_hit_px;
 std::vector<double> *sim_hit_py = nullptr;  std::vector<double> **sim_hit_py_address = &sim_hit_py;
 std::vector<double> *sim_hit_pz = nullptr;  std::vector<double> **sim_hit_pz_address = &sim_hit_pz;
 std::vector<double> *sim_hit_weight = nullptr;  std::vector<double> **sim_hit_weight_address = &sim_hit_weight;
 Double_t sim_NumGenParticles;
 std::vector<double> *sim_GenParticle_index = nullptr;  std::vector<double> **sim_GenParticle_index_address = &sim_GenParticle_index;
 std::vector<double> *sim_GenParticle_G4index = nullptr;  std::vector<double> **sim_GenParticle_G4index_address = &sim_GenParticle_G4index;
 std::vector<double> *sim_GenParticle_pdgid = nullptr;  std::vector<double> **sim_GenParticle_pdgid_address = &sim_GenParticle_pdgid;
 std::vector<double> *sim_GenParticle_status = nullptr;  std::vector<double> **sim_GenParticle_status_address = &sim_GenParticle_status;
 std::vector<double> *sim_GenParticle_time = nullptr;  std::vector<double> **sim_GenParticle_time_address = &sim_GenParticle_time;
 std::vector<double> *sim_GenParticle_x = nullptr;  std::vector<double> **sim_GenParticle_x_address = &sim_GenParticle_x;
 std::vector<double> *sim_GenParticle_y = nullptr;  std::vector<double> **sim_GenParticle_y_address = &sim_GenParticle_y;
 std::vector<double> *sim_GenParticle_z = nullptr;  std::vector<double> **sim_GenParticle_z_address = &sim_GenParticle_z;
 std::vector<double> *sim_GenParticle_energy = nullptr;  std::vector<double> **sim_GenParticle_energy_address = &sim_GenParticle_energy;
 std::vector<double> *sim_GenParticle_px = nullptr;  std::vector<double> **sim_GenParticle_px_address = &sim_GenParticle_px;
 std::vector<double> *sim_GenParticle_py = nullptr;  std::vector<double> **sim_GenParticle_py_address = &sim_GenParticle_py;
 std::vector<double> *sim_GenParticle_pz = nullptr;  std::vector<double> **sim_GenParticle_pz_address = &sim_GenParticle_pz;
 std::vector<double> *sim_GenParticle_mo1 = nullptr;  std::vector<double> **sim_GenParticle_mo1_address = &sim_GenParticle_mo1;
 std::vector<double> *sim_GenParticle_mo2 = nullptr;  std::vector<double> **sim_GenParticle_mo2_address = &sim_GenParticle_mo2;
 std::vector<double> *sim_GenParticle_dau1 = nullptr;  std::vector<double> **sim_GenParticle_dau1_address = &sim_GenParticle_dau1;
 std::vector<double> *sim_GenParticle_dau2 = nullptr;  std::vector<double> **sim_GenParticle_dau2_address = &sim_GenParticle_dau2;
 std::vector<double> *sim_GenParticle_mass = nullptr;  std::vector<double> **sim_GenParticle_mass_address = &sim_GenParticle_mass;
 std::vector<double> *sim_GenParticle_pt = nullptr;  std::vector<double> **sim_GenParticle_pt_address = &sim_GenParticle_pt;
 std::vector<double> *sim_GenParticle_eta = nullptr;  std::vector<double> **sim_GenParticle_eta_address = &sim_GenParticle_eta;
 std::vector<double> *sim_GenParticle_phi = nullptr;  std::vector<double> **sim_GenParticle_phi_address = &sim_GenParticle_phi;
 std::vector<double> *sim_COSMIC_EVENT_ID = nullptr;  std::vector<double> **sim_COSMIC_EVENT_ID_address = &sim_COSMIC_EVENT_ID;
 std::vector<double> *sim_COSMIC_CORE_X = nullptr;  std::vector<double> **sim_COSMIC_CORE_X_address = &sim_COSMIC_CORE_X;
 std::vector<double> *sim_COSMIC_CORE_Y = nullptr;  std::vector<double> **sim_COSMIC_CORE_Y_address = &sim_COSMIC_CORE_Y;
 std::vector<double> *sim_COSMIC_GEN_PRIMARY_ENERGY = nullptr;  std::vector<double> **sim_COSMIC_GEN_PRIMARY_ENERGY_address = &sim_COSMIC_GEN_PRIMARY_ENERGY;
 std::vector<double> *sim_COSMIC_GEN_THETA = nullptr;  std::vector<double> **sim_COSMIC_GEN_THETA_address = &sim_COSMIC_GEN_THETA;
 std::vector<double> *sim_COSMIC_GEN_PHI = nullptr;  std::vector<double> **sim_COSMIC_GEN_PHI_address = &sim_COSMIC_GEN_PHI;
 std::vector<double> *sim_COSMIC_GEN_FIRST_HEIGHT = nullptr;  std::vector<double> **sim_COSMIC_GEN_FIRST_HEIGHT_address = &sim_COSMIC_GEN_FIRST_HEIGHT;
 std::vector<double> *sim_COSMIC_GEN_ELECTRON_COUNT = nullptr;  std::vector<double> **sim_COSMIC_GEN_ELECTRON_COUNT_address = &sim_COSMIC_GEN_ELECTRON_COUNT;
 std::vector<double> *sim_COSMIC_GEN_MUON_COUNT = nullptr;  std::vector<double> **sim_COSMIC_GEN_MUON_COUNT_address = &sim_COSMIC_GEN_MUON_COUNT;
 std::vector<double> *sim_COSMIC_GEN_HADRON_COUNT = nullptr;  std::vector<double> **sim_COSMIC_GEN_HADRON_COUNT_address = &sim_COSMIC_GEN_HADRON_COUNT;
 std::vector<double> *sim_COSMIC_GEN_PRIMARY_ID = nullptr; std::vector<double> **sim_COSMIC_GEN_PRIMARY_ID_address = &sim_COSMIC_GEN_PRIMARY_ID;
 std::vector<double> *sim_EXTRA_11 = nullptr; std::vector<double> **sim_EXTRA_11_address = &sim_EXTRA_11;
 std::vector<double> *sim_EXTRA_12 = nullptr; std::vector<double> **sim_EXTRA_12_address = &sim_EXTRA_12;
 std::vector<double> *sim_EXTRA_13 = nullptr; std::vector<double> **sim_EXTRA_13_address = &sim_EXTRA_13;
 std::vector<double> *sim_EXTRA_14 = nullptr; std::vector<double> **sim_EXTRA_14_address = &sim_EXTRA_14;
 std::vector<double> *sim_EXTRA_15 = nullptr; std::vector<double> **sim_EXTRA_15_address = &sim_EXTRA_15;

 integral_tree.Branch("NumHits", &sim_numhits);
 integral_tree.Branch("Hit_energy", sim_hit_e_address);
 integral_tree.Branch("Hit_time", sim_hit_t_address);
 integral_tree.Branch("Hit_detId", sim_hit_detId_address);
 integral_tree.Branch("Hit_particlePdgId", sim_hit_particlePdgId_address);
 integral_tree.Branch("Hit_G4TrackId", sim_hit_G4TrackId_address);
 integral_tree.Branch("Hit_G4ParentTrackId", sim_hit_G4ParentTrackId_address);
 integral_tree.Branch("Hit_x",sim_hit_x_address);
 integral_tree.Branch("Hit_y",sim_hit_y_address);
 integral_tree.Branch("Hit_z",sim_hit_z_address);
 integral_tree.Branch("Hit_particleEnergy",sim_hit_particleEnergy_address);
 integral_tree.Branch("Hit_particlePx",sim_hit_px_address);
 integral_tree.Branch("Hit_particlePy",sim_hit_py_address);
 integral_tree.Branch("Hit_particlePz",sim_hit_pz_address);
 integral_tree.Branch("Hit_weight",sim_hit_weight_address);
 integral_tree.Branch("NumGenParticles", &sim_NumGenParticles);
 integral_tree.Branch("GenParticle_index",sim_GenParticle_index_address);
 integral_tree.Branch("GenParticle_G4index",sim_GenParticle_G4index_address);
 integral_tree.Branch("GenParticle_pdgid",sim_GenParticle_pdgid_address);
 integral_tree.Branch("GenParticle_status",sim_GenParticle_status_address);
 integral_tree.Branch("GenParticle_time",sim_GenParticle_time_address);
 integral_tree.Branch("GenParticle_x",sim_GenParticle_x_address);
 integral_tree.Branch("GenParticle_y",sim_GenParticle_y_address);
 integral_tree.Branch("GenParticle_z",sim_GenParticle_z_address);
 integral_tree.Branch("GenParticle_energy",sim_GenParticle_energy_address);
 integral_tree.Branch("GenParticle_px", sim_GenParticle_px_address);
 integral_tree.Branch("GenParticle_py", sim_GenParticle_py_address);
 integral_tree.Branch("GenParticle_pz", sim_GenParticle_pz_address);
 integral_tree.Branch("GenParticle_mo1", sim_GenParticle_mo1_address);
 integral_tree.Branch("GenParticle_mo2", sim_GenParticle_mo2_address);
 integral_tree.Branch("GenParticle_dau1", sim_GenParticle_dau1_address);
 integral_tree.Branch("GenParticle_dau2", sim_GenParticle_dau2_address);
 integral_tree.Branch("GenParticle_mass", sim_GenParticle_mass_address);
 integral_tree.Branch("GenParticle_pt", sim_GenParticle_pt_address);
 integral_tree.Branch("GenParticle_eta", sim_GenParticle_eta_address);
 integral_tree.Branch("GenParticle_phi", sim_GenParticle_phi_address);
 integral_tree.Branch("COSMIC_EVENT_ID", sim_COSMIC_EVENT_ID_address);
 integral_tree.Branch("COSMIC_CORE_X", sim_COSMIC_CORE_X_address);
 integral_tree.Branch("COSMIC_CORE_Y", sim_COSMIC_CORE_Y_address);
 integral_tree.Branch("COSMIC_GEN_PRIMARY_ENERGY", sim_COSMIC_GEN_PRIMARY_ENERGY_address);
 integral_tree.Branch("COSMIC_GEN_THETA", sim_COSMIC_GEN_THETA_address);
 integral_tree.Branch("COSMIC_GEN_PHI", sim_COSMIC_GEN_PHI_address);
 integral_tree.Branch("COSMIC_GEN_FIRST_HEIGHT", sim_COSMIC_GEN_FIRST_HEIGHT_address);
 integral_tree.Branch("COSMIC_GEN_ELECTRON_COUNT", sim_COSMIC_GEN_ELECTRON_COUNT_address);
 integral_tree.Branch("COSMIC_GEN_MUON_COUNT", sim_COSMIC_GEN_MUON_COUNT_address);
 integral_tree.Branch("COSMIC_GEN_HADRON_COUNT", sim_COSMIC_GEN_HADRON_COUNT_address);
 integral_tree.Branch("COSMIC_GEN_PRIMARY_ID", sim_COSMIC_GEN_PRIMARY_ID_address);
 integral_tree.Branch("EXTRA_11", sim_EXTRA_11_address);
 integral_tree.Branch("EXTRA_12", sim_EXTRA_12_address);
 integral_tree.Branch("EXTRA_13", sim_EXTRA_13_address);
 integral_tree.Branch("EXTRA_14", sim_EXTRA_14_address);
 integral_tree.Branch("EXTRA_15", sim_EXTRA_15_address);
 //_____________________________________________________________________________________________

  for (std::size_t event_counter{}; event_counter < import_size; ++event_counter) {

    for (const auto path : paths) {
        auto sim_file = TFile::Open(path.c_str(), "READ");
        if (sim_file && !sim_file->IsZombie()) {
            sim_file->cd();
            TTree* sim_tree;
            sim_file->GetObject("box_run", sim_tree);
            sim_tree->SetBranchAddress("NumHits", &sim_numhits);
            sim_tree->SetBranchAddress("Hit_energy", &sim_hit_e);
            sim_tree->SetBranchAddress("Hit_time", &sim_hit_t);
            sim_tree->SetBranchAddress("Hit_detId", &sim_hit_detId);
            sim_tree->SetBranchAddress("Hit_particlePdgId", &sim_hit_particlePdgId);
            sim_tree->SetBranchAddress("Hit_G4TrackId", &sim_hit_G4TrackId);
            sim_tree->SetBranchAddress("Hit_G4ParentTrackId", &sim_hit_G4ParentTrackId);
            sim_tree->SetBranchAddress("Hit_x", &sim_hit_x);
            sim_tree->SetBranchAddress("Hit_y", &sim_hit_y);
            sim_tree->SetBranchAddress("Hit_z", &sim_hit_z);
            sim_tree->SetBranchAddress("Hit_particleEnergy", &sim_hit_particleEnergy);
            sim_tree->SetBranchAddress("Hit_particlePx", &sim_hit_px);
            sim_tree->SetBranchAddress("Hit_particlePy", &sim_hit_py);
            sim_tree->SetBranchAddress("Hit_particlePz", &sim_hit_pz);
            sim_tree->SetBranchAddress("Hit_weight", &sim_hit_weight);
            sim_tree->SetBranchAddress("NumGenParticles", &sim_NumGenParticles);
            sim_tree->SetBranchAddress("GenParticle_index", &sim_GenParticle_index);
            sim_tree->SetBranchAddress("GenParticle_G4index", &sim_GenParticle_G4index);
            sim_tree->SetBranchAddress("GenParticle_pdgid", &sim_GenParticle_pdgid);
            sim_tree->SetBranchAddress("GenParticle_status", &sim_GenParticle_status);
            sim_tree->SetBranchAddress("GenParticle_time", &sim_GenParticle_time);
            sim_tree->SetBranchAddress("GenParticle_x", &sim_GenParticle_x);
            sim_tree->SetBranchAddress("GenParticle_y", &sim_GenParticle_y);
            sim_tree->SetBranchAddress("GenParticle_z", &sim_GenParticle_z);
            sim_tree->SetBranchAddress("GenParticle_energy", &sim_GenParticle_energy);
            sim_tree->SetBranchAddress("GenParticle_px", &sim_GenParticle_px);
            sim_tree->SetBranchAddress("GenParticle_py", &sim_GenParticle_py);
            sim_tree->SetBranchAddress("GenParticle_pz", &sim_GenParticle_pz);
            sim_tree->SetBranchAddress("GenParticle_mo1", &sim_GenParticle_mo1);
            sim_tree->SetBranchAddress("GenParticle_mo2", &sim_GenParticle_mo2);
            sim_tree->SetBranchAddress("GenParticle_dau1", &sim_GenParticle_dau1);
            sim_tree->SetBranchAddress("GenParticle_dau2", &sim_GenParticle_dau2);
            sim_tree->SetBranchAddress("GenParticle_mass", &sim_GenParticle_mass);
            sim_tree->SetBranchAddress("GenParticle_pt", &sim_GenParticle_pt);
            sim_tree->SetBranchAddress("GenParticle_eta", &sim_GenParticle_eta);
            sim_tree->SetBranchAddress("GenParticle_phi", &sim_GenParticle_phi);
            sim_tree->SetBranchAddress("COSMIC_EVENT_ID", &sim_COSMIC_EVENT_ID);
            sim_tree->SetBranchAddress("COSMIC_CORE_X", &sim_COSMIC_CORE_X);
            sim_tree->SetBranchAddress("COSMIC_CORE_Y", &sim_COSMIC_CORE_Y);
            sim_tree->SetBranchAddress("COSMIC_GEN_PRIMARY_ENERGY", &sim_COSMIC_GEN_PRIMARY_ENERGY);
            sim_tree->SetBranchAddress("COSMIC_GEN_THETA", &sim_COSMIC_GEN_THETA);
            sim_tree->SetBranchAddress("COSMIC_GEN_PHI", &sim_COSMIC_GEN_PHI);
            sim_tree->SetBranchAddress("COSMIC_GEN_FIRST_HEIGHT", &sim_COSMIC_GEN_FIRST_HEIGHT);
            sim_tree->SetBranchAddress("COSMIC_GEN_ELECTRON_COUNT", &sim_COSMIC_GEN_ELECTRON_COUNT);
            sim_tree->SetBranchAddress("COSMIC_GEN_MUON_COUNT", &sim_COSMIC_GEN_MUON_COUNT);
            sim_tree->SetBranchAddress("COSMIC_GEN_HADRON_COUNT", &sim_COSMIC_GEN_HADRON_COUNT);
            sim_tree->SetBranchAddress("COSMIC_GEN_PRIMARY_ID", &sim_COSMIC_GEN_PRIMARY_ID);
            sim_tree->SetBranchAddress("EXTRA_11", &sim_EXTRA_11);
            sim_tree->SetBranchAddress("EXTRA_12", &sim_EXTRA_12);
            sim_tree->SetBranchAddress("EXTRA_13", &sim_EXTRA_13);
            sim_tree->SetBranchAddress("EXTRA_14", &sim_EXTRA_14);
            sim_tree->SetBranchAddress("EXTRA_15", &sim_EXTRA_15);
            sim_tree->GetEntry(event_counter);
        }
		sim_file->Close();
    }

    const auto digitized_full_event = analysis::full_digi_event<box::geometry>(imported_events[event_counter], energy_events[event_counter], complete_events[event_counter]);

    const auto digi_event = analysis::add_digi_event<box::geometry>(imported_events[event_counter], energy_events[event_counter], complete_events[event_counter]);

	const auto event = analysis::add_width<box::geometry>(digi_event);
	const auto event_size = event.size();
	const auto event_counter_string = std::to_string(event_counter);

	const auto compressed_event_t = options.time_smearing ? mc::time_smear<box::geometry>(mc::compress<box::geometry>(event))
	   	                                                  : mc::compress<box::geometry>(event);

	const auto compressed_event = options.positionz_smearing ? mc::positionz_smear<box::geometry>(compressed_event_t)
		                                                     : compressed_event_t;

    const auto compression_size = event_size / static_cast<type::real>(compressed_event.size());

    const auto altered_event = alter_event(compressed_event, options, extension);

    const auto event_density = box::geometry::event_density(altered_event);
    box::io::print_event_summary(event_counter, altered_event.size(), compression_size, event_density);

    std::string canvas_title{"event"};
    for (const auto& path : paths) {
      canvas_title.push_back(' ');
      canvas_title.insert(canvas_title.cend(), path.begin(), path.end());
    }

    plot::canvas canvas("event" + event_counter_string,
                        canvas_title + " | " + event_counter_string,
                        !options.verbose_output);
    if (options.draw_events) {
      box::io::draw_detector(canvas, extension.layer_count);
      box::io::draw_mc_tracks(canvas, mc::convert_events(mc_imported_events[event_counter]));
      for (const auto& hit : altered_event)
        canvas.add_point(type::reduce_to_r3(hit), 0.8, plot::color::BLACK);
    }

    analysis::full_event non_primary_track_points, non_secondary_track_points;
    auto tracks = find_tracks(altered_event,
                              options,
                              20.0L,
                              1UL,
                              canvas,
                              non_primary_track_points);
    auto secondary_tracks = find_tracks(non_primary_track_points,
                                        options,
                                        100.0L,
                                        0UL,
                                        canvas,
                                        non_secondary_track_points);
    tracks.reserve(tracks.size() + secondary_tracks.size());
    tracks.insert(tracks.cend(),
                  std::make_move_iterator(secondary_tracks.cbegin()),
                  std::make_move_iterator(secondary_tracks.cend()));


//___Fill Digi varibales_________________________________________________________________________
    digi_hit_t.clear();
    digi_hit_x.clear();
    digi_hit_y.clear();
    digi_hit_z.clear();
    digi_hit_e.clear();
    digi_hit_px.clear();
    digi_hit_py.clear();
    digi_hit_pz.clear();
    digi_hit_indices.clear();

    for (const auto& h : digitized_full_event) {
        digi_hit_e.push_back(h.e / units::energy);
        digi_hit_px.push_back(h.px / units::momentum);
        digi_hit_py.push_back(h.py / units::momentum);
        digi_hit_pz.push_back(h.pz / units::momentum);
        digi_hit_indices.push_back(h.sim_indices);
    }


    Digi_numHits = digi_hit_e.size();

    for (const auto& h : compressed_event) {
        digi_hit_t.push_back(h.t / units::time);
        digi_hit_x.push_back(h.x / units::length);
        digi_hit_y.push_back(h.y / units::length);
        digi_hit_z.push_back(h.z / units::length);
    }
//______________________________________________________________________________________________

//___Fill Track varibales_______________________________________________________________________
// box::io::save_tracks(tracks, canvas, track_tree, options);
    for (const auto& track : tracks) {
        if (options.verbose_output)
            std::cout << track << "\n";
        if (options.draw_events)
            box::io::draw_track(canvas, track);
    }

    track_t.clear();
    track_x.clear();
    track_y.clear();
    track_z.clear();
    track_vx.clear();
    track_vy.clear();
    track_vz.clear();
    track_t_error.clear();
    track_x_error.clear();
    track_y_error.clear();
    track_z_error.clear();
    track_vx_error.clear();
    track_vy_error.clear();
    track_vz_error.clear();
    track_numHits.clear();
    track_chi2.clear();
    track_chi2_per_dof.clear();
    track_chi2_p_value.clear();
    track_beta.clear();
    track_beta_error.clear();
    track_angle.clear();
    track_angle_error.clear();
    unique_detector_count.clear();
	track_expected_hit_layer.clear();
	track_missing_hit_layer.clear();
  track_indices.clear();

    track_t.reserve(tracks.size());
    track_x.reserve(tracks.size());
    track_y.reserve(tracks.size());
    track_z.reserve(tracks.size());
    track_vx.reserve(tracks.size());
    track_vy.reserve(tracks.size());
    track_vz.reserve(tracks.size());
    track_t_error.reserve(tracks.size());
    track_x_error.reserve(tracks.size());
    track_y_error.reserve(tracks.size());
    track_z_error.reserve(tracks.size());
    track_vx_error.reserve(tracks.size());
    track_vy_error.reserve(tracks.size());
    track_vz_error.reserve(tracks.size());
    track_chi2.reserve(tracks.size());
    track_numHits.reserve(tracks.size());
    track_chi2_per_dof.reserve(tracks.size());
    track_chi2_p_value.reserve(tracks.size());
    track_beta.reserve(tracks.size());
    track_beta_error.reserve(tracks.size());
    track_angle.reserve(tracks.size());
    track_angle_error.reserve(tracks.size());
    unique_detector_count.reserve(tracks.size());
	track_expected_hit_layer.reserve(tracks.size());
	track_missing_hit_layer.reserve(tracks.size());
  track_indices.reserve(tracks.size());

    if (event_density >= options.event_density_limit) {
        for (const auto& track : tracks) {
            if (track.fit_converged()) {
                track_numHits.push_back(track.size());
                track_chi2.push_back(track.chi_squared());
                track_chi2_per_dof.push_back(track.chi_squared_per_dof());
                track_chi2_p_value.push_back(track.chi_squared_p_value());
                track_beta.push_back(track.beta());
                track_beta_error.push_back(track.beta_error());
                track_angle.push_back(track.angle());
                track_angle_error.push_back(track.angle_error());
                // for (const auto& point : track) {
                //         event_t.get().push_back(point.t / units::time);
                //         event_x.get().push_back(point.x / units::length);
                //         event_y.get().push_back(point.y / units::length);
                //         event_z.get().push_back(point.z / units::length);
                //         event_detector.get().push_back(geometry::volume(reduce_to_r3(point)));
                // }

                geometry::structure_vector unique_geometry;
                const auto detectors = track.detectors();
                std::unique_copy(detectors.cbegin(), detectors.cend(), std::back_inserter(unique_geometry));
                unique_detector_count.push_back(unique_geometry.size());

                track_t.push_back(track.final_fit().t0.value / units::time);
                track_x.push_back(track.final_fit().x0.value / units::length);
                track_y.push_back(track.final_fit().y0.value / units::length);
                track_z.push_back(track.final_fit().z0.value / units::length);
                track_vx.push_back(track.final_fit().vx.value / units::velocity);
                track_vy.push_back(track.final_fit().vy.value / units::velocity);
                track_vz.push_back(track.final_fit().vz.value / units::velocity);
                track_t_error.push_back(track.final_fit().t0.error / units::time);
                track_x_error.push_back(track.final_fit().x0.error / units::length);
                track_y_error.push_back(track.final_fit().y0.error / units::length);
                track_z_error.push_back(track.final_fit().z0.error / units::length);
                track_vx_error.push_back(track.final_fit().vx.error / units::velocity);
                track_vy_error.push_back(track.final_fit().vy.error / units::velocity);
                track_vz_error.push_back(track.final_fit().vz.error / units::velocity);
                track_indices.push_back(track.final_fit().track_indices);
                //td::cout << "hello worldddddddddddddddddddddd:" << track_indices << '\n';

				std::vector<double> hits_y;
				hits_y.clear();
				hits_y.reserve(track.event().size());
				for (const auto& hits : track.event()) {
					hits_y.push_back(hits.y / units::length);
				}

				std::vector<double> missing_hit_layer;
                missing_hit_layer.clear();
				int missing_hit_counter = 0;
				missing_hit_layer.reserve(missing_hit_counter);

				std::vector<double> expected_hit_layer;
                expected_hit_layer.clear();
                int expected_hit_counter = 0;
                expected_hit_layer.reserve(expected_hit_counter);

				// std::cout << "//////////////////////////////////////////////" << std::endl;
				//				std::cout << "digi_event : " << digi_event << std::endl;
				// std::cout << "compressed_event : " << compressed_event << std::endl;
				// std::cout << "track.event(): " << track.event() << std::endl;

                for (const auto& determined_y : layer_center) { //per layer
					if (track.final_fit().vy.value > 0) {
						const auto point_dt = ((determined_y - (track.final_fit().y0.value / units::length)) / (track.final_fit().vy.value / units::velocity)) + (track.final_fit().t0.value / units::time);
                        const auto point_t = ((determined_y - (track.final_fit().y0.value / units::length)) / (track.final_fit().vy.value / units::velocity));
						const auto point_x = (track.final_fit().x0.value / units::length) + (point_t*(track.final_fit().vx.value / units::velocity));
						const auto point_y = determined_y;
						const auto point_z = (track.final_fit().z0.value / units::length) + (point_t*(track.final_fit().vz.value / units::velocity));
						//						std::cout << "determined_y: " << determined_y << " track_y: " << track.final_fit().y0.value/ units::length << " point_t: " << point_t << " point_x: " << point_x << " point_y: " << point_y << " point_z: " << point_z << "\n";


						for (auto x_limit : detector_x_limits) {
							for (auto z_limit : detector_z_limits) {
								if (point_x < x_limit[1] && point_x > x_limit[0] && point_z < z_limit[1] && point_z > z_limit[0]) { // is there a detector for this point
									//expected hit in this layer
                                    ++expected_hit_counter;
									expected_hit_layer.push_back(determined_y);
									// std::cout << "point_t: " << point_dt << " point_x: " << point_x << " point_y: " << point_y << " point_z: " << point_z << std::endl;
									if (std::count(hits_y.begin(), hits_y.end(), determined_y)==0) { // if there is no hit in that layer for this track
										//missing hit
                                        ++missing_hit_counter;
										missing_hit_layer.push_back(determined_y);
										// std::cout << "missing in layer: " << determined_y << std::endl;
									}
								}
							}
						}
					}
                }

				track_missing_hit_layer.push_back(missing_hit_layer);
                track_expected_hit_layer.push_back(expected_hit_layer);

			}
		}
		numtracks = tracks.size();

	}
 //______________________________________________________________________________________________

 //___Fill Vertex varibales______________________________________________________________________

 ////////////////////////////////////////////////////////////////////////////////////////////////
    box::io::print_tracking_summary(event, tracks);

    analysis::track_vector converged_tracks;
    converged_tracks.reserve(tracks.size());
    std::copy_if(std::make_move_iterator(tracks.begin()),
                 std::make_move_iterator(tracks.end()),
                 std::back_inserter(converged_tracks),
                 [](const auto& track) { return track.fit_converged(); });

////////////////////////////////////////////////////////////////////////////////////////////////
    // box::io::save_vertices(analysis::pairwise_fit_tracks(converged_tracks), canvas, vertex_tree, options);
    const analysis::vertex_vector& vertices = analysis::pairwise_fit_tracks(converged_tracks);

      for (const auto& vertex : vertices) {
        if (options.verbose_output)
          std::cout << vertex << "\n";
        if (options.draw_events)
          box::io::draw_vertex_and_guess(canvas, vertex);
      }

      vertex_t.clear();
      vertex_x.clear();
      vertex_y.clear();
      vertex_z.clear();
      vertex_t_error.clear();
      vertex_x_error.clear();
      vertex_y_error.clear();
      vertex_z_error.clear();
      vertex_numTracks.clear();
      vertex_chi2.clear();
      vertex_chi2_per_dof.clear();
      vertex_chi2_p_value.clear();
      vertex_indices.clear();

      vertex_t.reserve(vertices.size());
      vertex_x.reserve(vertices.size());
      vertex_y.reserve(vertices.size());
      vertex_z.reserve(vertices.size());
      vertex_t_error.reserve(vertices.size());
      vertex_x_error.reserve(vertices.size());
      vertex_y_error.reserve(vertices.size());
      vertex_z_error.reserve(vertices.size());
      vertex_numTracks.reserve(vertices.size());
      vertex_chi2.reserve(vertices.size());
      vertex_chi2_per_dof.reserve(vertices.size());
      vertex_chi2_p_value.reserve(vertices.size());
      vertex_indices.reserve(vertices.size());

	  if (event_density >= options.event_density_limit) {
		  for (const auto& vertex : vertices) {
			  if (vertex.fit_converged()) {
				  vertex_numTracks.push_back(vertex.size());
				  vertex_chi2.push_back(vertex.chi_squared());
				  vertex_chi2_per_dof.push_back(vertex.chi_squared_per_dof());
				  vertex_chi2_p_value.push_back(vertex.chi_squared_p_value());
				  vertex_t.push_back(vertex.final_fit().t.value / units::time);
				  vertex_x.push_back(vertex.final_fit().x.value / units::length);
				  vertex_y.push_back(vertex.final_fit().y.value / units::length);
				  vertex_z.push_back(vertex.final_fit().z.value / units::length);
				  vertex_t_error.push_back(vertex.final_fit().t.error / units::time);
				  vertex_x_error.push_back(vertex.final_fit().x.error / units::length);
				  vertex_y_error.push_back(vertex.final_fit().y.error / units::length);
				  vertex_z_error.push_back(vertex.final_fit().z.error / units::length);
          vertex_indices.push_back(vertex.final_fit().vertex_indices);
			  }
		  }
		  numvertices = vertices.size();
	  }
//______________________________________________________________________________________________
    integral_tree.Fill();
    canvas.draw();


} //close event for loop


  plot::save_all(save_path,
    filetype_tag,
    plot::value_tag{"TIMESTAMP", util::time::GetString("%c %Z")},
    project_tag,
    plot::value_tag{"EVENTS", std::to_string(import_size)},
    plot::value_tag{"EFFICIENCY", std::to_string(options.simulated_efficiency)},
    plot::value_tag{"NOISE", std::to_string(options.simulated_noise_rate * units::time) + " / " + units::time_string},
    box::geometry::value_tags(),
    box::io::data_paths_value_tags(paths, options.data_timing_offsets));



  const auto path_count = paths.size();
  std::vector<std::string> prefixes;
  prefixes.reserve(path_count);
  for (std::size_t i{}; i < path_count; ++i)
      prefixes.push_back("SIM_" + std::to_string(i) + "_");



   auto file = TFile::Open(save_path.c_str(), "UPDATE");
   if (file && !file->IsZombie()) {
       file->cd();
       integral_tree.Write();
       file->Close();
   }


}
//----------------------------------------------------------------------------------------------

//__Box Tracking Algorithm______________________________________________________________________
int box_tracking(int argc,
                 char* argv[]) {
  box::io::extension_parser extension{};
  const auto options = script::parse_command_line(argc, argv, extension);

  plot::init(options.draw_events);
  geometry::open(options.geometry_file, options.default_time_error);
  box::geometry::import(extension);

  box::io::print_tracking_directories(options.data_directories);

  std::size_t path_counter{};
  const auto statistics_path_prefix = box::io::add_statistics_path(options);
  for (const auto& paths : reader::root::transpose_search_directories(options.data_directories, options.data_file_extension)) {
    const auto path_counter_string = std::to_string(path_counter++);
    const auto statistics_save_path = statistics_path_prefix
                                    + path_counter_string
                                    + "."
                                    + options.statistics_file_extension;

    box::io::print_bar();
    box::io::print_tracking_paths(paths);

    track_event_bundle(paths,
      reader::root::import_full_event_mc_bundle(paths,
          options.data_timing_offsets,
          options.data_track_id_key,
          options.data_t_key,
          options.data_x_key,
          options.data_y_key,
          options.data_z_key,
          options.data_e_key,
          options.data_px_key,
          options.data_py_key,
          options.data_pz_key,
          options.data_detector_key),
      options,
      extension,
      statistics_save_path);

    std::cout << "\n\nSaved: " << statistics_save_path << "\n";
  }

  box::io::print_bar();
  geometry::close();
  plot::end();
  std::cout << "Completed.\n";
  return 0;
}
//----------------------------------------------------------------------------------------------

//__Silent Box Tracking Algorithm_______________________________________________________________
int silent_box_tracking(int argc,
                        char* argv[]) {
  util::io::remove_buffer(std::cout, std::cerr, std::clog);
  return box_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------



} /* namespace MATHUSLA */

//__Main Function: Box Tracker__________________________________________________________________
int main(int argc,
         char* argv[]) {
  return MATHUSLA::box_tracking(argc, argv);
}
//----------------------------------------------------------------------------------------------
