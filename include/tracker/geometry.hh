/*
 * include/util/geometry.hh
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

#ifndef TRACKER__GEOMETRY_HH
#define TRACKER__GEOMETRY_HH
#pragma once

#include <tracker/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace geometry { ///////////////////////////////////////////////////////////////////////////

using namespace type;

//__Geometry Navigation System__________________________________________________________________
void open(const std::string& path);
void close();
//----------------------------------------------------------------------------------------------

//__Current Geometry File_______________________________________________________________________
const std::string& current_geometry_path();
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes________________________________________________________________
const std::vector<std::string> full_structure();
//----------------------------------------------------------------------------------------------

//__List of all Geometry Volumes Except those in the List_______________________________________
const std::vector<std::string> full_structure_except(const std::vector<std::string>& names);
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy of a Volume________________________________________________________________
// TODO: finish
const std::vector<std::string> volume_hierarchy(const std::string& name);
//----------------------------------------------------------------------------------------------

//__Volume Containment Check____________________________________________________________________
bool is_inside_volume(const r3_point& point,
                      const std::string& name);
bool is_inside_volume(const r4_point& point,
                      const std::string& name);
//----------------------------------------------------------------------------------------------

//__Volume Hierarchy Search_____________________________________________________________________
const std::vector<std::string> volume_hierarchy(const r3_point& point);
const std::vector<std::string> volume_hierarchy(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Volume Search_______________________________________________________________________________
const std::string volume(const r3_point& point);
const std::string volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Box Volume__________________________________________________________________________________
struct box_volume { r3_point center, min, max; };
using box_volume_vector = std::vector<box_volume>;
//----------------------------------------------------------------------------------------------

//__Limit Box of a Volume_______________________________________________________________________
const box_volume limits_of(const std::string& name);
//----------------------------------------------------------------------------------------------

//__Compute Intersection of Volumes then Find Bounding Box______________________________________
// TODO: finish
const box_volume intersection_volume(const std::string& first,
                                     const std::string& second);
//----------------------------------------------------------------------------------------------

//__Compute Union of Volumes then Find Bounding Box_____________________________________________
// TODO: finish
const box_volume union_volume(const std::string& first,
                              const std::string& second);
//----------------------------------------------------------------------------------------------

//__Compute Intersection of Box Volumes_________________________________________________________
const box_volume intersection_volume(const box_volume& first,
                                     const box_volume& second);
//----------------------------------------------------------------------------------------------

//__Compute Union of Box Volumes________________________________________________________________
const box_volume union_volume(const box_volume& first,
                              const box_volume& second);
//----------------------------------------------------------------------------------------------

//__Limit Box of the Volume with Point__________________________________________________________
const box_volume limits_of_volume(const r3_point& point);
const box_volume limits_of_volume(const r4_point& point);
//----------------------------------------------------------------------------------------------

//__Box Volume Containment Check________________________________________________________________
constexpr bool is_inside_volume(const r3_point& point,
                                const box_volume& box);
constexpr bool is_inside_volume(const r4_point& point,
                                const box_volume& box);
//----------------------------------------------------------------------------------------------

//__Find Center of Geometry around Point________________________________________________________
const r3_point find_center(const r3_point& point);
const r4_point find_center(const r4_point& point);
//----------------------------------------------------------------------------------------------

} /* namespace geometry */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__GEOMETRY_HH */
