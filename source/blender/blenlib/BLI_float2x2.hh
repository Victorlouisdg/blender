/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#pragma once

#include "BLI_float2.hh"
#include "BLI_math_matrix.h"

namespace blender {

struct float2x2 {
  float values[4][4];

  float2x2() = default;

  friend float2x2 operator*(const float2x2 &a, const float2 &b)
  {
    float2x2 result;
    // mul_m4_m4m4(result.values, a.values, b.values);
    return result;
  }

  float2x2 inverted() const
  {
    float2x2 result;
    // invert_m4_m4(result.values, values);
    return result;
  }
};

}  // namespace blender
