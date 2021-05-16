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

#include "BLI_float3.hh"
#include "BLI_math_matrix.h"

namespace blender {

struct float3x3 {
  float values[3][3];

  static float3x3 identity()
  {
    float3x3 mat;
    unit_m3(mat.values);
    return mat;
  }

  friend float3 operator*(const float3x3 &m, const float3 &v)
  {
    float3 result;
    mul_v3_m3v3(result, m.values, v);
    return result;
  }

  friend float3 operator*(const float3x3 &m, const float (*v)[3])
  {
    return m * float3(v);
  }

  float3x3 inverted() const
  {
    float3x3 result;
    invert_m3_m3(result.values, values);
    return result;
  }

  float3x3 transposed() const
  {
    float3x3 result;
    transpose_m3_m3(result.values, values);
    return result;
  }
};

}  // namespace blender
