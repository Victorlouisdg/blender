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
#include "BLI_index_range.hh"

namespace blender {

struct float3x3 {
  /* Note on initialization: blender interprets 2D-arrays in column-major format.
     so in the example below, {1, 2, 3} is the first column of the matrix.

    float array[3][3] = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };

  float3x3 block = float3x3(array);
  */


  float values[3][3];

  float3x3() = default;

  float3x3(const float *matrix)
  {
    memcpy(values, matrix, sizeof(float) * 9);
  }

  float3x3(const float matrix[3][3]) : float3x3(static_cast<const float *>(matrix[0]))
  {
  }

  static float3x3 identity()
  {
    float3x3 mat;
    unit_m3(mat.values);
    return mat;
  }

  static float3x3 outer(const float3 &a, const float3 &b)
  {
    float3x3 mat;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        mat.values[j][i] = a[i] * b[j];
      }
    }
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

  friend float3x3 operator-(const float3x3 &a, const float3x3 &b)
  {
    float3x3 result;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        result.values[i][j] = a.values[i][j] - b.values[i][j];
      }
    }
    return result;
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

std::ostream &operator<<(std::ostream &os, float3x3 const &m)
{
  /* Prints [[row0], [row1], [row2]]. */
  os << "[";
  for (int i : IndexRange(3)) {
    os << "[";
    for (int j : IndexRange(3)) {
      os << m.values[j][i];
      if (j != 2) {
        os << ", ";
      }
    }
    os << "]";
  }
  os << "]";
  return os;
};

}  // namespace blender