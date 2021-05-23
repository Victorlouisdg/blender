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

/* IMPORTANT TODO: store values in column-major format, like the other 3x3 adn 4x4 matrices in
 * blender. Also change all the places where a float2x2 is created. */
struct float2x2 {
  float values[2][2];

  float2x2() = default;

  float2x2(const float *matrix)
  {
    memcpy(values, matrix, sizeof(float) * 4);
  }

  float2x2(const float matrix[2][2]) : float2x2(static_cast<const float *>(matrix[0]))
  {
  }

  // friend float2x2 operator*(const float2x2 &a, const float2 &b)
  // {
  //   float2x2 result;
  //   mul_m2_v2(result.values, a.values, &b);
  //   return result;
  // }

  friend float2x2 operator*(const float s, const float2x2 &a)
  {
    float2x2 result;

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        result.values[i][j] = s * a.values[i][j];
      }
    }
    return result;
  }

  float2x2 inverted() const
  {
    /* TODO: maybe move this to matrix_math.c? And check for division by zero. */
    float a = values[0][0];
    float b = values[0][1];
    float c = values[1][0];
    float d = values[1][1];

    float result_array[2][2] = {{d, -b}, {-c, a}};
    float determinant = a * d - b * c;

    /* For Debug builds, assert that determinant is not too close to zero? */
    return (1.0f / determinant) * float2x2(result_array);
  }
};

}  // namespace blender
