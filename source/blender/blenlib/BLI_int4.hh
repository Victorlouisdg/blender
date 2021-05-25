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

namespace blender {

struct int4 {
  int i, j, k, l;

  int4() = default;

  int4(const int *ptr) : i{ptr[0]}, j{ptr[1]}, k{ptr[2]}, l{ptr[3]}
  {
  }

  int4(int i, int j, int k, int l) : i{i}, j{j}, k{k}, l{l}
  {
  }

  operator int *()
  {
    return &i;
  }
};

}  // namespace blender
