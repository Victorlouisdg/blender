#include "BLI_array.hh"
#include "BLI_float2.hh"
#include "BLI_float2x2.hh"
#include "BLI_float3.hh"
#include "BLI_float3x3.hh"
#include "BLI_index_range.hh"
#include "BLI_sparse_matrix.hh"
#include "BLI_vector.hh"

#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"

#include "BKE_customdata.h"
#include "BKE_deform.h"
#include "BKE_lib_query.h"
#include "BKE_object.h"

typedef std::array<int, 2> int2;
typedef std::array<int, 3> int3;
typedef std::array<int, 4> int4;

using blender::Array;
using blender::float2;
using blender::float2x2;
using blender::float3;
using blender::float3x3;
using blender::IndexRange;
using blender::MutableSpan;
using blender::Span;
using blender::Vector;

class ClothBWAttributeProvider {
 private:
  int amount_of_vertices;
  int amount_of_triangles;
  int amount_of_bend_edges;

  /* State */
  Array<float3> vertex_positions;
  Array<float3> vertex_velocities;

  /* Indices */
  Array<int2> spring_vertex_indices;
  Array<int3> triangle_vertex_indices;
  Array<int4> bend_vertex_indices;

  /* Parameters / Configuration */
  Array<float> vertex_masses;

  Array<float> triangle_stretch_stiffness_u;
  Array<float> triangle_stretch_stiffness_v;
  Array<float> triangle_shear_stiffness;
  Array<float> bend_stiffness;

  Array<float> triangle_stretch_damping_u;
  Array<float> triangle_stretch_damping_v;
  Array<float> triangle_shear_damping;
  Array<float> bend_damping;

  Array<float3> triangle_normals;
  Array<float> triangle_area_factors;
  Array<float2x2> triangle_inverted_delta_u_matrices;
  Array<float3> triangle_wu_derivatives;
  Array<float3> triangle_wv_derivatives;

  Array<float> bend_rest_lengths;

  Vector<int> pinned_vertices;

 public:
  ClothBWAttributeProvider()
  {
  }

  ClothBWAttributeProvider(const Mesh &mesh,
                           const ClothBWModifierData &modifier_data,
                           const Object &cloth_object)
  {
    amount_of_vertices = mesh.totvert;
    initialize_vertex_positions(mesh);
    initialize_vertex_velocities();
    initialize_vertex_masses_uniform(1.0);
    initialize_triangle_attributes(mesh, modifier_data);
    initialize_bend_attributes(mesh, modifier_data);
    initialize_pinned_vertices(mesh, modifier_data, cloth_object);
  }

  void initialize_vertex_positions(const Mesh &mesh)
  {
    vertex_positions = Array<float3>(amount_of_vertices);
    for (const int i : IndexRange(amount_of_vertices)) {
      vertex_positions[i] = mesh.mvert[i].co;
    }
  }

  void initialize_vertex_velocities()
  {
    vertex_velocities = Array<float3>(amount_of_vertices, float3(0.0f));
  }

  void initialize_vertex_masses_uniform(float total_mass)
  {
    float vertex_mass = total_mass / amount_of_vertices;
    vertex_masses = Array<float>(amount_of_vertices, vertex_mass);
  }

  void initialize_triangle_attributes(const Mesh &mesh, const ClothBWModifierData &md)
  {
    Span<MLoopTri> looptris = get_mesh_looptris(mesh);
    amount_of_triangles = looptris.size();

    triangle_stretch_stiffness_u = Array<float>(amount_of_triangles, md.stretch_stiffness);
    triangle_stretch_stiffness_v = Array<float>(amount_of_triangles, md.stretch_stiffness);
    triangle_shear_stiffness = Array<float>(amount_of_triangles, md.shear_stiffness);

    float stretch_damping = md.stretch_stiffness * md.stretch_damping_factor;
    float shear_damping = md.shear_stiffness * md.shear_damping_factor;

    triangle_stretch_damping_u = Array<float>(amount_of_triangles, stretch_damping);
    triangle_stretch_damping_v = Array<float>(amount_of_triangles, stretch_damping);
    triangle_shear_damping = Array<float>(amount_of_triangles, shear_damping);

    triangle_vertex_indices = Array<int3>(amount_of_triangles);
    triangle_normals = Array<float3>(amount_of_triangles);
    triangle_area_factors = Array<float>(amount_of_triangles);
    triangle_inverted_delta_u_matrices = Array<float2x2>(amount_of_triangles);
    triangle_wu_derivatives = Array<float3>(amount_of_triangles);
    triangle_wv_derivatives = Array<float3>(amount_of_triangles);

    /* Reusing the texture UVs for the rest state here is what is
     * also done in BW98, but feel a bit hacky.
     */
    const float2 v0_uv = mesh.mloopuv[0].uv;
    const float2 v1_uv = mesh.mloopuv[1].uv;
    const float3 v0_pos = vertex_positions[mesh.mloop[0].v];
    const float3 v1_pos = vertex_positions[mesh.mloop[1].v];
    float uv_distance = (v0_uv - v1_uv).length();
    float world_distance = (v0_pos - v1_pos).length();
    float uv_to_world_scale_factor = world_distance / uv_distance;

    Array<float2> vertex_positions_uv = Array<float2>(amount_of_vertices);
    for (const int i : IndexRange(mesh.totloop)) {
      const MLoop &loop = mesh.mloop[i];
      const float2 uv = mesh.mloopuv[i].uv;
      vertex_positions_uv[loop.v] = uv * uv_to_world_scale_factor;
    }

    for (const int looptri_index : looptris.index_range()) {
      const MLoopTri &looptri = looptris[looptri_index];
      const int v0_loop = looptri.tri[0];
      const int v1_loop = looptri.tri[1];
      const int v2_loop = looptri.tri[2];
      const int v0_index = mesh.mloop[v0_loop].v;
      const int v1_index = mesh.mloop[v1_loop].v;
      const int v2_index = mesh.mloop[v2_loop].v;
      const float3 v0_pos = vertex_positions[v0_index];
      const float3 v1_pos = vertex_positions[v1_index];
      const float3 v2_pos = vertex_positions[v2_index];

      triangle_vertex_indices[looptri_index] = int3({v0_index, v1_index, v2_index});

      const float3 normal = float3::cross(v1_pos - v0_pos, v2_pos - v0_pos);
      triangle_normals[looptri_index] = normal;

      const float triangle_area_uv = normal.length() / 2.0f;

      /* Sqrt such that energy is linear in area (see Bhat and Pritchard)*/
      triangle_area_factors[looptri_index] = sqrtf(triangle_area_uv);

      float2 uv0 = vertex_positions_uv[v0_index];
      float2 uv1 = vertex_positions_uv[v1_index];
      float2 uv2 = vertex_positions_uv[v2_index];

      auto [inverted_delta_u, dwu_dx, dwv_dx] = calculate_uv_derived_values(uv0, uv1, uv2);
      triangle_inverted_delta_u_matrices[looptri_index] = inverted_delta_u;
      triangle_wu_derivatives[looptri_index] = dwu_dx;
      triangle_wv_derivatives[looptri_index] = dwv_dx;
    }
  }

  static std::tuple<float2x2, float3, float3> calculate_uv_derived_values(float2 uv0,
                                                                          float2 uv1,
                                                                          float2 uv2)
  {
    float u0 = uv0[0];
    float u1 = uv1[0];
    float u2 = uv2[0];
    float v0 = uv0[1];
    float v1 = uv1[1];
    float v2 = uv2[1];

    float delta_u1 = u1 - u0;
    float delta_u2 = u2 - u0;
    float delta_v1 = v1 - v0;
    float delta_v2 = v2 - v0;

    float array[2][2] = {{delta_u1, delta_u2}, {delta_v1, delta_v2}};
    float2x2 delta_u = float2x2(array);
    float2x2 inverted_delta_u = delta_u.inverted();

    float dw_denominator = 1.0f / (delta_u1 * delta_v2 - delta_u2 * delta_v1);

    float3 dwu_dx;
    dwu_dx[0] = (delta_v1 - delta_v2) * dw_denominator;
    dwu_dx[1] = delta_v2 * dw_denominator;
    dwu_dx[2] = -delta_v1 * dw_denominator;

    float3 dwv_dx;
    dwv_dx[0] = (delta_u2 - delta_u1) * dw_denominator;
    dwv_dx[1] = -delta_u2 * dw_denominator;
    dwv_dx[2] = delta_u1 * dw_denominator;

    return {inverted_delta_u, dwu_dx, dwv_dx};
  }

  void initialize_bend_attributes(const Mesh &mesh, const ClothBWModifierData &md)
  {
    /* Maps from the 2 vertices of a looptri edge to an oppossing vertex. */
    std::multimap<std::pair<int, int>, int> looptri_edge_opposing_vertices =
        std::multimap<std::pair<int, int>, int>();

    Span<MLoopTri> looptris = get_mesh_looptris(mesh);

    for (const int looptri_index : looptris.index_range()) {
      const MLoopTri &looptri = looptris[looptri_index];
      const int v0_loop = looptri.tri[0];
      const int v1_loop = looptri.tri[1];
      const int v2_loop = looptri.tri[2];
      const int v0_index = mesh.mloop[v0_loop].v;
      const int v1_index = mesh.mloop[v1_loop].v;
      const int v2_index = mesh.mloop[v2_loop].v;

      /* bend edges. */
      looptri_edge_opposing_vertices.insert(
          std::pair<std::pair<int, int>, int>(std::minmax(v0_index, v1_index), v2_index));
      looptri_edge_opposing_vertices.insert(
          std::pair<std::pair<int, int>, int>(std::minmax(v0_index, v2_index), v1_index));
      looptri_edge_opposing_vertices.insert(
          std::pair<std::pair<int, int>, int>(std::minmax(v1_index, v2_index), v0_index));
    }

    Vector<int4> bend_indices_vector = Vector<int4>();
    Vector<float> bend_rest_lengths_vector = Vector<float>();

    /* Here we need to check which looptri edges have 2 opposing vertices, these edges will
     * become bend edges. */
    for (std::multimap<std::pair<int, int>, int>::iterator it =
             looptri_edge_opposing_vertices.begin();
         it != looptri_edge_opposing_vertices.end();
         it++) {
      std::pair<int, int> edge = it->first;

      /* Note that indices of the shared edge are stored in position 1 and 2, this is the
       * convention used in the "Implementing Baraff Wikin" by David Pritchard. */
      int v0 = it->second;
      int v1 = edge.first;
      int v2 = edge.second;

      std::multimap<std::pair<int, int>, int>::iterator it_next = std::next(it, 1);

      if (it_next != looptri_edge_opposing_vertices.end()) {
        std::pair<int, int> edge_next = it_next->first;
        if (edge == edge_next) { /* This means two looptris share this edge. */
          int v3 = it_next->second;
          bend_indices_vector.append(int4({v0, v1, v2, v3}));
          float3 x0 = vertex_positions[v0];
          float3 x3 = vertex_positions[v3];
          float rest_length = (x3 - x0).length();
          bend_rest_lengths_vector.append(rest_length);
        }
      }
    }

    amount_of_bend_edges = bend_indices_vector.size();

    bend_vertex_indices = Array<int4>(bend_indices_vector.as_span());
    bend_rest_lengths = Array<float>(bend_rest_lengths_vector.as_span());
    bend_stiffness = Array<float>(amount_of_bend_edges, md.bend_stiffness);

    float bend_damping_value = md.bend_stiffness * md.bend_damping_factor;
    bend_damping = Array<float>(amount_of_bend_edges, bend_damping_value);
  }

  void initialize_pinned_vertices(const Mesh &mesh,
                                  const ClothBWModifierData &md,
                                  const Object &cloth_object)
  {
    pinned_vertices = Vector<int>();

    const int defgrp_index = BKE_object_defgroup_name_index(&cloth_object, md.defgrp_name);
    if (defgrp_index != -1) {
      MDeformVert *dvert, *dv;
      dvert = (MDeformVert *)CustomData_get_layer(&(mesh.vdata), CD_MDEFORMVERT);
      if (dvert) {
        dv = &dvert[0];
        for (uint i = 0; i < mesh.totvert; i++, dv++) {
          const bool found = BKE_defvert_find_weight(dv, defgrp_index) > 0.0f;
          if (found) {
            pinned_vertices.append(i);
          }
        }
      }
    }
  }

  int get_amount_of_vertices()
  {
    return amount_of_vertices;
  }

  int get_amount_of_triangles()
  {
    return amount_of_triangles;
  }

  int get_amount_of_bend_edges()
  {
    return amount_of_bend_edges;
  }

  MutableSpan<float3> get_vertex_positions()
  {
    return vertex_positions;
  }

  MutableSpan<float3> get_vertex_velocities()
  {
    return vertex_velocities;
  }

  Span<float> get_vertex_masses()
  {
    return vertex_masses;
  }

  Span<float> get_triangle_stretch_stiffness_u()
  {
    return triangle_stretch_stiffness_u;
  }

  Span<float> get_triangle_stretch_stiffness_v()
  {
    return triangle_stretch_stiffness_v;
  }

  Span<float> get_triangle_shear_stiffness()
  {
    return triangle_shear_stiffness;
  }

  Span<float> get_bend_stiffness()
  {
    return bend_stiffness;
  }

  Span<float> get_triangle_stretch_damping_u()
  {
    return triangle_stretch_damping_u;
  }

  Span<float> get_triangle_stretch_damping_v()
  {
    return triangle_stretch_damping_v;
  }

  Span<float> get_triangle_shear_damping()
  {
    return triangle_shear_damping;
  }

  Span<float> get_bend_damping()
  {
    return bend_damping;
  }

  Span<float> get_triangle_area_factors()
  {
    return triangle_area_factors;
  }

  Span<float2x2> get_triangle_inverted_delta_u_matrices()
  {
    return triangle_inverted_delta_u_matrices;
  }

  Span<float3> get_triangle_wu_derivatives()
  {
    return triangle_wu_derivatives;
  }

  Span<float3> get_triangle_wv_derivatives()
  {
    return triangle_wv_derivatives;
  }

  Span<int3> get_triangle_vertex_indices()
  {
    return triangle_vertex_indices;
  }

  Span<int4> get_bend_vertex_indices()
  {
    return bend_vertex_indices;
  }

  Span<float> get_bend_rest_lengths()
  {
    return bend_rest_lengths;
  }

  Vector<int> get_pinned_vertices()
  {
    return pinned_vertices;
  }

  /* Copied from point_distirbute. */
  static Span<MLoopTri> get_mesh_looptris(const Mesh &mesh)
  {
    /* This only updates a cache and can be considered to be logically const. */
    const MLoopTri *looptris = BKE_mesh_runtime_looptri_ensure(const_cast<Mesh *>(&mesh));
    const int looptris_len = BKE_mesh_runtime_looptri_len(&mesh);
    return {looptris, looptris_len};
  }
};
