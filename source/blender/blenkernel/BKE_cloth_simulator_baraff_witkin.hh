#include <iostream>

#include "BLI_array.hh"
#include "BLI_float2.hh"
#include "BLI_float2x2.hh"
#include "BLI_float3.hh"

#include "BLI_index_range.hh"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"

#include <tuple>

// #include "bmesh.h"

using blender::Array;
using blender::float2;
using blender::float2x2;
using blender::float3;
using blender::IndexRange;
using blender::Span;

/* Utilities */

/* Copied from point_distirbute. */
static Span<MLoopTri> get_mesh_looptris(const Mesh &mesh)
{
  /* This only updates a cache and can be considered to be logically const. */
  const MLoopTri *looptris = BKE_mesh_runtime_looptri_ensure(const_cast<Mesh *>(&mesh));
  const int looptris_len = BKE_mesh_runtime_looptri_len(&mesh);
  return {looptris, looptris_len};
}

/* Adapted from multires_unsubdivide. */
// static BMesh *get_bmesh_from_mesh(const Mesh *mesh)
// {
//   const BMAllocTemplate allocsize = BMALLOC_TEMPLATE_FROM_ME(mesh);
//   BMesh *bmesh = BM_mesh_create(&allocsize,
//                                 &((struct BMeshCreateParams){
//                                     .use_toolflags = true,
//                                 }));

//   BM_mesh_bm_from_me(bmesh,
//                      mesh,
//                      (&(struct BMeshFromMeshParams){
//                          .calc_face_normal = true,
//                      }));

//   return bmesh;
// }

// BMIter viter, eiter, fiter;
// BMFace *f = NULL;

// // BMesh *bmesh = get_bmesh_from_mesh(&mesh);

// BM_ITER_MESH (f, &fiter, bmesh, BM_FACES_OF_MESH) {
// }

/* Cloth Simulator based off of the popular paper "Large steps in cloth simulation."
 * by Baraff and Witkin, hence the name. In comments below Ithis paper is referred to
 * as BW98.
 *
 * Each timestep, the simulator builds a linear system of equations that's filled with
 * several conditions that encode how we want the cloth to behave. (e.g. how it responds to
 * stretch, or to collisions.) sidenote: collisions response can be done in several ways.
 */
class ClothSimulatorBaraffWitkin {
 public:
  int n_vertices;
  int n_triangles;

  /* State */
  Array<float3> vertex_positions;
  Array<float3> vertex_velocities;

  /* Indices */
  Array<std::tuple<int, int, int>> triangle_vertex_indices;
  Array<std::tuple<int, int, int, int>> bending_vertex_indices;

  /* Parameters / Configuration */
  Array<float2> vertex_positions_uv;
  Array<float> vertex_masses;
  Array<float3> triangle_normals;
  Array<float> triangle_stretch_stiffness_u;
  Array<float> triangle_stretch_stiffness_v;
  Array<float> triangle_areas_uv;
  Array<float3> edges_normals;
  /* Note about the "uv" in the names above. */

  /* Datastructures for intermediate computations. */
  Array<float3> vertex_forces;

  /* Precomputed quantities. */
  Array<float2x2> triangle_inverted_delta_u_matrices;
  Array<float> triangle_wu_derivatives;
  Array<float> triangle_wv_derivatives;

  /* I decided that the simulator should have it's own local copy of the necessary
   * mesh attributes. I'll make setters() for certain attributes later that do
   * the required recalculations.
   */
  void initialize(const Mesh &mesh)
  {
    std::cout << "Cloth simulation initialisation" << std::endl;

    initialize_vertex_attributes(mesh);
    initialize_triangle_attributes(mesh);
  };

  void initialize_vertex_attributes(const Mesh &mesh)
  {
    n_vertices = mesh.totvert;

    vertex_positions = Array<float3>(n_vertices);
    vertex_velocities = Array<float3>(n_vertices, float3(0.0f));
    vertex_forces = Array<float3>(n_vertices, float3(0.0f));

    for (const int i : IndexRange(mesh.totvert)) {
      MVert vertex = mesh.mvert[i];
      vertex_positions[i] = vertex.co;
    }
  }

  void initialize_triangle_attributes(const Mesh &mesh)
  {
    Span<MLoopTri> looptris = get_mesh_looptris(mesh);
    for (const int looptri_index : looptris.index_range()) {
      const MLoopTri &looptri = looptris[looptri_index];
      const int v0_loop = looptri.tri[0];
      const int v1_loop = looptri.tri[1];
      const int v2_loop = looptri.tri[2];
      const int v0_index = mesh.mloop[v0_loop].v;
      const int v1_index = mesh.mloop[v1_loop].v;
      const int v2_index = mesh.mloop[v2_loop].v;
      const float3 v0 = vertex_positions[v0_index];
      const float3 v1 = vertex_positions[v1_index];
      const float3 v2 = vertex_positions[v2_index];

      const float3 normal = float3::cross(v1 - v0, v2 - v0);
      const float triangle_area_uv = normal.length() / 2;

      std::cout << "normal" << normal << std::endl;
      std::cout << "triangle_area_uv" << triangle_area_uv << std::endl;
    }
  };

  void step()
  {
    std::cout << "Cloth simulation step" << std::endl;
  };

  ClothSimulatorBaraffWitkin()
  {
    std::cout << "ClothSimulatorBaraffWitkin constructed" << std::endl;
  }

  ~ClothSimulatorBaraffWitkin()
  {
    std::cout << "ClothSimulatorBaraffWitkin destructed" << std::endl;
  }
};