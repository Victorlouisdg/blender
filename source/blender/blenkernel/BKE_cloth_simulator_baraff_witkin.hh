#include <iostream>

#include "BLI_array.hh"
#include "BLI_float2.hh"
#include "BLI_float3.hh"

using blender::Array;
using blender::float2;
using blender::float3;

class ClothSimulatorBaraffWitkin {
 public:
  int n_vertices;
  int n_triangles;

  /* State */
  Array<float3> vertex_positions;
  Array<float3> vertex_velocities;

  /* Indices */
  Array<int[3]> triangle_vertex_indices;
  Array<int[4]> bending_vertex_indices;

  /* Parameters / Configuration */
  Array<float2> vertex_positions_uv; /* = rest state */
  Array<float> vertex_masses;
  Array<float3> triangle_normals;
  Array<float> triangle_stretch_stiffness_u;
  Array<float> triangle_stretch_stiffness_v;
  Array<float> triangle_areas_uv;
  Array<float3> edges_normals;

  /* Datastructures for intermediate computations. */
  Array<float3> vertex_forces;

  /* Precomputed quantities. */
  Array<float[2][2]> triangle_inverted_delta_u_matrices;
  Array<float> triangle_wu_derivatives;
  Array<float> triangle_wv_derivatives;

  void initialize()
  {
    std::cout << "Cloth simulation initialisation" << std::endl;
    vertex_positions = Array<float3>(100000);
  };

  void step()
  {
    std::cout << "Cloth simulation step" << std::endl;
  };
};