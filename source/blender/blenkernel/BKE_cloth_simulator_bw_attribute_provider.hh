#include "BLI_index_range.hh"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_cloth_bw_general.hh"

using blender::IndexRange;

class ClothBWAttributeProvider {
 private:
  int amount_of_vertices;
  int system_size;

  Array<float3> vertex_positions;
  Array<float3> vertex_velocities;

  Array<float> vertex_masses;

 public:
  ClothBWAttributeProvider(const Mesh &mesh)
  {
    amount_of_vertices = mesh.totvert;
    system_size = 3 * amount_of_vertices;
    initialize_vertex_positions(mesh);
    initialize_vertex_velocities();
    initialize_vertex_masses_uniform(1.0);
  }

  //   void initialize_vertex_positions(const Mesh &mesh)
  //   {
  //     vertex_positions = Vector3N(amount_of_vertices);
  //     for (const int i : IndexRange(amount_of_vertices)) {
  //       vertex_positions[i] = Vector3(mesh.mvert[i].co);
  //     }
  //   }

  //   void initialize_vertex_velocities()
  //   {
  //     vertex_velocities = Vector3N(amount_of_vertices);
  //     vertex_velocities.setZero();
  //   }

  //   void initialize_vertex_masses_uniform(double total_mass)
  //   {
  //     double vertex_mass = total_mass / amount_of_vertices;
  //     vertex_masses = VectorN(amount_of_vertices);
  //     vertex_masses.fill(vertex_mass);
  //   }

  int get_amount_of_vertices()
  {
    return amount_of_vertices;
  }

  Vector3N &get_vertex_positions()
  {
    return vertex_positions;
  }

  Vector3N &get_vertex_velocities()
  {
    return vertex_velocities;
  }

  VectorN &get_vertex_masses()
  {
    return vertex_masses;
  }
};
