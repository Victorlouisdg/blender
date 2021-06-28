/* Apache License, Version 2.0 */

#include "testing/testing.h"

#include "BKE_cloth_simulator_bw_attribute_provider.hh"
#include "BKE_cloth_simulator_bw_force_elements.hh"

namespace blender::tests {

TEST(clothbw, DeformationGradient)
{
  float2 uv0 = float2(0.0f, 0.0f);
  float2 uv1 = float2(1.0f, 0.0f);
  float2 uv2 = float2(0.0f, 1.0f);
  auto [inverted_delta_u, dwu_dx, dwv_dx] = ClothBWAttributeProvider::calculate_uv_derived_values(
      uv0, uv1, uv2);

  float3 x0 = float3(0.0f, 0.0f, 0.0f);
  float3 x1 = float3(1.0f, 0.0f, 0.0f);
  float3 x2 = float3(0.0f, 1.0f, 0.0f);
  //   float2x2 delta_u = triangle_inverted_delta_u_matrices[ti];
  //   deformation_gradients[ti].calculate(x0, x1, x2, delta_u);

  DeformationGradient F = DeformationGradient();
  F.calculate(x0, x1, x2, inverted_delta_u);

  float3 Fdu = F.du;
  float3 Fdv = F.dv;

  float h = 0.0001f;

  std::cout << "dwu_dx" << dwu_dx << std::endl;

  for (int m : IndexRange(3)) {    // vertex0, vertex1, vertex2
    for (int s : IndexRange(3)) {  // x, y, z
      std::array<float3, 3> x = {x0, x1, x2};
      x[m][s] += h;
      F.calculate(x[0], x[1], x[2], inverted_delta_u);
      float3 Fdu_h = F.du;
      std::cout << "u:" << m << " " << s << (Fdu_h - Fdu) / h << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "dwv_dx" << dwv_dx << std::endl;

  for (int m : IndexRange(3)) {    // vertex0, vertex1, vertex2
    for (int s : IndexRange(3)) {  // x, y, z
      std::array<float3, 3> x = {x0, x1, x2};
      x[m][s] += h;
      F.calculate(x[0], x[1], x[2], inverted_delta_u);
      float3 Fdv_h = F.dv;
      std::cout << "v:" << m << " " << s << (Fdv_h - Fdv) / h << std::endl;
    }
  }

  // std::cout << F.du << std::endl;
  // std::cout << F.dv << std::endl;
}

TEST(clothbw, Stretch)
{

  /*  Simple 2D triangle that we scale 2x and place in the z=0 plane in 3D.
   *
   *   (0, 1)
   *   *\
   *   |  \                 v  ^
   *   |    \                  |
   *   |      \                |
   *   |        \   (1, 0)     |
   *   *----------*            -------> u
   *   (0, 0)
   *
   */

  float2 uv0 = float2(0.0f, 0.0f);
  float2 uv1 = float2(1.0f, 0.0f);
  float2 uv2 = float2(0.0f, 1.0f);
  auto [inverted_delta_u, dwu_dx, dwv_dx] = ClothBWAttributeProvider::calculate_uv_derived_values(
      uv0, uv1, uv2);

  float3 x0 = float3(0.0f, 0.0f, 0.0f);
  float3 x1 = float3(2.0f, 0.0f, 0.0f);
  float3 x2 = float3(0.0f, 2.0f, 0.0f);

  DeformationGradient F = DeformationGradient();
  F.calculate(x0, x1, x2, inverted_delta_u);

  StretchForceElementBW stretch = StretchForceElementBW();
  stretch.calculate(1.0f, 1.0f, 1.0f, F, dwu_dx, dwv_dx);

  float3 f0 = stretch.forces[0];
  float3 f1 = stretch.forces[1];
  float3 f2 = stretch.forces[2];

  float3 sum_of_internal_forces = f0 + f1 + f2;

  EXPECT_EQ(sum_of_internal_forces[0], 0.0f);
  EXPECT_EQ(sum_of_internal_forces[1], 0.0f);
  EXPECT_EQ(sum_of_internal_forces[2], 0.0f);

  EXPECT_GE(float3::dot(f0, float3(1.0f, 1.0f, 0.0f)), 0.0f);
  EXPECT_GE(float3::dot(f1, float3(-1.0f, 0.0f, 0.0f)), 0.0f);
  EXPECT_GE(float3::dot(f2, float3(0.0f, -1.0f, 0.0f)), 0.0f);

  // std::cout << "Deformation Gradient" << std::endl;
  // std::cout << F.wu << std::endl;
  // std::cout << F.wv << std::endl;

  // std::cout << "Forces" << std::endl;
  // std::cout << f0 << std::endl;
  // std::cout << f1 << std::endl;
  // std::cout << f2 << std::endl;
  // std::cout << sum_of_internal_forces << std::endl;
}

}  // namespace blender::tests
