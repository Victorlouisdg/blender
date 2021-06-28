#include "BLI_float3.hh"
#include "BLI_float3x3.hh"

#include "BLI_index_range.hh"

using blender::float2x2;
using blender::float3;
using blender::float3x3;
using blender::IndexRange;
using std::array;

/* The deformation gradient for a 2D triangle mapped to 3D. */
class DeformationGradient {
 public:
  float3 du;
  float3 dv;

  void calculate(float3 x0, float3 x1, float3 x2, float2x2 delta_u)
  {
    float3 delta_x1 = x1 - x0;
    float3 delta_x2 = x2 - x0;

    float du0 = delta_x1[0] * delta_u.values[0][0] + delta_x2[0] * delta_u.values[1][0];
    float du1 = delta_x1[1] * delta_u.values[0][0] + delta_x2[1] * delta_u.values[1][0];
    float du2 = delta_x1[2] * delta_u.values[0][0] + delta_x2[2] * delta_u.values[1][0];

    float dv0 = delta_x1[0] * delta_u.values[0][1] + delta_x2[0] * delta_u.values[1][1];
    float dv1 = delta_x1[1] * delta_u.values[0][1] + delta_x2[1] * delta_u.values[1][1];
    float dv2 = delta_x1[2] * delta_u.values[0][1] + delta_x2[2] * delta_u.values[1][1];

    du = float3(du0, du1, du2);
    dv = float3(dv0, dv1, dv2);
  }
};

class GravityForceElement {
 public:
  float3 force;

  void calculate(float mass, float standard_gravity)
  {
    force = float3(0.0f, 0.0f, mass * standard_gravity);
  }
};

class StretchForceElementBW {
 public:
  array<float3, 3> forces;
  array<array<float3x3, 3>, 3> force_derivatives;  // float3x3 force_derivatives[3][3];

  /* TODO add second calculate() that calculates the deformation gradient. */
  void calculate(
      float ku, float kv, float area_factor, DeformationGradient F, float3 dwu_dx, float3 dwv_dx)
  {
    set_zero();
    calculate_single_stretch_direction(ku, area_factor, F.du, dwu_dx);
    calculate_single_stretch_direction(kv, area_factor, F.dv, dwv_dx);
  }

  void set_zero()
  {
    for (int m : IndexRange(3)) {
      forces[m] = float3(0.0f);
      for (int n : IndexRange(3)) {
        force_derivatives[m][n] = float3x3(0.0f);
      }
    }
  }

  void calculate_single_stretch_direction(float k, float area_factor, float3 w, float3 dw_dx)
  {
    float w_norm = w.normalize_and_get_length();
    float C = area_factor * (w_norm - 1.0f);

    array<float3, 3> dC_dx = array<float3, 3>();

    for (int m : IndexRange(3)) {
      dC_dx[m] = area_factor * dw_dx[m] * w;
    }

    for (int m : IndexRange(3)) {
      forces[m] += -k * C * dC_dx[m];
    }

    /* Derivatives */
    float3x3 I_w_wT = float3x3::identity() - float3x3::outer(w, w);

    for (int m : IndexRange(3)) {
      for (int n : IndexRange(3)) {
        float3x3 dC_dx_mn = area_factor / w_norm * dw_dx[m] * dw_dx[n] * I_w_wT;
        float3x3 dC_outer = float3x3::outer(dC_dx[m], dC_dx[n]);
        float3x3 df_dx_mn = -k * (dC_outer + C * dC_dx_mn);

        force_derivatives[m][n] += df_dx_mn;
      }
    }
  }
};