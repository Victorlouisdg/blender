#include "BLI_float3.hh"
#include "BLI_float3x3.hh"

#include "BLI_index_range.hh"

#include "BKE_cloth_simulator_constants.hh"

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
  array<array<float3x3, 3>, 3> force_derivatives;

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

    /* Force derivatives */
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

class ShearForceElementBW {
 public:
  array<float3, 3> forces;
  array<array<float3x3, 3>, 3> force_derivatives;

  void calculate(float k, float area_factor, DeformationGradient F, float3 dwu_dx, float3 dwv_dx)
  {
    float3 wu = F.du;
    float3 wv = F.dv;

    /* Shear condition: section 4.3 in [BW98]. */
    float C = area_factor * float3::dot(wu, wv);

    array<float3, 3> dC_dx = array<float3, 3>();

    for (int m : IndexRange(3)) {
      dC_dx[m] = area_factor * (dwu_dx[m] * wv + dwv_dx[m] * wu);
    }

    for (int m : IndexRange(3)) {
      forces[m] = -k * C * dC_dx[m];
    }

    /* Force derivatives */
    for (int m : IndexRange(3)) {
      for (int n : IndexRange(3)) {
        float3x3 dC_dx_mn = area_factor * (dwu_dx[m] * dwv_dx[n] + dwu_dx[n] * dwv_dx[m]) *
                            float3x3::identity();

        float3x3 dC_outer = float3x3::outer(dC_dx[m], dC_dx[n]);
        float3x3 df_dx_mn = -k * (dC_outer + C * dC_dx_mn);

        force_derivatives[m][n] = df_dx_mn;
      }
    }
  }
};

class BendingForceElementBW {
 public:
  array<float3, 4> forces;
  array<array<float3x3, 4>, 4> force_derivatives;

  void calculate(float k, float3 x0, float3 x1, float3 x2, float3 x3)
  {
    float3 e = x1 - x2;
    float3 nA = float3::cross(x2 - x0, x1 - x0);
    float3 nB = float3::cross(x1 - x3, x2 - x3);

    float e_norm = e.normalize_and_get_length();
    float nA_norm = nA.normalize_and_get_length();
    float nB_norm = nB.normalize_and_get_length();

    float cos = float3::dot(nA, nB);
    float sin = float3::dot(float3::cross(nA, nB), e);

    float C = atan2f(sin, cos);

    float3 qA[4] = {x2 - x1, x0 - x2, x1 - x0, float3(0.0f)};
    float3 qB[4] = {float3(0.0f), x2 - x3, x3 - x1, x1 - x2};
    float3x3 qIe[4] = {float3x3(0.0f),
                       1.0f / e_norm * float3x3::diagonal(1.0f),
                       1.0f / e_norm * float3x3::diagonal(-1.0f),
                       float3x3(0.0f)};

    float3 dC_dx[4];
    float3 de_dx[4][3];
    float3 dnA_dx[4][3];
    float3 dnB_dx[4][3];
    float dcos_dx[4][3];
    float dsin_dx[4][3];

    float C_dot = 0.0f;

    /* Forces */
    for (int m : IndexRange(4)) {
      float3x3 SA = 1.0f / nA_norm * float3x3::skew(qA[m]);
      float3x3 SB = 1.0f / nB_norm * float3x3::skew(qB[m]);

      for (int s : IndexRange(3)) {
        de_dx[m][s] = qIe[m].row(s);
        dnA_dx[m][s] = SA.row(s);
        dnB_dx[m][s] = SB.row(s);
        dcos_dx[m][s] = float3::dot(dnA_dx[m][s], nB) + float3::dot(nA, dnB_dx[m][s]);
        dsin_dx[m][s] = float3::dot(
                            float3::cross(dnA_dx[m][s], nB) + float3::cross(nA, dnB_dx[m][s]), e) +
                        float3::dot(float3::cross(nA, nB), de_dx[m][s]);
        dC_dx[m][s] = cos * dsin_dx[m][s] - sin * dcos_dx[m][s];
      }
    }

    for (int m : IndexRange(4)) {
      float3 force_m = -k * C * dC_dx[m];
      forces[m] = force_m;
    }

    /* Force derivatives */

    float3x3 dC_dx_mn[4][4];

    for (int m : IndexRange(4)) {
      for (int n : IndexRange(4)) {
        for (int s : IndexRange(3)) {
          for (int t : IndexRange(3)) {
            /* TODO think over whether s and t might need to be swapped. */
            float3 dnA_dx_mnst = normalA_second_derivatives[m][n][s][t] / nA_norm;
            float3 dnB_dx_mnst = normalB_second_derivatives[m][n][s][t] / nB_norm;

            float3 dnA_dx_ms = dnA_dx[m][s];
            float3 dnA_dx_nt = dnA_dx[n][t];
            float3 dnB_dx_ms = dnB_dx[m][s];
            float3 dnB_dx_nt = dnB_dx[n][t];

            float dcos_dx_mnst = float3::dot(dnA_dx_mnst, nB) + float3::dot(dnB_dx_nt, dnA_dx_ms) +
                                 float3::dot(dnA_dx_nt, dnB_dx_ms) + float3::dot(nA, dnB_dx_mnst);

            float dsin_dx_mnst =
                float3::dot(float3::cross(dnA_dx_mnst, nB) + float3::cross(dnA_dx_ms, dnB_dx_nt) +
                                float3::cross(dnA_dx_nt, dnB_dx_ms) +
                                float3::cross(nA, dnB_dx_mnst),
                            e) +
                float3::dot(float3::cross(dnA_dx_ms, nB) + float3::cross(nA, dnB_dx_ms),
                            de_dx[n][t]) +
                float3::dot(float3::cross(dnA_dx_nt, nB) + float3::cross(nA, dnB_dx_nt),
                            de_dx[m][s]);

            /* t ans s are swapped here due to the column wise storage of small matrices in
             * blender. */
            dC_dx_mn[m][n].values[t][s] = dcos_dx[n][t] * dsin_dx[m][s] + cos * dsin_dx_mnst -
                                          dsin_dx[n][t] * dcos_dx[m][s] - sin * dcos_dx_mnst;
          }
        }

        float3x3 dC_outer = float3x3::outer(dC_dx[m], dC_dx[n]);
        float3x3 df_dx_mn = -k * (dC_outer + dC_dx_mn[m][n] * C);

        force_derivatives[m][n] = df_dx_mn;
      }
    }
  }
};