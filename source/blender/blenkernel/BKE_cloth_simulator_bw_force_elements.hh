#include "BLI_float3.hh"
#include "BLI_float3x3.hh"

#include "BLI_index_range.hh"

using blender::float3;
using blender::float3x3;
using blender::IndexRange;

class StretchForceElementBW {
 public:
  float3 forces[3];
  float3x3 force_derivatives[3][3];

  void calculate(
      float ku, float kv, float area_factor, float3 wu, float3 wv, float3 dwu_dx, float3 dwv_dx)
  {
    set_zero();
    calculate_single_stretch_direction(ku, area_factor, wu, dwu_dx);
    calculate_single_stretch_direction(kv, area_factor, wv, dwv_dx);
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

    float3 dC_dx[3];

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