#include "BLI_float3.hh"

using blender::float3;

/* These values where generated with a python script, based on the procedure described in:
 * Implementing Baraff Witkin by David Pritchard.
 *
 * These normal second derivatives are needed for the second derivatives of the bending condtion.
 * They are used for dot and vector product with other vectors whose values are only known at
 * runtime. It is unfortunate that about 90% of these values are 0.0, which means many floating
 * points operations are being wasted.
 */
const float3 normalA_second_derivatives[4][4][3][3] = {
    {
        {
            {
                /* m = 0, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 0, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 0, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 0, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, -1.0f), /* t = 1 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 0, n = 1, s = 1 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 1, s = 2 */
                float3(0.0f, -1.0f, 0.0f), /* t = 0 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 0, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 1 */
                float3(0.0f, -1.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 2, s = 1 */
                float3(0.0f, 0.0f, -1.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 0, n = 2, s = 2 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 0 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 0, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 3, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 3, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
    },
    {
        {
            {
                /* m = 1, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 1 */
                float3(0.0f, -1.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 0, s = 1 */
                float3(0.0f, 0.0f, -1.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 1, n = 0, s = 2 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 0 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 1, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 1, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 1, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 1, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, -1.0f), /* t = 1 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 1, n = 2, s = 1 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 2, s = 2 */
                float3(0.0f, -1.0f, 0.0f), /* t = 0 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 1, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 3, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 3, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
    },
    {
        {
            {
                /* m = 2, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, -1.0f), /* t = 1 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 2, n = 0, s = 1 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 0, s = 2 */
                float3(0.0f, -1.0f, 0.0f), /* t = 0 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 2, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 1 */
                float3(0.0f, -1.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 1, s = 1 */
                float3(0.0f, 0.0f, -1.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 2, n = 1, s = 2 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 0 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 2, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 2, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 2, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 2, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 3, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 3, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
    },
    {
        {
            {
                /* m = 3, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 0, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 0, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 3, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 1, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 1, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 3, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 2, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 2, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 3, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 3, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 3, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
    },
};

const float3 normalB_second_derivatives[4][4][3][3] = {
    {
        {
            {
                /* m = 0, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 0, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 0, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 0, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 1, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 1, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 0, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 2, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 2, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 0, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 3, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 0, n = 3, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
    },
    {
        {
            {
                /* m = 1, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 0, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 0, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 1, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 1, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 1, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 1, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 1 */
                float3(0.0f, -1.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 2, s = 1 */
                float3(0.0f, 0.0f, -1.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 1, n = 2, s = 2 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 0 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 1, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, -1.0f), /* t = 1 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 1, n = 3, s = 1 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 1, n = 3, s = 2 */
                float3(0.0f, -1.0f, 0.0f), /* t = 0 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
    },
    {
        {
            {
                /* m = 2, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 0, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 0, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 2, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, -1.0f), /* t = 1 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 2, n = 1, s = 1 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 1, s = 2 */
                float3(0.0f, -1.0f, 0.0f), /* t = 0 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 2, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 2, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 2, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 2, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 1 */
                float3(0.0f, -1.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 2, n = 3, s = 1 */
                float3(0.0f, 0.0f, -1.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 2, n = 3, s = 2 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 0 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
    },
    {
        {
            {
                /* m = 3, n = 0, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 0, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 0, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
        {
            {
                /* m = 3, n = 1, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 1 */
                float3(0.0f, -1.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 1, s = 1 */
                float3(0.0f, 0.0f, -1.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 3, n = 1, s = 2 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 0 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 3, n = 2, s = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 0 */
                float3(0.0f, 0.0f, -1.0f), /* t = 1 */
                float3(0.0f, 1.0f, 0.0f),  /* t = 2 */
            },
            {
                /* m = 3, n = 2, s = 1 */
                float3(0.0f, 0.0f, 1.0f),  /* t = 0 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(-1.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 2, s = 2 */
                float3(0.0f, -1.0f, 0.0f), /* t = 0 */
                float3(1.0f, 0.0f, 0.0f),  /* t = 1 */
                float3(0.0f, 0.0f, 0.0f),  /* t = 2 */
            },
        },
        {
            {
                /* m = 3, n = 3, s = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 3, s = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
            {
                /* m = 3, n = 3, s = 2 */
                float3(0.0f, 0.0f, 0.0f), /* t = 0 */
                float3(0.0f, 0.0f, 0.0f), /* t = 1 */
                float3(0.0f, 0.0f, 0.0f), /* t = 2 */
            },
        },
    },
};