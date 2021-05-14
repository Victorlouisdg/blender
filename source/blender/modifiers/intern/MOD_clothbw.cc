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

/** \file
 * \ingroup modifiers
 */

#include <vector>

#include "BKE_lib_query.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_volume.h"
#include "BKE_volume_to_mesh.hh"

#include "MOD_modifiertypes.h"
#include "MOD_ui_common.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_screen_types.h"
#include "DNA_volume_types.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "RNA_access.h"

#include "BLI_float4x4.hh"
#include "BLI_math_vector.h"
#include "BLI_span.hh"
#include "BLI_timeit.hh"

#include "DEG_depsgraph_query.h"

static void initData(ModifierData *md)
{
  ClothBWModifierData *cbwmd = reinterpret_cast<ClothBWModifierData *>(cbwmd);
  cbwmd->object = nullptr;
}

static void updateDepsgraph(ModifierData *md, const ModifierUpdateDepsgraphContext *ctx)
{
  //   VolumeToMeshModifierData *vmmd = reinterpret_cast<VolumeToMeshModifierData *>(md);
  //   DEG_add_modifier_to_transform_relation(ctx->node, "Volume to Mesh Modifier");
  //   if (vmmd->object) {
  //     DEG_add_object_relation(
  //         ctx->node, vmmd->object, DEG_OB_COMP_GEOMETRY, "Volume to Mesh Modifier");
  //     DEG_add_object_relation(
  //         ctx->node, vmmd->object, DEG_OB_COMP_TRANSFORM, "Volume to Mesh Modifier");
  //   }
}

static void foreachIDLink(ModifierData *md, Object *ob, IDWalkFunc walk, void *userData)
{
  ClothBWModifierData *cbwmd = reinterpret_cast<ClothBWModifierData *>(cbwmd);
  walk(userData, ob, (ID **)&cbwmd->object, IDWALK_CB_NOP);
}

static void panel_draw(const bContext *UNUSED(C), Panel *panel)
{
  uiLayout *layout = panel->layout;

  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, nullptr);
  ClothBWModifierData *cbwmd = static_cast<ClothBWModifierData *>(ptr->data);

  uiLayoutSetPropSep(layout, true);

  {
    uiLayout *col = uiLayoutColumn(layout, false);
    uiItemR(col, ptr, "object", 0, nullptr, ICON_NONE);
  }

  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  modifier_panel_register(region_type, eModifierType_ClothBW, panel_draw);
}

static Mesh *create_empty_mesh(const Mesh *input_mesh)
{
  Mesh *new_mesh = BKE_mesh_new_nomain(0, 0, 0, 0, 0);
  BKE_mesh_copy_settings(new_mesh, input_mesh);
  return new_mesh;
}

static Mesh *modifyMesh(ModifierData *md, const ModifierEvalContext *ctx, Mesh *input_mesh)
{
  UNUSED_VARS(md);
  BKE_modifier_set_error(ctx->object, md, "Compiled without OpenVDB");
  return create_empty_mesh(input_mesh);
}

ModifierTypeInfo modifierType_VolumeToMesh = {
    /* name */ "ClothBW",
    /* structName */ "ClothBWModifierData",
    /* structSize */ sizeof(ClothBWModifierData),
    /* srna */ &RNA_ClothBWModifier,
    /* type */ eModifierTypeType_Constructive,
    /* flags */ eModifierTypeFlag_AcceptsMesh,
    /* icon */ ICON_MOD_CLOTH,

    /* copyData */ BKE_modifier_copydata_generic,

    /* deformVerts */ nullptr,
    /* deformMatrices */ nullptr,
    /* deformVertsEM */ nullptr,
    /* deformMatricesEM */ nullptr,
    /* modifyMesh */ modifyMesh,
    /* modifyHair */ nullptr,
    /* modifyGeometrySet */ nullptr,

    /* initData */ initData,
    /* requiredDataMask */ nullptr,
    /* freeData */ nullptr,
    /* isDisabled */ nullptr,
    /* updateDepsgraph */ updateDepsgraph,
    /* dependsOnTime */ nullptr,
    /* dependsOnNormals */ nullptr,
    /* foreachIDLink */ foreachIDLink,
    /* foreachTexLink */ nullptr,
    /* freeRuntimeData */ nullptr,
    /* panelRegister */ panelRegister,
    /* blendWrite */ nullptr,
    /* blendRead */ nullptr,
};
