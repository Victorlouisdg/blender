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

#include "BKE_cloth_simulator_baraff_witkin.hh"
#include "BKE_lib_query.h"
#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"
#include "BKE_modifier.h"
#include "BKE_pointcache.h"

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

#include "BLI_math_vector.h"
#include "BLI_span.hh"
#include "BLI_timeit.hh"

#include "DEG_depsgraph_query.h"

using blender::Span;

/* The static simulator is a temporary solution. I tried instancing the simulator in initData() and
 * storing a pointer to it in the ModifierData, but this lead to segfaults because the simulator
 * was destructed too early. */
static ClothSimulatorBaraffWitkin simulator = ClothSimulatorBaraffWitkin();

// UNUSED_VARS(modifier_data);

static void initData(ModifierData *modifier_data)
{
  UNUSED_VARS(modifier_data);

  // ClothBWModifierData *clothbw_modifier_data = reinterpret_cast<ClothBWModifierData *>(
  //     modifier_data);
  // ClothSimulatorBaraffWitkin simulator = ClothSimulatorBaraffWitkin();
  // clothbw_modifier_data->simulator_object = &simulator;
  // simulator = ClothSimulatorBaraffWitkin();
}

static bool dependsOnTime(ModifierData *UNUSED(modifier_data))
{
  return true;
}

static void updateDepsgraph(ModifierData *modifier_data, const ModifierUpdateDepsgraphContext *ctx)
{
  std::cout << __func__ << std::endl;
  ClothBWModifierData *clothbw_modifier_data = reinterpret_cast<ClothBWModifierData *>(
      modifier_data);
  DEG_add_modifier_to_transform_relation(ctx->node, "ClothBW Modifier");
  if (clothbw_modifier_data->object) {
    DEG_add_object_relation(
        ctx->node, clothbw_modifier_data->object, DEG_OB_COMP_GEOMETRY, "ClothBW Modifier");
    DEG_add_object_relation(
        ctx->node, clothbw_modifier_data->object, DEG_OB_COMP_TRANSFORM, "ClothBW Modifier");
  }
}

static void panel_draw(const bContext *UNUSED(C), Panel *panel)
{
  uiLayout *layout = panel->layout;
  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, nullptr);
  uiLayoutSetPropSep(layout, true);
  {
    uiLayout *col = uiLayoutColumn(layout, false);
    uiItemR(col, ptr, "object", 0, nullptr, ICON_NONE);
  }
  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  std::cout << __func__ << std::endl;
  modifier_panel_register(region_type, eModifierType_ClothBW, panel_draw);
}

static Mesh *modifyMesh(ModifierData *modifier_data, const ModifierEvalContext *ctx, Mesh *mesh)
{
  ClothBWModifierData *clothbw_modifier_data = reinterpret_cast<ClothBWModifierData *>(
      modifier_data);
  UNUSED_VARS(clothbw_modifier_data);

  // ClothSimulatorBaraffWitkin *simulator = reinterpret_cast<ClothSimulatorBaraffWitkin *>(
  //     clothbw_modifier_data->simulator_object);

  /* TODO: figure out how the caching system works. */
  int framenr = DEG_get_ctime(ctx->depsgraph);

  if (framenr == 1) {
    simulator.initialize(*mesh);
    return mesh;
  }

  simulator.step();

  for (int i = 0; i < mesh->totvert; i++) {
    mesh->mvert[i].co[2] += 0.1;
  }

  return mesh;
}

ModifierTypeInfo modifierType_ClothBW = {
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
    /* dependsOnTime */ dependsOnTime,
    /* dependsOnNormals */ nullptr,
    /* foreachIDLink */ nullptr,
    /* foreachTexLink */ nullptr,
    /* freeRuntimeData */ nullptr,
    /* panelRegister */ panelRegister,
    /* blendWrite */ nullptr,
    /* blendRead */ nullptr,
};
