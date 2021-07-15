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

#include "BKE_cloth_simulator_bw.hh"

#include "BKE_deform.h"
#include "BKE_lib_query.h"
#include "BKE_mesh.h"
#include "BKE_mesh_runtime.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_pointcache.h"

#include "MOD_modifiertypes.h"
#include "MOD_ui_common.h"
#include "MOD_util.h"

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

static void initData(ModifierData *md)
{
  std::cout << "initializing Cloth BW data" << std::endl;
  ClothBWModifierData *d = reinterpret_cast<ClothBWModifierData *>(md);
  d->n_substeps = 5;

  d->stretch_stiffness = 200.0f;
  d->shear_stiffness = 50.0f;
  d->bend_stiffness = 0.01f;
  d->spring_stiffness = 2000.0f;

  d->stretch_damping_factor = 0.1f;
  d->shear_damping_factor = 0.1f;
  d->bend_damping_factor = 0.1f;
  d->spring_damping_factor = 0.1f;

  d->enable_shear = true;
  d->enable_bend = true;

  d->damp_stretch = true;
  d->damp_shear = true;
  d->damp_bend = true;
  d->damp_springs = true;

  d->use_explicit_integration = false;
}

static void freeData(ModifierData *md)
{
  std::cout << "freeing Cloth BW data" << std::endl;
  ClothBWModifierData *d = reinterpret_cast<ClothBWModifierData *>(md);
  if (d->simulator_object) {
    delete reinterpret_cast<ClothSimulatorBW *>(d->simulator_object);
  }
}

static bool dependsOnTime(ModifierData *UNUSED(modifier_data))
{
  return true;
}

static void updateDepsgraph(ModifierData *modifier_data, const ModifierUpdateDepsgraphContext *ctx)
{
  std::cout << __func__ << std::endl;
  ClothBWModifierData *cbw_modifier_data = reinterpret_cast<ClothBWModifierData *>(modifier_data);
  DEG_add_modifier_to_transform_relation(ctx->node, "ClothBW Modifier");
  if (cbw_modifier_data->collision_object) {
    DEG_add_object_relation(
        ctx->node, cbw_modifier_data->collision_object, DEG_OB_COMP_GEOMETRY, "ClothBW Modifier");
    DEG_add_object_relation(
        ctx->node, cbw_modifier_data->collision_object, DEG_OB_COMP_TRANSFORM, "ClothBW Modifier");
  }
}

static void panel_draw(const bContext *UNUSED(C), Panel *panel)
{
  uiLayout *layout = panel->layout;

  PointerRNA ob_ptr;
  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, &ob_ptr);

  uiLayoutSetPropSep(layout, true);
  {
    uiLayout *col = uiLayoutColumn(layout, false);
    uiItemR(col, ptr, "collision_object", 0, nullptr, ICON_NONE);
    modifier_vgroup_ui(layout, ptr, &ob_ptr, "vertex_group", NULL, NULL);
    uiItemR(layout, ptr, "n_substeps", 0, NULL, ICON_NONE);

    uiItemR(layout, ptr, "stretch_stiffness", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "shear_stiffness", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "bend_stiffness", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "spring_stiffness", 0, NULL, ICON_NONE);

    uiItemR(layout, ptr, "stretch_damping_factor", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "shear_damping_factor", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "bend_damping_factor", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "spring_damping_factor", 0, NULL, ICON_NONE);

    uiItemR(layout, ptr, "enable_shear", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "enable_bend", 0, NULL, ICON_NONE);

    uiItemR(layout, ptr, "damp_stretch", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "damp_shear", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "damp_bend", 0, NULL, ICON_NONE);
    uiItemR(layout, ptr, "damp_springs", 0, NULL, ICON_NONE);

    uiItemR(layout, ptr, "use_explicit_integration", 0, NULL, ICON_NONE);
  }
  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  std::cout << __func__ << std::endl;
  modifier_panel_register(region_type, eModifierType_ClothBW, panel_draw);
}

static Mesh *modifyMesh(ModifierData *md, const ModifierEvalContext *ctx, Mesh *mesh)
{
  ClothBWModifierData *modifier_data = reinterpret_cast<ClothBWModifierData *>(md);

  if (!modifier_data->simulator_object) {
    modifier_data->simulator_object = new ClothSimulatorBW();
  }

  ClothSimulatorBW *simulator = reinterpret_cast<ClothSimulatorBW *>(
      modifier_data->simulator_object);

  /* TODO: figure out how the caching system works. */
  int framenr = DEG_get_ctime(ctx->depsgraph);

  /* Currently added the modifier on a frame the is not 1 results in a crash because of this. */
  if (framenr == 1) {
    simulator->initialize(*mesh, *modifier_data, *(ctx->object));

    // Object *collision_object = modifier_data->collision_object;
    // if (collision_object) {
    //   Mesh *collision_mesh = BKE_object_get_pre_modified_mesh(collision_object);
    //   simulator->set_collision_mesh(*collision_mesh);
    // }
    return mesh;
  }

  simulator->step();

  for (int i : IndexRange(mesh->totvert)) {
    // Make a cleaner interface to get the vertex positions.
    mesh->mvert[i].co[0] = simulator->vertex_positions_eigen[3 * i];
    mesh->mvert[i].co[1] = simulator->vertex_positions_eigen[3 * i + 1];
    mesh->mvert[i].co[2] = simulator->vertex_positions_eigen[3 * i + 2];
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
    /* freeData */ freeData,
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
