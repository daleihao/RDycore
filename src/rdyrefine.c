#include <muParserDLL.h>
#include <petscdmceed.h>
#include <petscdmplex.h>
#include <petscsys.h>
#include <private/rdycoreimpl.h>
#include <private/rdydmimpl.h>
#include <private/rdymathimpl.h>

#define MAX_COMP_NAME_LENGTH 20

/// This function creates one Section with 3 DOFs for SWE.
static PetscErrorCode CreateSectionForSWE(RDy rdy, PetscSection *sec) {
  PetscInt n_field                             = 1;
  PetscInt n_field_comps[1]                    = {3};
  char     comp_names[3][MAX_COMP_NAME_LENGTH] = {
      "Height",
      "MomentumX",
      "MomentumY",
  };

  PetscFunctionBeginUser;
  PetscCall(PetscSectionCreate(rdy->comm, sec));
  PetscCall(PetscSectionSetNumFields(*sec, n_field));
  PetscInt n_field_dof_tot = 0;
  for (PetscInt f = 0; f < n_field; ++f) {
    PetscCall(PetscSectionSetFieldComponents(*sec, f, n_field_comps[f]));
    for (PetscInt c = 0; c < n_field_comps[f]; ++c, ++n_field_dof_tot) {
      PetscCall(PetscSectionSetComponentName(*sec, f, c, comp_names[c]));
    }
  }

  // set the number of degrees of freedom in each cell
  PetscInt c_start, c_end;  // starting and ending cell points
  PetscCall(DMPlexGetHeightStratum(rdy->dm, 0, &c_start, &c_end));
  PetscCall(PetscSectionSetChart(*sec, c_start, c_end));
  for (PetscInt c = c_start; c < c_end; ++c) {
    for (PetscInt f = 0; f < n_field; ++f) {
      PetscCall(PetscSectionSetFieldDof(*sec, c, f, n_field_comps[f]));
    }
    PetscCall(PetscSectionSetDof(*sec, c, n_field_dof_tot));
  }
  PetscCall(PetscSectionSetUp(*sec));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode InitRegionsFromCoarseDM(RDy rdy_coarse, RDy rdy_fine) {
  PetscFunctionBegin;

  rdy_fine->num_regions = rdy_coarse->num_regions;

  PetscCall(PetscCalloc1(rdy_fine->num_regions, &rdy_fine->regions));

  for (PetscInt r = 0; r < rdy_coarse->num_regions; ++r) {
    RDyRegion *region_coarse = &rdy_coarse->regions[r];
    RDyRegion *region_fine   = &rdy_fine->regions[r];

    region_fine->id = region_coarse->id;
    strcpy(region_fine->name, region_coarse->name);

    PetscInt num_refined_cells = 4; // assuming homogeneous refinement over the entire domain
    region_fine->num_cells = region_coarse->num_cells * num_refined_cells;

    if (region_fine->num_cells > 0) {
      PetscCall(PetscCalloc1(region_fine->num_cells, &region_fine->cell_ids));

      PetscInt count = 0;
      for (PetscInt i = 0; i < region_coarse->num_cells; ++i) {
        for (PetscInt j = 0; j < num_refined_cells; j++) {
          region_fine->cell_ids[count] = region_coarse->cell_ids[i] * num_refined_cells + j;
          count++;
        }
      }
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyRefine (RDy rdy_coarse) {
  PetscFunctionBegin;

  printf("++++++++++++++++++++++++++++++++++\n");
  printf("In RDyRefine\n");
  printf("++++++++++++++++++++++++++++++++++\n");

  DM dm_fine;
  PetscCall(DMRefine(rdy_coarse->dm, PETSC_COMM_WORLD, &dm_fine));
  PetscCall(DMCopyDisc(rdy_coarse->dm, dm_fine));
  PetscCall(DMViewFromOptions(dm_fine, NULL, "-dm_fine_view"));
  printf("\n");

  RDy rdy_fine;
  PetscCall(RDyCreate(rdy_coarse->global_comm, rdy_coarse->config_file, &rdy_fine));

  // note: default config values are specified in the YAML input schema!
  PetscCall(ReadConfigFile(rdy_fine));

  rdy_fine->log = stdout;

  // override parameters using command line arguments
  PetscCall(OverrideParameters(rdy_fine));

  // print configuration info
  PetscCall(PrintConfig(rdy_fine));

  // ++++++++++++++++++++++++++++++++++++++++++++++++++
  // Create DM
  RDyLogDebug(rdy_coarse, "FINE: Creating DMs...");
  rdy_fine->dm = dm_fine;
  PetscSection sec;
  PetscCall(CreateSectionForSWE(rdy_fine, &sec));
  // embed the section's data in our grid and toss the section
  PetscCall(DMSetLocalSection(rdy_fine->dm, sec));
  PetscCall(DMViewFromOptions(rdy_fine->dm, NULL, "-dm_fine_view"));
  // ++++++++++++++++++++++++++++++++++++++++++++++++++

  PetscCall(CreateAuxiliaryDM(rdy_fine));  // for diagnostics
  PetscCall(CreateVectors(rdy_fine));      // global and local vectors, residuals

  RDyLogDebug(rdy_coarse, "FINE: Initializing regions...");
  PetscCall(InitRegionsFromCoarseDM(rdy_coarse, rdy_fine));

  RDyLogDebug(rdy_coarse, "FINE: Initializing initial conditions and sources...");
  PetscCall(InitInitialConditions(rdy_fine));
  PetscCall(InitSources(rdy_fine));

  RDyLogDebug(rdy_coarse, "FINE: Creating FV mesh...");
  PetscCall(RDyMeshCreateFromDM(rdy_fine->dm, &rdy_fine->mesh));

  RDyLogDebug(rdy_coarse, "FINE Initializing boundaries and boundary conditions...");
  PetscCall(InitBoundaries(rdy_fine));
  //PetscCall(InitBoundaryConditions(rdy_fine));

  PetscCall(RDyMeshCreateFromDM(rdy_fine->dm, &rdy_fine->mesh));

  PetscFunctionReturn(PETSC_SUCCESS);
}