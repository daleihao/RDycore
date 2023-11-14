#include <petscdmceed.h>
#include <petscdmplex.h>
#include <private/rdycoreimpl.h>
#include <private/rdysweimpl.h>
#include <rdycore.h>

// Maximum length of the name of a prognostic or diagnostic field component
#define MAX_COMP_NAME_LENGTH 20

// time conversion factors
static const PetscReal secs_in_min = 60.0;
static const PetscReal mins_in_hr  = 60.0;
static const PetscReal hrs_in_day  = 24.0;
static const PetscReal days_in_mon = 30.0;
static const PetscReal days_in_yr  = 365.0;

/// Returns a string corresponding to the given time unit
const char *TimeUnitAsString(RDyTimeUnit time_unit) {
  static const char *time_unit_strings[6] = {
      "sec",  // seconds
      "min",  // minutes
      "hr",   // hours
      "day",  // days
      "mon",  // months
      "yr",   // years
  };
  PetscFunctionBegin;
  PetscFunctionReturn(time_unit_strings[time_unit]);
}

/// Converts the given time (expressed in the given units) to seconds.
/// @param [in] time the time as expressed in the given units
/// @param [in] time_unit the units in which the time is expressed
PetscReal ConvertTimeToSeconds(PetscReal time, RDyTimeUnit time_unit) {
  PetscFunctionBegin;

  PetscReal time_in_sec;
  switch (time_unit) {
    case TIME_SECONDS:
      time_in_sec = time;
      break;
    case TIME_MINUTES:
      time_in_sec = time * secs_in_min;
      break;
    case TIME_HOURS:
      time_in_sec = time * mins_in_hr * secs_in_min;
      break;
    case TIME_DAYS:
      time_in_sec = time * hrs_in_day * mins_in_hr * secs_in_min;
      break;
    case TIME_MONTHS:
      time_in_sec = time * days_in_mon * hrs_in_day * mins_in_hr * secs_in_min;
      break;
    case TIME_YEARS:
      time_in_sec = time * days_in_yr * hrs_in_day * mins_in_hr * secs_in_min;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unsupported time unit");
      break;
  }

  PetscFunctionReturn(time_in_sec);
}

/// Converts the given time (expressed in seconds) to the given units.
/// @param [in] time the time as expressed in seconds
/// @param [in] time_unit the units to which the time is to be converted
PetscReal ConvertTimeFromSeconds(PetscReal time, RDyTimeUnit time_unit) {
  PetscFunctionBegin;
  PetscReal time_in_units;
  switch (time_unit) {
    case TIME_SECONDS:
      time_in_units = time;
      break;
    case TIME_MINUTES:
      time_in_units = time / secs_in_min;
      break;
    case TIME_HOURS:
      time_in_units = time / mins_in_hr / secs_in_min;
      break;
    case TIME_DAYS:
      time_in_units = time / hrs_in_day / mins_in_hr / secs_in_min;
      break;
    case TIME_MONTHS:
      time_in_units = time / days_in_mon / hrs_in_day / mins_in_hr / secs_in_min;
      break;
    case TIME_YEARS:
      time_in_units = time / days_in_yr / hrs_in_day / mins_in_hr / secs_in_min;
      break;
    default:
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unsupported time unit");
      break;
  }

  PetscFunctionReturn(time_in_units);
}

// overrides parameters with command line arguments
static PetscErrorCode OverrideParameters(RDy rdy) {
  PetscFunctionBegin;

  if (rdy->dt <= 0.0) {
    // ѕet a default timestep if needed
    rdy->dt = ConvertTimeToSeconds(rdy->config.time.final_time, rdy->config.time.unit) / rdy->config.time.max_step;
  } else {
    // convert dt to seconds in any case
    rdy->dt = ConvertTimeToSeconds(rdy->dt, rdy->config.time.unit);
  }

  PetscOptionsBegin(rdy->comm, NULL, "RDycore options", "");
  {
    PetscCall(PetscOptionsReal("-dt", "dt (seconds)", "", rdy->dt, &rdy->dt, NULL));
    PetscCall(PetscOptionsString("-ceed", "Ceed resource (/cpu/self, /gpu/cuda, /gpu/hip, ...)", "", rdy->ceed_resource, rdy->ceed_resource,
                                 sizeof rdy->ceed_resource, NULL));
  }
  PetscOptionsEnd();

  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CreateSection(RDy rdy, PetscSection *sec) {
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

static PetscErrorCode CreateDM(RDy rdy) {
  PetscSection sec;
  PetscMPIInt  size;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(rdy->comm, &size));

  PetscCall(DMCreate(rdy->comm, &rdy->dm));
  PetscCall(DMSetType(rdy->dm, DMPLEX));

  // if we're using CEED, set Vec and Mat types based on the selected backend
  if (rdy->ceed_resource[0]) {
    CeedMemType mem_type_backend;
    PetscCallCEED(CeedGetPreferredMemType(rdy->ceed, &mem_type_backend));
    VecType vec_type = NULL;
    switch (mem_type_backend) {
      case CEED_MEM_HOST:
        vec_type = VECSTANDARD;
        break;
      case CEED_MEM_DEVICE: {
        const char *resolved;
        PetscCallCEED(CeedGetResource(rdy->ceed, &resolved));
        if (strstr(resolved, "/gpu/cuda")) vec_type = VECCUDA;
        else if (strstr(resolved, "/gpu/hip")) vec_type = VECKOKKOS;
        else if (strstr(resolved, "/gpu/sycl")) vec_type = VECKOKKOS;
        else vec_type = VECSTANDARD;
      }
    }
    PetscCall(DMSetVecType(rdy->dm, vec_type));

    MatType mat_type = NULL;
    if (strstr(vec_type, VECCUDA)) mat_type = MATAIJCUSPARSE;
    else if (strstr(vec_type, VECKOKKOS)) mat_type = MATAIJKOKKOS;
    else mat_type = MATAIJ;
    PetscCall(DMSetMatType(rdy->dm, mat_type));
  }

  PetscCall(DMPlexDistributeSetDefault(rdy->dm, PETSC_FALSE));
  PetscCall(DMSetFromOptions(rdy->dm));

  // name the grid
  PetscCall(PetscObjectSetName((PetscObject)rdy->dm, "grid"));

  // NOTE Need to create section before distribution, so that natural map can be created
  // create a section with (h, hu, hv) as degrees of freedom
  if (!rdy->refine) {
    PetscCall(CreateSection(rdy, &sec));
    // embed the section's data in our grid and toss the section
    PetscCall(DMSetLocalSection(rdy->dm, sec));
  }

  // distribution phase
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)rdy->dm, "dist_"));
  PetscCall(DMPlexDistributeSetDefault(rdy->dm, PETSC_TRUE));
  PetscCall(DMSetFromOptions(rdy->dm));
  PetscCall(DMViewFromOptions(rdy->dm, NULL, "-dm_view"));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)rdy->dm, NULL));

  // parallel refinement phase
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)rdy->dm, "ref_"));
  PetscCall(DMPlexDistributeSetDefault(rdy->dm, PETSC_FALSE));
  PetscCall(DMSetFromOptions(rdy->dm));
  PetscCall(DMViewFromOptions(rdy->dm, NULL, "-dm_view"));
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)rdy->dm, NULL));

  // Overlap meshes after refinement
  if (size > 1) {
    DM      dmOverlap;
    PetscSF sfOverlap, sfMigration, sfMigrationNew;

    PetscCall(DMPlexGetMigrationSF(rdy->dm, &sfMigration));
    PetscCall(DMPlexDistributeOverlap(rdy->dm, 1, &sfOverlap, &dmOverlap));
    PetscCall(DMPlexRemapMigrationSF(sfOverlap, sfMigration, &sfMigrationNew));
    PetscCall(PetscSFDestroy(&sfOverlap));
    PetscCall(DMPlexSetMigrationSF(dmOverlap, sfMigrationNew));
    PetscCall(PetscSFDestroy(&sfMigrationNew));
    PetscCall(DMDestroy(&rdy->dm));
    rdy->dm = dmOverlap;
  }

  // mark boundary edges so we can enforce reflecting BCs on them if needed
  {
    DMLabel boundary_edges;
    PetscCall(DMCreateLabel(rdy->dm, "boundary_edges"));
    PetscCall(DMGetLabel(rdy->dm, "boundary_edges", &boundary_edges));
    PetscCall(DMPlexMarkBoundaryFaces(rdy->dm, 1, boundary_edges));
  }

  // create parallel section and global-to-natural mapping
  if (rdy->refine) {
    PetscCall(CreateSection(rdy, &sec));
    PetscCall(DMSetLocalSection(rdy->dm, sec));
  } else if (size > 1) {
    PetscSF      sfMigration, sfNatural;
    PetscSection psec;
    PetscInt    *remoteOffsets;

    PetscCall(DMPlexGetMigrationSF(rdy->dm, &sfMigration));
    PetscCall(DMPlexCreateGlobalToNaturalSF(rdy->dm, sec, sfMigration, &sfNatural));
    PetscCall(DMPlexSetGlobalToNaturalSF(rdy->dm, sfNatural));
    PetscCall(PetscSFDestroy(&sfNatural));

    PetscCall(PetscSectionCreate(rdy->comm, &psec));
    PetscCall(PetscSFDistributeSection(sfMigration, sec, &remoteOffsets, psec));
    PetscCall(DMSetLocalSection(rdy->dm, psec));
    PetscCall(PetscFree(remoteOffsets));
    PetscCall(PetscSectionDestroy(&sec));
    PetscCall(PetscSectionDestroy(&psec));
  }
  PetscCall(PetscSectionDestroy(&sec));

  // set grid adacency
  PetscCall(DMSetBasicAdjacency(rdy->dm, PETSC_TRUE, PETSC_TRUE));

  PetscCall(DMViewFromOptions(rdy->dm, NULL, "-dm_view"));

  PetscFunctionReturn(PETSC_SUCCESS);
}

// retrieves the index of a flow condition using its name
static PetscErrorCode FindFlowCondition(RDy rdy, const char *name, PetscInt *index) {
  PetscFunctionBegin;

  // Currently, we do a linear search on the name of the condition, which is O(N)
  // for N regions. If this is too slow, we can sort the conditions by name and
  // use binary search, which is O(log2 N).
  *index = -1;
  for (PetscInt i = 0; i < rdy->config.num_flow_conditions; ++i) {
    if (!strcmp(rdy->config.flow_conditions[i].name, name)) {
      *index = i;
      break;
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CreateAuxiliaryDM(RDy rdy) {
  PetscFunctionBegin;

  PetscCall(DMClone(rdy->dm, &rdy->aux_dm));

  // create an auxiliary section with a diagnostic parameter.
  PetscInt     n_aux_field            = 1;
  PetscInt     n_aux_field_dof[1]     = {1};
  char         aux_field_names[1][20] = {"Parameter"};
  PetscSection aux_sec;
  PetscCall(PetscSectionCreate(rdy->comm, &aux_sec));
  PetscCall(PetscSectionSetNumFields(aux_sec, n_aux_field));
  PetscInt n_aux_field_dof_tot = 0;
  for (PetscInt f = 0; f < n_aux_field; ++f) {
    PetscCall(PetscSectionSetFieldName(aux_sec, f, &aux_field_names[f][0]));
    PetscCall(PetscSectionSetFieldComponents(aux_sec, f, n_aux_field_dof[f]));
    n_aux_field_dof_tot += n_aux_field_dof[f];
  }

  // set the number of auxiliary degrees of freedom in each cell
  PetscInt c_start, c_end;  // starting and ending cell points
  DMPlexGetHeightStratum(rdy->dm, 0, &c_start, &c_end);
  PetscCall(PetscSectionSetChart(aux_sec, c_start, c_end));
  for (PetscInt c = c_start; c < c_end; ++c) {
    for (PetscInt f = 0; f < n_aux_field; ++f) {
      PetscCall(PetscSectionSetFieldDof(aux_sec, c, f, n_aux_field_dof[f]));
    }
    PetscCall(PetscSectionSetDof(aux_sec, c, n_aux_field_dof_tot));
  }

  // embed the section's data in the auxiliary DM and toss the section
  PetscCall(PetscSectionSetUp(aux_sec));
  PetscCall(DMSetLocalSection(rdy->aux_dm, aux_sec));
  PetscCall(PetscSectionViewFromOptions(aux_sec, NULL, "-aux_layout_view"));
  PetscCall(PetscSectionDestroy(&aux_sec));

  if (!rdy->refine) {
    // copy adjacency info from the primary DM
    PetscSF sf_migration, sf_natural;
    PetscCall(DMPlexGetMigrationSF(rdy->dm, &sf_migration));
    PetscCall(DMPlexCreateGlobalToNaturalSF(rdy->aux_dm, aux_sec, sf_migration, &sf_natural));
    PetscCall(DMPlexSetGlobalToNaturalSF(rdy->aux_dm, sf_natural));
    PetscCall(PetscSFDestroy(&sf_natural));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// retrieves the index of a sediment condition using its name
static PetscErrorCode FindSedimentCondition(RDy rdy, const char *name, PetscInt *index) {
  PetscFunctionBegin;

  // Currently, we do a linear search on the name of the condition, which is O(N)
  // for N regions. If this is too slow, we can sort the conditions by name and
  // use binary search, which is O(log2 N).
  *index = -1;
  for (PetscInt i = 0; i < rdy->config.num_sediment_conditions; ++i) {
    if (!strcmp(rdy->config.sediment_conditions[i].name, name)) {
      *index = i;
      break;
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// retrieves the index of a salinity condition using its name
static PetscErrorCode FindSalinityCondition(RDy rdy, const char *name, PetscInt *index) {
  PetscFunctionBegin;

  // Currently, we do a linear search on the name of the condition, which is O(N)
  // for N regions. If this is too slow, we can sort the conditions by name and
  // use binary search, which is O(log2 N).
  *index = -1;
  for (PetscInt i = 0; i < rdy->config.num_salinity_conditions; ++i) {
    if (!strcmp(rdy->config.salinity_conditions[i].name, name)) {
      *index = i;
      break;
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// initializes mesh region data
//   can be run after refinement
static PetscErrorCode InitRegions(RDy rdy) {
  PetscFunctionBegin;

  // Count and fetch regions.
  PetscInt c_start, c_end;  // starting and ending cell points
  PetscCall(DMPlexGetHeightStratum(rdy->dm, 0, &c_start, &c_end));
  DMLabel label;
  PetscCall(DMGetLabel(rdy->dm, "Cell Sets", &label));
  // If we didn't find any regions, we can't perform the simulation.
  PetscCheck(label, rdy->comm, PETSC_ERR_USER, "No regions (cell sets) found in grid! Cannot assign initial conditions.");
  PetscCall(DMLabelGetNumValues(label, &rdy->num_regions));
  PetscCheck(rdy->num_regions <= MAX_NUM_REGIONS, rdy->comm, PETSC_ERR_USER, "Number of regions in mesh (%" PetscInt_FMT ") exceeds MAX_NUM_REGIONS (%d)",
             rdy->num_regions, MAX_NUM_REGIONS);

  // fetch region IDs
  IS region_id_is;
  PetscCall(DMLabelGetValueIS(label, &region_id_is));
  const PetscInt *region_ids;
  PetscCall(ISGetIndices(region_id_is, &region_ids));

  // allocate and set region data
  PetscCall(PetscCalloc1(rdy->num_regions, &rdy->regions));
  for (PetscInt r = 0; r < rdy->num_regions; ++r) {
    PetscInt region_id = region_ids[r];
    IS       cell_is;  // cell index space
    PetscCall(DMLabelGetStratumIS(label, region_id, &cell_is));
    if (cell_is) {
      RDyRegion *region = &rdy->regions[r];
      region->id        = region_id;

      // fish the region's name out of our region config data
      for (PetscInt r1 = 0; r1 < rdy->config.num_regions; ++r1) {
        if (rdy->config.regions[r1].grid_region_id == region_id) {
          strcpy(region->name, rdy->config.regions[r1].name);
          break;
        }
      }

      PetscInt num_cells;
      PetscCall(ISGetLocalSize(cell_is, &num_cells));
      if (num_cells > 0) {
        RDyLogDebug(rdy, "  Found region %d (%d cells)", region_id, num_cells);
        region->num_cells = num_cells;
        PetscCall(PetscCalloc1(region->num_cells, &region->cell_ids));
      }
      const PetscInt *cell_ids;
      PetscCall(ISGetIndices(cell_is, &cell_ids));
      for (PetscInt i = 0; i < num_cells; ++i) {
        region->cell_ids[i] = cell_ids[i] - c_start;
      }
      PetscCall(ISRestoreIndices(cell_is, &cell_ids));
      PetscCall(ISDestroy(&cell_is));
    }
  }
  PetscCall(ISRestoreIndices(region_id_is, &region_ids));
  PetscCall(ISDestroy(&region_id_is));

  PetscFunctionReturn(PETSC_SUCCESS);
}

// initializes mesh boundary data
//   can be run after refinement
static PetscErrorCode InitBoundaries(RDy rdy) {
  PetscFunctionBegin;

  // Extract edges on the domain boundary.
  DMLabel boundary_edge_label;
  PetscCall(DMGetLabel(rdy->dm, "boundary_edges", &boundary_edge_label));
  IS boundary_edge_is;
  PetscCall(DMLabelGetStratumIS(boundary_edge_label, 1, &boundary_edge_is));
  PetscBool boundary_edge_present = (boundary_edge_is != NULL);

  // Keep track of whether edges on the domain boundary have been assigned to
  // any boundaries.
  IS unassigned_edges_is;
  if (boundary_edge_present) {
    ISDuplicate(boundary_edge_is, &unassigned_edges_is);
  }
  PetscInt unassigned_edge_boundary_id = 0;  // boundary ID for unassigned edges

  // Count boundaries. We rely on face sets in our grids to express
  // boundary conditions. All edges on the domain boundary not assigned to other
  // boundaries are assigned to a special boundary to which we apply reflecting
  // boundary conditions.
  PetscInt e_start, e_end;  // starting and ending edge points
  DMPlexGetHeightStratum(rdy->dm, 1, &e_start, &e_end);
  DMLabel label;
  PetscCall(DMGetLabel(rdy->dm, "Face Sets", &label));
  PetscInt        num_boundaries_in_file = 0;
  IS              boundary_id_is;
  const PetscInt *boundary_ids;
  if (label) {  // found face sets!
    PetscCall(DMLabelGetNumValues(label, &num_boundaries_in_file));
    PetscCheck(num_boundaries_in_file <= MAX_NUM_BOUNDARIES, rdy->comm, PETSC_ERR_USER,
               "Number of boundaries in mesh (%" PetscInt_FMT ") exceeds MAX_NUM_BOUNDARIES (%d)", num_boundaries_in_file, MAX_NUM_BOUNDARIES);

    // fetch boundary IDs
    PetscCall(DMLabelGetValueIS(label, &boundary_id_is));
    PetscCall(ISGetIndices(boundary_id_is, &boundary_ids));

    for (PetscInt b = 0; b < num_boundaries_in_file; ++b) {
      PetscInt boundary_id = boundary_ids[b];
      IS       edge_is;
      PetscCall(DMLabelGetStratumIS(label, boundary_id, &edge_is));
      if (edge_is) {
        PetscInt num_edges;
        PetscCall(ISGetLocalSize(edge_is, &num_edges));
        RDyLogDebug(rdy, "  Found boundary %" PetscInt_FMT " (%" PetscInt_FMT " edges)", boundary_id, num_edges);

        // we can't use this boundary ID for unassigned edges
        if (unassigned_edge_boundary_id == boundary_id) ++unassigned_edge_boundary_id;

        // intersect this IS with our domain boundary IS to produce the edges
        // to subtract from our unassigned edge IS
        IS assigned_edges_is, new_unassigned_edges_is;
        PetscCall(ISIntersect(edge_is, boundary_edge_is, &assigned_edges_is));
        PetscCall(ISDifference(unassigned_edges_is, assigned_edges_is, &new_unassigned_edges_is));
        ISDestroy(&assigned_edges_is);
        ISDestroy(&unassigned_edges_is);
        unassigned_edges_is = new_unassigned_edges_is;

        ++rdy->num_boundaries;
      }
      PetscCall(ISDestroy(&edge_is));
    }
    PetscCall(ISRestoreIndices(boundary_id_is, &boundary_ids));
    PetscCall(ISDestroy(&boundary_id_is));
  }

  // add an additional boundary for unassigned boundary edges if needed
  PetscInt num_unassigned_edges = 0;
  if (boundary_edge_present) {
    PetscCall(ISGetLocalSize(unassigned_edges_is, &num_unassigned_edges));
  }
  if (num_unassigned_edges > 0) {
    RDyLogDebug(rdy, "Adding boundary %" PetscInt_FMT " for %" PetscInt_FMT " unassigned boundary edges", unassigned_edge_boundary_id,
                num_unassigned_edges);
    if (!label) {
      // create a "Face Sets" label if one doesn't already exist
      PetscCall(DMCreateLabel(rdy->dm, "Face Sets"));
      PetscCall(DMGetLabel(rdy->dm, "Face Sets", &label));
    }
    // add these edges to a new boundary with the given ID
    PetscCall(DMLabelSetStratumIS(label, unassigned_edge_boundary_id, unassigned_edges_is));
    ++rdy->num_boundaries;
  }
  if (boundary_edge_present) {
    PetscCall(ISDestroy(&boundary_edge_is));
    PetscCall(ISDestroy(&unassigned_edges_is));
  }

  // allocate resources for boundaries
  PetscCall(PetscCalloc1(rdy->num_boundaries, &rdy->boundaries));

  // now fetch boundary edge IDs
  if (label) {
    // fetch boundary IDs once again
    PetscCall(DMLabelGetValueIS(label, &boundary_id_is));
    PetscCall(ISGetIndices(boundary_id_is, &boundary_ids));

    for (PetscInt b = 0; b < rdy->num_boundaries; ++b) {
      PetscInt boundary_id = (b < num_boundaries_in_file) ? boundary_ids[b] : unassigned_edge_boundary_id;
      IS       edge_is;  // edge index space//
      PetscCall(DMLabelGetStratumIS(label, boundary_id, &edge_is));
      if (edge_is) {
        RDyBoundary *boundary = &rdy->boundaries[b];
        boundary->index       = b;
        boundary->id          = boundary_id;

        // fish the boundary's name out of our boundary config data
        for (PetscInt b1 = 0; b1 < rdy->config.num_boundaries; ++b1) {
          if (rdy->config.boundaries[b1].grid_boundary_id == boundary_id) {
            strcpy(boundary->name, rdy->config.boundaries[b1].name);
            break;
          }
        }

        // find the number of edges for this boundary
        PetscInt num_edges;
        PetscCall(ISGetLocalSize(edge_is, &num_edges));
        if (num_edges > 0) {
          boundary->num_edges = num_edges;
          PetscCall(PetscCalloc1(boundary->num_edges, &boundary->edge_ids));
        }

        // extract edge IDs
        const PetscInt *edge_ids;
        PetscCall(ISGetIndices(edge_is, &edge_ids));
        for (PetscInt i = 0; i < num_edges; ++i) {
          PetscCheck((edge_ids[i] >= e_start) && (edge_ids[i] < e_end), PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                     "Mesh point %" PetscInt_FMT " is not an edge. Likely the option -dm_plex_transform_label_match_strata is missing", edge_ids[i]);
          boundary->edge_ids[i] = edge_ids[i] - e_start;
        }
        PetscCall(ISRestoreIndices(edge_is, &edge_ids));
        PetscCall(ISDestroy(&edge_is));
      }
    }
    PetscCall(ISRestoreIndices(boundary_id_is, &boundary_ids));
    PetscCall(ISDestroy(&boundary_id_is));
  }

  // make sure we have at least one region and boundary across all mpi ranks
  PetscInt num_global_boundaries = 0;
  MPI_Allreduce(&rdy->num_boundaries, &num_global_boundaries, 1, MPI_INT, MPI_SUM, rdy->comm);
  PetscCheck(num_global_boundaries > 0, rdy->comm, PETSC_ERR_USER, "No boundaries were found in the grid!");

  PetscFunctionReturn(PETSC_SUCCESS);
}

// reads data for a single DOF from a binary file into a Vec
// TODO: currently, we only support binary PETSc files
static PetscErrorCode ReadOneDOFVecFromFile(RDy rdy, const char filename[], Vec *local) {
  PetscFunctionBegin;

  PetscViewer viewer;
  PetscCall(PetscViewerBinaryOpen(rdy->comm, filename, FILE_MODE_READ, &viewer));

  // create a naturally-ordered vector with a stride equal to the number of
  Vec natural, global;

  PetscCall(DMPlexCreateNaturalVector(rdy->aux_dm, &natural));
  PetscCall(DMCreateGlobalVector(rdy->aux_dm, &global));
  PetscCall(DMCreateLocalVector(rdy->aux_dm, local));

  // load the properties into the vector and copy them into place
  PetscCall(VecLoad(natural, viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  // scatter natural-to-global
  PetscCall(DMPlexNaturalToGlobalBegin(rdy->aux_dm, natural, global));
  PetscCall(DMPlexNaturalToGlobalEnd(rdy->aux_dm, natural, global));

  // scatter global-to-local
  PetscCall(DMGlobalToLocalBegin(rdy->aux_dm, global, INSERT_VALUES, *local));
  PetscCall(DMGlobalToLocalEnd(rdy->aux_dm, global, INSERT_VALUES, *local));

  PetscCall(VecDestroy(&natural));
  PetscCall(VecDestroy(&global));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#define READ_MATERIAL_PROPERTY(property, mat_props_spec, region, materials_by_cell) \
  if (mat_props_spec.property.file[0]) { /* read property from a file */            \
    Vec local;                                                                      \
    PetscCall(ReadOneDOFVecFromFile(rdy, mat_props_spec.property.file, &local));    \
    PetscScalar *x_ptr;                                                             \
    PetscCall(VecGetArray(local, &x_ptr));                                          \
    for (PetscInt c = 0; c < region.num_cells; ++c) {                               \
      PetscInt cell                    = region.cell_ids[c];                        \
      materials_by_cell[cell].property = x_ptr[c];                                  \
    }                                                                               \
    PetscCall(VecRestoreArray(local, &x_ptr));                                      \
    PetscCall(VecDestroy(&local));                                                  \
  } else {                                                                          \
    /* set this material property for all cells in the region */                    \
    for (PetscInt c = 0; c < region.num_cells; ++c) {                               \
      PetscInt cell                    = region.cell_ids[c];                        \
      materials_by_cell[cell].property = mat_props_spec.property.value;             \
    }                                                                               \
  }

// sets up materials
//   unsafe for refinement if file is given for surface composition
static PetscErrorCode InitMaterials(RDy rdy) {
  PetscFunctionBegin;

  // allocate storage for materials
  PetscCall(PetscCalloc1(rdy->mesh.num_cells, &rdy->materials_by_cell));

  // assign materials to each region as needed
  for (PetscInt r = 0; r < rdy->num_regions; ++r) {
    RDyRegion region = rdy->regions[r];

    // find the material specification corresponding to this region
    PetscInt region_mat_index = -1;
    for (PetscInt isurf_comp = 0; isurf_comp < rdy->config.num_material_assignments; ++isurf_comp) {
      if (!strcmp(rdy->config.surface_composition[isurf_comp].region, region.name)) {
        RDySurfaceCompositionSpec surf_comp = rdy->config.surface_composition[isurf_comp];
        for (PetscInt imat = 0; imat < rdy->config.num_materials; ++imat) {
          if (!strcmp(rdy->config.materials[imat].name, surf_comp.material)) {
            region_mat_index = imat;
            break;
          }
        }
        if (region_mat_index != -1) break;
      }
    }
    PetscCheck(region_mat_index != -1, rdy->comm, PETSC_ERR_USER, "Region '%s' has no assigned material!", region.name);
    RDyMaterialPropertiesSpec mat_props_spec = rdy->config.materials[region_mat_index].properties;

    // set the region's material properties
    READ_MATERIAL_PROPERTY(manning, mat_props_spec, region, rdy->materials_by_cell);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// sets up initial conditions
//   can be run after refinement
static PetscErrorCode InitInitialConditions(RDy rdy) {
  PetscFunctionBegin;

  // allocate storage for by-region initial conditions
  PetscCall(PetscCalloc1(rdy->num_regions, &rdy->initial_conditions));

  // assign іnitial conditions to each region as indicated in our config
  for (PetscInt r = 0; r < rdy->num_regions; ++r) {
    RDyCondition *ic              = &rdy->initial_conditions[r];
    RDyRegion     region          = rdy->regions[r];
    PetscInt      region_ic_index = -1;
    for (PetscInt ic = 0; ic < rdy->config.num_initial_conditions; ++ic) {
      if (!strcmp(rdy->config.initial_conditions[ic].region, region.name)) {
        region_ic_index = ic;
        break;
      }
    }
    PetscCheck(region_ic_index != -1, rdy->comm, PETSC_ERR_USER, "Region '%s' has no initial conditions!", region.name);

    RDyRegionConditionSpec ic_spec = rdy->config.initial_conditions[region_ic_index];
    PetscCheck(strlen(ic_spec.flow), rdy->comm, PETSC_ERR_USER, "Region '%s' has no initial flow condition!", region.name);
    PetscInt flow_index;
    PetscCall(FindFlowCondition(rdy, ic_spec.flow, &flow_index));
    RDyFlowCondition *flow_cond = &rdy->config.flow_conditions[flow_index];
    PetscCheck(flow_cond->type == CONDITION_DIRICHLET, rdy->comm, PETSC_ERR_USER,
               "initial flow condition '%s' for region '%s' is not of dirichlet type!", flow_cond->name, region.name);
    ic->flow = flow_cond;

    if (rdy->config.physics.sediment) {
      PetscCheck(strlen(ic_spec.sediment), rdy->comm, PETSC_ERR_USER, "Region '%s' has no initial sediment condition!", region.name);
      PetscInt sed_index;
      PetscCall(FindSedimentCondition(rdy, ic_spec.sediment, &sed_index));
      RDySedimentCondition *sed_cond = &rdy->config.sediment_conditions[sed_index];
      PetscCheck(sed_cond->type == CONDITION_DIRICHLET, rdy->comm, PETSC_ERR_USER,
                 "initial sediment condition '%s' for region '%s' is not of dirichlet type!", sed_cond->name, region.name);
      ic->sediment = sed_cond;
    }
    if (rdy->config.physics.salinity) {
      PetscCheck(strlen(ic_spec.salinity), rdy->comm, PETSC_ERR_USER, "Region '%s' has no initial salinity condition!", region.name);
      PetscInt sal_index;
      PetscCall(FindSalinityCondition(rdy, ic_spec.salinity, &sal_index));
      RDySalinityCondition *sal_cond = &rdy->config.salinity_conditions[sal_index];
      PetscCheck(sal_cond->type == CONDITION_DIRICHLET, rdy->comm, PETSC_ERR_USER,
                 "initial salinity condition '%s' for region '%s' is not of dirichlet type!", sal_cond->name, region.name);
      ic->salinity = sal_cond;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// sets up sources
//   can be run after refinement
static PetscErrorCode InitSources(RDy rdy) {
  PetscFunctionBegin;
  if (rdy->config.num_sources > 0) {
    // allocate storage for sources
    PetscCall(PetscCalloc1(rdy->num_regions, &rdy->sources));

    // assign sources to each region as needed
    for (PetscInt r = 0; r < rdy->num_regions; ++r) {
      RDyCondition *src              = &rdy->sources[r];
      RDyRegion     region           = rdy->regions[r];
      PetscInt      region_src_index = -1;
      for (PetscInt isrc = 0; isrc < rdy->config.num_sources; ++isrc) {
        if (!strcmp(rdy->config.sources[isrc].region, region.name)) {
          region_src_index = isrc;
          break;
        }
      }
      if (region_src_index != -1) {
        RDyRegionConditionSpec src_spec = rdy->config.sources[region_src_index];
        if (strlen(src_spec.flow)) {
          PetscInt flow_index;
          PetscCall(FindFlowCondition(rdy, src_spec.flow, &flow_index));
          RDyFlowCondition *flow_cond = &rdy->config.flow_conditions[flow_index];
          PetscCheck(flow_cond->type == CONDITION_DIRICHLET, rdy->comm, PETSC_ERR_USER, "flow source '%s' for region '%s' is not of dirichlet type!",
                     flow_cond->name, region.name);
          src->flow = flow_cond;
        }

        if (rdy->config.physics.sediment && strlen(src_spec.sediment)) {
          PetscInt sed_index;
          PetscCall(FindSedimentCondition(rdy, src_spec.sediment, &sed_index));
          RDySedimentCondition *sed_cond = &rdy->config.sediment_conditions[sed_index];
          PetscCheck(sed_cond->type == CONDITION_DIRICHLET, rdy->comm, PETSC_ERR_USER,
                     "sediment source '%s' for region '%s' is not of dirichlet type!", sed_cond->name, region.name);
          src->sediment = sed_cond;
        }
        if (rdy->config.physics.salinity && strlen(src_spec.salinity)) {
          PetscInt sal_index;
          PetscCall(FindSalinityCondition(rdy, src_spec.salinity, &sal_index));
          RDySalinityCondition *sal_cond = &rdy->config.salinity_conditions[sal_index];
          PetscCheck(sal_cond->type == CONDITION_DIRICHLET, rdy->comm, PETSC_ERR_USER,
                     "initial salinity condition '%s' for region '%s' is not of dirichlet type!", sal_cond->name, region.name);
          src->salinity = sal_cond;
        }
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// sets up boundary conditions
//   can be run after refinement
static PetscErrorCode InitBoundaryConditions(RDy rdy) {
  PetscFunctionBegin;
  // Set up a reflecting flow boundary condition.
  RDyFlowCondition *reflecting_flow = NULL;
  for (PetscInt c = 0; c < MAX_NUM_CONDITIONS; ++c) {
    if (!strlen(rdy->config.flow_conditions[c].name)) {
      reflecting_flow = &rdy->config.flow_conditions[c];
      strcpy((char *)reflecting_flow->name, "reflecting");
      reflecting_flow->type = CONDITION_REFLECTING;
      break;
    }
  }
  PetscCheck(reflecting_flow, rdy->comm, PETSC_ERR_USER, "Could not allocate a reflecting flow condition! Please increase MAX_BOUNDARY_ID.");

  // Allocate storage for boundary conditions.
  PetscCall(PetscCalloc1(rdy->num_boundaries, &rdy->boundary_conditions));

  // Assign a boundary condition to each boundary.
  for (PetscInt b = 0; b < rdy->num_boundaries; ++b) {
    RDyCondition *bc       = &rdy->boundary_conditions[b];
    RDyBoundary   boundary = rdy->boundaries[b];

    // identify the index of the boundary condition assigned to this boundary.
    PetscInt bc_index = -1;
    for (PetscInt ib = 0; ib < rdy->config.num_boundary_conditions; ++ib) {
      RDyBoundaryConditionSpec bc_spec = rdy->config.boundary_conditions[ib];
      for (PetscInt ib1 = 0; ib1 < bc_spec.num_boundaries; ++ib1) {
        if (!strcmp(bc_spec.boundaries[ib1], boundary.name)) {
          PetscCheck(bc_index == -1, rdy->comm, PETSC_ERR_USER, "Boundary '%s' is assigned to more than one boundary condition!", boundary.name);
          bc_index = ib;
          break;
        }
      }
    }
    if (bc_index != -1) {
      RDyBoundaryConditionSpec bc_spec = rdy->config.boundary_conditions[bc_index];

      // If no flow condition was specified for a boundary, we set it to our
      // reflecting flow condition.
      if (!strlen(bc_spec.flow)) {
        bc->flow = reflecting_flow;
      } else {
        PetscInt flow_index;
        PetscCall(FindFlowCondition(rdy, bc_spec.flow, &flow_index));
        bc->flow = &rdy->config.flow_conditions[flow_index];
      }

      if (rdy->config.physics.sediment) {
        PetscCheck(strlen(bc_spec.sediment), rdy->comm, PETSC_ERR_USER, "Boundary '%s' has no sediment boundary condition!", boundary.name);
        PetscInt sed_index;
        PetscCall(FindSedimentCondition(rdy, bc_spec.sediment, &sed_index));
        bc->sediment = &rdy->config.sediment_conditions[sed_index];
      }
      if (rdy->config.physics.salinity) {
        PetscCheck(strlen(bc_spec.salinity), rdy->comm, PETSC_ERR_USER, "Boundary '%s' has no salinity boundary condition!", boundary.name);
        PetscInt sal_index;
        PetscCall(FindSalinityCondition(rdy, bc_spec.salinity, &sal_index));
        bc->salinity = &rdy->config.salinity_conditions[sal_index];
      }
    } else {
      // this boundary wasn't explicitly requested, so set up an auto-generated
      // reflecting BC
      bc->flow           = reflecting_flow;
      bc->auto_generated = PETSC_TRUE;
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// create solvers and vectors
static PetscErrorCode CreateSolvers(RDy rdy) {
  PetscFunctionBegin;

  // set up vectors
  PetscCall(DMCreateGlobalVector(rdy->dm, &rdy->X));
  PetscCall(VecDuplicate(rdy->X, &rdy->R));
  PetscCall(VecDuplicate(rdy->X, &rdy->Soln));
  PetscCall(VecViewFromOptions(rdy->X, NULL, "-vec_view"));
  PetscCall(DMCreateLocalVector(rdy->dm, &rdy->X_local));

  if (!rdy->ceed_resource[0]) {
    // water_src is only needed for PETSc source operator
    PetscCall(DMCreateGlobalVector(rdy->aux_dm, &rdy->water_src));
    PetscCall(VecZeroEntries(rdy->water_src));
  }

  PetscInt n_dof;
  PetscCall(VecGetSize(rdy->X, &n_dof));

  // set up a TS solver
  PetscCall(TSCreate(rdy->comm, &rdy->ts));
  PetscCall(TSSetProblemType(rdy->ts, TS_NONLINEAR));
  switch (rdy->config.numerics.temporal) {
    case TEMPORAL_EULER:
      PetscCall(TSSetType(rdy->ts, TSEULER));
      break;
    case TEMPORAL_RK4:
      PetscCall(TSSetType(rdy->ts, TSRK));
      PetscCall(TSRKSetType(rdy->ts, TSRK4));
      break;
    case TEMPORAL_BEULER:
      PetscCall(TSSetType(rdy->ts, TSBEULER));
      break;
  }
  PetscCall(TSSetDM(rdy->ts, rdy->dm));

  PetscCheck(rdy->config.physics.flow.mode == FLOW_SWE, rdy->comm, PETSC_ERR_USER, "Only the 'swe' flow mode is currently supported.");
  PetscCall(InitSWE(rdy));  // initialize SWE physics
  PetscCall(TSSetRHSFunction(rdy->ts, rdy->R, RHSFunctionSWE, rdy));

  PetscCall(TSSetMaxSteps(rdy->ts, rdy->config.time.max_step));
  PetscCall(TSSetExactFinalTime(rdy->ts, TS_EXACTFINALTIME_MATCHSTEP));
  PetscCall(TSSetSolution(rdy->ts, rdy->X));
  PetscCall(TSSetTime(rdy->ts, 0.0));
  PetscCall(TSSetTimeStep(rdy->ts, rdy->dt));

  // apply any solver-related options supplied on the command line
  PetscCall(TSSetFromOptions(rdy->ts));
  PetscCall(TSGetTimeStep(rdy->ts, &rdy->dt));  // just in case!

  PetscFunctionReturn(PETSC_SUCCESS);
}

// initializes solution vector data
//   unsafe for refinement if file is given with initial conditions
static PetscErrorCode InitSolution(RDy rdy) {
  PetscFunctionBegin;

  PetscCall(VecZeroEntries(rdy->X));

  // now initialize or override initial conditions by looping over regions
  // and writing values for corresponding cells
  PetscInt n_local;
  PetscCall(VecGetLocalSize(rdy->X, &n_local));
  PetscScalar *x_ptr;
  PetscCall(VecGetArray(rdy->X, &x_ptr));
  for (PetscInt r = 0; r < rdy->num_regions; ++r) {
    RDyRegion    region = rdy->regions[r];
    RDyCondition ic     = rdy->initial_conditions[r];
    PetscCheck(ic.flow, rdy->comm, PETSC_ERR_USER, "No initial condition specified for region '%s'", region.name);
    if (strlen(ic.flow->file)) {  // read regional data from file
      // FIXME: Figure this out!
      PetscViewer viewer;
      PetscCall(PetscViewerBinaryOpen(rdy->comm, ic.flow->file, FILE_MODE_READ, &viewer));
      Vec natural;
      PetscCall(DMPlexCreateNaturalVector(rdy->dm, &natural));
      PetscCall(VecLoad(natural, viewer));
      PetscCall(DMPlexNaturalToGlobalBegin(rdy->dm, natural, rdy->X));
      PetscCall(DMPlexNaturalToGlobalEnd(rdy->dm, natural, rdy->X));
      PetscCall(PetscViewerDestroy(&viewer));
      PetscCall(VecDestroy(&natural));
    } else {  // set components to specified values
      for (PetscInt c = 0; c < region.num_cells; ++c) {
        PetscInt cell_id = region.cell_ids[c];
        if (3 * cell_id < n_local) {  // skip ghost cells
          x_ptr[3 * cell_id]     = ic.flow->height;
          x_ptr[3 * cell_id + 1] = ic.flow->momentum[0];
          x_ptr[3 * cell_id + 2] = ic.flow->momentum[1];
        }
      }
    }
    // TODO: salinity and sediment initial conditions go here.
  }
  PetscCall(VecRestoreArray(rdy->X, &x_ptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}

// initialize the data on the right hand side of the boundary edges
static PetscErrorCode SetInitialBoundaryConditions(RDy rdy) {
  PetscFunctionBegin;

  for (PetscInt b = 0; b < rdy->num_boundaries; ++b) {
    RDyBoundary       boundary      = rdy->boundaries[b];
    RDyCondition      boundary_cond = rdy->boundary_conditions[b];
    RDyFlowCondition *flow_bc       = boundary_cond.flow;

    PetscReal boundary_values[3 * boundary.num_edges];
    switch (boundary_cond.flow->type) {
      case CONDITION_DIRICHLET:

        // initialize the relevant boundary values
        for (PetscInt e = 0; e < boundary.num_edges; ++e) {
          boundary_values[3 * e]     = flow_bc->height;
          boundary_values[3 * e + 1] = flow_bc->momentum[0];
          boundary_values[3 * e + 2] = flow_bc->momentum[1];
        }
        PetscCall(RDySetDirichletBoundaryValues(rdy, boundary.index, boundary.num_edges, 3, boundary_values));
        break;
      case CONDITION_REFLECTING:
        break;
      case CONDITION_CRITICAL_OUTFLOW:
        break;
      default:
        PetscCheck(PETSC_FALSE, PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid boundary condition encountered for boundary %" PetscInt_FMT "\n",
                   boundary.id);
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// the name of the log file set by RDySetLogFile below
static char overridden_logfile_[PETSC_MAX_PATH_LEN] = {0};

/// Sets the name of the log file for RDycore. If this function is called
/// before RDySetup, the name passed to it overrides any log filename set in
/// the YAML config file.
/// @param log_file [in] the name of the log file written by RDycore
PetscErrorCode RDySetLogFile(RDy rdy, const char *filename) {
  PetscFunctionBegin;
  strncpy(overridden_logfile_, filename, PETSC_MAX_PATH_LEN - 1);
  overridden_logfile_[PETSC_MAX_PATH_LEN - 1] = '\0';
  PetscFunctionReturn(PETSC_SUCCESS);
}

/// Performs any setup needed by RDy, reading from the specified configuration
/// file.
PetscErrorCode RDySetup(RDy rdy) {
  PetscFunctionBegin;

  // note: default config values are specified in the YAML input schema!
  PetscCall(ReadConfigFile(rdy));

  // override the log file name if necessary
  if (overridden_logfile_[0]) {
    strcpy(rdy->config.logging.file, overridden_logfile_);
  }

  // open the primary log file
  if (strlen(rdy->config.logging.file)) {
    PetscCall(PetscFOpen(rdy->comm, rdy->config.logging.file, "w", &rdy->log));
  } else {
    rdy->log = stdout;
  }

  // override parameters using command line arguments
  PetscCall(OverrideParameters(rdy));

  // print configuration info
  PetscCall(PrintConfig(rdy));

  // initialize CEED if needed
  if (rdy->ceed_resource[0]) {
    CeedInit(rdy->ceed_resource, &rdy->ceed);
  }

  RDyLogDebug(rdy, "Creating DMs...");
  PetscCall(CreateDM(rdy));           // for mesh and solution vector
  PetscCall(CreateAuxiliaryDM(rdy));  // for diagnostics

  RDyLogDebug(rdy, "Initializing regions and boundaries...");
  PetscCall(InitRegions(rdy));
  PetscCall(InitBoundaries(rdy));

  RDyLogDebug(rdy, "Initializing initial/boundary conditions and sources...");
  PetscCall(InitInitialConditions(rdy));
  PetscCall(InitSources(rdy));
  PetscCall(InitBoundaryConditions(rdy));

  RDyLogDebug(rdy, "Creating solvers and vectors...");
  PetscCall(CreateSolvers(rdy));

  RDyLogDebug(rdy, "Creating FV mesh...");
  // note: this must be done after global vectors are created so a global
  // note: section exists for the DM
  PetscCall(RDyMeshCreateFromDM(rdy->dm, &rdy->mesh));

  RDyLogDebug(rdy, "Initializing materials...");
  PetscCall(InitMaterials(rdy));

  RDyLogDebug(rdy, "Initializing solution data...");
  PetscCall(InitSolution(rdy));

  if (rdy->ceed_resource[0]) {
    RDyLogDebug(rdy, "Setting up CEED Operators...");

    // create the operators themselves
    PetscCall(CreateSWEFluxOperator(rdy->ceed, &rdy->mesh, rdy->num_boundaries, rdy->boundaries, rdy->boundary_conditions,
                                    rdy->config.physics.flow.tiny_h, &rdy->ceed_rhs.op_edges));

    PetscCall(CreateSWESourceOperator(rdy->ceed, &rdy->mesh, rdy->mesh.num_cells, rdy->materials_by_cell, rdy->config.physics.flow.tiny_h,
                                      &rdy->ceed_rhs.op_src));

    // create associated vectors for storage
    int num_comp = 3;
    CeedVectorCreate(rdy->ceed, rdy->mesh.num_cells * num_comp, &rdy->ceed_rhs.u_local_ceed);
    CeedVectorCreate(rdy->ceed, rdy->mesh.num_cells * num_comp, &rdy->ceed_rhs.f_ceed);
    CeedVectorCreate(rdy->ceed, rdy->mesh.num_cells_local * num_comp, &rdy->ceed_rhs.s_ceed);
    CeedVectorCreate(rdy->ceed, rdy->mesh.num_cells_local * num_comp, &rdy->ceed_rhs.u_ceed);

    // reset the time step size
    rdy->ceed_rhs.dt = 0.0;
  } else {
    // allocate storage for our PETSc implementation of the  flux and
    // source terms
    RDyLogDebug(rdy, "Allocating PETSc data structures for fluxes and sources...");
    PetscCall(CreatePetscSWEFlux(rdy->mesh.num_internal_edges, rdy->num_boundaries, rdy->boundaries, &rdy->petsc_rhs));
    PetscCall(CreatePetscSWESource(&rdy->mesh, rdy->petsc_rhs));
  }

  SetInitialBoundaryConditions(rdy);
  PetscFunctionReturn(PETSC_SUCCESS);
}
