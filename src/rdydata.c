#include <private/rdycoreimpl.h>
#include <private/rdysweimpl.h>
#include <rdycore.h>

PetscErrorCode RDyGetNumLocalCells(RDy rdy, PetscInt *num_cells) {
  PetscFunctionBegin;
  *num_cells = rdy->mesh.num_cells_local;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetNumBoundaryConditions(RDy rdy, PetscInt *num_bnd_conds) {
  PetscFunctionBegin;
  *num_bnd_conds = rdy->num_boundaries;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CheckBoundaryConditionIndex(RDy rdy, PetscInt boundary_index) {
  PetscFunctionBegin;
  PetscCheck(boundary_index < rdy->num_boundaries, rdy->comm, PETSC_ERR_USER,
             "Boundary condition index (%" PetscInt_FMT ") exceeds the max number of boundary conditions (%" PetscInt_FMT ")", boundary_index,
             rdy->num_boundaries);
  PetscCheck(boundary_index >= 0, rdy->comm, PETSC_ERR_USER, "Boundary condition index (%" PetscInt_FMT ") cannot be less than zero.",
             boundary_index);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CheckBoundaryNumEdges(RDy rdy, PetscInt boundary_index, PetscInt num_edges) {
  PetscFunctionBegin;

  RDyBoundary boundary = rdy->boundaries[boundary_index];

  PetscCheck(boundary.num_edges == num_edges, rdy->comm, PETSC_ERR_USER,
             "The given number of edges (%" PetscInt_FMT ") for boundary with index %" PetscInt_FMT " is incorrect (should be %" PetscInt_FMT ")",
             num_edges, boundary_index, boundary.num_edges);

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetNumBoundaryEdges(RDy rdy, PetscInt boundary_index, PetscInt *num_edges) {
  PetscFunctionBegin;
  PetscCall(CheckBoundaryConditionIndex(rdy, boundary_index));
  RDyBoundary *boundary = &rdy->boundaries[boundary_index];
  *num_edges            = boundary->num_edges;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetBoundaryConditionFlowType(RDy rdy, PetscInt boundary_index, PetscInt *bc_type) {
  PetscFunctionBegin;
  PetscCall(CheckBoundaryConditionIndex(rdy, boundary_index));
  RDyCondition *boundary_cond = &rdy->boundary_conditions[boundary_index];
  *bc_type                    = boundary_cond->flow->type;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDySetDirichletBoundaryValues(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscInt ndof, PetscReal *boundary_values) {
  PetscFunctionBegin;

  PetscCall(CheckBoundaryConditionIndex(rdy, boundary_index));

  PetscCheck(ndof == 3, rdy->comm, PETSC_ERR_USER, "The number of DOFs (%" PetscInt_FMT ") for the boundary condition need to be three.", ndof);

  RDyBoundary boundary = rdy->boundaries[boundary_index];
  PetscCheck(boundary.num_edges == num_edges, rdy->comm, PETSC_ERR_USER,
             "The given number of edges (%" PetscInt_FMT ") for boundary with index %" PetscInt_FMT " is incorrect (should be %" PetscInt_FMT ")",
             num_edges, boundary_index, boundary.num_edges);

  RDyCondition boundary_cond = rdy->boundary_conditions[boundary_index];
  PetscCheck(boundary_cond.flow->type == CONDITION_DIRICHLET, rdy->comm, PETSC_ERR_USER,
             "Trying to set dirichlet values for boundary with index %" PetscInt_FMT ", but it has a different type (%d)", boundary_index,
             boundary_cond.flow->type);

  // dispatch this call to CEED or PETSc
  PetscReal tiny_h = rdy->config.physics.flow.tiny_h;
  if (rdy->ceed_resource[0]) {  // ceed
    PetscInt size = 3 * rdy->boundaries[boundary_index].num_edges;
    PetscCall(SWEFluxOperatorSetDirichletBoundaryValues(rdy->ceed_rhs.op_edges, &rdy->mesh, rdy->boundaries[boundary_index], size, boundary_values));
  } else {  // petsc
    // fetch the boundary data
    RiemannDataSWE bdata;
    PetscCall(GetPetscSWEDirichletBoundaryValues(rdy->petsc_rhs, boundary_index, &bdata));

    // set the boundary values
    RDyCells *cells = &rdy->mesh.cells;
    RDyEdges *edges = &rdy->mesh.edges;
    for (PetscInt e = 0; e < boundary.num_edges; ++e) {
      PetscInt iedge = boundary.edge_ids[e];
      PetscInt icell = edges->cell_ids[2 * iedge];
      if (cells->is_local[icell]) {
        bdata.h[e]  = boundary_values[3 * e];
        bdata.hu[e] = boundary_values[3 * e + 1];
        bdata.hv[e] = boundary_values[3 * e + 2];

        if (bdata.h[e] > tiny_h) {
          bdata.u[e] = boundary_values[3 * e + 1] / bdata.h[e];
          bdata.v[e] = boundary_values[3 * e + 2] / bdata.h[e];
        } else {
          bdata.u[e] = 0.0;
          bdata.v[e] = 0.0;
        }
      }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetHeight(RDy rdy, PetscReal *h) {
  PetscFunctionBegin;

  PetscReal *x;
  PetscCall(VecGetArray(rdy->X, &x));
  for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
    h[i] = x[3 * i];
  }
  PetscCall(VecRestoreArray(rdy->X, &x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetXVelocity(RDy rdy, PetscReal *vx) {
  PetscFunctionBegin;

  PetscReal *x;
  PetscCall(VecGetArray(rdy->X, &x));
  for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
    vx[i] = x[3 * i + 1];
  }
  PetscCall(VecRestoreArray(rdy->X, &x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetYVelocity(RDy rdy, PetscReal *vy) {
  PetscFunctionBegin;

  PetscReal *x;
  PetscCall(VecGetArray(rdy->X, &x));
  for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
    vy[i] = x[3 * i + 2];
  }
  PetscCall(VecRestoreArray(rdy->X, &x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDySetWaterSource(RDy rdy, PetscReal *watsrc) {
  PetscFunctionBegin;

  if (rdy->ceed_resource[0]) {  // ceed
    PetscCall(SWESourceOperatorSetWaterSource(rdy->ceed_rhs.op_src, watsrc));
  } else {  // petsc
    PetscReal *s;
    PetscCall(VecGetArray(rdy->water_src, &s));
    for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
      s[i] = watsrc[i];
    }
    PetscCall(VecRestoreArray(rdy->water_src, &s));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDySetMomentumrSource(RDy rdy, Vec momentum_vec, PetscReal *momentum_value) {
  PetscFunctionBegin;

  if (rdy->ceed_resource[0]) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Extend RDySetMomentumrSource for Ceed");
  } else {
    PetscReal *m;
    PetscCall(VecGetArray(momentum_vec, &m));
    for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
      m[i] = momentum_value[i];
    }
    PetscCall(VecRestoreArray(momentum_vec, &m));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDySetXMomentumSource(RDy rdy, PetscReal *x_momentum) {
  PetscFunctionBegin;
  PetscCall(RDySetMomentumrSource(rdy, rdy->x_momentum_src, x_momentum));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDySetYMomentumSource(RDy rdy, PetscReal *y_momentum) {
  PetscFunctionBegin;
  PetscCall(RDySetMomentumrSource(rdy, rdy->y_momentum_src, y_momentum));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RDyGetIDimCentroidOfLocalCell(RDy rdy, PetscInt idim, PetscReal *x) {
  PetscFunctionBegin;

  RDyCells *cells = &rdy->mesh.cells;

  PetscInt count = 0;
  for (PetscInt icell = 0; icell < rdy->mesh.num_cells; ++icell) {
    if (cells->is_local[icell]) {
      x[count++] = cells->centroids[icell].X[idim];
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetXCentroidOfLocalCell(RDy rdy, PetscReal *x) {
  PetscFunctionBegin;
  PetscInt idim = 0;  // x-dim
  PetscCall(RDyGetIDimCentroidOfLocalCell(rdy, idim, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetYCentroidOfLocalCell(RDy rdy, PetscReal *y) {
  PetscFunctionBegin;
  PetscInt idim = 1;  // y-dim
  PetscCall(RDyGetIDimCentroidOfLocalCell(rdy, idim, y));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetZCentroidOfLocalCell(RDy rdy, PetscReal *z) {
  PetscFunctionBegin;
  PetscInt idim = 2;  // z-dim
  PetscCall(RDyGetIDimCentroidOfLocalCell(rdy, idim, z));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetNatIDOfLocalCell(RDy rdy, PetscInt *nat_id) {
  PetscFunctionBegin;
  RDyCells *cells = &rdy->mesh.cells;

  PetscInt count = 0;
  for (PetscInt icell = 0; icell < rdy->mesh.num_cells; ++icell) {
    if (cells->is_local[icell]) {
      nat_id[count++] = cells->natural_ids[icell];
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RDyGetIDimCentroidOfBoundaryEdgeOrCell(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscBool data_for_edge,
                                                             PetscInt idim, PetscReal *x) {
  PetscFunctionBegin;

  PetscCall(CheckBoundaryConditionIndex(rdy, boundary_index));
  PetscCall(CheckBoundaryNumEdges(rdy, boundary_index, num_edges));

  RDyBoundary boundary = rdy->boundaries[boundary_index];
  RDyEdges   *edges    = &rdy->mesh.edges;

  if (data_for_edge) {
    for (PetscInt e = 0; e < boundary.num_edges; ++e) {
      PetscInt iedge = boundary.edge_ids[e];
      x[e]           = edges->centroids[iedge].X[idim];
    }
  } else {
    RDyCells *cells = &rdy->mesh.cells;
    for (PetscInt e = 0; e < boundary.num_edges; ++e) {
      PetscInt iedge = boundary.edge_ids[e];
      PetscInt icell = edges->cell_ids[2 * iedge];
      x[e]           = cells->centroids[icell].X[idim];
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetXCentroidOfBoundaryEdge(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscReal *x) {
  PetscFunctionBegin;
  PetscBool data_for_edge = PETSC_TRUE;
  PetscInt  idim          = 0;  // x-dim
  PetscCall(RDyGetIDimCentroidOfBoundaryEdgeOrCell(rdy, boundary_index, num_edges, data_for_edge, idim, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetYCentroidOfBoundaryEdge(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscReal *x) {
  PetscFunctionBegin;
  PetscBool data_for_edge = PETSC_TRUE;
  PetscInt  idim          = 1;  // y-dim
  PetscCall(RDyGetIDimCentroidOfBoundaryEdgeOrCell(rdy, boundary_index, num_edges, data_for_edge, idim, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetZCentroidOfBoundaryEdge(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscReal *x) {
  PetscFunctionBegin;
  PetscBool data_for_edge = PETSC_TRUE;
  PetscInt  idim          = 2;  // z-dim
  PetscCall(RDyGetIDimCentroidOfBoundaryEdgeOrCell(rdy, boundary_index, num_edges, data_for_edge, idim, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetXCentroidOfBoundaryCell(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscReal *x) {
  PetscFunctionBegin;
  PetscBool data_for_edge = PETSC_FALSE;
  PetscInt  idim          = 0;  // x-dim
  PetscCall(RDyGetIDimCentroidOfBoundaryEdgeOrCell(rdy, boundary_index, num_edges, data_for_edge, idim, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetYCentroidOfBoundaryCell(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscReal *x) {
  PetscFunctionBegin;
  PetscBool data_for_edge = PETSC_FALSE;
  PetscInt  idim          = 1;  // y-dim
  PetscCall(RDyGetIDimCentroidOfBoundaryEdgeOrCell(rdy, boundary_index, num_edges, data_for_edge, idim, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetZCentroidOfBoundaryCellCell(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscReal *x) {
  PetscFunctionBegin;
  PetscBool data_for_edge = PETSC_FALSE;
  PetscInt  idim          = 2;  // z-dim
  PetscCall(RDyGetIDimCentroidOfBoundaryEdgeOrCell(rdy, boundary_index, num_edges, data_for_edge, idim, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDyGetNaturalIDOfBoundaryCell(RDy rdy, PetscInt boundary_index, PetscInt num_edges, PetscInt *nat_ids) {
  PetscFunctionBegin;
  PetscCall(CheckBoundaryConditionIndex(rdy, boundary_index));
  PetscCall(CheckBoundaryNumEdges(rdy, boundary_index, num_edges));

  RDyBoundary boundary = rdy->boundaries[boundary_index];
  RDyCells   *cells    = &rdy->mesh.cells;
  RDyEdges   *edges    = &rdy->mesh.edges;

  for (PetscInt e = 0; e < boundary.num_edges; ++e) {
    PetscInt iedge = boundary.edge_ids[e];
    PetscInt icell = edges->cell_ids[2 * iedge];
    nat_ids[e]     = cells->natural_ids[icell];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RDySetManningN(RDy rdy, PetscReal *n_values) {
  PetscFunctionBegin;

  if (rdy->ceed_resource[0]) {  // ceed
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Extend RDySetManningN for Ceed");
  } else {  // petsc
    for (PetscInt icell = 0; icell < rdy->mesh.num_cells_local; ++icell) {
      rdy->materials_by_cell[icell].manning = n_values[icell];
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}
