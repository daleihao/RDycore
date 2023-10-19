#include <private/rdycoreimpl.h>
#include <rdycore.h>

PetscErrorCode RDyGetNumLocalCells(RDy rdy, PetscInt *num_cells) {
  PetscFunctionBegin;
  *num_cells = rdy->mesh.num_cells_local;
  PetscFunctionReturn(0);
}

PetscErrorCode RDyGetNumBoundaryConditions(RDy rdy, PetscInt *num_bnd_conds) {
  PetscFunctionBegin;
  *num_bnd_conds = rdy->num_boundaries;
  PetscFunctionReturn(0);
}

PetscErrorCode RDyGetNumEdgesInABoundaryConditions(RDy rdy, PetscInt bnd_cond_id, PetscInt *num_edges) {
  PetscFunctionBegin;
  PetscCheck(bnd_cond_id < rdy->num_boundaries, rdy->comm, PETSC_ERR_USER,
             "Boundary condition ID (%d) exceeds the max number of boundary conditions (%d)", bnd_cond_id, rdy->num_boundaries);
  PetscCheck(bnd_cond_id >= 0, rdy->comm, PETSC_ERR_USER, "Boundary condition ID (%d) cannot be less than zero.", bnd_cond_id);
  RDyBoundary *boundary = &rdy->boundaries[bnd_cond_id];
  *num_edges            = boundary->num_edges;
  PetscFunctionReturn(0);
}

PetscErrorCode RDyGetHeight(RDy rdy, PetscReal *h) {
  PetscFunctionBegin;

  PetscReal *x;
  PetscCall(VecGetArray(rdy->X, &x));
  for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
    h[i] = x[3 * i];
  }
  PetscCall(VecRestoreArray(rdy->X, &x));
  PetscFunctionReturn(0);
}

PetscErrorCode RDyGetXVelocity(RDy rdy, PetscReal *vx) {
  PetscFunctionBegin;

  PetscReal *x;
  PetscCall(VecGetArray(rdy->X, &x));
  for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
    vx[i] = x[3 * i + 1];
  }
  PetscCall(VecRestoreArray(rdy->X, &x));
  PetscFunctionReturn(0);
}

PetscErrorCode RDyGetYVelocity(RDy rdy, PetscReal *vy) {
  PetscFunctionBegin;

  PetscReal *x;
  PetscCall(VecGetArray(rdy->X, &x));
  for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
    vy[i] = x[3 * i + 2];
  }
  PetscCall(VecRestoreArray(rdy->X, &x));
  PetscFunctionReturn(0);
}

PetscErrorCode RDySetWaterSource(RDy rdy, PetscReal *watsrc) {
  PetscFunctionBegin;

  PetscReal *s;
  PetscCall(VecGetArray(rdy->water_src, &s));
  for (PetscInt i = 0; i < rdy->mesh.num_cells_local; ++i) {
    s[i] = watsrc[i];
  }
  PetscCall(VecRestoreArray(rdy->water_src, &s));
  rdy->ceed_rhs.water_src_updated = PETSC_FALSE;
  PetscFunctionReturn(0);
}
