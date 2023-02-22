#include <private/rdycoreimpl.h>
#include <private/rdymathimpl.h>

// gravitational acceleration [m/s/s]
static const PetscReal GRAVITY = 9.806;

// computes velocities in x and y-dir based on momentum in x and y-dir
// N - Size of the array
// tiny_h - Threshold value for height
// h - Height
// hu - Momentum in x-dir
// hv - Momentum in y-dir
// u - Velocity in x-dir
// v - Velocity in y-dir
static PetscErrorCode GetVelocityFromMomentum(PetscInt N, PetscReal tiny_h, const PetscReal h[N], const PetscReal hu[N], const PetscReal hv[N],
                                              PetscReal u[N], PetscReal v[N]) {
  PetscFunctionBeginUser;

  for (PetscInt n = 0; n < N; n++) {
    if (h[n] < tiny_h) {
      u[n] = 0.0;
      v[n] = 0.0;
    } else {
      u[n] = hu[n] / h[n];
      v[n] = hv[n] / h[n];
    }
  }

  PetscFunctionReturn(0);
}

// Computes flux based on Roe solver
// N    - Size of the array
// hl   - Height left of the edge
// hr   - Height right of the edge
// ul   - Velocity in x-dir left of the edge
// ur   - Velocity in x-dir right of the edge
// vl   - Velocity in y-dir left of the edge
// vr   - Velocity in y-dir right of the edge
// sn   - sine of the angle between edge and y-axis
// cn   - cosine of the angle between edge and y-axis
// fij  - flux
// amax - maximum wave speed
static PetscErrorCode ComputeRoeFlux(PetscInt N, const PetscReal hl[N], const PetscReal hr[N], const PetscReal ul[N], const PetscReal ur[N],
                                     const PetscReal vl[N], const PetscReal vr[N], const PetscReal sn[N], const PetscReal cn[N], PetscReal fij[N][3],
                                     PetscReal amax[N]) {
  PetscFunctionBeginUser;

  for (PetscInt n = 0; n < N; n++) {
    // Compute Roe averages
    PetscReal duml  = pow(hl[n], 0.5);
    PetscReal dumr  = pow(hr[n], 0.5);
    PetscReal cl    = pow(GRAVITY * hl[n], 0.5);
    PetscReal cr    = pow(GRAVITY * hr[n], 0.5);
    PetscReal hhat  = duml * dumr;
    PetscReal uhat  = (duml * ul[n] + dumr * ur[n]) / (duml + dumr);
    PetscReal vhat  = (duml * vl[n] + dumr * vr[n]) / (duml + dumr);
    PetscReal chat  = pow(0.5 * GRAVITY * (hl[n] + hr[n]), 0.5);
    PetscReal uperp = uhat * cn[n] + vhat * sn[n];

    PetscReal dh     = hr[n] - hl[n];
    PetscReal du     = ur[n] - ul[n];
    PetscReal dv     = vr[n] - vl[n];
    PetscReal dupar  = -du * sn[n] + dv * cn[n];
    PetscReal duperp = du * cn[n] + dv * sn[n];

    PetscReal dW[3];
    dW[0] = 0.5 * (dh - hhat * duperp / chat);
    dW[1] = hhat * dupar;
    dW[2] = 0.5 * (dh + hhat * duperp / chat);

    PetscReal uperpl = ul[n] * cn[n] + vl[n] * sn[n];
    PetscReal uperpr = ur[n] * cn[n] + vr[n] * sn[n];
    PetscReal al1    = uperpl - cl;
    PetscReal al3    = uperpl + cl;
    PetscReal ar1    = uperpr - cr;
    PetscReal ar3    = uperpr + cr;

    PetscReal R[3][3];
    R[0][0] = 1.0;
    R[0][1] = 0.0;
    R[0][2] = 1.0;
    R[1][0] = uhat - chat * cn[n];
    R[1][1] = -sn[n];
    R[1][2] = uhat + chat * cn[n];
    R[2][0] = vhat - chat * sn[n];
    R[2][1] = cn[n];
    R[2][2] = vhat + chat * sn[n];

    PetscReal da1 = fmax(0.0, 2.0 * (ar1 - al1));
    PetscReal da3 = fmax(0.0, 2.0 * (ar3 - al3));
    PetscReal a1  = fabs(uperp - chat);
    PetscReal a2  = fabs(uperp);
    PetscReal a3  = fabs(uperp + chat);

    // Critical flow fix
    if (a1 < da1) {
      a1 = 0.5 * (a1 * a1 / da1 + da1);
    }
    if (a3 < da3) {
      a3 = 0.5 * (a3 * a3 / da3 + da3);
    }

    // Compute interface flux
    PetscReal A[3][3];
    for (PetscInt i = 0; i < 3; i++) {
      for (PetscInt j = 0; j < 3; j++) {
        A[i][j] = 0.0;
      }
    }
    A[0][0] = a1;
    A[1][1] = a2;
    A[2][2] = a3;

    PetscReal FL[3], FR[3];
    FL[0] = uperpl * hl[n];
    FL[1] = ul[n] * uperpl * hl[n] + 0.5 * GRAVITY * hl[n] * hl[n] * cn[n];
    FL[2] = vl[n] * uperpl * hl[n] + 0.5 * GRAVITY * hl[n] * hl[n] * sn[n];

    FR[0] = uperpr * hr[n];
    FR[1] = ur[n] * uperpr * hr[n] + 0.5 * GRAVITY * hr[n] * hr[n] * cn[n];
    FR[2] = vr[n] * uperpr * hr[n] + 0.5 * GRAVITY * hr[n] * hr[n] * sn[n];

    // fij = 0.5*(FL + FR - matmul(R,matmul(A,dW))
    fij[n][0] = 0.5 * (FL[0] + FR[0] - R[0][0] * A[0][0] * dW[0] - R[0][1] * A[1][1] * dW[1] - R[0][2] * A[2][2] * dW[2]);
    fij[n][1] = 0.5 * (FL[1] + FR[1] - R[1][0] * A[0][0] * dW[0] - R[1][1] * A[1][1] * dW[1] - R[1][2] * A[2][2] * dW[2]);
    fij[n][2] = 0.5 * (FL[2] + FR[2] - R[2][0] * A[0][0] * dW[0] - R[2][1] * A[1][1] * dW[1] - R[2][2] * A[2][2] * dW[2]);

    amax[n] = chat + fabs(uperp);
  }

  PetscFunctionReturn(0);
}

// computes RHS on internal edges
static PetscErrorCode RHSFunctionForInternalEdges(RDy rdy, Vec F, PetscReal *amax_value) {
  PetscFunctionBeginUser;

  RDyMesh  *mesh  = &rdy->mesh;
  RDyCells *cells = &mesh->cells;
  RDyEdges *edges = &mesh->edges;

  // Get pointers to vector data
  PetscScalar *x_ptr, *f_ptr, *b_ptr;
  PetscCall(VecGetArray(rdy->X_local, &x_ptr));
  PetscCall(VecGetArray(F, &f_ptr));
  PetscCall(VecGetArray(rdy->B_local, &b_ptr));

  const PetscInt ndof = 3;
  PetscInt       num  = mesh->num_internal_edges;
  PetscReal      hl_vec_int[num], hul_vec_int[num], hvl_vec_int[num], ul_vec_int[num], vl_vec_int[num];
  PetscReal      hr_vec_int[num], hur_vec_int[num], hvr_vec_int[num], ur_vec_int[num], vr_vec_int[num];
  PetscReal      sn_vec_int[num], cn_vec_int[num];
  PetscReal      flux_vec_int[num][3], amax_vec_int[num];

  // Collect the h/hu/hv for left and right cells to compute u/v
  for (PetscInt ii = 0; ii < mesh->num_internal_edges; ii++) {
    PetscInt iedge      = edges->internal_edge_ids[ii];
    PetscInt cellOffset = edges->cell_offsets[iedge];
    PetscInt l          = edges->cell_ids[cellOffset];
    PetscInt r          = edges->cell_ids[cellOffset + 1];

    hl_vec_int[ii]  = x_ptr[l * ndof + 0];
    hul_vec_int[ii] = x_ptr[l * ndof + 1];
    hvl_vec_int[ii] = x_ptr[l * ndof + 2];

    hr_vec_int[ii]  = x_ptr[r * ndof + 0];
    hur_vec_int[ii] = x_ptr[r * ndof + 1];
    hvr_vec_int[ii] = x_ptr[r * ndof + 2];
  }

  // Compute u/v for left and right cells
  const PetscReal tiny_h = rdy->config.tiny_h;
  PetscCall(GetVelocityFromMomentum(num, tiny_h, hl_vec_int, hul_vec_int, hvl_vec_int, ul_vec_int, vl_vec_int));
  PetscCall(GetVelocityFromMomentum(num, tiny_h, hr_vec_int, hur_vec_int, hvr_vec_int, ur_vec_int, vr_vec_int));

  // Update u/v for reflective internal edges
  for (PetscInt ii = 0; ii < mesh->num_internal_edges; ii++) {
    PetscInt  iedge      = edges->internal_edge_ids[ii];
    PetscInt  cellOffset = edges->cell_offsets[iedge];
    PetscInt  l          = edges->cell_ids[cellOffset];
    PetscInt  r          = edges->cell_ids[cellOffset + 1];
    PetscReal bl         = b_ptr[l];
    PetscReal br         = b_ptr[r];

    cn_vec_int[ii] = edges->cn[iedge];
    sn_vec_int[ii] = edges->sn[iedge];

    if (bl == 1 && br == 0) {
      // Update left values as it is a reflective boundary wall
      hl_vec_int[ii] = hr_vec_int[ii];

      PetscReal dum1 = Square(sn_vec_int[ii]) - Square(cn_vec_int[ii]);
      PetscReal dum2 = 2.0 * sn_vec_int[ii] * cn_vec_int[ii];

      ul_vec_int[ii] = ur_vec_int[ii] * dum1 - vr_vec_int[ii] * dum2;
      vl_vec_int[ii] = -ur_vec_int[ii] * dum2 - vr_vec_int[ii] * dum1;

    } else if (bl == 0 && br == 1) {
      // Update right values as it is a reflective boundary wall
      hr_vec_int[ii] = hl_vec_int[ii];

      PetscReal dum1 = Square(sn_vec_int[ii]) - Square(cn_vec_int[ii]);
      PetscReal dum2 = 2.0 * sn_vec_int[ii] * cn_vec_int[ii];

      ur_vec_int[ii] = ul_vec_int[ii] * dum1 - vl_vec_int[ii] * dum2;
      vr_vec_int[ii] = -ul_vec_int[ii] * dum2 - vl_vec_int[ii] * dum1;
    }
  }

  // Call Riemann solver (only Roe currently supported)
  PetscCheck(rdy->config.riemann == RIEMANN_ROE, rdy->comm, PETSC_ERR_USER, "Invalid Riemann solver selected! (Only roe is supported)");
  PetscCall(ComputeRoeFlux(num, hl_vec_int, hr_vec_int, ul_vec_int, ur_vec_int, vl_vec_int, vr_vec_int, sn_vec_int, cn_vec_int, flux_vec_int,
                           amax_vec_int));

  // Save the flux values in the Vec based by TS
  for (PetscInt ii = 0; ii < mesh->num_internal_edges; ii++) {
    PetscInt  iedge      = edges->internal_edge_ids[ii];
    PetscInt  cellOffset = edges->cell_offsets[iedge];
    PetscInt  l          = edges->cell_ids[cellOffset];
    PetscInt  r          = edges->cell_ids[cellOffset + 1];
    PetscReal edgeLen    = edges->lengths[iedge];

    PetscReal hl = x_ptr[l * ndof + 0];
    PetscReal hr = x_ptr[r * ndof + 0];
    PetscReal bl = b_ptr[l];
    PetscReal br = b_ptr[r];

    *amax_value = fmax(*amax_value, amax_vec_int[ii]);

    if (bl == 0 && br == 0) {
      // Both, left and right cells are not reflective boundary walls
      if (!(hr < tiny_h && hl < tiny_h)) {
        PetscReal areal = cells->areas[l];
        PetscReal arear = cells->areas[r];

        for (PetscInt idof = 0; idof < ndof; idof++) {
          if (cells->is_local[l]) f_ptr[l * ndof + idof] -= flux_vec_int[ii][idof] * edgeLen / areal;
          if (cells->is_local[r]) f_ptr[r * ndof + idof] += flux_vec_int[ii][idof] * edgeLen / arear;
        }
      }

    } else if (bl == 1 && br == 0) {
      // Left cell is a reflective boundary wall and right cell is an internal cell

      PetscReal arear = cells->areas[r];
      for (PetscInt idof = 0; idof < ndof; idof++) {
        if (cells->is_local[r]) f_ptr[r * ndof + idof] += flux_vec_int[ii][idof] * edgeLen / arear;
      }

    } else if (bl == 0 && br == 1) {
      // Left cell is an internal cell and right cell is a reflective boundary wall

      PetscReal areal = cells->areas[l];
      for (PetscInt idof = 0; idof < ndof; idof++) {
        if (cells->is_local[l]) f_ptr[l * ndof + idof] -= flux_vec_int[ii][idof] * edgeLen / areal;
      }
    }
  }

  // Restore vectors
  PetscCall(VecRestoreArray(rdy->X_local, &x_ptr));
  PetscCall(VecRestoreArray(F, &f_ptr));
  PetscCall(VecRestoreArray(rdy->B_local, &b_ptr));

  PetscFunctionReturn(0);
}

// computes RHS on boundary edges
static PetscErrorCode RHSFunctionForBoundaryEdges(RDy rdy, Vec F, PetscReal *amax_value) {
  PetscFunctionBeginUser;

  RDyMesh  *mesh  = &rdy->mesh;
  RDyCells *cells = &mesh->cells;
  RDyEdges *edges = &mesh->edges;

  // Get pointers to vector data
  PetscScalar *x_ptr, *f_ptr, *b_ptr;
  PetscCall(VecGetArray(rdy->X_local, &x_ptr));
  PetscCall(VecGetArray(F, &f_ptr));
  PetscCall(VecGetArray(rdy->B_local, &b_ptr));

  const PetscInt ndof = 3;
  PetscInt       num  = mesh->num_boundary_edges;
  PetscReal      hl_vec_bnd[num], hul_vec_bnd[num], hvl_vec_bnd[num], ul_vec_bnd[num], vl_vec_bnd[num];
  PetscReal      hr_vec_bnd[num], ur_vec_bnd[num], vr_vec_bnd[num];
  PetscReal      sn_vec_bnd[num], cn_vec_bnd[num];
  PetscReal      flux_vec_bnd[num][3], amax_vec_bnd[num];

  const PetscReal tiny_h = rdy->config.tiny_h;

  // Collect the h/hu/hv for left cells to compute u/v
  for (PetscInt ii = 0; ii < mesh->num_boundary_edges; ii++) {
    PetscInt iedge      = edges->boundary_edge_ids[ii];
    PetscInt cellOffset = edges->cell_offsets[iedge];
    PetscInt l          = edges->cell_ids[cellOffset];

    hl_vec_bnd[ii]  = x_ptr[l * ndof + 0];
    hul_vec_bnd[ii] = x_ptr[l * ndof + 1];
    hvl_vec_bnd[ii] = x_ptr[l * ndof + 2];
  }

  // Compute u/v for left cells
  PetscCall(GetVelocityFromMomentum(num, tiny_h, hl_vec_bnd, hul_vec_bnd, hvl_vec_bnd, ul_vec_bnd, vl_vec_bnd));

  // Compute h/u/v for right cells
  for (PetscInt ii = 0; ii < mesh->num_boundary_edges; ii++) {
    PetscInt iedge      = edges->boundary_edge_ids[ii];
    PetscInt cellOffset = edges->cell_offsets[iedge];
    PetscInt l          = edges->cell_ids[cellOffset];

    cn_vec_bnd[ii] = edges->cn[iedge];
    sn_vec_bnd[ii] = edges->sn[iedge];

    if (cells->is_local[l] && b_ptr[l] == 0) {
      // Perform computation for a boundary edge

      if (cells->is_local[l] && b_ptr[l] == 0) {
        hr_vec_bnd[ii] = hl_vec_bnd[ii];

        PetscReal dum1 = Square(sn_vec_bnd[ii]) - Square(cn_vec_bnd[ii]);
        PetscReal dum2 = 2.0 * sn_vec_bnd[ii] * cn_vec_bnd[ii];

        ur_vec_bnd[ii] = ul_vec_bnd[ii] * dum1 - vl_vec_bnd[ii] * dum2;
        vr_vec_bnd[ii] = -ul_vec_bnd[ii] * dum2 - vl_vec_bnd[ii] * dum1;
      }
    }
  }

  // Call Riemann solver (only Roe is currently supported)
  PetscCheck(rdy->config.riemann == RIEMANN_ROE, rdy->comm, PETSC_ERR_USER, "Invalid Riemann solver selected! (Only roe is supported)");
  PetscCall(ComputeRoeFlux(num, hl_vec_bnd, hr_vec_bnd, ul_vec_bnd, ur_vec_bnd, vl_vec_bnd, vr_vec_bnd, sn_vec_bnd, cn_vec_bnd, flux_vec_bnd,
                           amax_vec_bnd));

  // Save the flux values in the Vec based by TS
  for (PetscInt ii = 0; ii < mesh->num_boundary_edges; ii++) {
    PetscInt  iedge      = edges->boundary_edge_ids[ii];
    PetscInt  cellOffset = edges->cell_offsets[iedge];
    PetscInt  l          = edges->cell_ids[cellOffset];
    PetscReal edgeLen    = edges->lengths[iedge];
    PetscReal areal      = cells->areas[l];

    if (cells->is_local[l] && b_ptr[l] == 0) {
      // Perform computation for a boundary edge

      PetscReal hl = x_ptr[l * ndof + 0];

      if (!(hl < tiny_h)) {
        *amax_value = fmax(*amax_value, amax_vec_bnd[ii]);
        for (PetscInt idof = 0; idof < ndof; idof++) {
          f_ptr[l * ndof + idof] -= flux_vec_bnd[ii][idof] * edgeLen / areal;
        }
      }
    }
  }

  // Restore vectors
  PetscCall(VecRestoreArray(rdy->X_local, &x_ptr));
  PetscCall(VecRestoreArray(F, &f_ptr));
  PetscCall(VecRestoreArray(rdy->B_local, &b_ptr));

  PetscFunctionReturn(0);
}

// adds source terms to the right hand side vector F
static PetscErrorCode AddSourceTerm(RDy rdy, Vec F) {
  PetscFunctionBeginUser;

  RDyMesh  *mesh  = &rdy->mesh;
  RDyCells *cells = &mesh->cells;

  // Get access to Vec
  PetscScalar *x_ptr, *f_ptr;
  PetscCall(VecGetArray(rdy->X_local, &x_ptr));
  PetscCall(VecGetArray(F, &f_ptr));

  const PetscInt ndof = 3;

  PetscInt  N = mesh->num_cells;
  PetscReal h_vec[N], hu_vec[N], hv_vec[N], u_vec[N], v_vec[N];

  // Collect the h/hu/hv for cells to compute u/v
  for (PetscInt icell = 0; icell < mesh->num_cells; icell++) {
    h_vec[icell]  = x_ptr[icell * ndof + 0];
    hu_vec[icell] = x_ptr[icell * ndof + 1];
    hv_vec[icell] = x_ptr[icell * ndof + 2];
  }

  // Compute u/v for cells
  PetscCall(GetVelocityFromMomentum(N, rdy->config.tiny_h, h_vec, hu_vec, hv_vec, u_vec, v_vec));

  for (PetscInt icell = 0; icell < mesh->num_cells; icell++) {
    if (cells->is_local[icell]) {
      PetscReal h  = h_vec[icell];
      PetscReal hu = hu_vec[icell];
      PetscReal hv = hv_vec[icell];

      PetscReal dz_dx = cells->dz_dx[icell];
      PetscReal dz_dy = cells->dz_dy[icell];

      PetscReal bedx = dz_dx * GRAVITY * h;
      PetscReal bedy = dz_dy * GRAVITY * h;

      PetscReal u = u_vec[icell];
      PetscReal v = u_vec[icell];

      PetscReal Fsum_x = f_ptr[icell * ndof + 1];
      PetscReal Fsum_y = f_ptr[icell * ndof + 2];

      PetscReal tbx = 0.0, tby = 0.0;

      if (h >= rdy->config.tiny_h) {
        // Manning's coefficient
        PetscReal Uniform_roughness = 0.015;
        PetscReal N_mannings        = GRAVITY * Uniform_roughness * Uniform_roughness;

        // Cd = g n^2 h^{-1/3}, where n is Manning's coefficient
        PetscReal Cd = GRAVITY * Square(N_mannings) * PetscPowReal(h, -1.0 / 3.0);

        PetscReal velocity = PetscSqrtReal(Square(u) + Square(v));

        PetscReal tb = Cd * velocity / h;

        PetscReal dt     = rdy->dt;
        PetscReal factor = tb / (1.0 + dt * tb);

        tbx = (hu + dt * Fsum_x - dt * bedx) * factor;
        tby = (hv + dt * Fsum_y - dt * bedy) * factor;
      }

      f_ptr[icell * ndof + 0] += 0.0;
      f_ptr[icell * ndof + 1] += -bedx - tbx;
      f_ptr[icell * ndof + 2] += -bedy - tby;
    }
  }

  // Restore vectors
  PetscCall(VecRestoreArray(rdy->X_local, &x_ptr));
  PetscCall(VecRestoreArray(F, &f_ptr));

  PetscFunctionReturn(0);
}

// This is the right-hand-side function used by our timestepping solver.
// Parameters:
//  ts  - the solver
//  t   - the simulation time [minutes]
//  X   - the solution vector at time t
//  F   - the right hand side vector to be evaluated at time t
//  ctx - a generic pointer to our RDy object
PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec X, Vec F, void *ctx) {
  PetscFunctionBegin;

  RDy rdy = ctx;
  DM  dm  = rdy->dm;

  PetscCall(VecZeroEntries(F));

  // populate the local X vector
  PetscCall(DMGlobalToLocalBegin(dm, X, INSERT_VALUES, rdy->X_local));
  PetscCall(DMGlobalToLocalEnd(dm, X, INSERT_VALUES, rdy->X_local));

  // compute the right hand side
  PetscReal amax_value = 0.0;
  PetscCall(RHSFunctionForInternalEdges(rdy, F, &amax_value));
  PetscCall(RHSFunctionForBoundaryEdges(rdy, F, &amax_value));
  PetscCall(AddSourceTerm(rdy, F));

  /* disable this for now
  if (rdy->save) {
    char fname[PETSC_MAX_PATH_LEN];
    sprintf(fname, "outputs/ex2b_Nx_%d_Ny_%d_dt_%f_%d_np%d.dat", rdy->Nx, rdy->Ny, rdy->dt, rdy->tstep - 1, rdy->comm_size);
    PetscViewer viewer;
    PetscCall(PetscViewerBinaryOpen(rdy->comm, fname, FILE_MODE_WRITE, &viewer));

    Vec natural;
    PetscCall(DMPlexCreateNaturalVector(rdy->dm, &natural));
    PetscCall(DMPlexGlobalToNaturalBegin(rdy->dm, X, natural));
    PetscCall(DMPlexGlobalToNaturalEnd(rdy->dm, X, natural));
    PetscCall(VecView(natural, viewer));
    PetscCall(PetscViewerDestroy(&viewer));

    sprintf(fname, "outputs/ex2b_flux_Nx_%d_Ny_%d_dt_%f_%d_np%d.dat", rdy->Nx, rdy->Ny, rdy->dt, rdy->tstep - 1, rdy->comm_size);
    PetscCall(PetscViewerBinaryOpen(rdy->comm, fname, FILE_MODE_WRITE, &viewer));
    PetscCall(DMPlexGlobalToNaturalBegin(rdy->dm, F, natural));
    PetscCall(DMPlexGlobalToNaturalEnd(rdy->dm, F, natural));
    PetscCall(VecView(natural, viewer));
    PetscCall(PetscViewerDestroy(&viewer));

    PetscCall(VecDestroy(&natural));
  }
  */

  RDyLogInfo(rdy, "Time step %d: t = %f, CFL = %f", rdy->step, t, amax_value * rdy->dt * 2);
  rdy->step++;

  PetscFunctionReturn(0);
}
