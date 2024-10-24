#include <petscdmceed.h>
#include <private/rdysweimpl.h>

#include "swe_ceed_impl.h"

// CEED uses C99 VLA features for shaping multidimensional
// arrays, which don't have the same drawbacks as VLA allocations.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wvla"

// frees a data context allocated using PETSc, returning a libCEED error code
static int FreeContextPetsc(void *data) {
  if (PetscFree(data)) return CeedError(NULL, CEED_ERROR_ACCESS, "PetscFree failed");
  return CEED_ERROR_SUCCESS;
}

// creates a QFunction context for a flux or source operator with the given
// minimum water height threshold
static PetscErrorCode CreateQFunctionContext(Ceed ceed, PetscReal tiny_h, CeedQFunctionContext *qf_context) {
  PetscFunctionBeginUser;

  SWEContext swe_ctx;
  PetscCall(PetscCalloc1(1, &swe_ctx));

  swe_ctx->dtime   = 0.0;
  swe_ctx->tiny_h  = tiny_h;
  swe_ctx->gravity = 9.806;

  PetscCallCEED(CeedQFunctionContextCreate(ceed, qf_context));
  PetscCallCEED(CeedQFunctionContextSetData(*qf_context, CEED_MEM_HOST, CEED_USE_POINTER, sizeof(*swe_ctx), swe_ctx));
  PetscCallCEED(CeedQFunctionContextSetDataDestroy(*qf_context, CEED_MEM_HOST, FreeContextPetsc));
  PetscCallCEED(CeedQFunctionContextRegisterDouble(*qf_context, "time step", offsetof(struct SWEContext_, dtime), 1, "Time step of TS"));
  PetscCallCEED(CeedQFunctionContextRegisterDouble(*qf_context, "small h value", offsetof(struct SWEContext_, tiny_h), 1,
                                                   "Height threshold below which dry condition is assumed"));
  PetscCallCEED(CeedQFunctionContextRegisterDouble(*qf_context, "gravity", offsetof(struct SWEContext_, gravity), 1, "Accelaration due to gravity"));

  PetscFunctionReturn(CEED_ERROR_SUCCESS);
}

static PetscErrorCode CreateSWEInteriorFluxSubOperator(RDyMesh *mesh, PetscReal tiny_h, CeedOperator *flux_op) {
  PetscFunctionBeginUser;

  CeedInt   num_comp = 3;
  RDyCells *cells    = &mesh->cells;
  RDyEdges *edges    = &mesh->edges;

  Ceed ceed = CeedContext();

  CeedQFunction qf;
  CeedInt       num_comp_geom = 4, num_comp_cnum = 2;
  PetscCallCEED(CeedQFunctionCreateInterior(ceed, 1, SWEFlux_Roe, SWEFlux_Roe_loc, &qf));
  PetscCallCEED(CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddInput(qf, "q_left", num_comp, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddInput(qf, "q_right", num_comp, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddOutput(qf, "cell_left", num_comp, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddOutput(qf, "cell_right", num_comp, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddOutput(qf, "flux", num_comp, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddOutput(qf, "courant_number", num_comp_cnum, CEED_EVAL_NONE));

  CeedQFunctionContext qf_context;
  PetscCallCEED(CreateQFunctionContext(ceed, tiny_h, &qf_context));
  if (0) PetscCallCEED(CeedQFunctionContextView(qf_context, stdout));
  PetscCallCEED(CeedQFunctionSetContext(qf, qf_context));
  PetscCallCEED(CeedQFunctionContextDestroy(&qf_context));

  CeedElemRestriction q_restrict_l, q_restrict_r, c_restrict_l, c_restrict_r, restrict_geom, restrict_flux, restrict_cnum;
  CeedVector          geom, flux, cnum;
  {
    CeedInt num_edges = mesh->num_owned_internal_edges;

    // create an element restriction for geometric factors that convert
    // fluxes to cell states
    CeedInt g_strides[] = {num_comp_geom, 1, num_comp_geom};
    PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_edges, 1, num_comp_geom, num_edges * num_comp_geom, g_strides, &restrict_geom));
    PetscCallCEED(CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL));
    PetscCallCEED(CeedVectorSetValue(geom, 0.0));

    // create an element restriction for accumulated fluxes
    CeedInt f_strides[] = {num_comp, 1, num_comp};
    PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_edges, 1, num_comp, num_edges * num_comp, f_strides, &restrict_flux));
    PetscCallCEED(CeedElemRestrictionCreateVector(restrict_flux, &flux, NULL));
    PetscCallCEED(CeedVectorSetValue(flux, 0.0));

    // create an element restriction for courant number
    CeedInt cnum_strides[] = {num_comp_cnum, 1, num_comp_cnum};
    PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_edges, 1, num_comp_cnum, num_edges * num_comp_cnum, cnum_strides, &restrict_cnum));
    PetscCallCEED(CeedElemRestrictionCreateVector(restrict_cnum, &cnum, NULL));
    PetscCallCEED(CeedVectorSetValue(cnum, 0.0));

    // create element restrictions for left and right input/output states,
    // populate offsets for these states, and set the (invariant)
    // geometric parameters
    CeedInt *q_offset_l, *q_offset_r, *c_offset_l, *c_offset_r;
    PetscCall(PetscMalloc2(num_edges, &q_offset_l, num_edges, &q_offset_r));
    PetscCall(PetscMalloc2(num_edges, &c_offset_l, num_edges, &c_offset_r));
    CeedScalar(*g)[4];
    PetscCallCEED(CeedVectorGetArray(geom, CEED_MEM_HOST, (CeedScalar **)&g));
    for (CeedInt e = 0, oe = 0; e < mesh->num_internal_edges; e++) {
      CeedInt iedge = edges->internal_edge_ids[e];
      if (!edges->is_owned[iedge]) continue;
      CeedInt l      = edges->cell_ids[2 * iedge];
      CeedInt r      = edges->cell_ids[2 * iedge + 1];
      q_offset_l[oe] = l * num_comp;
      q_offset_r[oe] = r * num_comp;
      c_offset_l[oe] = cells->local_to_owned[l] * num_comp;
      c_offset_r[oe] = cells->local_to_owned[r] * num_comp;

      g[oe][0] = edges->sn[iedge];
      g[oe][1] = edges->cn[iedge];
      g[oe][2] = -edges->lengths[iedge] / cells->areas[l];
      g[oe][3] = edges->lengths[iedge] / cells->areas[r];
      oe++;
    }
    PetscCallCEED(CeedVectorRestoreArray(geom, (CeedScalar **)&g));

    // create element restrictions for left and right cell states
    PetscCallCEED(CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, q_offset_l,
                                            &q_restrict_l));
    PetscCallCEED(CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, q_offset_r,
                                            &q_restrict_r));
    PetscCallCEED(CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, c_offset_l,
                                            &c_restrict_l));
    PetscCallCEED(CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, c_offset_r,
                                            &c_restrict_r));
    PetscCall(PetscFree2(q_offset_l, q_offset_r));
    PetscCall(PetscFree2(c_offset_l, c_offset_r));
    if (0) {
      PetscCallCEED(CeedElemRestrictionView(q_restrict_l, stdout));
      PetscCallCEED(CeedElemRestrictionView(q_restrict_r, stdout));
      PetscCallCEED(CeedElemRestrictionView(c_restrict_l, stdout));
      PetscCallCEED(CeedElemRestrictionView(c_restrict_r, stdout));
    }
  }

  PetscCallCEED(CeedOperatorCreate(ceed, qf, NULL, NULL, flux_op));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "geom", restrict_geom, CEED_BASIS_COLLOCATED, geom));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "q_left", q_restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "q_right", q_restrict_r, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "cell_left", c_restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "cell_right", c_restrict_r, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "flux", restrict_flux, CEED_BASIS_COLLOCATED, flux));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "courant_number", restrict_cnum, CEED_BASIS_COLLOCATED, cnum));

  PetscCallCEED(CeedElemRestrictionDestroy(&restrict_geom));
  PetscCallCEED(CeedElemRestrictionDestroy(&restrict_flux));
  PetscCallCEED(CeedElemRestrictionDestroy(&restrict_cnum));
  PetscCallCEED(CeedElemRestrictionDestroy(&q_restrict_l));
  PetscCallCEED(CeedElemRestrictionDestroy(&q_restrict_r));
  PetscCallCEED(CeedElemRestrictionDestroy(&c_restrict_l));
  PetscCallCEED(CeedElemRestrictionDestroy(&c_restrict_r));
  PetscCallCEED(CeedVectorDestroy(&geom));
  PetscCallCEED(CeedVectorDestroy(&flux));
  PetscCallCEED(CeedVectorDestroy(&cnum));

  PetscFunctionReturn(CEED_ERROR_SUCCESS);
}

static PetscErrorCode CreateSWEBoundaryFluxSubOperator(RDyMesh *mesh, RDyBoundary boundary, RDyCondition boundary_condition, PetscReal tiny_h,
                                                       CeedOperator *flux_op) {
  PetscFunctionBeginUser;

  Ceed ceed = CeedContext();

  CeedInt   num_comp = 3;
  RDyCells *cells    = &mesh->cells;
  RDyEdges *edges    = &mesh->edges;

  CeedQFunctionUser func;
  const char       *func_loc;
  switch (boundary_condition.flow->type) {
    case CONDITION_DIRICHLET:
      func     = SWEBoundaryFlux_Dirichlet_Roe;
      func_loc = SWEBoundaryFlux_Dirichlet_Roe_loc;
      break;
    case CONDITION_REFLECTING:
      func     = SWEBoundaryFlux_Reflecting_Roe;
      func_loc = SWEBoundaryFlux_Reflecting_Roe_loc;
      break;
    case CONDITION_CRITICAL_OUTFLOW:
      func     = SWEBoundaryFlux_Outflow_Roe;
      func_loc = SWEBoundaryFlux_Outflow_Roe_loc;
      break;
    default:
      PetscCheck(PETSC_FALSE, PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid boundary condition encountered for boundary %" PetscInt_FMT "\n",
                 boundary.id);
  }
  CeedQFunction qf;
  CeedInt       num_comp_geom = 3, num_comp_cnum = 1;
  PetscCallCEED(CeedQFunctionCreateInterior(ceed, 1, func, func_loc, &qf));
  PetscCallCEED(CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddInput(qf, "q_left", num_comp, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddOutput(qf, "cell_left", num_comp, CEED_EVAL_NONE));
  if (boundary_condition.flow->type == CONDITION_DIRICHLET) {
    PetscCallCEED(CeedQFunctionAddInput(qf, "q_dirichlet", num_comp, CEED_EVAL_NONE));
  }
  PetscCallCEED(CeedQFunctionAddOutput(qf, "flux", num_comp, CEED_EVAL_NONE));
  PetscCallCEED(CeedQFunctionAddOutput(qf, "courant_number", num_comp_cnum, CEED_EVAL_NONE));

  CeedQFunctionContext qf_context;
  PetscCallCEED(CreateQFunctionContext(ceed, tiny_h, &qf_context));
  if (0) PetscCallCEED(CeedQFunctionContextView(qf_context, stdout));
  PetscCallCEED(CeedQFunctionSetContext(qf, qf_context));
  PetscCallCEED(CeedQFunctionContextDestroy(&qf_context));

  CeedElemRestriction q_restrict_l, c_restrict_l, restrict_dirichlet, restrict_geom, restrict_flux, restrict_cnum;
  CeedVector          geom, flux, dirichlet, cnum;
  {
    CeedInt num_edges = boundary.num_edges;

    // create element restrictions for left and right input/output states
    CeedInt *q_offset_l, *c_offset_l, *offset_dirichlet = NULL;
    PetscCall(PetscMalloc1(num_edges, &q_offset_l));
    PetscCall(PetscMalloc1(num_edges, &c_offset_l));
    if (boundary_condition.flow->type == CONDITION_DIRICHLET) {
      PetscCall(PetscMalloc1(num_edges, &offset_dirichlet));
    }

    // create an element restriction for geometric factors that convert
    // fluxes to cell states
    CeedInt num_owned_edges = 0;
    for (CeedInt e = 0; e < boundary.num_edges; e++) {
      CeedInt iedge = boundary.edge_ids[e];
      if (edges->is_owned[iedge]) num_owned_edges++;
    }
    CeedInt g_strides[] = {num_comp_geom, 1, num_comp_geom};
    PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_owned_edges, 1, num_comp_geom, num_edges * num_comp_geom, g_strides, &restrict_geom));
    PetscCallCEED(CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL));
    PetscCallCEED(CeedVectorSetValue(geom, 0.0));

    // create an element restriction for accumulated fluxes
    CeedInt f_strides[] = {num_comp, 1, num_comp};
    PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_owned_edges, 1, num_comp, num_edges * num_comp, f_strides, &restrict_flux));
    PetscCallCEED(CeedElemRestrictionCreateVector(restrict_flux, &flux, NULL));
    PetscCallCEED(CeedVectorSetValue(flux, 0.0));

    // create an element restriction for courant number
    CeedInt cnum_strides[] = {num_comp_cnum, 1, num_comp_cnum};
    PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_owned_edges, 1, num_comp_cnum, num_edges * num_comp_cnum, cnum_strides, &restrict_cnum));
    PetscCallCEED(CeedElemRestrictionCreateVector(restrict_cnum, &cnum, NULL));
    PetscCallCEED(CeedVectorSetValue(cnum, 0.0));

    // create an element restriction for the "left" (interior) input/output
    // states, populate offsets for these states, and set the (invariant)
    // geometric parameters
    CeedScalar(*g)[3];
    PetscCallCEED(CeedVectorGetArray(geom, CEED_MEM_HOST, (CeedScalar **)&g));
    for (CeedInt e = 0, oe = 0; e < num_edges; e++) {
      CeedInt iedge = boundary.edge_ids[e];
      if (!edges->is_owned[iedge]) continue;
      CeedInt l      = edges->cell_ids[2 * iedge];
      q_offset_l[oe] = l * num_comp;
      c_offset_l[oe] = cells->local_to_owned[l] * num_comp;
      if (offset_dirichlet) {  // Dirichlet boundary values
        offset_dirichlet[oe] = e * num_comp;
      }

      g[oe][0] = edges->sn[iedge];
      g[oe][1] = edges->cn[iedge];
      g[oe][2] = -edges->lengths[iedge] / cells->areas[l];
      oe++;
    }
    PetscCallCEED(CeedVectorRestoreArray(geom, (CeedScalar **)&g));

    // create the element restriction for the left cell states
    PetscCallCEED(CeedElemRestrictionCreate(ceed, num_owned_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES,
                                            q_offset_l, &q_restrict_l));
    PetscCallCEED(CeedElemRestrictionCreate(ceed, num_owned_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES,
                                            c_offset_l, &c_restrict_l));
    PetscCall(PetscFree(q_offset_l));
    PetscCall(PetscFree(c_offset_l));
    if (0) {
      PetscCallCEED(CeedElemRestrictionView(q_restrict_l, stdout));
      PetscCallCEED(CeedElemRestrictionView(c_restrict_l, stdout));
    }

    // if we have Dirichlet boundary values, create a restriction and passive
    // input vector for them
    if (offset_dirichlet) {
      PetscCallCEED(CeedElemRestrictionCreate(ceed, num_owned_edges, 1, num_comp, 1, num_edges * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES,
                                              offset_dirichlet, &restrict_dirichlet));
      PetscCall(PetscFree(offset_dirichlet));
      if (0) PetscCallCEED(CeedElemRestrictionView(restrict_dirichlet, stdout));
      PetscCallCEED(CeedElemRestrictionCreateVector(restrict_dirichlet, &dirichlet, NULL));
      PetscCallCEED(CeedVectorSetValue(dirichlet, 0.0));
    }
  }

  PetscCallCEED(CeedOperatorCreate(ceed, qf, NULL, NULL, flux_op));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "geom", restrict_geom, CEED_BASIS_COLLOCATED, geom));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "q_left", q_restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
  if (boundary_condition.flow->type == CONDITION_DIRICHLET) {
    PetscCallCEED(CeedOperatorSetField(*flux_op, "q_dirichlet", restrict_dirichlet, CEED_BASIS_COLLOCATED, dirichlet));
  }
  PetscCallCEED(CeedOperatorSetField(*flux_op, "cell_left", c_restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "flux", restrict_flux, CEED_BASIS_COLLOCATED, flux));
  PetscCallCEED(CeedOperatorSetField(*flux_op, "courant_number", restrict_cnum, CEED_BASIS_COLLOCATED, cnum));

  PetscCallCEED(CeedElemRestrictionDestroy(&restrict_geom));
  PetscCallCEED(CeedElemRestrictionDestroy(&restrict_flux));
  PetscCallCEED(CeedElemRestrictionDestroy(&restrict_cnum));
  PetscCallCEED(CeedElemRestrictionDestroy(&q_restrict_l));
  PetscCallCEED(CeedElemRestrictionDestroy(&c_restrict_l));
  PetscCallCEED(CeedVectorDestroy(&geom));
  PetscCallCEED(CeedVectorDestroy(&flux));
  PetscCallCEED(CeedVectorDestroy(&cnum));

  PetscFunctionReturn(CEED_ERROR_SUCCESS);
}

// Given a computational mesh, creates a source operator for the shallow water
// equations that computes source terms. The resulting operator can be
// manipulated by libCEED calls.
// @param [in]  ceed The Ceed context used to create the operator
// @param [in]  mesh The computational mesh for which the operator is created
// @param [in]  num_cells Number of cells
// @param [in]  materials_by_cell An array of RDyMaterials defining cellwise material properties
// @param [in]  tiny_h the minimum height threshold for water flow
// @param [out] flux_op A pointer to the flux operator to be created
PetscErrorCode CreateSWESourceOperator(RDyMesh *mesh, PetscInt num_cells, RDyMaterial *materials_by_cell, PetscReal tiny_h, CeedOperator *source_op) {
  PetscFunctionBeginUser;

  Ceed ceed = CeedContext();

  CeedInt   num_comp = 3;
  RDyCells *cells    = &mesh->cells;

  {
    // source term
    CeedQFunction qf;
    CeedInt       num_comp_geom = 2, num_comp_swe_src = 3, num_comp_mannings_n = 1;
    PetscCallCEED(CeedQFunctionCreateInterior(ceed, 1, SWESourceTerm, SWESourceTerm_loc, &qf));
    PetscCallCEED(CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE));
    PetscCallCEED(CeedQFunctionAddInput(qf, "swe_src", num_comp_swe_src, CEED_EVAL_NONE));
    PetscCallCEED(CeedQFunctionAddInput(qf, "mannings_n", num_comp_mannings_n, CEED_EVAL_NONE));
    PetscCallCEED(CeedQFunctionAddInput(qf, "riemannf", num_comp, CEED_EVAL_NONE));
    PetscCallCEED(CeedQFunctionAddInput(qf, "q", num_comp, CEED_EVAL_NONE));
    PetscCallCEED(CeedQFunctionAddOutput(qf, "cell", num_comp, CEED_EVAL_NONE));

    CeedQFunctionContext qf_context;
    PetscCallCEED(CreateQFunctionContext(ceed, tiny_h, &qf_context));
    if (0) PetscCallCEED(CeedQFunctionContextView(qf_context, stdout));
    PetscCallCEED(CeedQFunctionSetContext(qf, qf_context));
    PetscCallCEED(CeedQFunctionContextDestroy(&qf_context));

    CeedElemRestriction restrict_c, restrict_q, restrict_geom, restrict_swe, restrict_mannings_n, restrict_riemannf;
    CeedVector          geom;
    CeedVector          swe_src;
    CeedVector          mannings_n;
    CeedVector          riemannf;
    {  // Create element restrictions for state
      CeedInt *offset_c, *offset_q;
      CeedScalar(*g)[num_comp_geom];
      CeedScalar(*n)[num_comp_mannings_n];
      CeedInt num_owned_cells = mesh->num_owned_cells;
      CeedInt num_cells       = mesh->num_cells;

      CeedInt strides_geom[] = {num_comp_geom, 1, num_comp_geom};
      PetscCallCEED(
          CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp_geom, num_owned_cells * num_comp_geom, strides_geom, &restrict_geom));
      PetscCallCEED(CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL));
      PetscCallCEED(CeedVectorSetValue(geom, 0.0));

      CeedInt strides_swe_src[] = {num_comp_swe_src, 1, num_comp_swe_src};
      PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp_swe_src, num_owned_cells * num_comp_swe_src, strides_swe_src,
                                                     &restrict_swe));
      PetscCallCEED(CeedElemRestrictionCreateVector(restrict_swe, &swe_src, NULL));
      PetscCallCEED(CeedVectorSetValue(swe_src, 0.0));

      CeedInt strides_mannings_n[] = {num_comp_mannings_n, 1, num_comp_mannings_n};
      PetscCallCEED(CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp_mannings_n, num_owned_cells * num_comp_mannings_n,
                                                     strides_mannings_n, &restrict_mannings_n));
      PetscCallCEED(CeedElemRestrictionCreateVector(restrict_mannings_n, &mannings_n, NULL));
      PetscCallCEED(CeedVectorSetValue(mannings_n, 0.0));

      CeedInt strides_riemannf[] = {num_comp, 1, num_comp};
      PetscCallCEED(
          CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp, num_owned_cells * num_comp, strides_riemannf, &restrict_riemannf));
      PetscCallCEED(CeedElemRestrictionCreateVector(restrict_riemannf, &riemannf, NULL));
      PetscCallCEED(CeedVectorSetValue(riemannf, 0.0));

      PetscCall(PetscMalloc1(num_owned_cells, &offset_q));
      PetscCall(PetscMalloc1(num_owned_cells, &offset_c));
      PetscCallCEED(CeedVectorGetArray(geom, CEED_MEM_HOST, (CeedScalar **)&g));
      PetscCallCEED(CeedVectorGetArray(mannings_n, CEED_MEM_HOST, (CeedScalar **)&n));
      for (CeedInt c = 0, oc = 0; c < mesh->num_cells; c++) {
        if (!cells->is_local[c]) continue;

        offset_q[oc] = c * num_comp;
        offset_c[oc] = cells->local_to_owned[c] * num_comp;

        g[oc][0] = cells->dz_dx[c];
        g[oc][1] = cells->dz_dy[c];

        n[oc][0] = materials_by_cell[c].manning;

        oc++;
      }
      PetscCallCEED(CeedVectorRestoreArray(geom, (CeedScalar **)&g));
      PetscCallCEED(CeedVectorRestoreArray(mannings_n, (CeedScalar **)&n));
      PetscCallCEED(CeedElemRestrictionCreate(ceed, num_owned_cells, 1, num_comp, 1, num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, offset_q,
                                              &restrict_q));
      PetscCallCEED(CeedElemRestrictionCreate(ceed, num_owned_cells, 1, num_comp, 1, num_owned_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES,
                                              offset_c, &restrict_c));
      PetscCall(PetscFree(offset_c));
      PetscCall(PetscFree(offset_q));
      if (0) {
        PetscCallCEED(CeedElemRestrictionView(restrict_q, stdout));
        PetscCallCEED(CeedElemRestrictionView(restrict_c, stdout));
      }
    }

    {
      PetscCallCEED(CeedOperatorCreate(ceed, qf, NULL, NULL, source_op));
      PetscCallCEED(CeedOperatorSetField(*source_op, "geom", restrict_geom, CEED_BASIS_COLLOCATED, geom));
      PetscCallCEED(CeedOperatorSetField(*source_op, "swe_src", restrict_swe, CEED_BASIS_COLLOCATED, swe_src));
      PetscCallCEED(CeedOperatorSetField(*source_op, "mannings_n", restrict_mannings_n, CEED_BASIS_COLLOCATED, mannings_n));
      PetscCallCEED(CeedOperatorSetField(*source_op, "riemannf", restrict_riemannf, CEED_BASIS_COLLOCATED, riemannf));
      PetscCallCEED(CeedOperatorSetField(*source_op, "q", restrict_q, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
      PetscCallCEED(CeedOperatorSetField(*source_op, "cell", restrict_c, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE));
    }
    PetscCallCEED(CeedElemRestrictionDestroy(&restrict_geom));
    PetscCallCEED(CeedElemRestrictionDestroy(&restrict_mannings_n));
    PetscCallCEED(CeedElemRestrictionDestroy(&restrict_c));
    PetscCallCEED(CeedElemRestrictionDestroy(&restrict_q));
    PetscCallCEED(CeedVectorDestroy(&geom));
  }

  PetscFunctionReturn(CEED_ERROR_SUCCESS);
}

// Given a shallow water equations source operator created by
// CreateSWESourceOperator, fetches the field representing the Riemann flux.
PetscErrorCode SWESourceOperatorGetRiemannFlux(Operator *op, CeedOperatorField *riemann_flux_field) {
  PetscFunctionBeginUser;

  // get the source sub-operator responsible for the water source (the first one)
  CeedOperator *sub_ops;
  PetscCallCEED(CeedCompositeOperatorGetSubList(op->ceed.source_operator, &sub_ops));
  CeedOperator riemannf_source_op = sub_ops[0];

  // fetch the field
  PetscCallCEED(CeedOperatorGetFieldByName(riemannf_source_op, "riemannf", riemann_flux_field));
  PetscFunctionReturn(CEED_ERROR_SUCCESS);
}

static PetscErrorCode UpdateLocalCourantNumberForInternalEdges(Operator *op, CourantNumberDiagnostics *diags) {
  PetscFunctionBegin;

  CeedOperator interior_flux_op;
  PetscCall(GetCeedInteriorFluxSubOperator(op, &interior_flux_op));

  // fetch the field
  CeedOperatorField courant_num;
  PetscCallCEED(CeedOperatorGetFieldByName(interior_flux_op, "courant_number", &courant_num));

  CeedVector courant_num_vec;
  PetscCallCEED(CeedOperatorFieldGetVector(courant_num, &courant_num_vec));

  CeedScalar(*courant_num_data)[2];  // values to the left/right of an edge
  PetscCallCEED(CeedVectorGetArray(courant_num_vec, CEED_MEM_HOST, (CeedScalar **)&courant_num_data));

  RDyEdges *edges       = &op->mesh->edges;
  RDyCells *cells       = &op->mesh->cells;
  PetscInt  iedge_owned = 0;
  for (PetscInt ii = 0; ii < op->mesh->num_internal_edges; ii++) {
    if (edges->is_owned[ii]) {
      CeedScalar local_max = fmax(courant_num_data[iedge_owned][0], courant_num_data[iedge_owned][1]);
      if (local_max > diags->max_courant_num) {
        diags->max_courant_num = local_max;
        diags->global_edge_id  = edges->global_ids[ii];
        diags->global_cell_id  = cells->global_ids[edges->cell_ids[2 * ii]];
      }
      iedge_owned++;
    }
  }
  PetscCallCEED(CeedVectorRestoreArray(courant_num_vec, (CeedScalar **)&courant_num_data));

  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode UpdateLocalCourantNumberForBoundaryEdges(Operator *op, CourantNumberDiagnostics *diags) {
  PetscFunctionBegin;

  // loop over all boundaries
  for (PetscInt b = 0; b < op->num_boundaries; ++b) {
    RDyBoundary boundary = op->boundaries[b];

    CeedOperator boundary_flux_op;
    PetscCall(GetCeedBoundaryFluxSubOperator(op, boundary, &boundary_flux_op));

    // fetch the field
    CeedOperatorField courant_num;
    PetscCallCEED(CeedOperatorGetFieldByName(boundary_flux_op, "courant_number", &courant_num));

    // get access to the data
    CeedVector courant_num_vec;
    PetscCallCEED(CeedOperatorFieldGetVector(courant_num, &courant_num_vec));
    CeedScalar(*courant_num_data)[1];
    PetscCallCEED(CeedVectorGetArray(courant_num_vec, CEED_MEM_HOST, (CeedScalar **)&courant_num_data));

    // find the maximum value
    RDyEdges *edges = &op->mesh->edges;
    RDyCells *cells = &op->mesh->cells;
    for (PetscInt e = 0; e < boundary.num_edges; ++e) {
      CeedScalar local_max = courant_num_data[e][0];
      if (local_max > diags->max_courant_num) {
        diags->max_courant_num = local_max;
        PetscInt edge_id       = boundary.edge_ids[e];
        diags->global_edge_id  = edges->global_ids[edge_id];
        diags->global_cell_id  = cells->global_ids[edges->cell_ids[2 * edge_id]];
      }
    }

    // restores the pointer
    PetscCallCEED(CeedVectorRestoreArray(courant_num_vec, (CeedScalar **)&courant_num_data));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

/// @brief Finds the global maximum Courant number across all internal and boundary edges.
/// @param [in] op_edges A CeedOperator object for edges
/// @param [in] num_boundaries Total number of boundaries
/// @param [in] boundaries A RDyBoundary object
/// @param [in] comm A MPI_Comm object
/// @param [out] *max_courant_number Global maximum value of courant number
/// @return 0 on sucess, or a non-zero error code on failure
static PetscErrorCode UpdateLocalCourantDiagnostics(Operator *op, CourantNumberDiagnostics *diags) {
  PetscFunctionBegin;

  PetscCall(UpdateLocalCourantNumberForInternalEdges(op, diags));
  PetscCall(UpdateLocalCourantNumberForBoundaryEdges(op, diags));

  PetscFunctionReturn(PETSC_SUCCESS);
}

// FIXME: this is here to support our janky E3SM boundary edge flux coupling;
// FIXME: there's probably a better way to do this
extern PetscErrorCode CreatePetscSWEFluxForBoundaryEdges(RDyEdges *, PetscInt, PetscInt, RDyBoundary *, PetscBool, PetscReal, Operator *);

PetscErrorCode CreateCeedSWEOperator(RDy rdy, Operator *operator) {
  PetscFunctionBegin;

  const PetscReal tiny_h = rdy->config.physics.flow.tiny_h;

  PetscInt num_comp = 3;
  PetscCall(CreateCeedOperator(rdy->dm, &rdy->mesh, num_comp, rdy->num_boundaries, rdy->boundaries, operator));

  CeedOperator interior_flux_op;
  PetscCall(CreateSWEInteriorFluxSubOperator(&rdy->mesh, tiny_h, &interior_flux_op));
  if (0) PetscCallCEED(CeedOperatorView(interior_flux_op, stdout));
  PetscCall(AddCeedInteriorFluxSubOperator(operator, interior_flux_op));

  for (PetscInt b = 0; b < rdy->num_boundaries; ++b) {
    CeedOperator boundary_flux_op;
    RDyBoundary  boundary           = rdy->boundaries[b];
    RDyCondition boundary_condition = rdy->boundary_conditions[b];
    PetscCall(CreateSWEBoundaryFluxSubOperator(&rdy->mesh, boundary, boundary_condition, tiny_h, &boundary_flux_op));
    if (0) PetscCallCEED(CeedOperatorView(boundary_flux_op, stdout));
    PetscCall(AddCeedBoundaryFluxSubOperator(operator, boundary, boundary_flux_op));
  }

  CeedOperator source_op;
  PetscCall(CreateSWESourceOperator(&rdy->mesh, rdy->mesh.num_cells, rdy->materials_by_cell, tiny_h, &source_op));
  PetscCall(AddCeedSourceSubOperator(operator, source_op));

  // create associated vectors for storage
  Ceed ceed = CeedContext();
  PetscCallCEED(CeedVectorCreate(ceed, rdy->mesh.num_cells * num_comp, &operator->ceed.u_local));
  PetscCallCEED(CeedVectorCreate(ceed, rdy->mesh.num_cells * num_comp, &operator->ceed.rhs));
  PetscCallCEED(CeedVectorCreate(ceed, rdy->mesh.num_owned_cells * num_comp, &operator->ceed.sources));

  // FIXME: I'm not sure what to thing of this vvv
  PetscBool ceed_enabled = PETSC_TRUE;
  PetscCall(CreatePetscSWEFluxForBoundaryEdges(&rdy->mesh.edges, num_comp, rdy->num_boundaries, rdy->boundaries, ceed_enabled, tiny_h, operator));

  operator->update_local_courant_diags = UpdateLocalCourantDiagnostics;

  // reset the operator's time step size
  operator->ceed.dt = 0.0;

  PetscFunctionReturn(PETSC_SUCCESS);
}

#pragma GCC diagnostic   pop
#pragma clang diagnostic pop
