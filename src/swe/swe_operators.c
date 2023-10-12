#include "swe_operators.h"

#include "swe_operators_impl.h"

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

  CeedQFunctionContextCreate(ceed, qf_context);
  CeedQFunctionContextSetData(*qf_context, CEED_MEM_HOST, CEED_USE_POINTER, sizeof(*swe_ctx), swe_ctx);
  CeedQFunctionContextSetDataDestroy(*qf_context, CEED_MEM_HOST, FreeContextPetsc);
  CeedQFunctionContextRegisterDouble(*qf_context, "time step", offsetof(struct SWEContext_, dtime), 1, "Time step of TS");
  CeedQFunctionContextRegisterDouble(*qf_context, "small h value", offsetof(struct SWEContext_, tiny_h), 1,
                                     "Height threshold below which dry condition is assumed");
  CeedQFunctionContextRegisterDouble(*qf_context, "gravity", offsetof(struct SWEContext_, gravity), 1, "Accelaration due to gravity");

  PetscFunctionReturn(0);
}

PetscErrorCode CreateSWEFluxOperator(Ceed ceed, RDyMesh *mesh, int num_boundaries, RDyBoundary boundaries[num_boundaries],
                                     RDyCondition boundary_conditions[num_boundaries], PetscReal tiny_h, CeedOperator *flux_op) {
  PetscFunctionBeginUser;

  CeedCompositeOperatorCreate(ceed, flux_op);
  CeedInt   num_comp = 3;
  RDyCells *cells    = &mesh->cells;
  RDyEdges *edges    = &mesh->edges;
  // interior sub-operator
  {
    CeedQFunction qf;
    CeedInt       num_comp_geom = 4;
    CeedQFunctionCreateInterior(ceed, 1, SWEFlux_Roe, SWEFlux_Roe_loc, &qf);
    CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE);
    CeedQFunctionAddInput(qf, "q_left", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddInput(qf, "q_right", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf, "cell_left", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf, "cell_right", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf, "flux", num_comp, CEED_EVAL_NONE);

    CeedQFunctionContext qf_context;
    CreateQFunctionContext(ceed, tiny_h, &qf_context);
    if (0) CeedQFunctionContextView(qf_context, stdout);
    CeedQFunctionSetContext(qf, qf_context);
    CeedQFunctionContextDestroy(&qf_context);

    CeedElemRestriction restrict_l, restrict_r, restrict_geom, restrict_flux;
    CeedVector          geom, flux;
    {
      CeedInt num_edges = mesh->num_owned_internal_edges;

      // create an element restriction for geometric factors that convert
      // fluxes to cell states
      CeedInt g_strides[] = {num_comp_geom, 1, num_comp_geom};
      CeedElemRestrictionCreateStrided(ceed, num_edges, 1, num_comp_geom, num_edges * num_comp_geom, g_strides, &restrict_geom);
      CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL);
      CeedVectorSetValue(geom, 0.0);

      // create an element restriction for accumulated fluxes
      CeedInt f_strides[] = {num_comp, 1, num_comp};
      CeedElemRestrictionCreateStrided(ceed, num_edges, 1, num_comp, num_edges * num_comp, f_strides, &restrict_flux);
      CeedElemRestrictionCreateVector(restrict_flux, &flux, NULL);
      CeedVectorSetValue(flux, 0.0);

      // create element restrictions for left and right input/output states,
      // populate offsets for these states, and set the (invariant)
      // geometric parameters
      CeedInt *offset_l, *offset_r;
      PetscCall(PetscMalloc2(num_edges, &offset_l, num_edges, &offset_r));
      CeedScalar(*g)[4];
      CeedVectorGetArray(geom, CEED_MEM_HOST, (CeedScalar **)&g);
      for (CeedInt e = 0, oe = 0; e < mesh->num_internal_edges; e++) {
        PetscInt iedge = edges->internal_edge_ids[e];
        if (!edges->is_owned[iedge]) continue;
        PetscInt l   = edges->cell_ids[2 * iedge];
        PetscInt r   = edges->cell_ids[2 * iedge + 1];
        offset_l[oe] = l * num_comp;
        offset_r[oe] = r * num_comp;

        g[oe][0] = edges->sn[iedge];
        g[oe][1] = edges->cn[iedge];
        g[oe][2] = -edges->lengths[iedge] / cells->areas[l];
        g[oe][3] = edges->lengths[iedge] / cells->areas[r];
        oe++;
      }
      CeedVectorRestoreArray(geom, (CeedScalar **)&g);

      // create element restrictions for left and right cell states
      CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, offset_l, &restrict_l);
      CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, offset_r, &restrict_r);
      PetscCall(PetscFree2(offset_l, offset_r));
      if (0) {
        CeedElemRestrictionView(restrict_l, stdout);
        CeedElemRestrictionView(restrict_r, stdout);
      }
    }

    {
      CeedOperator op;
      CeedOperatorCreate(ceed, qf, NULL, NULL, &op);
      CeedOperatorSetField(op, "geom", restrict_geom, CEED_BASIS_COLLOCATED, geom);
      CeedOperatorSetField(op, "q_left", restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedOperatorSetField(op, "q_right", restrict_r, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedOperatorSetField(op, "cell_left", restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedOperatorSetField(op, "cell_right", restrict_r, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedOperatorSetField(op, "flux", restrict_flux, CEED_BASIS_COLLOCATED, flux);
      CeedCompositeOperatorAddSub(*flux_op, op);
      CeedOperatorDestroy(&op);
    }
    CeedElemRestrictionDestroy(&restrict_geom);
    CeedElemRestrictionDestroy(&restrict_flux);
    CeedElemRestrictionDestroy(&restrict_l);
    CeedElemRestrictionDestroy(&restrict_r);
    CeedVectorDestroy(&geom);
    CeedVectorDestroy(&flux);
  }

  // boundary sub-operators
  for (PetscInt b = 0; b < num_boundaries; b++) {
    RDyBoundary      *boundary = &boundaries[b];
    CeedQFunctionUser func;
    const char       *func_loc;
    switch (boundary_conditions[b].flow->type) {
      case CONDITION_REFLECTING:
        func     = SWEBoundaryFlux_Reflecting_Roe;
        func_loc = SWEBoundaryFlux_Reflecting_Roe_loc;
        break;
      case CONDITION_CRITICAL_OUTFLOW:
        func     = SWEBoundaryFlux_Outflow_Roe;
        func_loc = SWEBoundaryFlux_Outflow_Roe_loc;
        break;
      default:
        PetscCheck(PETSC_FALSE, PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid boundary condition encountered for boundary %d\n", boundaries[b].id);
    }
    CeedQFunction qf;
    CeedInt       num_comp_geom = 3;
    CeedQFunctionCreateInterior(ceed, 1, func, func_loc, &qf);
    CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE);
    CeedQFunctionAddInput(qf, "q_left", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf, "cell_left", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf, "flux", num_comp, CEED_EVAL_NONE);

    CeedQFunctionContext qf_context;
    CreateQFunctionContext(ceed, tiny_h, &qf_context);
    if (0) CeedQFunctionContextView(qf_context, stdout);
    CeedQFunctionSetContext(qf, qf_context);
    CeedQFunctionContextDestroy(&qf_context);

    CeedElemRestriction restrict_l, restrict_geom, restrict_flux;
    CeedVector          geom, flux;
    {
      CeedInt num_edges = boundary->num_edges;

      // create element restrictions for left and right input/output states
      CeedInt *offset_l;
      PetscCall(PetscMalloc1(num_edges, &offset_l));

      // create an element restriction for geometric factors that convert
      // fluxes to cell states
      CeedInt num_owned_edges = 0;
      for (CeedInt e = 0; e < boundary->num_edges; e++) {
        PetscInt iedge = boundary->edge_ids[e];
        if (edges->is_owned[iedge]) num_owned_edges++;
      }
      CeedInt g_strides[] = {num_comp_geom, 1, num_comp_geom};
      CeedElemRestrictionCreateStrided(ceed, num_owned_edges, 1, num_comp_geom, num_edges * num_comp_geom, g_strides, &restrict_geom);
      CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL);
      CeedVectorSetValue(geom, 0.0);  // initialize to ensure the arrays is allocated
      //
      // create an element restriction for accumulated fluxes
      CeedInt f_strides[] = {num_comp, 1, num_comp};
      CeedElemRestrictionCreateStrided(ceed, num_owned_edges, 1, num_comp, num_edges * num_comp, f_strides, &restrict_flux);
      CeedElemRestrictionCreateVector(restrict_flux, &flux, NULL);
      CeedVectorSetValue(flux, 0.0);

      // create an element restrictions for the "left" (interior) input/output
      // states, populate offsets for these states, and set the (invariant)
      // geometric parameters
      CeedScalar(*g)[3];
      CeedVectorGetArray(geom, CEED_MEM_HOST, (CeedScalar **)&g);
      for (CeedInt e = 0, oe = 0; e < num_edges; e++) {
        PetscInt iedge = boundary->edge_ids[e];
        if (!edges->is_owned[iedge]) continue;
        PetscInt l   = edges->cell_ids[2 * iedge];
        offset_l[oe] = l * num_comp;

        g[oe][0] = edges->sn[iedge];
        g[oe][1] = edges->cn[iedge];
        g[oe][2] = -edges->lengths[iedge] / cells->areas[l];
        oe++;
      }
      CeedVectorRestoreArray(geom, (CeedScalar **)&g);

      // create the element restriction for the left cell states
      CeedElemRestrictionCreate(ceed, num_owned_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, offset_l,
                                &restrict_l);
      PetscCall(PetscFree(offset_l));
      if (0) CeedElemRestrictionView(restrict_l, stdout);
    }

    {
      CeedOperator op;
      CeedOperatorCreate(ceed, qf, NULL, NULL, &op);
      CeedOperatorSetField(op, "geom", restrict_geom, CEED_BASIS_COLLOCATED, geom);
      CeedOperatorSetField(op, "q_left", restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedOperatorSetField(op, "cell_left", restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedOperatorSetField(op, "flux", restrict_flux, CEED_BASIS_COLLOCATED, flux);
      CeedCompositeOperatorAddSub(*flux_op, op);
      CeedOperatorDestroy(&op);
    }
    CeedElemRestrictionDestroy(&restrict_geom);
    CeedElemRestrictionDestroy(&restrict_flux);
    CeedElemRestrictionDestroy(&restrict_l);
    CeedVectorDestroy(&geom);
    CeedVectorDestroy(&flux);
  }

  if (0) CeedOperatorView(*flux_op, stdout);
  PetscFunctionReturn(0);
}

PetscErrorCode CreateSWESourceOperator(Ceed ceed, RDyMesh *mesh, RDyMaterial materials_by_cell[mesh->num_cells], PetscReal tiny_h,
                                       CeedOperator *source_op) {
  PetscFunctionBeginUser;

  CeedCompositeOperatorCreate(ceed, source_op);
  CeedInt   num_comp = 3;
  RDyCells *cells    = &mesh->cells;

  {
    // source term
    CeedQFunction qf;
    CeedInt       num_comp_geom = 2, num_comp_water_src = 1, num_comp_mannings_n = 1;
    CeedQFunctionCreateInterior(ceed, 1, SWESourceTerm, SWESourceTerm_loc, &qf);
    CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE);
    CeedQFunctionAddInput(qf, "water_src", num_comp_water_src, CEED_EVAL_NONE);
    CeedQFunctionAddInput(qf, "mannings_n", num_comp_mannings_n, CEED_EVAL_NONE);
    CeedQFunctionAddInput(qf, "riemannf", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddInput(qf, "q", num_comp, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf, "cell", num_comp, CEED_EVAL_NONE);

    CeedQFunctionContext qf_context;
    CreateQFunctionContext(ceed, tiny_h, &qf_context);
    if (0) CeedQFunctionContextView(qf_context, stdout);
    CeedQFunctionSetContext(qf, qf_context);
    CeedQFunctionContextDestroy(&qf_context);

    CeedElemRestriction restrict_c, restrict_geom, restrict_water_src, restrict_mannings_n, restrict_riemannf;
    CeedVector          geom;
    CeedVector          water_src;
    CeedVector          mannings_n;
    CeedVector          riemannf;
    {  // Create element restrictions for state
      CeedInt *offset_c;
      CeedScalar(*g)[num_comp_geom];
      CeedScalar(*n)[num_comp_mannings_n];
      CeedInt num_owned_cells = mesh->num_cells_local;

      CeedInt strides_geom[] = {num_comp_geom, 1, num_comp_geom};
      CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp_geom, num_owned_cells * num_comp_geom, strides_geom, &restrict_geom);
      CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL);
      CeedVectorSetValue(geom, 0.0);  // initialize to ensure the arrays is allocated

      CeedInt strides_water_src[] = {num_comp_water_src, 1, num_comp_water_src};
      CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp_water_src, num_owned_cells * num_comp_water_src, strides_water_src,
                                       &restrict_water_src);
      CeedElemRestrictionCreateVector(restrict_water_src, &water_src, NULL);
      CeedVectorSetValue(water_src, 0.0);  // initialize to ensure the arrays is allocated

      CeedInt strides_mannings_n[] = {num_comp_mannings_n, 1, num_comp_mannings_n};
      CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp_mannings_n, num_owned_cells * num_comp_mannings_n, strides_mannings_n,
                                       &restrict_mannings_n);
      CeedElemRestrictionCreateVector(restrict_mannings_n, &mannings_n, NULL);
      CeedVectorSetValue(mannings_n, 0.0);  // initialize to ensure the arrays is allocated

      CeedInt strides_riemannf[] = {num_comp, 1, num_comp};
      CeedElemRestrictionCreateStrided(ceed, num_owned_cells, 1, num_comp, num_owned_cells * num_comp, strides_riemannf, &restrict_riemannf);
      CeedElemRestrictionCreateVector(restrict_riemannf, &riemannf, NULL);
      CeedVectorSetValue(riemannf, 0.0);  // initialize to ensure the arrays is allocated

      PetscCall(PetscMalloc1(num_owned_cells, &offset_c));
      CeedVectorGetArray(geom, CEED_MEM_HOST, (CeedScalar **)&g);
      CeedVectorGetArray(mannings_n, CEED_MEM_HOST, (CeedScalar **)&n);
      for (CeedInt c = 0, oc = 0; c < mesh->num_cells; c++) {
        if (!cells->is_local[c]) continue;

        offset_c[oc] = c * num_comp;

        g[oc][0] = cells->dz_dx[c];
        g[oc][1] = cells->dz_dy[c];

        n[oc][0] = materials_by_cell[c].manning;

        oc++;
      }
      CeedVectorRestoreArray(geom, (CeedScalar **)&g);
      CeedVectorRestoreArray(mannings_n, (CeedScalar **)&n);
      CeedElemRestrictionCreate(ceed, num_owned_cells, 1, num_comp, 1, num_owned_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, offset_c,
                                &restrict_c);
      PetscCall(PetscFree(offset_c));
      if (0) CeedElemRestrictionView(restrict_c, stdout);
    }

    {
      CeedOperator op;
      CeedOperatorCreate(ceed, qf, NULL, NULL, &op);
      CeedOperatorSetField(op, "geom", restrict_geom, CEED_BASIS_COLLOCATED, geom);
      CeedOperatorSetField(op, "water_src", restrict_water_src, CEED_BASIS_COLLOCATED, water_src);
      CeedOperatorSetField(op, "mannings_n", restrict_water_src, CEED_BASIS_COLLOCATED, mannings_n);
      CeedOperatorSetField(op, "riemannf", restrict_riemannf, CEED_BASIS_COLLOCATED, riemannf);
      CeedOperatorSetField(op, "q", restrict_c, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedOperatorSetField(op, "cell", restrict_c, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
      CeedCompositeOperatorAddSub(*source_op, op);
      CeedOperatorDestroy(&op);
    }
    CeedElemRestrictionDestroy(&restrict_geom);
    CeedElemRestrictionDestroy(&restrict_mannings_n);
    CeedElemRestrictionDestroy(&restrict_c);
    CeedVectorDestroy(&geom);
  }

  if (0) CeedOperatorView(*source_op, stdout);

  PetscFunctionReturn(0);
}

PetscErrorCode GetWaterSourceFromSWESourceOperator(CeedOperator source_op, CeedOperatorField *water_source_field) {
  PetscFunctionBeginUser;

  // get the source sub-operator responsible for the water source (the first one)
  CeedOperator *sub_ops;
  CeedCompositeOperatorGetSubList(source_op, &sub_ops);
  CeedOperator water_source_op = sub_ops[0];

  // fetch the field
  CeedOperatorGetFieldByName(water_source_op, "water_src", water_source_field);
  PetscFunctionReturn(0);
}

PetscErrorCode GetRiemannFluxFromSWESourceOperator(CeedOperator source_op, CeedOperatorField *riemann_flux_field) {
  PetscFunctionBeginUser;

  // get the source sub-operator responsible for the water source (the first one)
  CeedOperator *sub_ops;
  CeedCompositeOperatorGetSubList(source_op, &sub_ops);
  CeedOperator riemannf_source_op = sub_ops[0];

  // fetch the field
  CeedOperatorGetFieldByName(riemannf_source_op, "riemannf", riemann_flux_field);
  PetscFunctionReturn(0);
}
