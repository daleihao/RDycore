#include <private/rdycoreimpl.h>
#include <private/rdymathimpl.h>
#include <private/rdymemoryimpl.h>
#include <stddef.h>  // for offsetof

#include "swe_flux_ceed.h"
#include "swe_flux_petsc.h"

PetscClassId  RDY_CLASSID;
PetscLogEvent RDY_CeedOperatorApply;

//-----------------------
// Debugging diagnostics
//-----------------------

// MPI datatype corresponding to CourantNumberDiagnostics. Created during InitSWE.
static MPI_Datatype courant_num_diags_type;

// MPI operator used to determine the prevailing diagnostics for the maximum
// courant number on all processes. Created during InitSWE.
static MPI_Op courant_num_diags_op;

// function backing the above MPI operator
static void FindCourantNumberDiagnostics(void *in_vec, void *result_vec, int *len, MPI_Datatype *type) {
  CourantNumberDiagnostics *in_diags     = in_vec;
  CourantNumberDiagnostics *result_diags = result_vec;

  // select the item with the maximum courant number
  for (int i = 0; i < *len; ++i) {
    if (in_diags[i].max_courant_num > result_diags[i].max_courant_num) {
      result_diags[i] = in_diags[i];
    }
  }
}

// this function destroys the MPI type and operator associated with CourantNumberDiagnostics
static void DestroyCourantNumberDiagnostics(void) {
  MPI_Op_free(&courant_num_diags_op);
  MPI_Type_free(&courant_num_diags_type);
}

// this function is called by InitSWE to initialize the above type(s) and op(s).
static PetscErrorCode InitMPITypesAndOps(void) {
  PetscFunctionBegin;

  // create an MPI data type for the CourantNumberDiagnostics struct
  const int      num_blocks             = 3;
  const int      block_lengths[3]       = {1, 1, 1};
  const MPI_Aint block_displacements[3] = {
      offsetof(CourantNumberDiagnostics, max_courant_num),
      offsetof(CourantNumberDiagnostics, global_edge_id),
      offsetof(CourantNumberDiagnostics, global_cell_id),
  };
  MPI_Datatype block_types[3] = {MPI_DOUBLE, MPI_INT, MPI_INT};
  MPI_Type_create_struct(num_blocks, block_lengths, block_displacements, block_types, &courant_num_diags_type);
  MPI_Type_commit(&courant_num_diags_type);

  // create a corresponding reduction operator for the new type
  MPI_Op_create(FindCourantNumberDiagnostics, 1, &courant_num_diags_op);

  // make sure the operator and the type are destroyed upon exit
  PetscCall(RDyOnFinalize(DestroyCourantNumberDiagnostics));

  PetscFunctionReturn(0);
}

//---------------------------
// End debugging diagnostics
//---------------------------

// This function initializes SWE physics for the given dycore.
PetscErrorCode InitSWE(RDy rdy) {
  PetscFunctionBeginUser;

  // set up MPI types and operators used by SWE physics
  PetscCall(InitMPITypesAndOps());

  PetscCall(PetscClassIdRegister("RDycore", &RDY_CLASSID));
  PetscCall(PetscLogEventRegister("CeedOperatorApp", RDY_CLASSID, &RDY_CeedOperatorApply));
  PetscFunctionReturn(0);
}

static PetscErrorCode RDyCeedOperatorSetUp(RDy rdy) {
  PetscFunctionBeginUser;
  if (rdy->ceed_resource[0] && !rdy->ceed_rhs.op) {
    Ceed ceed;
    CeedInit(rdy->ceed_resource, &ceed);
    CeedCompositeOperatorCreate(ceed, &rdy->ceed_rhs.op);
    CeedInt   num_comp = 3;
    RDyMesh  *mesh     = &rdy->mesh;
    RDyCells *cells    = &mesh->cells;
    RDyEdges *edges    = &mesh->edges;
    {  // interior operator
      CeedQFunction qf;
      CeedInt       num_comp_geom = 4;
      CeedQFunctionCreateInterior(ceed, 1, SWEFlux_Roe, SWEFlux_Roe_loc, &qf);
      CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE);
      CeedQFunctionAddInput(qf, "q_left", num_comp, CEED_EVAL_NONE);
      CeedQFunctionAddInput(qf, "q_right", num_comp, CEED_EVAL_NONE);
      CeedQFunctionAddOutput(qf, "cell_left", num_comp, CEED_EVAL_NONE);
      CeedQFunctionAddOutput(qf, "cell_right", num_comp, CEED_EVAL_NONE);

      CeedElemRestriction restrict_l, restrict_r, restrict_geom;
      CeedVector          geom;
      {  // Create element restrictions for state
        CeedInt *offset_l, *offset_r;
        CeedScalar(*g)[4];
        CeedInt num_edges = mesh->num_owned_internal_edges;
        CeedInt strides[] = {num_comp_geom, 1, num_comp_geom};
        CeedElemRestrictionCreateStrided(ceed, num_edges, 1, num_comp_geom, num_edges * num_comp_geom, strides, &restrict_geom);
        CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL);
        CeedVectorSetValue(geom, 0.0);  // initialize to ensure the arrays is allocated
        PetscCall(PetscMalloc2(num_edges, &offset_l, num_edges, &offset_r));
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
        CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, offset_l,
                                  &restrict_l);
        CeedElemRestrictionCreate(ceed, num_edges, 1, num_comp, 1, mesh->num_cells * num_comp, CEED_MEM_HOST, CEED_COPY_VALUES, offset_r,
                                  &restrict_r);
        PetscCall(PetscFree2(offset_l, offset_r));
        if (0) {
          CeedElemRestrictionView(restrict_l, stdout);
          CeedElemRestrictionView(restrict_r, stdout);
        }

        CeedVectorCreate(ceed, mesh->num_cells * num_comp, &rdy->ceed_rhs.x_ceed);
        CeedVectorCreate(ceed, mesh->num_cells * num_comp, &rdy->ceed_rhs.y_ceed);
      }

      {
        CeedOperator op;
        CeedOperatorCreate(ceed, qf, NULL, NULL, &op);
        CeedOperatorSetField(op, "geom", restrict_geom, CEED_BASIS_COLLOCATED, geom);
        CeedOperatorSetField(op, "q_left", restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
        CeedOperatorSetField(op, "q_right", restrict_r, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
        CeedOperatorSetField(op, "cell_left", restrict_l, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
        CeedOperatorSetField(op, "cell_right", restrict_r, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
        CeedOperatorSetNumQuadraturePoints(op, 1);
        CeedCompositeOperatorAddSub(rdy->ceed_rhs.op, op);
        CeedOperatorDestroy(&op);
      }
      CeedElemRestrictionDestroy(&restrict_geom);
      CeedElemRestrictionDestroy(&restrict_l);
      CeedElemRestrictionDestroy(&restrict_r);
      CeedVectorDestroy(&geom);
    }
    for (PetscInt b = 0; b < rdy->num_boundaries; b++) {
      RDyBoundary      *boundary = &rdy->boundaries[b];
      CeedQFunctionUser func;
      const char       *func_loc;
      switch (rdy->boundary_conditions[b].flow->type) {
        case CONDITION_REFLECTING:
          func     = SWEBoundaryFlux_Reflecting_Roe;
          func_loc = SWEBoundaryFlux_Reflecting_Roe_loc;
          break;
        case CONDITION_CRITICAL_OUTFLOW:
          func     = SWEBoundaryFlux_Outflow_Roe;
          func_loc = SWEBoundaryFlux_Outflow_Roe_loc;
          break;
        default:
          PetscCheck(PETSC_FALSE, rdy->comm, PETSC_ERR_USER, "Invalid boundary condition encountered for boundary %d\n", rdy->boundary_ids[b]);
      }
      CeedQFunction qf;
      CeedInt       num_comp_geom = 3;
      CeedQFunctionCreateInterior(ceed, 1, func, func_loc, &qf);
      CeedQFunctionAddInput(qf, "geom", num_comp_geom, CEED_EVAL_NONE);
      CeedQFunctionAddInput(qf, "q_left", num_comp, CEED_EVAL_NONE);
      CeedQFunctionAddOutput(qf, "cell_left", num_comp, CEED_EVAL_NONE);

      CeedElemRestriction restrict_l, restrict_geom;
      CeedVector          geom;
      {  // Create element restrictions for state
        CeedInt *offset_l;
        CeedScalar(*g)[3];
        CeedInt num_edges = boundary->num_edges, num_owned_edges = 0;
        CeedInt strides[] = {num_comp_geom, 1, num_comp_geom};
        for (CeedInt e = 0; e < boundary->num_edges; e++) {
          PetscInt iedge = boundary->edge_ids[e];
          if (edges->is_owned[iedge]) num_owned_edges++;
        }
        CeedElemRestrictionCreateStrided(ceed, num_owned_edges, 1, num_comp_geom, num_edges * num_comp_geom, strides, &restrict_geom);
        CeedElemRestrictionCreateVector(restrict_geom, &geom, NULL);
        CeedVectorSetValue(geom, 0.0);  // initialize to ensure the arrays is allocated
        PetscCall(PetscMalloc1(num_edges, &offset_l));
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
        CeedOperatorSetNumQuadraturePoints(op, 1);
        CeedCompositeOperatorAddSub(rdy->ceed_rhs.op, op);
        CeedOperatorDestroy(&op);
      }
      CeedElemRestrictionDestroy(&restrict_geom);
      CeedElemRestrictionDestroy(&restrict_l);
      CeedVectorDestroy(&geom);
    }
    if (0) CeedOperatorView(rdy->ceed_rhs.op, stdout);
  }
  PetscFunctionReturn(0);
}

static inline CeedMemType MemTypeP2C(PetscMemType mem_type) { return PetscMemTypeDevice(mem_type) ? CEED_MEM_DEVICE : CEED_MEM_HOST; }

static PetscErrorCode RDyCeedOperatorApply(RDy rdy, Vec U_local, Vec F) {
  PetscScalar *u, *f;
  PetscMemType mem_type;
  Vec          F_local;
  PetscFunctionBeginUser;
  PetscCall(RDyCeedOperatorSetUp(rdy));
  CeedVector u_ceed = rdy->ceed_rhs.x_ceed;
  CeedVector f_ceed = rdy->ceed_rhs.y_ceed;
  PetscCall(DMGetLocalVector(rdy->dm, &F_local));
  PetscCall(VecGetArrayAndMemType(U_local, &u, &mem_type));
  CeedVectorSetArray(u_ceed, MemTypeP2C(mem_type), CEED_USE_POINTER, u);
  PetscCall(VecGetArrayAndMemType(F_local, &f, &mem_type));
  CeedVectorSetArray(f_ceed, MemTypeP2C(mem_type), CEED_USE_POINTER, f);

  PetscCall(PetscLogEventBegin(RDY_CeedOperatorApply, U_local, F, 0, 0));
  PetscCall(PetscLogGpuTimeBegin());
  CeedOperatorApply(rdy->ceed_rhs.op, u_ceed, f_ceed, CEED_REQUEST_IMMEDIATE);
  PetscCall(PetscLogGpuTimeEnd());
  PetscCall(PetscLogEventEnd(RDY_CeedOperatorApply, U_local, F, 0, 0));

  CeedVectorTakeArray(u_ceed, MemTypeP2C(mem_type), &u);
  PetscCall(VecRestoreArrayAndMemType(U_local, &u));
  CeedVectorTakeArray(f_ceed, MemTypeP2C(mem_type), &f);
  PetscCall(VecRestoreArrayAndMemType(F_local, &f));
  PetscCall(VecZeroEntries(F));
  PetscCall(DMLocalToGlobal(rdy->dm, F_local, ADD_VALUES, F));
  PetscCall(DMRestoreLocalVector(rdy->dm, &F_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// This is the right-hand-side function used by our timestepping solver for
// the shallow water equations.
// Parameters:
//  ts  - the solver
//  t   - the simulation time [seconds]
//  X   - the solution vector at time t
//  F   - the right hand side vector to be evaluated at time t
//  ctx - a generic pointer to our RDy object
PetscErrorCode RHSFunctionSWE(TS ts, PetscReal t, Vec X, Vec F, void *ctx) {
  PetscFunctionBegin;

  RDy rdy = ctx;
  DM  dm  = rdy->dm;

  PetscCall(VecZeroEntries(F));

  // populate the local X vector
  PetscCall(DMGlobalToLocalBegin(dm, X, INSERT_VALUES, rdy->X_local));
  PetscCall(DMGlobalToLocalEnd(dm, X, INSERT_VALUES, rdy->X_local));

  // compute the right hand side
  CourantNumberDiagnostics courant_num_diags = {
      .max_courant_num = 0.0,
      .global_edge_id  = -1,
      .global_cell_id  = -1,
  };
  if (rdy->ceed_resource[0]) {
    PetscCall(RDyCeedOperatorApply(rdy, rdy->X_local, F));
  } else {
    PetscCall(RHSFunctionForInternalEdges(rdy, F, &courant_num_diags));
    PetscCall(RHSFunctionForBoundaryEdges(rdy, F, &courant_num_diags));
    PetscCall(AddSourceTerm(rdy, F));  // TODO: move source term to use libCEED
  }

  // write out debugging info for maximum courant number
  if (rdy->config.logging.level >= LOG_DEBUG) {
    MPI_Allreduce(MPI_IN_PLACE, &courant_num_diags, 1, courant_num_diags_type, courant_num_diags_op, rdy->comm);
    PetscReal time;
    PetscInt  stepnum;
    PetscCall(TSGetTime(ts, &time));
    PetscCall(TSGetStepNumber(ts, &stepnum));
    RDyLogDebug(rdy, "[%d] Time = %f Max courant number %g encountered at edge %d of cell %d is %f", stepnum, time, courant_num_diags.max_courant_num,
                courant_num_diags.global_edge_id, courant_num_diags.global_cell_id, courant_num_diags.max_courant_num);
  }

  PetscFunctionReturn(0);
}