#include <petscdmceed.h>
#include <private/rdycoreimpl.h>

// global CEED context
static Ceed ceed_context;

/// returns true iff CEED is enabled
PetscBool CeedEnabled(void) { return (ceed_context != NULL); }

/// returns the global CEED context, which is only valid if CeedEnabled()
/// returns PETSC_TRUE
Ceed CeedContext(void) { return ceed_context; }

/// retrieves the appropriate PETSc Vec type for the selected CEED backend
PetscErrorCode GetCeedVecType(VecType *vec_type) {
  PetscFunctionBegin;

  CeedMemType mem_type_backend;
  PetscCallCEED(CeedGetPreferredMemType(ceed_context, &mem_type_backend));
  switch (mem_type_backend) {
    case CEED_MEM_HOST:
      *vec_type = VECSTANDARD;
      break;
    case CEED_MEM_DEVICE: {
      const char *resolved;
      PetscCallCEED(CeedGetResource(ceed_context, &resolved));
      if (strstr(resolved, "/gpu/cuda")) *vec_type = VECCUDA;
      else if (strstr(resolved, "/gpu/hip")) *vec_type = VECKOKKOS;
      else if (strstr(resolved, "/gpu/sycl")) *vec_type = VECKOKKOS;
      else *vec_type = VECSTANDARD;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/// Sets the CEED resource string to the given string, initializing the global
/// CEED context if the argument is specified. If CEED has already been enabled
/// with a different resource, a call to this function deletes the global
/// context and recreates it. If the resource string is empty or NULL, CEED is
/// disabled.
/// @param resource a CEED resource string, possibly empty or NULL
PetscErrorCode SetCeedResource(char *resource) {
  PetscFunctionBegin;
  if (ceed_context) {  // we already have a context
    CeedDestroy(&ceed_context);
    ceed_context = NULL;
  }
  if (resource && resource[0]) {
    PetscCallCEED(CeedInit(resource, &ceed_context));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
