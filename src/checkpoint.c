#include <petscviewerhdf5.h>
#include <private/rdycoreimpl.h>

// this struct contains metadata that is written to checkpoint files
typedef struct {
  PetscInt  nproc;  // number of MPI processes / tasks
  PetscReal dt;     // timestep
} CheckpointMetadata;

// creates a new PetscBag that can store checkpoint metadata, populated with
// data from rdy
static PetscErrorCode CreateMetadataBag(RDy rdy, PetscBag *bag) {
  PetscFunctionBegin;

  CheckpointMetadata *metadata;
  PetscCall(PetscBagCreate(rdy->comm, sizeof(CheckpointMetadata), bag));
  PetscCall(PetscBagGetData(*bag, (void **)&metadata));
  PetscCall(PetscBagSetName(*bag, "metadata", "Checkpoint metadata"));
  PetscCall(PetscBagRegisterInt(*bag, &metadata->nproc, rdy->nproc, "nproc", "Number of MPI tasks"));
  PetscCall(PetscBagRegisterReal(*bag, &metadata->dt, rdy->dt, "dt", "Simulation timestep"));

  PetscFunctionReturn(PETSC_SUCCESS);
}

// reads checkpoint metadata into rdy from the given PetscBag and destroys it
static PetscErrorCode ConsumeMetadataBag(RDy rdy, PetscBag bag) {
  PetscFunctionBegin;

  CheckpointMetadata *metadata;
  PetscCall(PetscBagGetData(bag, (void **)&metadata));
  rdy->nproc = metadata->nproc;
  rdy->dt    = metadata->dt;
  PetscCall(PetscBagDestroy(&bag));

  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode WriteCheckpoint(TS ts, PetscInt step, PetscReal time, Vec X, void *ctx) {
  PetscFunctionBegin;
  RDy rdy = ctx;

  PetscViewer viewer;
  char        prefix[PETSC_MAX_PATH_LEN], filename[PETSC_MAX_PATH_LEN];
  PetscCall(DetermineConfigPrefix(rdy, prefix));

  const PetscViewerFormat format = rdy->config.checkpoint.format;
  if (format == PETSC_VIEWER_NATIVE) {  // binary
    PetscCall(GenerateIndexedFilename(prefix, step, rdy->config.time.max_step, "bin", filename));
    PetscCall(PetscViewerBinaryOpen(rdy->comm, filename, FILE_MODE_WRITE, &viewer));
  } else {  // HDF5
    PetscCall(GenerateIndexedFilename(prefix, step, rdy->config.time.max_step, "h5", filename));
    PetscCall(PetscViewerHDF5Open(rdy->comm, filename, FILE_MODE_WRITE, &viewer));
  }

  RDyLogInfo(rdy, "Writing checkpoint file %s...", filename);

  PetscBag bag;
  PetscCall(CreateMetadataBag(rdy, &bag));
  PetscCall(PetscBagView(bag, viewer));
  PetscCall(VecView(X, viewer));
  PetscCall(PetscBagDestroy(&bag));

  PetscCall(PetscViewerDestroy(&viewer));
  RDyLogInfo(rdy, "Finished writing checkpoint file.");
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Initializes the periodic writing of checkpoint files as requested in the
// YAML config file.
PetscErrorCode InitCheckpoints(RDy rdy) {
  PetscFunctionBegin;
  if (rdy->config.checkpoint.interval) {
    const PetscViewerFormat format = rdy->config.checkpoint.format;
    PetscCheck((format == PETSC_VIEWER_NATIVE) || (format == PETSC_VIEWER_HDF5_PETSC), rdy->comm, PETSC_ERR_USER, "Invalid checkpoint format!");
    PetscCall(TSMonitorSet(rdy->ts, WriteCheckpoint, rdy, NULL));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// reads the checkpoint file with the given name, populating the state of the
// given RDy object
PetscErrorCode ReadCheckpointFile(RDy rdy, const char *filename) {
  PetscFunctionBegin;

  RDyLogInfo(rdy, "Reading checkpoint file %s...", filename);

  // determine the file's format from its suffix (it need not match the
  // format specified in the checkpoint section), and open an appropriate
  // viewer.
  PetscViewer viewer;
  if (strstr(filename, ".bin")) {  // binary
    PetscCall(PetscViewerBinaryOpen(rdy->comm, filename, FILE_MODE_READ, &viewer));
  } else if (strstr(filename, ".h5")) {  // HDF5
    PetscCall(PetscViewerHDF5Open(rdy->comm, filename, FILE_MODE_READ, &viewer));
  } else {
    PetscCheck(PETSC_FALSE, rdy->comm, PETSC_ERR_USER, "Invalid checkpoint file: %s", filename);
  }

  PetscBag bag;
  PetscCall(CreateMetadataBag(rdy, &bag));
  PetscCall(PetscBagLoad(viewer, bag));
  PetscCall(ConsumeMetadataBag(rdy, bag));
  PetscCall(VecLoad(rdy->X, viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  RDyLogInfo(rdy, "Finished reading checkpoint file.");

  PetscFunctionReturn(PETSC_SUCCESS);
}
