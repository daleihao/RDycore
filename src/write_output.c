#include <errno.h>
#include <petscviewerhdf5.h>
#include <private/rdycoreimpl.h>
#include <rdycore.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

// output directory name (relative to current working directory)
static const char *output_dir = "output";

PetscErrorCode CreateOutputDir(RDy rdy) {
  PetscFunctionBegin;

  RDyLogDebug(rdy, "Creating output directory %s...", output_dir);

  int result_and_errno[2];
  if (rdy->rank == 0) {
    result_and_errno[0] = mkdir(output_dir, 0755);
    result_and_errno[1] = errno;
  }
  MPI_Bcast(&result_and_errno, 2, MPI_INT, 0, rdy->comm);
  int result = result_and_errno[0];
  int err_no = result_and_errno[1];
  PetscCheck((result == 0) || (err_no == EEXIST), rdy->comm, PETSC_ERR_USER, "Could not create output directory: %s (errno = %d)", output_dir,
             err_no);

  PetscFunctionReturn(0);
}

static PetscErrorCode DetermineOutputFile(RDy rdy, PetscInt step, PetscReal time, const char *suffix, char *filename) {
  PetscFunctionBegin;

  char *p = strstr(rdy->config_file, ".yaml");
  if (!p) {  // could be .yml, I suppose (Windows habits die hard!)
    p = strstr(rdy->config_file, ".yml");
  }
  if (p) {
    size_t prefix_len = p - rdy->config_file;
    size_t n          = strlen(output_dir) + 1 + prefix_len;
    snprintf(filename, n + 1, "%s/%s", output_dir, rdy->config_file);
  } else {
    snprintf(filename, PETSC_MAX_PATH_LEN - 1, "%s/%s", output_dir, rdy->config_file);
  }

  // encode specific information into the filename based on its format
  char ending[PETSC_MAX_PATH_LEN];
  if (rdy->config.output_format == OUTPUT_BINARY) {  // PETSc native binary format
    snprintf(ending, PETSC_MAX_PATH_LEN - 1, "_dt_%f_%d_%d_np%d.%s", rdy->dt, rdy->config.max_step, step, rdy->nproc, suffix);
  } else {
    snprintf(ending, PETSC_MAX_PATH_LEN - 1, "_dt_%f_%d_np%d.%s", rdy->dt, rdy->config.max_step, rdy->nproc, suffix);
  }

  // concatenate some config parameters
  strncat(filename, ending, PETSC_MAX_PATH_LEN - 1 - strlen(filename));
  RDyLogDetail(rdy, "Step %d: writing output to %s", step, filename);

  PetscFunctionReturn(0);
}

// writes output in Petsc native binary format
static PetscErrorCode WriteBinaryOutput(RDy rdy, PetscInt step, PetscReal time) {
  PetscFunctionBegin;

  // Determine the output file name.
  char fname[PETSC_MAX_PATH_LEN];
  PetscCall(DetermineOutputFile(rdy, step, time, "dat", fname));

  PetscViewer viewer;
  PetscCall(PetscViewerBinaryOpen(rdy->comm, fname, FILE_MODE_WRITE, &viewer));
  PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_NATIVE));

  // dump the solution vector in natural ordering
  Vec natural;
  PetscCall(DMPlexCreateNaturalVector(rdy->dm, &natural));
  PetscCall(DMPlexGlobalToNaturalBegin(rdy->dm, rdy->X, natural));
  PetscCall(DMPlexGlobalToNaturalEnd(rdy->dm, rdy->X, natural));
  PetscCall(VecView(natural, viewer));
  PetscCall(VecDestroy(&natural));

  PetscCall(PetscViewerPopFormat(viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscFunctionReturn(0);
}

// writes a XDMF "heavy data" to an HDF5 file
static PetscErrorCode WriteXDMFHeavyData(RDy rdy, PetscInt step, PetscReal time) {
  PetscFunctionBegin;

  PetscViewer viewer;

  // Determine the output file name.
  char fname[PETSC_MAX_PATH_LEN];
  PetscCall(DetermineOutputFile(rdy, step, time, "h5", fname));

  // write the grid if we're on the first step
  if (step == 0) {
    PetscCall(PetscViewerHDF5Open(rdy->comm, fname, FILE_MODE_WRITE, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_XDMF));
    PetscCall(DMView(rdy->dm, viewer));
  } else {
    PetscCall(PetscViewerHDF5Open(rdy->comm, fname, FILE_MODE_APPEND, &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_XDMF));
  }

  // write solution data to a new GROUP with components in separate datasets
  PetscReal time_in_days = time / (24 * 3600);  // seconds -> days
  char      groupName[1025];
  snprintf(groupName, 1024, "%d %E d", step, time_in_days);
  PetscCall(PetscViewerHDF5PushGroup(viewer, groupName));

  // create and populate a multi-component natural vector
  Vec natural;
  PetscCall(DMPlexCreateNaturalVector(rdy->dm, &natural));
  PetscCall(DMPlexGlobalToNaturalBegin(rdy->dm, rdy->X, natural));
  PetscCall(DMPlexGlobalToNaturalEnd(rdy->dm, rdy->X, natural));

  // extract each component into a separate vector and write it to the group
  // FIXME: This setup is specific to the shallow water equations. We can
  // FIXME: generalize it later.
  const char *comp_names[3] = {
      "Water_Height",
      "X_Momentum",
      "Y_Momentum",
  };
  Vec      comp;  // single-component natural vector
  PetscInt n, N, bs;
  PetscCall(VecGetLocalSize(natural, &n));
  PetscCall(VecGetSize(natural, &N));
  PetscCall(VecGetBlockSize(natural, &bs));
  PetscCall(VecCreateMPI(rdy->comm, n / bs, N / bs, &comp));
  PetscReal *Xi;  // multi-component natural vector data
  PetscCall(VecGetArray(natural, &Xi));
  for (PetscInt c = 0; c < bs; ++c) {
    PetscObjectSetName((PetscObject)comp, comp_names[c]);
    PetscReal *Xci;  // single-component natural vector data
    PetscCall(VecGetArray(comp, &Xci));
    for (PetscInt i = 0; i < n / bs; ++i) Xci[i] = Xi[bs * i + c];
    PetscCall(VecRestoreArray(comp, &Xci));
    PetscCall(VecView(comp, viewer));
  }
  PetscCall(VecRestoreArray(natural, &Xi));

  // clean up
  PetscCall(VecDestroy(&comp));
  PetscCall(VecDestroy(&natural));
  PetscCall(PetscViewerHDF5PopGroup(viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscFunctionReturn(0);
}

// generates an XMDF "light data" file (.xmf)
static PetscErrorCode WriteXDMFLightData(RDy rdy, PetscInt num_files) {
  PetscFunctionBegin;

  // data (HDF5) paths
  const char *geom_path = "/geometry";
  const char *topo_path = "/viz/topology";

  // mesh metadata
  PetscInt num_cells = rdy->mesh.num_cells;
  PetscInt num_corners;
  PetscInt num_vertices = rdy->mesh.num_vertices;
  PetscInt cell_dim = 2, space_dim = 2;

  // time-stepping metadata
  PetscInt num_times = num_files;

  char h5_name[PETSC_MAX_PATH_LEN], xmf_name[PETSC_MAX_PATH_LEN];
  PetscCall(DetermineOutputFile(rdy, 0, 0.0, "h5", h5_name));
  PetscCall(DetermineOutputFile(rdy, 0, 0.0, "xmf", xmf_name));
  RDyLogDetail(rdy, "Writing XDMF light data to %s...", xmf_name);
  FILE *fp;
  PetscCall(PetscFOpen(rdy->comm, xmf_name, "w", &fp));

  // write header
  PetscCall(PetscFPrintf(rdy->comm, fp,
                         "<?xml version=\"1.0\" ?>\n"
                         "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
                         "<!ENTITY HeavyData \"%s\">\n"
                         "]>\n\n"
                         "<Xdmf>\n  <Domain Name=\"domain\">\n",
                         h5_name));

  // write cell metadata
  PetscCall(PetscFPrintf(rdy->comm, fp,
                         "    <DataItem Name=\"cells\"\n"
                         "              ItemType=\"Uniform\"\n"
                         "              Format=\"HDF\"\n"
                         "              NumberType=\"Float\" Precision=\"8\"\n"
                         "              Dimensions=\"%d %d\">\n"
                         "      &HeavyData;:/%s/cells\n"
                         "    </DataItem>\n",
                         num_cells, num_corners, topo_path));

  // write vertex metadata
  PetscCall(PetscFPrintf(rdy->comm, fp,
                         "    <DataItem Name=\"vertices\"\n"
                         "              Format=\"HDF\"\n"
                         "              Dimensions=\"%d %d\">\n"
                         "      &HeavyData;:/%s/vertices\n"
                         "    </DataItem>\n",
                         num_vertices, space_dim, geom_path));

  // write time grid header
  PetscCall(PetscFPrintf(rdy->comm, fp,
                         "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
                         "      <Time TimeType=\"List\">\n"
                         "        <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"%d\">\n"
                         "          ",
                         num_times));
  for (int n = 0; n < num_times; ++n) {
    PetscReal tn = n * rdy->dt;
    PetscCall(PetscFPrintf(rdy->comm, fp, "%g ", tn));
  }
  PetscCall(PetscFPrintf(rdy->comm, fp,
                         "        </DataItem>\n"
                         "      </Time>\n"));

  // write field data for each time
  for (int n = 0; n < num_times; ++n) {
    // write space grid header
    const char *topo_type[5] = {"Invalid", "Invalid", "Invalid", "Triangle", "Quadrilateral"};
    PetscCall(PetscFPrintf(rdy->comm, fp,
                           "      <Grid Name=\"domain\" GridType=\"Uniform\">\n"
                           "        <Topology\n"
                           "           TopologyType=\"%s\"\n"
                           "           NumberOfElements=\"%d\">\n"
                           "          <DataItem Reference=\"XML\">\n"
                           "            /Xdmf/Domain/DataItem[@Name=\"cells\"]\n"
                           "          </DataItem>\n"
                           "        </Topology>\n"
                           "        <Geometry GeometryType=\"XY\">\n"
                           "          <DataItem Reference=\"XML\">\n"
                           "            /Xdmf/Domain/DataItem[@Name=\"vertices\"]\n"
                           "          </DataItem>\n"
                           "        </Geometry>\n",
                           topo_type[num_corners], num_cells));

    // write vertex field metadata
    // (none so far!)

    // write cell field metadata
    PetscReal   tn                  = n * rdy->dt;
    const char *cell_field_names[3] = {"Water_Height", "X_Momentum", "Y_Momentum"};
    for (int f = 0; f < 3; ++f) {
      PetscCall(PetscFPrintf(rdy->comm, fp,
                             "        <Attribute Name=\"%s\" Center=\"Cell\">\n"
                             "          <DataItem DataType=\"Float\" Precision=\"8\" Dimensions=\"%d\" Format=\"HDF\">\n"
                             "            &HeavyData;:/%d %.7E d/%s\n"
                             "          </DataItem>\n"
                             "        </Attribute>\n",
                             cell_field_names[f], num_cells, n, tn, cell_field_names[f]));

      // write space grid footer
      PetscCall(PetscFPrintf(rdy->comm, fp, "      </Grid>\n"));
    }
  }

  // write time grid footer
  PetscCall(PetscFPrintf(rdy->comm, fp, "    </Grid>\n"));

  // write footer
  PetscCall(PetscFPrintf(rdy->comm, fp, "  </Domain>\n</Xdmf>\n"));

  PetscCall(PetscFClose(rdy->comm, fp));

  PetscFunctionReturn(0);
}

// TS monitoring routine used to write output files
PetscErrorCode WriteOutputFiles(TS ts, PetscInt step, PetscReal time, Vec X, void *ctx) {
  PetscFunctionBegin;
  RDy rdy = ctx;

  PetscReal final_time = ConvertTimeToSeconds(rdy->config.final_time, rdy->config.time_unit);
  if ((step == -1) ||          // last step (interpolated)
      (time >= final_time) ||  // last step without interpolation
      ((step % rdy->config.output_frequency) == 0)) {
    if (rdy->config.output_format == OUTPUT_XDMF) {
      PetscCall(WriteXDMFHeavyData(rdy, step, time));
    } else {  // binary
      PetscCall(WriteBinaryOutput(rdy, step, time));
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PostprocessOutput(RDy rdy) {
  PetscFunctionBegin;

  PetscReal final_time = ConvertTimeToSeconds(rdy->config.final_time, rdy->config.time_unit);
  PetscInt  num_files  = (int)(final_time / rdy->dt);
  if (rdy->config.output_format == OUTPUT_XDMF) {
    PetscCall(WriteXDMFLightData(rdy, num_files));
  }

  PetscFunctionReturn(0);
}
