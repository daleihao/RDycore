#ifndef RDYCOREIMPL_H
#define RDYCOREIMPL_H

#include <rdycore.h>
#include <private/rdyconfig.h>
#include <private/rdylogimpl.h>
#include <private/rdymeshimpl.h>

#include <petsc/private/petscimpl.h>

// This type identifies available spatial discretization methods.
typedef enum {
  SPATIAL_FV = 0, // finite volume method
  SPATIAL_FE      // finite element method
} RDySpatial;

// This type identifies available temporal discretization methods.
typedef enum {
  TEMPORAL_EULER = 0, // forward euler method
  TEMPORAL_RK4,       // 4th-order Runge-Kutta method
  TEMPORAL_BEULER     // backward euler method
} RDyTemporal;

// This type identifies available Riemann solvers for horizontal flow.
typedef enum {
  RIEMANN_ROE = 0, // Roe solver
  RIEMANN_HLLC     // Harten, Lax, van Leer Contact solver
} RDyRiemann;

// This type identifies available physics flow modes.
typedef enum {
  FLOW_SWE = 0,
  FLOW_DIFFUSION
} RDyFlowMode;

// This type identifies available bed friction models.
typedef enum {
  BED_FRICTION_NOT_SET = 0,
  BED_FRICTION_NONE,
  BED_FRICTION_CHEZY,
  BED_FRICTION_MANNING
} RDyBedFriction;

// This type identifies a time unit.
typedef enum {
  TIME_MINUTES = 0,
  TIME_HOURS,
  TIME_DAYS,
  TIME_MONTHS,
  TIME_YEARS
} RDyTimeUnit;

// This type specifies a "kind" of condition that indicates how that condition
// is to be enforced on a region or surface.
typedef enum {
  CONDITION_DIRICHLET = 0, // Dirichlet condition (value is specified)
  CONDITION_NEUMANN,       // Neumann condition (derivative is specified)
  CONDITION_REFLECTING     // Reflecting condition
} RDyConditionType;

// This type defines a flow-related condition.
typedef struct {
  const char       *name;
  RDyConditionType  type;
  PetscReal         height;
  PetscReal         momentum[2];
} RDyFlowCondition;

// This type defines a sediment-related condition.
typedef struct {
  const char       *name;
  RDyConditionType  type;
  PetscReal         concentration;
} RDySedimentCondition;

// This type defines a salinity-related condition.
typedef struct {
  const char       *name;
  RDyConditionType  type;
  PetscReal         concentration;
} RDySalinityCondition;

// This type defines a "condition" representing
// * an initial condition or source/sink associated with a region
// * a boundary condition associated with a surface
typedef struct {
  // flow, sediment, salinity conditions (NULL for none)
  RDyFlowCondition     *flow;
  RDySedimentCondition *sediment;
  RDySalinityCondition *salinity;

  // value(s) associated with the condition
  PetscReal value;
} RDyCondition;

// This type defines a region consisting of cells identified by their local
// indices.
typedef struct {
  PetscInt   *cell_ids;
  PetscInt    num_cells;
} RDyRegion;

// This type defines a surface consisting of edges identified by their local
// indices.
typedef struct {
  PetscInt   *edge_ids;
  PetscInt    num_edges;
} RDySurface;

// This type serves as a "virtual table" containing function pointers that
// define the behavior of the dycore.
typedef struct _RDyOps *RDyOps;
struct _RDyOps {
  // Called by RDyCreate to allocate implementation-specific resources, storing
  // the result in the given context pointer.
  PetscErrorCode (*create)(void**);
  // Called by RDyDestroy to free implementation-specific resources.
  PetscErrorCode (*destroy)(void*);
};

// an application context that stores data relevant to a simulation
struct _p_RDy {
  PETSCHEADER(struct _RDyOps);

  // Implementation-specific context pointer
  void *context;

  //------------------------------------------------------------------
  // TODO: The fields below are subject to change as we find our way!
  //------------------------------------------------------------------

  // MPI communicator used for the simulation
  MPI_Comm comm;
  // MPI rank of local process
  PetscInt rank;
  // Number of processes in the communicator
  PetscInt nproc;
  // file storing input data for the simulation
  char config_file[PETSC_MAX_PATH_LEN];

  //---------
  // Physics
  //---------

  // flow and related parameterization(s)
  RDyFlowMode    flow_mode;
  RDyBedFriction bed_friction;
  PetscReal      bed_friction_coef;

  // sediment
  PetscBool sediment;

  // salinity
  PetscBool salinity;


  //----------
  // Numerics
  //----------
  RDySpatial spatial;
  RDyTemporal temporal;
  RDyRiemann riemann;

  //------------------------
  // Spatial discretization
  //------------------------

  // mesh file
  char mesh_file[PETSC_MAX_PATH_LEN];

  // PETSc (DMPlex) grid
  DM dm;

  // mesh representing simulation domain
  RDyMesh mesh;

  // mesh regions
  PetscInt  num_regions;
  PetscInt  region_ids[1+MAX_REGION_ID];
  RDyRegion regions[1+MAX_REGION_ID];

  // mesh surfaces
  PetscInt   num_surfaces;
  PetscInt   surface_ids[1+MAX_SURFACE_ID];
  RDySurface surfaces[1+MAX_SURFACE_ID];

  // Table of named sets of flow, sediment, and salinity conditions. Used by
  // initial/source/boundary conditions below.
  PetscInt             num_flow_conditions;
  RDyFlowCondition     flow_conditions[1+MAX_CONDITION_ID];
  PetscInt             num_sediment_conditions;
  RDySedimentCondition sediment_conditions[1+MAX_CONDITION_ID];
  PetscInt             num_salinity_conditions;
  RDySalinityCondition salinity_conditions[1+MAX_CONDITION_ID];

  // initial conditions (either a filename or a set of conditions associated
  // with mesh regions (1 per region)
  char         initial_conditions_file[PETSC_MAX_PATH_LEN];
  RDyCondition initial_conditions[1+MAX_REGION_ID];

  // sources (and sinks) associated with mesh regions (1 per region)
  RDyCondition sources[1+MAX_REGION_ID];

  // boundary conditions associated with mesh surfaces (1 per surface)
  RDyCondition boundary_conditions[1+MAX_SURFACE_ID];

  //--------------
  // Timestepping
  //--------------

  // simulation time at which to end
  PetscReal final_time;
  // Units in which final_time is expressed
  RDyTimeUnit time_unit;
  // Maximum number of time steps
  PetscInt max_step;

  //----------
  // Restarts
  //----------

  // Restart file format
  char restart_format[4];
  // Restart frequency (in steps)
  PetscInt restart_frequency;

  //---------
  // Logging
  //---------

  // log filename
  char        log_file[PETSC_MAX_PATH_LEN];
  // selected log level
  RDyLogLevel log_level;
  // log file handle
  FILE       *log;

  //-----------------
  // Simulation data
  //-----------------

  // time step size
  PetscReal dt;
  // index of current timestep
  PetscInt tstep;

  PetscInt  dof;
  Vec       B, localB;
  Vec       localX;
  PetscBool debug, save, add_building;
  PetscBool interpolate;
};

#endif
