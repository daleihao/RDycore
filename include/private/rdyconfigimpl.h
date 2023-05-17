#ifndef RDYCONFIG_H
#define RDYCONFIG_H

#include <private/rdylogimpl.h>

// The types in this file ѕerve as an intermediate representation for our input
// configuration file:
//
// https://rdycore.atlassian.net/wiki/spaces/PD/pages/24576001/RDycore+configuration+file
//

// the maximum number of regions that can be defined on a domain
#define MAX_NUM_REGIONS 32

// the maximum number of boundaries that can be defined on a domain
#define MAX_NUM_BOUNDARIES 32

// the maximum number of flow/sediment/salinity conditions that can be defined for a
// simulation
#define MAX_NUM_CONDITIONS 32

// the maximum length of a string referring to a name in the config file
#define MAX_NAME_LEN 128

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
  TIME_SECONDS = 0,
  TIME_MINUTES,
  TIME_HOURS,
  TIME_DAYS,
  TIME_MONTHS,
  TIME_YEARS
} RDyTimeUnit;

// This type identifies a format for output files.
typedef enum {
  OUTPUT_BINARY = 0,
  OUTPUT_XDMF,
  OUTPUT_CGNS
} RDyOutputFormat;

// This type specifies a "kind" of condition that indicates how that condition
// is to be enforced on a region or boundary.
typedef enum {
  CONDITION_DIRICHLET = 0,    // Dirichlet condition (value is specified)
  CONDITION_NEUMANN,          // Neumann condition (derivative is specified)
  CONDITION_REFLECTING,       // Reflecting condition
  CONDITION_CRITICAL_OUTFLOW  // Critical flow
} RDyConditionType;

// This type defines an initial/boundary condition and/or source/sink with named
// flow, sediment, and salinity conditions.
typedef struct {
  char flow_name[MAX_NAME_LEN+1];     // name of related flow condition
  char sediment_name[MAX_NAME_LEN+1]; // name of related sediment condition
  char salinity_name[MAX_NAME_LEN+1]; // name of related salinity condition
} RDyConditionSpec;

// This type defines a flow-related condition.
typedef struct {
  char             name[MAX_NAME_LEN+1];
  RDyConditionType type;
  PetscReal        height;
  PetscReal        momentum[2];
} RDyFlowCondition;

// This type defines a sediment-related condition.
typedef struct {
  char             name[MAX_NAME_LEN+1];
  RDyConditionType type;
  PetscReal        concentration;
} RDySedimentCondition;

// This type defines a salinity-related condition.
typedef struct {
  char             name[MAX_NAME_LEN+1];
  RDyConditionType type;
  PetscReal        concentration;
} RDySalinityCondition;

// This type is a representation of the config file itself.
typedef struct {

  //---------
  // Physics
  //---------

  // flow and related parameterization(s)
  RDyFlowMode    flow_mode;
  PetscReal      tiny_h; // depth below which no flow occurs (hardwired)
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

  // initial conditions file (if given)
  char              initial_conditions_file[PETSC_MAX_PATH_LEN];
  PetscViewerFormat initial_conditions_format;

  // IDs of all regions mentioned in an input file
  PetscInt num_regions;
  PetscInt region_ids[MAX_NUM_REGIONS];

  // IDs of all boundaries mentioned in an input file
  PetscInt num_boundaries;
  PetscInt boundary_ids[MAX_NUM_BOUNDARIES];

  // table of initial conditions mapping regions to names of conditions
  // (unless file above is given)
  PetscInt         num_initial_conditions;
  RDyConditionSpec initial_conditions[MAX_NUM_REGIONS];

  // table of sources/sinks mapping regions to names of conditions
  PetscInt         num_sources;
  RDyConditionSpec sources[MAX_NUM_REGIONS];

  // table of boundary conditions mapping boundaries to names of conditions
  PetscInt         num_boundary_conditions;
  RDyConditionSpec boundary_conditions[MAX_NUM_BOUNDARIES];

  // tables of named sets of flow, sediment, and salinity conditions
  PetscInt             num_flow_conditions;
  RDyFlowCondition     flow_conditions[MAX_NUM_CONDITIONS];
  PetscInt             num_sediment_conditions;
  RDySedimentCondition sediment_conditions[MAX_NUM_CONDITIONS];
  PetscInt             num_salinity_conditions;
  RDySalinityCondition salinity_conditions[MAX_NUM_CONDITIONS];

  //--------------
  // Timestepping
  //--------------

  // simulation time at which to end
  PetscReal final_time;
  // Units in which final_time is expressed
  RDyTimeUnit time_unit;
  // Maximum number of time steps
  PetscInt max_step;
  // Time step
  PetscReal dtime;

  //----------
  // Restarts
  //----------

  // Restart file format
  PetscViewerFormat restart_format;
  // Restart frequency (in steps)
  PetscInt restart_frequency;

  //--------
  // Output
  //--------

  // Output file format
  RDyOutputFormat output_format;
  // Output frequency (in steps)
  PetscInt output_frequency;

  //---------
  // Logging
  //---------

  // log filename
  char        log_file[PETSC_MAX_PATH_LEN];
  // selected log level
  RDyLogLevel log_level;

} RDyConfig;

PETSC_INTERN PetscErrorCode RDyConfigFindRegion(RDyConfig*, PetscInt, PetscInt*);
PETSC_INTERN PetscErrorCode RDyConfigFindBoundary(RDyConfig*, PetscInt, PetscInt*);


#endif

