#include <cyaml/cyaml.h>
#include <float.h>
#include <muParserDLL.h>
#include <petscdmplex.h>
#include <private/rdycoreimpl.h>

// =============
//  YAML Parser
// =============
//
// The design of the options database in PETSc does not currently allow for one
// to traverse items in YAML mappings--one must know the exact name of every
// item one wants to retrieve. This puts a large restriction on the types of
// YAML files that can be parsed, so we've rolled our own.
//
// Here we've implemented a parser (using libcyaml) to handle the RDycore
// configuration file whose specification can be found at
//
// https://rdycore.github.io/RDycore/user/input.html
//
// In the style of libcyaml (https://github.com/tlsa/libcyaml/blob/main/docs/guide.md),
// this parser defines a schema for each section and populates the appropriate struct(s)
// within RDyConfig (include/private/rdyconfigimpl.h) accordingly. The schema for each
// section appears below and must remain consistent with the data structures in rdyconfigimpl.h.

// ====================
//  Schema definitions
// ====================

// Below, you'll see a bunch of schemas describing YAML sections and subsections
// with fields defined using macros like CYAML_FIELD_ENUM, CYAML_FIELD_INT,
// CYAML_FIELD_FLOAT, etc.
//
// clang-format off

// ---------------
// physics section
// ---------------
// physics:
//   flow:
//     mode: <swe|diffusion> # swe by default
//     tiny_h: <value> # 1e-7 by default
//   sediment: <true|false> # off by default
//   salinity: <true|false> # off by default

// mapping of strings to physics flow types
static const cyaml_strval_t physics_flow_modes[] = {
    {"swe",       FLOW_SWE      },
    {"diffusion", FLOW_DIFFUSION},
};

// mapping of physics.flow fields to members of RDyPhysicsFlow
static const cyaml_schema_field_t physics_flow_fields_schema[] = {
    CYAML_FIELD_ENUM("mode", CYAML_FLAG_DEFAULT, RDyPhysicsFlow, mode, physics_flow_modes, CYAML_ARRAY_LEN(physics_flow_modes)),
    CYAML_FIELD_FLOAT("tiny_h", CYAML_FLAG_OPTIONAL, RDyPhysicsFlow, tiny_h),
    CYAML_FIELD_END
};

// mapping of physics fields to members of RDyPhysicsSection
static const cyaml_schema_field_t physics_fields_schema[] = {
    CYAML_FIELD_MAPPING("flow", CYAML_FLAG_DEFAULT, RDyPhysicsSection, flow, physics_flow_fields_schema),
    CYAML_FIELD_BOOL("sediment", CYAML_FLAG_OPTIONAL, RDyPhysicsSection, sediment),
    CYAML_FIELD_BOOL("salinity", CYAML_FLAG_OPTIONAL, RDyPhysicsSection, salinity),
    CYAML_FIELD_END
};

// ----------------
// numerics section
// ----------------
// numerics:
//   spatial: <fv|fe>                 # To begin with, we will only have FV method
//   temporal: <euler,rk4,beuler,...> # These should map on PETSc -ts_type
//   riemann: <roe|hllc>              # To begin with, we will only have Roe

// mapping of strings to numerics spatial types
static const cyaml_strval_t numerics_spatial_types[] = {
    {"fv", SPATIAL_FV},
    {"fe", SPATIAL_FE},
};

// mapping of strings to numerics temporal types
static const cyaml_strval_t numerics_temporal_types[] = {
    {"euler",  TEMPORAL_EULER },
    {"rk4",    TEMPORAL_RK4   },
    {"beuler", TEMPORAL_BEULER},
};

// mapping of strings to numerics riemann solver types
static const cyaml_strval_t numerics_riemann_types[] = {
    {"roe",  RIEMANN_ROE },
    {"hllc", RIEMANN_HLLC},
};

// mapping of numerics fields to members of RDyNumericsSection
static const cyaml_schema_field_t numerics_fields_schema[] = {
    CYAML_FIELD_ENUM("spatial", CYAML_FLAG_DEFAULT, RDyNumericsSection, spatial, numerics_spatial_types, CYAML_ARRAY_LEN(numerics_spatial_types)),
    CYAML_FIELD_ENUM("temporal", CYAML_FLAG_DEFAULT, RDyNumericsSection, temporal, numerics_temporal_types, CYAML_ARRAY_LEN(numerics_temporal_types)),
    CYAML_FIELD_ENUM("riemann", CYAML_FLAG_DEFAULT, RDyNumericsSection, riemann, numerics_riemann_types, CYAML_ARRAY_LEN(numerics_riemann_types)),
    CYAML_FIELD_END
};

// ------------
// time section
// ------------
// time:
//   final_time: <value>
//   unit: <seconds|minutes|hours|days|months|years> # applies to final_time (stored internally in seconds)
//   max_step: <value>
//   time_step: <value> # optional
//   coupling_interval: <value> # optional

// mapping of strings to time units
static const cyaml_strval_t time_units[] = {
    {"seconds", RDY_TIME_SECONDS},
    {"minutes", RDY_TIME_MINUTES},
    {"hours",   RDY_TIME_HOURS  },
    {"days",    RDY_TIME_DAYS   },
    {"months",  RDY_TIME_MONTHS },
    {"years",   RDY_TIME_YEARS  },
};

static const cyaml_strval_t enable_units[] = {
    {"true", PETSC_TRUE},
    {"false", PETSC_FALSE},
};

static const cyaml_schema_field_t adaptive_fields_schema[] = {
    CYAML_FIELD_ENUM("enable", CYAML_FLAG_DEFAULT, RDyTimeAdaptiveSection, enable, enable_units, CYAML_ARRAY_LEN(enable_units)),
    CYAML_FIELD_FLOAT("target_courant_number", CYAML_FLAG_DEFAULT, RDyTimeAdaptiveSection, target_courant_number),
    CYAML_FIELD_FLOAT("max_increase_factor", CYAML_FLAG_DEFAULT, RDyTimeAdaptiveSection, max_increase_factor),
    CYAML_FIELD_FLOAT("initial_time_step", CYAML_FLAG_DEFAULT, RDyTimeAdaptiveSection, initial_time_step),
    CYAML_FIELD_END
};


// mapping of time fields to members of RDyTimeSection
static const cyaml_schema_field_t time_fields_schema[] = {
    CYAML_FIELD_FLOAT("final_time", CYAML_FLAG_OPTIONAL, RDyTimeSection, final_time),
    CYAML_FIELD_ENUM("unit", CYAML_FLAG_DEFAULT, RDyTimeSection, unit, time_units, CYAML_ARRAY_LEN(time_units)),
    CYAML_FIELD_INT("max_step", CYAML_FLAG_OPTIONAL, RDyTimeSection, max_step),
    CYAML_FIELD_FLOAT("time_step", CYAML_FLAG_OPTIONAL, RDyTimeSection, time_step),
    CYAML_FIELD_FLOAT( "coupling_interval", CYAML_FLAG_OPTIONAL, RDyTimeSection, coupling_interval),
    CYAML_FIELD_MAPPING("adaptive", CYAML_FLAG_OPTIONAL, RDyTimeSection, adaptive, adaptive_fields_schema),
    CYAML_FIELD_END
};

// ---------------
// logging section
// ---------------
// logging:
//   file: <path> # default: stdout
//   level: <none|warning|info|detail|debug> # <-- increasing levels of logging (default: none)

// mapping of strings to log levels
static const cyaml_strval_t logging_levels[] = {
    {"none",    LOG_NONE   },
    {"warning", LOG_WARNING},
    {"info",    LOG_INFO   },
    {"detail",  LOG_DETAIL },
    {"debug",   LOG_DEBUG  },
};

// mapping of logging fields to members of RDyLoggingSection
static const cyaml_schema_field_t logging_fields_schema[] = {
    CYAML_FIELD_STRING("file", CYAML_FLAG_OPTIONAL, RDyLoggingSection, file, 0),
    CYAML_FIELD_ENUM("level", CYAML_FLAG_OPTIONAL, RDyLoggingSection, level, logging_levels, CYAML_ARRAY_LEN(logging_levels)),
    CYAML_FIELD_END
};

// ------------------
// checkpoint section
// ------------------
// checkpoint:
//   format: <binary|hdf5>
//   interval: <value-in-steps>        # default: 0 (no checkpoints)
//   prefix: <checkpoint-file-prefix>  # default: same as input file prefix

// mapping of strings to file formats
static const cyaml_strval_t checkpoint_file_formats[] = {
    {"binary", PETSC_VIEWER_NATIVE    },
    {"hdf5",   PETSC_VIEWER_HDF5_PETSC},
};

// mapping of checkpoint fields to members of RDyCheckpointSection
static const cyaml_schema_field_t checkpoint_fields_schema[] = {
    CYAML_FIELD_ENUM("format", CYAML_FLAG_DEFAULT, RDyCheckpointSection, format, checkpoint_file_formats, CYAML_ARRAY_LEN(checkpoint_file_formats)),
    CYAML_FIELD_INT("interval", CYAML_FLAG_DEFAULT, RDyCheckpointSection, interval),
    CYAML_FIELD_STRING("prefix", CYAML_FLAG_OPTIONAL, RDyCheckpointSection, prefix, 0),
    CYAML_FIELD_END
};

// ---------------
// restart section
// ---------------
// restart:
//   file: <checkpoint-filename>
//   reinitialize: <true/false>  # default: false

// mapping of restart fields to members of RDyRestartSection
static const cyaml_schema_field_t restart_fields_schema[] = {
    CYAML_FIELD_STRING("file", CYAML_FLAG_OPTIONAL, RDyRestartSection, file, 0),
    CYAML_FIELD_BOOL("reinitialize", CYAML_FLAG_OPTIONAL, RDyRestartSection, reinitialize),
    CYAML_FIELD_END
};

// ---------------
// output section
// ---------------
// output:
//   format: <binary|xdmf|cgns>
//   interval: <number-of-steps-between-output-dumps> # default: 0 (no output)
//   batch_size: <number-of-steps-stored-in-each-output-file> # default: 1

// mapping of strings to file formats
static const cyaml_strval_t output_file_formats[] = {
    {"none",   OUTPUT_NONE},
    {"binary", OUTPUT_BINARY},
    {"xdmf",   OUTPUT_XDMF  },
    {"cgns",   OUTPUT_CGNS  },
};

// mapping of time_series fields to members of RDyTimeSeries
static const cyaml_schema_field_t output_time_series_fields_schema[] = {
    CYAML_FIELD_INT("boundary_fluxes", CYAML_FLAG_DEFAULT, RDyTimeSeries, boundary_fluxes),
    CYAML_FIELD_END
};

// mapping of output fields to members of RDyOutputSection
static const cyaml_schema_field_t output_fields_schema[] = {
    CYAML_FIELD_ENUM("format", CYAML_FLAG_OPTIONAL, RDyOutputSection, format, output_file_formats, CYAML_ARRAY_LEN(output_file_formats)),
    CYAML_FIELD_INT("step_interval", CYAML_FLAG_OPTIONAL, RDyOutputSection, step_interval),
    CYAML_FIELD_FLOAT("time_interval", CYAML_FLAG_OPTIONAL, RDyOutputSection, time_interval),
    CYAML_FIELD_ENUM("time_unit", CYAML_FLAG_OPTIONAL, RDyOutputSection, time_unit, time_units, CYAML_ARRAY_LEN(time_units)),
    CYAML_FIELD_INT("batch_size", CYAML_FLAG_OPTIONAL, RDyOutputSection, batch_size),
    CYAML_FIELD_MAPPING("time_series", CYAML_FLAG_OPTIONAL, RDyOutputSection, time_series, output_time_series_fields_schema),
    CYAML_FIELD_END
};

// ------------
// grid section
// ------------
// grid:
//   file: <path-to-file/grid.{msh,h5,exo}>

// mapping of grid fields to members of RDyGridSection
static const cyaml_schema_field_t grid_fields_schema[] = {
  CYAML_FIELD_STRING("file", CYAML_FLAG_DEFAULT, RDyGridSection, file, 1),
  CYAML_FIELD_END
};

// mapping of strings to input file formats
static const cyaml_strval_t input_file_formats[] = {
    {"",       PETSC_VIEWER_NOFORMAT  },
    {"binary", PETSC_VIEWER_NATIVE    },
    {"hdf5",   PETSC_VIEWER_HDF5_PETSC},
};

// ---------------------------
// surface_composition section
// ---------------------------

// mapping of material specification fields to members of RDyMaterialSpec
static const cyaml_schema_field_t surface_composition_fields_schema[] = {
    CYAML_FIELD_STRING("region", CYAML_FLAG_DEFAULT, RDySurfaceCompositionSpec, region, 1),
    CYAML_FIELD_STRING("material", CYAML_FLAG_DEFAULT, RDySurfaceCompositionSpec, material, 1),
    CYAML_FIELD_END
};

// a single surface composition entry
static const cyaml_schema_value_t surface_composition_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDySurfaceCompositionSpec, surface_composition_fields_schema),
};

// -----------------
// materials section
// -----------------

// mapping of material property fields to RDyMaterialPropertySpec
static const cyaml_schema_field_t material_property_fields_schema[] = {
    CYAML_FIELD_STRING("value", CYAML_FLAG_OPTIONAL, RDyMaterialPropertySpec, expression, 1),
    CYAML_FIELD_STRING("file", CYAML_FLAG_OPTIONAL, RDyMaterialPropertySpec, file, 1),
    CYAML_FIELD_ENUM("format", CYAML_FLAG_OPTIONAL, RDyMaterialPropertySpec, format, input_file_formats, CYAML_ARRAY_LEN(input_file_formats)),
    CYAML_FIELD_END
};

// mapping of material property fields to RDyMaterialPropertіesSpec
static const cyaml_schema_field_t material_properties_fields_schema[] = {
    CYAML_FIELD_MAPPING("manning", CYAML_FLAG_DEFAULT, RDyMaterialPropertiesSpec, manning, material_property_fields_schema),
    CYAML_FIELD_END
};

// mapping of material fields to RDyMaterialSpec
static const cyaml_schema_field_t material_fields_schema[] = {
    CYAML_FIELD_STRING("name", CYAML_FLAG_DEFAULT, RDyMaterialSpec, name, 1),
    CYAML_FIELD_MAPPING("properties", CYAML_FLAG_DEFAULT, RDyMaterialSpec, properties, material_properties_fields_schema),
    CYAML_FIELD_END
};

// a single material entry
static const cyaml_schema_value_t material_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDyMaterialSpec, material_fields_schema),
};

// ---------------
// regions section
// ---------------
// regions:
//  - name: downstream   # human-readable name for the region
//    grid_region_id: 2  # grid identifier for the region
//  - name: upstream
//    grid_region_id: 1

// schema for region fields
static const cyaml_schema_field_t region_spec_fields_schema[] = {
    CYAML_FIELD_STRING("name", CYAML_FLAG_DEFAULT, RDyRegionSpec, name, 1),
    CYAML_FIELD_INT("grid_region_id", CYAML_FLAG_DEFAULT, RDyRegionSpec, grid_region_id),
    CYAML_FIELD_END
};

// a single region entry
static const cyaml_schema_value_t region_spec_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDyRegionSpec, region_spec_fields_schema),
};

// ---------------------------------------
// initial_conditions and sources sections
// ---------------------------------------
// initial_conditions/sources:
//  - region: <region-name>
//    flow: <name-of-a-flow-condition>
//    sediment: <name-of-a-sediment-condition> # used if physics.sediment == true above
//    salinity: <name-of-a-salinity-condition> # used if physics.salinity == true above
//  - region: <region-name>
//    flow: <name-of-a-flow-condition>
//    sediment: <name-of-a-sediment-condition> # used only if physics.sediment == true above
//    salinity: <name-of-a-salinity-condition> # used only if physics.salinity == true above
// ...

// mapping of conditions fields to members of RDyRegionConditionSpec
static const cyaml_schema_field_t region_condition_spec_fields_schema[] = {
    CYAML_FIELD_STRING("region", CYAML_FLAG_DEFAULT, RDyRegionConditionSpec, region, 1),
    CYAML_FIELD_STRING("flow", CYAML_FLAG_DEFAULT, RDyRegionConditionSpec, flow, 1),
    CYAML_FIELD_STRING("sediment", CYAML_FLAG_OPTIONAL, RDyRegionConditionSpec, sediment, 0),
    CYAML_FIELD_STRING("salinity", CYAML_FLAG_OPTIONAL, RDyRegionConditionSpec, salinity, 0),
    CYAML_FIELD_END
};

// a single regional initial condition / source spec entry
static const cyaml_schema_value_t region_condition_spec_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDyRegionConditionSpec, region_condition_spec_fields_schema),
};

// ------------------
// boundaries section
// ------------------
// boundaries:
//   - name: bottom_wall
//     grid_boundary_id: 3  # grid identifier for the boundary
//   - name: top_wall
//     grid_boundary_id: 2
//   - name: exterior
//     grid_boundary_id: 1

// schema for boundary fields
static const cyaml_schema_field_t boundary_spec_fields_schema[] = {
    CYAML_FIELD_STRING("name", CYAML_FLAG_DEFAULT, RDyBoundarySpec, name, 1),
    CYAML_FIELD_INT("grid_boundary_id", CYAML_FLAG_DEFAULT, RDyBoundarySpec, grid_boundary_id),
    CYAML_FIELD_END
};

// a single boundary entry
static const cyaml_schema_value_t boundary_spec_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDyBoundarySpec, boundary_spec_fields_schema),
};

// ---------------------------
// boundary_conditions section
// ---------------------------
// boundary_conditions:
//   - boundaries: [<boundary-name1>, <boundary-name2>, ...]
//     flow: <name-of-a-flow-condition>
//     sediment: <name-of-a-sediment-condition> # used if physics.sediment = true above
//     salinity: <name-of-a-salinity-condition> # used if physics.salinity = true above
//   - boundaries: [<boundary-name1>, <boundary-name2>, ...]
//     flow: <name-of-a-flow-condition>
//     sediment: <name-of-a-sediment-condition> # used only if physics.sediment = true above
//     salinity: <name-of-a-salinity-condition> # used only if physics.salinity = true above
// ...

// schema for boundary name
static const cyaml_schema_value_t boundary_name_entry = {
    CYAML_VALUE_STRING(CYAML_FLAG_DEFAULT, char, 1, MAX_NAME_LEN),
};

// mapping of conditions fields to members of RDyBoundaryConditionSpec
static const cyaml_schema_field_t boundary_condition_spec_fields_schema[] = {
    CYAML_FIELD_SEQUENCE_COUNT("boundaries", CYAML_FLAG_DEFAULT, RDyBoundaryConditionSpec, boundaries, num_boundaries, &boundary_name_entry, 0, MAX_NUM_BOUNDARIES),
    CYAML_FIELD_STRING("flow", CYAML_FLAG_DEFAULT, RDyBoundaryConditionSpec, flow, 1),
    CYAML_FIELD_STRING("sediment", CYAML_FLAG_OPTIONAL, RDyBoundaryConditionSpec, sediment, 0),
    CYAML_FIELD_STRING("salinity", CYAML_FLAG_OPTIONAL, RDyBoundaryConditionSpec, salinity, 0),
    CYAML_FIELD_END
};

// a single boundary condition spec entry
static const cyaml_schema_value_t boundary_condition_spec_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDyBoundaryConditionSpec, boundary_condition_spec_fields_schema),
};

// -----------------------
// flow_conditions section
// -----------------------
// - name: <name-of-flow-condition-1>
//   type: <dirichlet|neumann|reflecting|critical-outflow>
//   height: <h>  # used only by initial conditions + dirichlet bcs
//   x_momentum: <px> # used only by initial conditions + dirichlet bcs
//   y_momentum: <py> # used only by initial conditions + dirichlet bcs
// - name: <name-of-flow-condition-2>
//   type: <dirichlet|neumann|reflecting|critical-outflow>
//   file: <filename>      # used only by initial conditions + dirichlet bcs
//   format: <binary|hdf5> # used only by initial conditions + dirichlet bcs
//   ...

// mapping of strings to types of conditions
static const cyaml_strval_t condition_types[] = {
    {"dirichlet",        CONDITION_DIRICHLET       },
    {"neumann",          CONDITION_NEUMANN         },
    {"reflecting",       CONDITION_REFLECTING      },
    {"critical-outflow", CONDITION_CRITICAL_OUTFLOW},
};

// schema for flow condition fields
static const cyaml_schema_field_t flow_condition_fields_schema[] = {
    CYAML_FIELD_STRING("name", CYAML_FLAG_DEFAULT, RDyFlowCondition, name, 1),
    CYAML_FIELD_ENUM("type", CYAML_FLAG_DEFAULT, RDyFlowCondition, type, condition_types, CYAML_ARRAY_LEN(condition_types)),
    CYAML_FIELD_STRING("height", CYAML_FLAG_OPTIONAL, RDyFlowCondition, height_expression, 1),
    CYAML_FIELD_STRING("x_momentum", CYAML_FLAG_OPTIONAL, RDyFlowCondition, x_momentum_expression, 1),
    CYAML_FIELD_STRING("y_momentum", CYAML_FLAG_OPTIONAL, RDyFlowCondition, y_momentum_expression, 1),
    CYAML_FIELD_STRING("file", CYAML_FLAG_OPTIONAL, RDyFlowCondition, file, 1),
    CYAML_FIELD_ENUM("format", CYAML_FLAG_OPTIONAL, RDyFlowCondition, format, input_file_formats, CYAML_ARRAY_LEN(input_file_formats)),
    CYAML_FIELD_END
};

// a single flow_conditions entry
static const cyaml_schema_value_t flow_condition_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDyFlowCondition, flow_condition_fields_schema),
};

// ---------------------------
// sediment_conditions section
// ---------------------------
// - name: <name-of-sediment-condition-1>
//   type: <dirichlet|neumann|reflecting|critical>
//   concentration: <value> # used only by dirichlet
// - name: <name-of-sediment-condition-2>
//   type: <dirichlet|neumann|reflecting|critical>
//   file: <filename>      # used only by dirichlet
//   format: <binary|hdf5> # used only by dirichlet
//   ...

// schema for sediment_condition fields
static const cyaml_schema_field_t sediment_condition_fields_schema[] = {
    CYAML_FIELD_STRING("name", CYAML_FLAG_DEFAULT, RDySedimentCondition, name, 1),
    CYAML_FIELD_ENUM("type", CYAML_FLAG_DEFAULT, RDySedimentCondition, type, condition_types, CYAML_ARRAY_LEN(condition_types)),
    CYAML_FIELD_STRING("concentration", CYAML_FLAG_OPTIONAL, RDySedimentCondition, expression, 1),
    CYAML_FIELD_STRING("file", CYAML_FLAG_OPTIONAL, RDySedimentCondition, file, 1),
    CYAML_FIELD_ENUM("format", CYAML_FLAG_OPTIONAL, RDySedimentCondition, format, input_file_formats, CYAML_ARRAY_LEN(input_file_formats)),
    CYAML_FIELD_END
};

// a single sediment_conditions entry
static const cyaml_schema_value_t sediment_condition_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDySedimentCondition, sediment_condition_fields_schema),
};

// ---------------------------
// salinity_conditions section
// ---------------------------
// - name: <name-of-salinity-condition-1>
//   type: <dirichlet|neumann|reflecting|critical>
//   concentration: <value> # used only by dirichlet
// - name: <name-of-salinity-condition-2>
//   type: <dirichlet|neumann|reflecting|critical>
//   file: <filename>      # used only by dirichlet
//   format: <binary|hdf5> # used only by dirichlet
//   ...

// schema for salinity fields
static const cyaml_schema_field_t salinity_condition_fields_schema[] = {
    CYAML_FIELD_STRING("name", CYAML_FLAG_DEFAULT, RDySalinityCondition, name, 1),
    CYAML_FIELD_ENUM("type", CYAML_FLAG_DEFAULT, RDySalinityCondition, type, condition_types, CYAML_ARRAY_LEN(condition_types)),
    CYAML_FIELD_STRING("concentration", CYAML_FLAG_OPTIONAL, RDySedimentCondition, expression, 1),
    CYAML_FIELD_STRING("file", CYAML_FLAG_OPTIONAL, RDySalinityCondition, file, 1),
    CYAML_FIELD_ENUM("format", CYAML_FLAG_OPTIONAL, RDySalinityCondition, format, input_file_formats, CYAML_ARRAY_LEN(input_file_formats)),
    CYAML_FIELD_END
};

// a single salinity_conditions entry
static const cyaml_schema_value_t salinity_condition_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDySalinityCondition, salinity_condition_fields_schema),
};

// ----------------
// ensemble section
// ----------------
//   size: <number-of-ensemble-members>
//   members:
//   - name: <name-of-first-member> # optional
//     <set-of-overridden-sections-for-first-member>
//   - name: <name-of-second-member> # optional
//     <set-of-overridden-sections-for-second-member>
//   ...
//   - name: <name-of-last-member> # optional
//     <set-of-overridden-sections-for-last-member>

// schema for individual ensemble members
static const cyaml_schema_field_t ensemble_member_fields_schema[] = {
    CYAML_FIELD_STRING("name", CYAML_FLAG_OPTIONAL, RDyEnsembleMember, name, 1),
    CYAML_FIELD_MAPPING("grid", CYAML_FLAG_OPTIONAL, RDyEnsembleMember, grid, grid_fields_schema),
    CYAML_FIELD_SEQUENCE("materials", CYAML_FLAG_OPTIONAL | CYAML_FLAG_POINTER, RDyEnsembleMember, materials, &material_entry, 0, MAX_NUM_MATERIALS),
    CYAML_FIELD_SEQUENCE("flow_conditions", CYAML_FLAG_OPTIONAL | CYAML_FLAG_POINTER, RDyEnsembleMember, flow_conditions, &flow_condition_entry, 0,
                               MAX_NUM_CONDITIONS),
    CYAML_FIELD_SEQUENCE("sediment_conditions", CYAML_FLAG_OPTIONAL | CYAML_FLAG_POINTER, RDyEnsembleMember, sediment_conditions,
                               &sediment_condition_entry, 0, MAX_NUM_CONDITIONS),
    CYAML_FIELD_SEQUENCE("salinity_conditions", CYAML_FLAG_OPTIONAL | CYAML_FLAG_POINTER, RDyEnsembleMember, salinity_conditions,
                               &salinity_condition_entry, 0, MAX_NUM_CONDITIONS),
    CYAML_FIELD_END
};

// a single ensemble member entry
static const cyaml_schema_value_t ensemble_member_entry = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_DEFAULT, RDyEnsembleMember, ensemble_member_fields_schema),
};

// ensemble specification
static const cyaml_schema_field_t ensemble_fields_schema[] = {
    CYAML_FIELD_INT("size", CYAML_FLAG_DEFAULT, RDyEnsembleSection, size),
    CYAML_FIELD_SEQUENCE("members", CYAML_FLAG_POINTER, RDyEnsembleSection, members, &ensemble_member_entry, 0, CYAML_UNLIMITED),
    CYAML_FIELD_END
};

//-------------------------------
// mms section (MMS driver only!)
//-------------------------------

static const cyaml_schema_field_t mms_constants_fields_schema[] = {
    CYAML_FIELD_FLOAT("A", CYAML_FLAG_OPTIONAL, RDyMMSConstants, A),
    CYAML_FIELD_FLOAT("B", CYAML_FLAG_OPTIONAL, RDyMMSConstants, B),
    CYAML_FIELD_FLOAT("C", CYAML_FLAG_OPTIONAL, RDyMMSConstants, C),
    CYAML_FIELD_FLOAT("D", CYAML_FLAG_OPTIONAL, RDyMMSConstants, D),
    CYAML_FIELD_FLOAT("E", CYAML_FLAG_OPTIONAL, RDyMMSConstants, E),
    CYAML_FIELD_FLOAT("F", CYAML_FLAG_OPTIONAL, RDyMMSConstants, F),
    CYAML_FIELD_FLOAT("G", CYAML_FLAG_OPTIONAL, RDyMMSConstants, G),
    CYAML_FIELD_FLOAT("H", CYAML_FLAG_OPTIONAL, RDyMMSConstants, H),
    CYAML_FIELD_FLOAT("I", CYAML_FLAG_OPTIONAL, RDyMMSConstants, I_),
    CYAML_FIELD_FLOAT("J", CYAML_FLAG_OPTIONAL, RDyMMSConstants, J),
    CYAML_FIELD_FLOAT("K", CYAML_FLAG_OPTIONAL, RDyMMSConstants, K),
    CYAML_FIELD_FLOAT("L", CYAML_FLAG_OPTIONAL, RDyMMSConstants, L),
    CYAML_FIELD_FLOAT("M", CYAML_FLAG_OPTIONAL, RDyMMSConstants, M),
    CYAML_FIELD_FLOAT("N", CYAML_FLAG_OPTIONAL, RDyMMSConstants, N),
    CYAML_FIELD_FLOAT("O", CYAML_FLAG_OPTIONAL, RDyMMSConstants, O),
    CYAML_FIELD_FLOAT("P", CYAML_FLAG_OPTIONAL, RDyMMSConstants, P),
    CYAML_FIELD_FLOAT("Q", CYAML_FLAG_OPTIONAL, RDyMMSConstants, Q),
    CYAML_FIELD_FLOAT("R", CYAML_FLAG_OPTIONAL, RDyMMSConstants, R),
    CYAML_FIELD_FLOAT("S", CYAML_FLAG_OPTIONAL, RDyMMSConstants, S),
    CYAML_FIELD_FLOAT("T", CYAML_FLAG_OPTIONAL, RDyMMSConstants, T),
    CYAML_FIELD_FLOAT("U", CYAML_FLAG_OPTIONAL, RDyMMSConstants, U),
    CYAML_FIELD_FLOAT("V", CYAML_FLAG_OPTIONAL, RDyMMSConstants, V),
    CYAML_FIELD_FLOAT("W", CYAML_FLAG_OPTIONAL, RDyMMSConstants, W),
    CYAML_FIELD_FLOAT("X", CYAML_FLAG_OPTIONAL, RDyMMSConstants, X),
    CYAML_FIELD_FLOAT("Y", CYAML_FLAG_OPTIONAL, RDyMMSConstants, Y),
    CYAML_FIELD_FLOAT("Z", CYAML_FLAG_OPTIONAL, RDyMMSConstants, Z),
    CYAML_FIELD_END
};

static const cyaml_schema_field_t mms_swe_error_norms_fields_schema[] = {
    CYAML_FIELD_FLOAT("L1", CYAML_FLAG_OPTIONAL, RDyMMSSWEErrorNorms, L1),
    CYAML_FIELD_FLOAT("L2", CYAML_FLAG_OPTIONAL, RDyMMSSWEErrorNorms, L2),
    CYAML_FIELD_FLOAT("Linf", CYAML_FLAG_OPTIONAL, RDyMMSSWEErrorNorms, Linf),
    CYAML_FIELD_END
};

static const cyaml_schema_field_t mms_swe_convergence_rates_fields_schema[] = {
    CYAML_FIELD_MAPPING("h", CYAML_FLAG_OPTIONAL, RDyMMSSWEConvergenceRates, h, mms_swe_error_norms_fields_schema),
    CYAML_FIELD_MAPPING("hu", CYAML_FLAG_OPTIONAL, RDyMMSSWEConvergenceRates, hu, mms_swe_error_norms_fields_schema),
    CYAML_FIELD_MAPPING("hv", CYAML_FLAG_OPTIONAL, RDyMMSSWEConvergenceRates, hv, mms_swe_error_norms_fields_schema),
    CYAML_FIELD_END
};

static const cyaml_schema_field_t mms_swe_convergence_fields_schema[] = {
    CYAML_FIELD_INT("num_refinements", CYAML_FLAG_DEFAULT, RDyMMSSWEConvergence, num_refinements),
    CYAML_FIELD_INT("base_refinement", CYAML_FLAG_OPTIONAL, RDyMMSSWEConvergence, base_refinement),
    CYAML_FIELD_MAPPING("expected_rates", CYAML_FLAG_DEFAULT, RDyMMSSWEConvergence,
                        expected_rates, mms_swe_convergence_rates_fields_schema),
    CYAML_FIELD_END
};

static const cyaml_schema_field_t mms_swe_fields_schema[] = {
    CYAML_FIELD_STRING("h", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.h, 1),
    CYAML_FIELD_STRING("dhdx", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dhdx, 1),
    CYAML_FIELD_STRING("dhdy", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dhdy, 1),
    CYAML_FIELD_STRING("dhdt", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dhdt, 1),
    CYAML_FIELD_STRING("u", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.u, 1),
    CYAML_FIELD_STRING("dudx", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dudx, 1),
    CYAML_FIELD_STRING("dudy", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dudy, 1),
    CYAML_FIELD_STRING("dudt", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dudt, 1),
    CYAML_FIELD_STRING("v", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.v, 1),
    CYAML_FIELD_STRING("dvdx", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dvdx, 1),
    CYAML_FIELD_STRING("dvdy", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dvdy, 1),
    CYAML_FIELD_STRING("dvdt", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dvdt, 1),
    CYAML_FIELD_STRING("z", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.z, 1),
    CYAML_FIELD_STRING("dzdx", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dzdx, 1),
    CYAML_FIELD_STRING("dzdy", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.dzdy, 1),
    CYAML_FIELD_STRING("n", CYAML_FLAG_DEFAULT, RDyMMSSWESolutions, expressions.n, 1),
    CYAML_FIELD_MAPPING("convergence", CYAML_FLAG_OPTIONAL, RDyMMSSWESolutions, convergence, mms_swe_convergence_fields_schema),
    CYAML_FIELD_END
};

static const cyaml_schema_field_t mms_fields_schema[] = {
    CYAML_FIELD_MAPPING("constants", CYAML_FLAG_OPTIONAL, RDyMMSSection, constants, mms_constants_fields_schema),
    CYAML_FIELD_MAPPING("swe", CYAML_FLAG_OPTIONAL, RDyMMSSection, swe, mms_swe_fields_schema),
    CYAML_FIELD_END
};

// -----------------
// top-level schemas
// -----------------

// RDycore's schema for top-level configuration fields
static const cyaml_schema_field_t config_fields_schema[] = {
    CYAML_FIELD_MAPPING("physics", CYAML_FLAG_DEFAULT, RDyConfig, physics, physics_fields_schema),
    CYAML_FIELD_MAPPING("numerics", CYAML_FLAG_DEFAULT, RDyConfig, numerics, numerics_fields_schema),
    CYAML_FIELD_MAPPING("time", CYAML_FLAG_DEFAULT, RDyConfig, time, time_fields_schema),
    CYAML_FIELD_MAPPING("logging", CYAML_FLAG_OPTIONAL, RDyConfig, logging, logging_fields_schema),
    CYAML_FIELD_MAPPING("checkpoint", CYAML_FLAG_OPTIONAL, RDyConfig, checkpoint, checkpoint_fields_schema),
    CYAML_FIELD_MAPPING("restart", CYAML_FLAG_OPTIONAL, RDyConfig, restart, restart_fields_schema),
    CYAML_FIELD_MAPPING("output", CYAML_FLAG_OPTIONAL, RDyConfig, output, output_fields_schema),
    CYAML_FIELD_MAPPING("grid", CYAML_FLAG_DEFAULT, RDyConfig, grid, grid_fields_schema),
    CYAML_FIELD_SEQUENCE_COUNT("surface_composition", CYAML_FLAG_DEFAULT, RDyConfig, surface_composition, num_material_assignments, &surface_composition_entry, 0, MAX_NUM_REGIONS),
    CYAML_FIELD_SEQUENCE_COUNT("materials", CYAML_FLAG_DEFAULT, RDyConfig, materials, num_materials, &material_entry, 0, MAX_NUM_MATERIALS),
    CYAML_FIELD_SEQUENCE_COUNT("regions", CYAML_FLAG_DEFAULT, RDyConfig, regions, num_regions,
                               &region_spec_entry, 0, MAX_NUM_REGIONS),
    CYAML_FIELD_SEQUENCE_COUNT("initial_conditions", CYAML_FLAG_DEFAULT, RDyConfig, initial_conditions, num_initial_conditions, &region_condition_spec_entry, 0, MAX_NUM_REGIONS),
    CYAML_FIELD_SEQUENCE_COUNT("boundaries", CYAML_FLAG_OPTIONAL, RDyConfig, boundaries, num_boundaries,
                               &boundary_spec_entry, 0, MAX_NUM_BOUNDARIES),
    CYAML_FIELD_SEQUENCE_COUNT("boundary_conditions", CYAML_FLAG_OPTIONAL, RDyConfig, boundary_conditions, num_boundary_conditions,
                               &boundary_condition_spec_entry, 0, MAX_NUM_BOUNDARIES),
    CYAML_FIELD_SEQUENCE_COUNT("sources", CYAML_FLAG_OPTIONAL, RDyConfig, sources, num_sources, &region_condition_spec_entry, 0, MAX_NUM_REGIONS),
    CYAML_FIELD_SEQUENCE_COUNT("flow_conditions", CYAML_FLAG_DEFAULT, RDyConfig, flow_conditions, num_flow_conditions, &flow_condition_entry, 0,
                               MAX_NUM_CONDITIONS),
    CYAML_FIELD_SEQUENCE_COUNT("sediment_conditions", CYAML_FLAG_OPTIONAL, RDyConfig, sediment_conditions, num_sediment_conditions,
                               &sediment_condition_entry, 0, MAX_NUM_CONDITIONS),
    CYAML_FIELD_SEQUENCE_COUNT("salinity_conditions", CYAML_FLAG_OPTIONAL, RDyConfig, salinity_conditions, num_salinity_conditions,
                               &salinity_condition_entry, 0, MAX_NUM_CONDITIONS),
    CYAML_FIELD_MAPPING("ensemble", CYAML_FLAG_OPTIONAL, RDyConfig, ensemble, ensemble_fields_schema),
    CYAML_FIELD_END
};

// schema for top-level configuration datum itself
static const cyaml_schema_value_t config_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_POINTER, RDyConfig, config_fields_schema),
};

// MMS driver's schema for top-level configuration fields
static const cyaml_schema_field_t mms_config_fields_schema[] = {
    CYAML_FIELD_MAPPING("physics", CYAML_FLAG_DEFAULT, RDyConfig, physics, physics_fields_schema),
    CYAML_FIELD_MAPPING("numerics", CYAML_FLAG_DEFAULT, RDyConfig, numerics, numerics_fields_schema),
    CYAML_FIELD_MAPPING("time", CYAML_FLAG_DEFAULT, RDyConfig, time, time_fields_schema),
    CYAML_FIELD_MAPPING("logging", CYAML_FLAG_OPTIONAL, RDyConfig, logging, logging_fields_schema),
    CYAML_FIELD_MAPPING("output", CYAML_FLAG_OPTIONAL, RDyConfig, output, output_fields_schema),
    CYAML_FIELD_MAPPING("grid", CYAML_FLAG_DEFAULT, RDyConfig, grid, grid_fields_schema),
    CYAML_FIELD_SEQUENCE_COUNT("regions", CYAML_FLAG_DEFAULT, RDyConfig, regions, num_regions,
                               &region_spec_entry, 0, MAX_NUM_REGIONS),
    CYAML_FIELD_SEQUENCE_COUNT("boundaries", CYAML_FLAG_OPTIONAL, RDyConfig, boundaries, num_boundaries,
                               &boundary_spec_entry, 0, MAX_NUM_BOUNDARIES),
    CYAML_FIELD_MAPPING("mms", CYAML_FLAG_DEFAULT, RDyConfig, mms, mms_fields_schema),
    CYAML_FIELD_END
};

// schema for top-level MMS configuration datum itself
static const cyaml_schema_value_t mms_config_schema = {
    CYAML_VALUE_MAPPING(CYAML_FLAG_POINTER, RDyConfig, mms_config_fields_schema),
};

// clang-format on

// CYAML log function (of type cyaml_log_fn_t)
static void YamlLog(cyaml_log_t level, void *ctx, const char *fmt, va_list args) {
  // render a log string
  char message[1024];
  vsnprintf(message, 1023, fmt, args);
  MPI_Comm *comm = ctx;
  PetscFPrintf(*comm, stdout, "%s", message);
}

// CYAML memory allocation function (of type cyaml_mem_fn_t)
static void *YamlAlloc(void *ctx, void *ptr, size_t size) {
  if (size) {
    PetscRealloc(size, &ptr);
    return ptr;
  } else {  // free
    PetscFree(ptr);
    return NULL;
  }
}

// parses the given YAML string into the given config representation using the
// given configuration schema
static PetscErrorCode ParseYaml(MPI_Comm comm, const char *yaml_str, const cyaml_schema_value_t *config_schema, RDyConfig **config) {
  PetscFunctionBegin;

  // configure our YAML parser
  cyaml_config_t yaml_config = {
      .log_fn    = YamlLog,
      .log_ctx   = &comm,
      .mem_fn    = YamlAlloc,
      .log_level = CYAML_LOG_WARNING,
  };

  const uint8_t *yaml_data     = (const uint8_t *)yaml_str;
  size_t         yaml_data_len = strlen(yaml_str);
  cyaml_err_t    err           = cyaml_load_data(yaml_data, yaml_data_len, &yaml_config, config_schema, (void **)config, NULL);
  PetscCheck(err == CYAML_OK, comm, PETSC_ERR_USER, "Error parsing config file: %s", cyaml_strerror(err));

  PetscFunctionReturn(PETSC_SUCCESS);
}

// sets missing parameters to their default values for the parameters
#define SET_MISSING_PARAMETER(param, value) \
  if (!param) param = value
static PetscErrorCode SetMissingValues(RDyConfig *config) {
  PetscFunctionBegin;

  SET_MISSING_PARAMETER(config->physics.flow.tiny_h, 1e-7);

  SET_MISSING_PARAMETER(config->time.final_time, INVALID_REAL);
  SET_MISSING_PARAMETER(config->time.max_step, INVALID_INT);
  SET_MISSING_PARAMETER(config->time.time_step, INVALID_REAL);
  SET_MISSING_PARAMETER(config->time.coupling_interval, INVALID_REAL);

  PetscFunctionReturn(PETSC_SUCCESS);
}
#undef SET_MISSING_PARAMETER

// checks config for any invalid or omitted parameters
static PetscErrorCode ValidateConfig(MPI_Comm comm, RDyConfig *config, PetscBool mms_mode) {
  PetscFunctionBegin;

  // check ensemble settings
  PetscCheck(config->ensemble.size == 0 || config->ensemble.size > 1, comm, PETSC_ERR_USER,
             "Ensemble size (%" PetscInt_FMT ") must be greater than 1", config->ensemble.size);
  PetscCheck(config->ensemble.size == config->ensemble.members_count, comm, PETSC_ERR_USER,
             "Declared ensemble size (%" PetscInt_FMT ") does not match the number of members (%" PetscInt_FMT ")", config->ensemble.size,
             config->ensemble.members_count);
  if (config->ensemble.size) {
    PetscMPIInt nproc;
    MPI_Comm_size(comm, &nproc);
    PetscCheck((nproc % config->ensemble.size) == 0, comm, PETSC_ERR_USER,
               "In ensemble mode, the ensemble size must evenly divide the number of processes.");
  }

  // check numerics settings
  if (config->numerics.spatial != SPATIAL_FV) {
    PetscCheck(PETSC_FALSE, comm, PETSC_ERR_USER, "Only the finite volume spatial method (FV) is currently implemented.");
  }
  if (config->numerics.temporal != TEMPORAL_EULER) {
    PetscCheck(PETSC_FALSE, comm, PETSC_ERR_USER, "Only the forward euler temporal method (EULER) is currently implemented.");
  }
  if (config->numerics.riemann != RIEMANN_ROE) {
    PetscCheck(PETSC_FALSE, comm, PETSC_ERR_USER, "Only the roe riemann solver (ROE) is currently implemented.");
  }

  PetscCheck(strlen(config->grid.file), comm, PETSC_ERR_USER, "grid.file not specified!");

  // check adaptive time step setting
  if (config->time.adaptive.enable) {
    PetscCheck(config->time.adaptive.target_courant_number > 0.0, comm, PETSC_ERR_USER,
               "time.adaptive.target_courant_number must be greater than 0.0");
    PetscCheck(config->time.adaptive.target_courant_number < 1.0, comm, PETSC_ERR_USER, "time.adaptive.target_courant_number must be less than 1.0");
    PetscCheck(config->time.adaptive.max_increase_factor > 1.0, comm, PETSC_ERR_USER, "time.adaptive.max_increase_factor must be greater than 1.0");
    PetscCheck(config->time.adaptive.initial_time_step > 0.0, comm, PETSC_ERR_USER, "time.adaptive.initial_time_step must be greater than 0.0");

    // ensure that max_step and time_step is not specified
    PetscCheck(config->time.max_step == INVALID_INT, comm, PETSC_ERR_USER, "max_step cannot be specified with adaptive time step enabled");
    PetscCheck(config->time.time_step == INVALID_REAL, comm, PETSC_ERR_USER, "time_step cannot be specified with adaptive time step enabled");

    // set time step
    config->time.time_step = config->time.adaptive.initial_time_step;
  }

  // check time settings
  // 'final_time', 'max_step', 'time_step': exactly two of these three can be specified in the .yaml file.
  PetscInt num_time_settings = 0;
  if (config->time.final_time != INVALID_REAL) ++num_time_settings;
  if (config->time.max_step != INVALID_REAL) ++num_time_settings;
  if (config->time.time_step != INVALID_REAL) ++num_time_settings;
  PetscCheck(num_time_settings, comm, PETSC_ERR_USER,
             "Exactly 2 of time.final_time, time.max_step, time.time_step must be specified (%" PetscInt_FMT " given)", num_time_settings);

  // set the third parameter based on the two that are given
  if (config->time.final_time == INVALID_REAL) {
    config->time.final_time = config->time.max_step * config->time.time_step;
  } else if (config->time.max_step == INVALID_INT) {
    config->time.max_step = (PetscInt)(config->time.final_time / config->time.time_step);
  } else {  // config->time.time_step == INVALID_REAL
    config->time.time_step = config->time.final_time / config->time.max_step;
  }

  // by default, the coupling interval is set to the final time (but in any case,
  // it shouldn't be greater!)
  if (config->time.coupling_interval == INVALID_REAL) {
    config->time.coupling_interval = config->time.final_time;
  } else {
    PetscCheck(config->time.coupling_interval > 0.0, comm, PETSC_ERR_USER, "time.coupling_interval must be positive");
    PetscCheck(config->time.coupling_interval <= config->time.final_time, comm, PETSC_ERR_USER,
               "time.coupling_interval must not exceed time.final_time");
  }

  // we need initial conditions and material properties specified for each
  // region for non-MMS runs
  if (!mms_mode) {
    PetscCheck(config->num_initial_conditions > 0, comm, PETSC_ERR_USER, "No initial conditions were specified!");
    PetscCheck(config->num_initial_conditions == config->num_regions, comm, PETSC_ERR_USER,
               "%" PetscInt_FMT " initial conditions were specified in initial_conditions (exactly %" PetscInt_FMT " needed)",
               config->num_initial_conditions, config->num_regions);

    PetscCheck(config->num_material_assignments == config->num_regions, comm, PETSC_ERR_USER,
               "Only %" PetscInt_FMT " material <-> region assignments were found in surface_composition (%" PetscInt_FMT " needed)",
               config->num_material_assignments, config->num_regions);

    // validate our materials
    PetscCheck(config->num_materials > 0, comm, PETSC_ERR_USER, "No materials specified!");
  } else {  // mms mode
    // check that we have expressions for our manufactured solutions
    if (config->physics.flow.mode == FLOW_SWE) {
      PetscCheck(config->mms.swe.expressions.h[0], comm, PETSC_ERR_USER, "No expression for h was specified!");
      PetscCheck(config->mms.swe.expressions.dhdx[0], comm, PETSC_ERR_USER, "No expression for dh/dx was specified!");
      PetscCheck(config->mms.swe.expressions.dhdy[0], comm, PETSC_ERR_USER, "No expression for dh/dy was specified!");
      PetscCheck(config->mms.swe.expressions.dhdt[0], comm, PETSC_ERR_USER, "No expression for dh/dt was specified!");
      PetscCheck(config->mms.swe.expressions.u[0], comm, PETSC_ERR_USER, "No expression for u was specified!");
      PetscCheck(config->mms.swe.expressions.dudx[0], comm, PETSC_ERR_USER, "No expression for du/dx was specified!");
      PetscCheck(config->mms.swe.expressions.dudy[0], comm, PETSC_ERR_USER, "No expression for du/dy was specified!");
      PetscCheck(config->mms.swe.expressions.dudt[0], comm, PETSC_ERR_USER, "No expression for du/dt was specified!");
      PetscCheck(config->mms.swe.expressions.v[0], comm, PETSC_ERR_USER, "No expression for v was specified!");
      PetscCheck(config->mms.swe.expressions.dvdx[0], comm, PETSC_ERR_USER, "No expression for dv/dx was specified!");
      PetscCheck(config->mms.swe.expressions.dvdy[0], comm, PETSC_ERR_USER, "No expression for dv/dy was specified!");
      PetscCheck(config->mms.swe.expressions.dvdt[0], comm, PETSC_ERR_USER, "No expression for dv/dt was specified!");
      PetscCheck(config->mms.swe.expressions.n[0], comm, PETSC_ERR_USER, "No expression for n was specified!");
    }
  }

  // validate our flow conditions
  for (PetscInt i = 0; i < config->num_flow_conditions; ++i) {
    const RDyFlowCondition *flow_cond = &config->flow_conditions[i];
    PetscCheck(flow_cond->type >= 0, comm, PETSC_ERR_USER, "Flow condition type not set in flow_conditions.%s", flow_cond->name);
    if (flow_cond->type != CONDITION_REFLECTING && flow_cond->type != CONDITION_CRITICAL_OUTFLOW) {
      PetscCheck(flow_cond->height_expression[0] || flow_cond->file[0], comm, PETSC_ERR_USER, "Missing height specification for flow_conditions.%s",
                 flow_cond->name);
      PetscCheck(flow_cond->file[0] || (flow_cond->x_momentum_expression[0] && flow_cond->y_momentum_expression[0]), comm, PETSC_ERR_USER,
                 "Missing or incomplete momentum specification for flow_conditions.%s", flow_cond->name);
    }
  }

  // validate sediment conditions
  for (PetscInt i = 0; i < config->num_sediment_conditions; ++i) {
    const RDySedimentCondition *sed_cond = &config->sediment_conditions[i];
    PetscCheck(sed_cond->type >= 0, comm, PETSC_ERR_USER, "Sediment condition type not set in sediment_conditions.%s", sed_cond->name);
    PetscCheck(sed_cond->expression[0], comm, PETSC_ERR_USER, "Missing sediment concentration for sediment_conditions.%s", sed_cond->name);
  }

  // validate salinity conditions
  for (PetscInt i = 0; i < config->num_salinity_conditions; ++i) {
    const RDySalinityCondition *sal_cond = &config->salinity_conditions[i];
    PetscCheck(sal_cond->type >= 0, comm, PETSC_ERR_USER, "Salinity condition type not set in salinity_conditions.%s", sal_cond->name);
    PetscCheck(sal_cond->expression[0], comm, PETSC_ERR_USER, "Missing salinity concentration for salinity_conditions.%s", sal_cond->name);
  }

  // validate output options
  if (config->output.format != OUTPUT_NONE || config->output.step_interval > 0 || config->output.time_interval > 0.0) {
    config->output.enable           = PETSC_TRUE;
    config->output.prev_output_time = -1.0;
  } else {
    config->output.enable = PETSC_FALSE;
  }

  if (config->output.enable) {
    PetscCheck((config->output.format != OUTPUT_NONE), comm, PETSC_ERR_USER, "Output requested, but the format is not specified.");
    PetscCheck(!(config->output.step_interval == 0 && config->output.time_interval == 0.0), comm, PETSC_ERR_USER,
               "Output requested, but neither step_interval nor time_interval specified.");
    PetscCheck((config->output.step_interval >= 0), comm, PETSC_ERR_USER, "Output step interval must be specified as a positive number of steps.");
    PetscCheck((config->output.time_interval >= 0.0), comm, PETSC_ERR_USER, "Output time interval must be specified as a positive number of steps.");
    PetscCheck((config->output.batch_size == 0) || (config->output.format != OUTPUT_BINARY), comm, PETSC_ERR_USER,
               "Binary output does not support output batching");
    if ((config->output.batch_size == 0) && (config->output.format != OUTPUT_NONE) && config->output.format != OUTPUT_BINARY) {
      config->output.batch_size = 1;
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#define DEFINE_CONSTANT(model, function, constants, name) mupDefineConst(model->solutions.function, #name, constants->name)

#define DEFINE_CONSTANT_I(model, function, constants) mupDefineConst(model->solutions.function, "I", constants->I_)

#define DEFINE_CONSTANTS_FOR_FUNCTION(model, function, constants) \
  DEFINE_CONSTANT(model, function, constants, A);                 \
  DEFINE_CONSTANT(model, function, constants, B);                 \
  DEFINE_CONSTANT(model, function, constants, C);                 \
  DEFINE_CONSTANT(model, function, constants, D);                 \
  DEFINE_CONSTANT(model, function, constants, E);                 \
  DEFINE_CONSTANT(model, function, constants, F);                 \
  DEFINE_CONSTANT(model, function, constants, G);                 \
  DEFINE_CONSTANT(model, function, constants, H);                 \
  DEFINE_CONSTANT_I(model, function, constants);                  \
  DEFINE_CONSTANT(model, function, constants, J);                 \
  DEFINE_CONSTANT(model, function, constants, K);                 \
  DEFINE_CONSTANT(model, function, constants, L);                 \
  DEFINE_CONSTANT(model, function, constants, M);                 \
  DEFINE_CONSTANT(model, function, constants, N);                 \
  DEFINE_CONSTANT(model, function, constants, O);                 \
  DEFINE_CONSTANT(model, function, constants, P);                 \
  DEFINE_CONSTANT(model, function, constants, Q);                 \
  DEFINE_CONSTANT(model, function, constants, R);                 \
  DEFINE_CONSTANT(model, function, constants, S);                 \
  DEFINE_CONSTANT(model, function, constants, T);                 \
  DEFINE_CONSTANT(model, function, constants, U);                 \
  DEFINE_CONSTANT(model, function, constants, V);                 \
  DEFINE_CONSTANT(model, function, constants, W);                 \
  DEFINE_CONSTANT(model, function, constants, X);                 \
  DEFINE_CONSTANT(model, function, constants, Y);                 \
  DEFINE_CONSTANT(model, function, constants, Z)

#define DEFINE_FUNCTION(model, constants, function)                   \
  model->solutions.function = mupCreate(muBASETYPE_FLOAT);            \
  mupSetExpr(model->solutions.function, model->expressions.function); \
  DEFINE_CONSTANTS_FOR_FUNCTION(model, function, constants)

static PetscErrorCode ParseSWEManufacturedSolutions(MPI_Comm comm, RDyMMSConstants *constants, RDyMMSSWESolutions *swe) {
  PetscFunctionBegin;

  // NOTE: you must define the relavent variables (e.g. x, y or x, y, t)
  // NOTE: at the time of evaluation using mupDefineVar or mupDefineBulkVar.

  DEFINE_FUNCTION(swe, constants, h);
  DEFINE_FUNCTION(swe, constants, dhdx);
  DEFINE_FUNCTION(swe, constants, dhdy);
  DEFINE_FUNCTION(swe, constants, dhdt);

  DEFINE_FUNCTION(swe, constants, u);
  DEFINE_FUNCTION(swe, constants, dudx);
  DEFINE_FUNCTION(swe, constants, dudy);
  DEFINE_FUNCTION(swe, constants, dudt);

  DEFINE_FUNCTION(swe, constants, v);
  DEFINE_FUNCTION(swe, constants, dvdx);
  DEFINE_FUNCTION(swe, constants, dvdy);
  DEFINE_FUNCTION(swe, constants, dvdt);

  DEFINE_FUNCTION(swe, constants, z);
  DEFINE_FUNCTION(swe, constants, dzdx);
  DEFINE_FUNCTION(swe, constants, dzdy);

  DEFINE_FUNCTION(swe, constants, n);

  PetscFunctionReturn(PETSC_SUCCESS);
}

// parses mathematical expressions given for manufactured solutions, material
// properties, initial/boundary conditions, etc
static PetscErrorCode ParseMathExpressions(MPI_Comm comm, RDyConfig *config) {
  PetscFunctionBegin;

  if (config->mms.swe.expressions.h[0]) {
    PetscCall(ParseSWEManufacturedSolutions(comm, &config->mms.constants, &config->mms.swe));
  }

  // material properties
  for (PetscInt m = 0; m < config->num_materials; ++m) {
    RDyMaterialPropertiesSpec *properties = &config->materials[m].properties;
    properties->manning.value             = mupCreate(muBASETYPE_FLOAT);
    mupSetExpr(properties->manning.value, properties->manning.expression);
  }

  // flow conditions
  for (PetscInt f = 0; f < config->num_flow_conditions; ++f) {
    RDyFlowCondition *flow_cond = &config->flow_conditions[f];
    flow_cond->height           = mupCreate(muBASETYPE_FLOAT);
    mupSetExpr(flow_cond->height, flow_cond->height_expression);
    flow_cond->x_momentum = mupCreate(muBASETYPE_FLOAT);
    flow_cond->y_momentum = mupCreate(muBASETYPE_FLOAT);
    mupSetExpr(flow_cond->x_momentum, flow_cond->x_momentum_expression);
    mupSetExpr(flow_cond->y_momentum, flow_cond->y_momentum_expression);
  }

  // sediment conditions
  for (PetscInt s = 0; s < config->num_sediment_conditions; ++s) {
    RDySedimentCondition *sed_cond = &config->sediment_conditions[s];
    sed_cond->concentration        = mupCreate(muBASETYPE_FLOAT);
    mupSetExpr(sed_cond->concentration, sed_cond->expression);
  }

  // salinity conditions
  for (PetscInt s = 0; s < config->num_salinity_conditions; ++s) {
    RDySalinityCondition *sal_cond = &config->salinity_conditions[s];
    sal_cond->concentration        = mupCreate(muBASETYPE_FLOAT);
    mupSetExpr(sal_cond->concentration, sal_cond->expression);
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

// adds entries to PETSc's options database based on parsed config info
// (necessary because not all functionality is exposed via PETSc's C API)
static PetscErrorCode SetAdditionalOptions(RDy rdy) {
#define VALUE_LEN 128
  PetscFunctionBegin;
  PetscBool has_param, loc_refine;
  char      value[VALUE_LEN + 1] = {0};

  PetscCall(PetscOptionsHasName(NULL, NULL, "-dm_refine", &loc_refine));
  PetscCall(PetscOptionsHasName(NULL, "ref_", "-dm_refine", &rdy->refine));
  if (rdy->config.grid.file[0] != '\0') PetscCall(PetscOptionsSetValue(NULL, "-dm_plex_filename", rdy->config.grid.file));
  PetscCall(PetscOptionsSetValue(NULL, "-dist_dm_distribute_save_sf", "true"));
  if (loc_refine) PetscCall(PetscOptionsSetValue(NULL, "-dm_plex_transform_label_match_strata", "true"));
  if (rdy->refine) PetscCall(PetscOptionsSetValue(NULL, "-ref_dm_plex_transform_label_match_strata", "true"));

  //--------
  // Output
  //--------

  // set the solution monitoring interval (except for XDMF, which does its own thing)
  if ((rdy->config.output.step_interval > 0) && (rdy->config.output.format != OUTPUT_XDMF)) {
    PetscCall(PetscOptionsHasName(NULL, NULL, "-ts_monitor_solution_interval", &has_param));
    if (!has_param) {
      snprintf(value, VALUE_LEN, "%" PetscInt_FMT "", rdy->config.output.step_interval);
      PetscOptionsSetValue(NULL, "-ts_monitor_solution_interval", value);
    }
  }

  // allow the user to override our output format
  PetscBool ts_monitor = PETSC_FALSE;
  char      ts_monitor_opt[128];
  PetscCall(PetscOptionsGetString(NULL, NULL, "-ts_monitor_solution", ts_monitor_opt, 128, &ts_monitor));
  if (ts_monitor) {
    // which is it?
    if (strstr(ts_monitor_opt, "cgns")) {
      rdy->config.output.format = OUTPUT_CGNS;
    } else if (strstr(ts_monitor_opt, "xdmf")) {
      rdy->config.output.format = OUTPUT_XDMF;
    } else if (strstr(ts_monitor_opt, "binary")) {
      rdy->config.output.format = OUTPUT_BINARY;
    }
  } else {  // TS monitoring not set on command line
    // the CGNS viewer's options aren't exposed at all by the C API, so we have
    // to set them here
    if (rdy->config.output.format == OUTPUT_CGNS) {
      // configure timestep monitoring parameters
      char file_pattern[PETSC_MAX_PATH_LEN];
      PetscCall(DetermineOutputFile(rdy, 0, 0.0, "cgns", file_pattern));
      snprintf(value, VALUE_LEN, "cgns:%s", file_pattern);
      PetscOptionsSetValue(NULL, "-ts_monitor_solution", value);
    }
  }

  // allow the user to override our checkpoint format
  if (rdy->config.checkpoint.interval) {
    PetscBool format_overridden = PETSC_FALSE;
    char      format[16];
    PetscCall(PetscOptionsGetString(NULL, NULL, "-checkpoint_format", format, 16, &format_overridden));
    if (format_overridden) {
      PetscCheck(!strcmp(format, "binary") || !strcmp(format, "hdf5"), rdy->comm, PETSC_ERR_USER, "Invalid checkpoint format given: %s", format);
      if (!strcmp(format, "binary")) {
        rdy->config.checkpoint.format = PETSC_VIEWER_NATIVE;
      } else {  // hdf5
        rdy->config.checkpoint.format = PETSC_VIEWER_HDF5_PETSC;
      }
    }
  }

  // set the solution monitoring interval (except for XDMF, which does its own thing)
  if ((rdy->config.output.step_interval > 0) && (rdy->config.output.format != OUTPUT_XDMF)) {
    PetscCall(PetscOptionsHasName(NULL, NULL, "-ts_monitor_solution_interval", &has_param));
    if (!has_param) {
      snprintf(value, VALUE_LEN, "%" PetscInt_FMT, rdy->config.output.step_interval);
      PetscOptionsSetValue(NULL, "-ts_monitor_solution_interval", value);
    }
  }

  // adjust the CGNS output batch size if needed
  if (rdy->config.output.format == OUTPUT_CGNS) {
    PetscCall(PetscOptionsHasName(NULL, NULL, "-viewer_cgns_batch_size", &has_param));
    if (!has_param) {
      snprintf(value, VALUE_LEN, "%" PetscInt_FMT "", rdy->config.output.batch_size);
      PetscOptionsSetValue(NULL, "-viewer_cgns_batch_size", value);
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
#undef VALUE_LEN
}

typedef struct {
  const char *pattern;
  const char *substitution;
} Substitution;

// supported string substitutions
static const Substitution substitutions[] = {
    {"${PETSC_ID_TYPE}", PETSC_ID_TYPE},
    {NULL,               NULL         }, // terminator
};

// ON RANK 0 ONLY, reads the given file and performs the given set of string
// substitutions, storing the resulting (newly allocated) string in content
// and its size in content_size
static PetscErrorCode ReadAndSubstitute(MPI_Comm comm, const char *filename, const Substitution substitutions[], char **content,
                                        PetscMPIInt *content_size) {
  PetscFunctionBegin;

  FILE *file = NULL;
  PetscCall(PetscFOpen(comm, filename, "r", &file));

  // determine the file's size and read it into a buffer
  fseek(file, 0, SEEK_END);
  PetscMPIInt raw_size = (PetscMPIInt)ftell(file);
  rewind(file);
  char *raw_content;
  PetscCall(PetscCalloc1(raw_size + 1, &raw_content));
  fread(raw_content, sizeof(char), raw_size, file);
  PetscCall(PetscFClose(comm, file));
  raw_content[raw_size] = 0;  // null termination for C string

  // determine the size of the content with all substitutions applied
  PetscInt num_substitutions = 0;
  *content_size              = raw_size;
  for (PetscInt s = 0; substitutions[s].pattern; ++s) {
    const Substitution sub         = substitutions[s];
    PetscInt           pattern_len = (PetscInt)strlen(sub.pattern);
    PetscInt           subst_len   = (PetscInt)strlen(sub.substitution);
    char              *p           = raw_content;
    while (p != NULL) {
      p = strstr(p, sub.pattern);
      if (p != NULL) {
        *content_size += subst_len - pattern_len;
        p += pattern_len;
        ++num_substitutions;
      }
    }
  }

  // perform any needed string substitutions or just use the raw input
  if (num_substitutions > 0) {
    PetscCall(PetscCalloc1(*content_size + 1, content));
    for (PetscInt s = 0; substitutions[s].pattern; ++s) {
      const Substitution sub         = substitutions[s];
      PetscInt           subst_len   = (PetscInt)strlen(sub.substitution);
      PetscInt           pattern_len = (PetscInt)strlen(sub.pattern);
      char              *p = raw_content, *q = *content;
      while (p != NULL) {
        char *new_p = strstr(p, sub.pattern);
        if (new_p != NULL) {
          memcpy(q, p, new_p - p);
          q += new_p - p;
          p = new_p + pattern_len;
          memcpy(q, sub.substitution, subst_len);
          q += subst_len;
        } else {
          memcpy(q, p, raw_size - (p - raw_content));
          q += raw_size - (p - raw_content);
          p = new_p;
        }
      }
      PetscCheck(q - *content == *content_size, comm, PETSC_ERR_USER, "error performing string substitutions in %s!", filename);
    }
    PetscFree(raw_content);
  } else {
    *content = raw_content;
  }
  (*content)[*content_size] = 0;

  PetscFunctionReturn(PETSC_SUCCESS);
}

/// Determines the prefix of the YAML configuration file.
PetscErrorCode DetermineConfigPrefix(RDy rdy, char *prefix) {
  PetscFunctionBegin;

  memset(prefix, 0, sizeof(char) * (strlen(rdy->config_file) + 1));
  char *p = strstr(rdy->config_file, ".yaml");
  if (!p) {  // could be .yml, I suppose (Windows habits die hard!)
    p = strstr(rdy->config_file, ".yml");
  }
  if (p) {
    size_t prefix_len = p - rdy->config_file;
    strncpy(prefix, rdy->config_file, prefix_len);
  } else {
    strcpy(prefix, rdy->config_file);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// reads the config file on process 0, broadcasts it as a string to all other
// processes, making it available as config_str
static PetscErrorCode ReadAndBroadcastConfigFile(RDy rdy, char **config_str) {
  PetscFunctionBegin;
  PetscMPIInt config_size;
  if (rdy->rank == 0) {
    // process 0: read the file and perform substitutions
    PetscCall(ReadAndSubstitute(rdy->comm, rdy->config_file, substitutions, config_str, &config_size));

    // broadcast the size of the content and then the content itself
    MPI_Bcast(&config_size, 1, MPI_INT, 0, rdy->comm);
    MPI_Bcast(*config_str, config_size, MPI_CHAR, 0, rdy->comm);
  } else {
    // other processes: read the size of the content
    MPI_Bcast(&config_size, 1, MPI_INT, 0, rdy->comm);

    // recreate the configuration string.
    PetscCall(PetscCalloc1(config_size + 1, config_str));
    MPI_Bcast(*config_str, config_size, MPI_CHAR, 0, rdy->comm);
    (*config_str)[config_size] = 0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// reads the config file into a string and parses the string into rdy->config
PetscErrorCode ReadConfigFile(RDy rdy) {
  PetscFunctionBegin;

  char *config_str;
  PetscCall(ReadAndBroadcastConfigFile(rdy, &config_str));

  // parse the YAML config file into a new config struct and validate it
  RDyConfig *config;
  PetscCall(ParseYaml(rdy->comm, config_str, &config_schema, &config));
  PetscCall(SetMissingValues(config));
  PetscCall(ValidateConfig(rdy->comm, config, PETSC_FALSE));
  PetscCall(ParseMathExpressions(rdy->comm, config));

  // copy the config into place and dispose of the original
  rdy->config = *config;
  PetscFree(config);

  // if this is an ensemble run, split our communicator, assign ranks to
  // ensemble members, and override parameters
  if (rdy->config.ensemble.size > 1) {
    PetscCall(ConfigureEnsembleMember(rdy));
  } else {
    rdy->ensemble_member_index = -1;  // not a member of an ensemble
  }

  // set any additional options needed in PETSc's options database
  PetscCall(SetAdditionalOptions(rdy));

  // clean up
  PetscFree(config_str);

  PetscFunctionReturn(PETSC_SUCCESS);
}

// reads the config file for the MMS driver into a string and parses the string
// into rdy->config
PetscErrorCode ReadMMSConfigFile(RDy rdy) {
  PetscFunctionBegin;

  char *config_str;
  PetscCall(ReadAndBroadcastConfigFile(rdy, &config_str));

  // parse the YAML config file into a new config struct and validate it
  RDyConfig *config;
  PetscCall(ParseYaml(rdy->comm, config_str, &mms_config_schema, &config));
  PetscCall(SetMissingValues(config));
  PetscCall(ValidateConfig(rdy->comm, config, PETSC_TRUE));
  PetscCall(ParseMathExpressions(rdy->comm, config));

  // copy the config into place and dispose of the original
  rdy->config = *config;
  PetscFree(config);

  // if this is an ensemble run, split our communicator, assign ranks to
  // ensemble members, and override parameters
  if (rdy->config.ensemble.size > 1) {
    PetscCall(ConfigureEnsembleMember(rdy));
  } else {
    rdy->ensemble_member_index = -1;  // not a member of an ensemble
  }

  // set any additional options needed in PETSc's options database
  PetscCall(SetAdditionalOptions(rdy));

  // clean up
  PetscFree(config_str);

  PetscFunctionReturn(PETSC_SUCCESS);
}

// =============
//  PrintConfig
// =============

static const char *FlagString(PetscBool flag) { return flag ? "enabled" : "disabled"; }

static PetscErrorCode PrintEnsemble(RDy rdy) {
  PetscFunctionBegin;
  if (rdy->config.ensemble.size) {
    RDyLogDetail(rdy, "Ensemble:");
    RDyLogDetail(rdy, "  Size: %" PetscInt_FMT, rdy->config.ensemble.size);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PrintPhysics(RDy rdy) {
  PetscFunctionBegin;
  RDyLogDetail(rdy, "Physics:");
  RDyLogDetail(rdy, "  Flow:");
  RDyLogDetail(rdy, "  Sediment model: %s", FlagString(rdy->config.physics.sediment));
  RDyLogDetail(rdy, "  Salinity model: %s", FlagString(rdy->config.physics.salinity));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static const char *SpatialString(RDyNumericsSpatial method) {
  static const char *strings[2] = {"finite volume (FV)", "finite element (FE)"};
  return strings[method];
}

static const char *TemporalString(RDyNumericsTemporal method) {
  static const char *strings[3] = {"forward euler", "4th-order Runge-Kutta", "backward euler"};
  return strings[method];
}

static const char *RiemannString(RDyNumericsRiemann solver) {
  static const char *strings[2] = {"roe", "hllc"};
  return strings[solver];
}

static PetscErrorCode PrintNumerics(RDy rdy) {
  PetscFunctionBegin;
  RDyLogDetail(rdy, "Numerics:");
  RDyLogDetail(rdy, "  Spatial discretization: %s", SpatialString(rdy->config.numerics.spatial));
  RDyLogDetail(rdy, "  Temporal discretization: %s", TemporalString(rdy->config.numerics.temporal));
  RDyLogDetail(rdy, "  Riemann solver: %s", RiemannString(rdy->config.numerics.riemann));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static const char *TimeUnitString(RDyTimeUnit unit) {
  static const char *strings[6] = {"seconds", "minutes", "hours", "days", "months", "years"};
  return strings[unit];
}

static PetscErrorCode PrintTime(RDy rdy) {
  PetscFunctionBegin;
  RDyLogDetail(rdy, "Time:");
  RDyLogDetail(rdy, "  Final time: %g %s", rdy->config.time.final_time, TimeUnitString(rdy->config.time.unit));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PrintCheckpoint(RDy rdy) {
  PetscFunctionBegin;
  RDyLogDetail(rdy, "Checkpoint:");
  if (rdy->config.checkpoint.interval > 0) {
    char format[12];
    if (rdy->config.checkpoint.format == PETSC_VIEWER_NATIVE) {
      strcpy(format, "binary");
    } else {
      strcpy(format, "hdf5");
    }
    RDyLogDetail(rdy, "  File format: %s", format);
    RDyLogDetail(rdy, "  interval: %" PetscInt_FMT, rdy->config.checkpoint.interval);
  } else {
    RDyLogDetail(rdy, "  (disabled)");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PrintRestart(RDy rdy) {
  PetscFunctionBegin;
  RDyLogDetail(rdy, "Restart:");
  if (rdy->config.restart.file[0]) {
    RDyLogDetail(rdy, "  File: %s", rdy->config.restart.file);
    if (rdy->config.restart.reinitialize) {
      RDyLogDetail(rdy, "  (time is reinitialized to zero)");
    }
  } else {
    RDyLogDetail(rdy, "  (disabled)");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PrintLogging(RDy rdy) {
  PetscFunctionBegin;
  RDyLogDetail(rdy, "Logging:");
  if (strlen(rdy->config.logging.file)) {
    RDyLogDetail(rdy, "  Primary log file: %s", rdy->config.logging.file);
  } else {
    RDyLogDetail(rdy, "  Primary log file: <stdout>");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

// prints config information at the requested log level (single-run mode only)
PetscErrorCode PrintConfig(RDy rdy) {
  PetscFunctionBegin;

  RDyLogDetail(rdy, "==========================================================");
  RDyLogDetail(rdy, "RDycore (input read from %s)", rdy->config_file);
  RDyLogDetail(rdy, "==========================================================");

  PetscCall(PrintEnsemble(rdy));
  PetscCall(PrintPhysics(rdy));
  PetscCall(PrintNumerics(rdy));
  PetscCall(PrintTime(rdy));
  PetscCall(PrintLogging(rdy));
  PetscCall(PrintCheckpoint(rdy));
  PetscCall(PrintRestart(rdy));

  PetscFunctionReturn(PETSC_SUCCESS);
}
