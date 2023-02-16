#ifndef RDYLOGIMPL_H
#define RDYLOGIMPL_H

#include <petscsys.h>

//---------
// Logging
//---------

// RDycore logs messages at varying levels of detail to a log file. The file
// and the log level are set in the configuration file.

// This type defines the various logging levels in order of ascending detail.
typedef enum {
  LOG_NONE = 0, // no messages are written to the log
  LOG_WARNING,  // only warnings are written
  LOG_INFO,     // warnings and basic information are written
  LOG_DETAIL,   // detailed information is added
  LOG_DEBUG,    // debugging information is added
} RDyLogLevel;

// These are labels that get prepended to log messages.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static const char *RDyLogLabels[5] = {
  "",
  "WARNING: ",
  "INFO: ",
  "DETAIL: ",
  "DEBUG: "
};
#pragma GCC diagnostic pop

// Writes a log message at the given level of detail. You don't need to wrap
// this call in PetscCall.
#define RDyLog(rdy, level, ...) \
  if (((rdy)->config.log_level > LOG_NONE) && (level >= (rdy)->config.log_level)) { \
    PetscCall(PetscFPrintf((rdy)->comm, (rdy)->log, "%s", RDyLogLabels[(rdy)->config.log_level])); \
    PetscCall(PetscFPrintf((rdy)->comm, (rdy)->log, __VA_ARGS__)); \
  }

// These are convenient macros for writing log messages at specific levels.
#define RDyLogWarning(rdy, ...) RDyLog(rdy, LOG_WARNING, __VA_ARGS__)
#define RDyLogInfo(rdy, ...) RDyLog(rdy, LOG_INFO, __VA_ARGS__)
#define RDyLogDetail(rdy, ...) RDyLog(rdy, LOG_DETAIL, __VA_ARGS__)
#define RDyLogDebug(rdy, ...) RDyLog(rdy, LOG_DEBUG, __VA_ARGS__)

#endif