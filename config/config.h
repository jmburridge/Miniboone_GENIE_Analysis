//////////////////////////////////////////////////////
// config/config.h
//   Configuration paths for MiniBooNE LEE analysis
//////////////////////////////////////////////////////
#ifndef CONFIG_H
#define CONFIG_H

// ======================================================
// Absolute base path
// ======================================================
#define BASE_PATH \
"/exp/uboone/app/users/jburridg/Geometry/Analysis"

// ======================================================
// Fortran-generated code
// ======================================================
#define FORTRAN_DIR \
BASE_PATH "/MiniBooNEDatasets2023/CombinedFunctions_from_Fortran"

#define COMBINED_FUNCTIONS_H \
FORTRAN_DIR "/CombinedFunctions.h"

#define COMBINED_TYPES_H \
FORTRAN_DIR "/CombinedTypes.h"

#define COMBINED_FUNCTIONS_CXX \
FORTRAN_DIR "/CombinedFunctions.cxx"

// ======================================================
// Output directories
// ======================================================
#define OUTPUT_BASE \
"/exp/uboone/app/users/jburridg/Geometry/MiniBooNE_GENIE_Analysis/outputs"

#define MATRIX_PNG_OUTPUT_DIR OUTPUT_BASE "/matrices/png"
#define MATRIX_ROOT_OUTPUT_DIR OUTPUT_BASE "/matrices/root"

#define TEST_CACHE_DIR OUTPUT_BASE "/test/cache"
#define TEST_PLOTS_DIR OUTPUT_BASE "/test/plots"

#define RUN_CACHE_DIR OUTPUT_BASE "/run/cache"
#define RUN_PLOTS_DIR OUTPUT_BASE "/run/plots"

// ======================================================
// ROOT input files
// ======================================================
#define PATH_FF_TREES_NUANCE "truth_trees_for_FF_TEST.root"

#define PATH_FF_TREES_GENIE "truth_trees_for_FF_RUN.root"

#define PATH_PI0_MATRIX \
MATRIX_ROOT_OUTPUT_DIR "/response_pi0_LEE.root"

#define PATH_DELTA_MATRIX \
MATRIX_ROOT_OUTPUT_DIR "/response_ncdelta_LEE.root"

#define PATH_NUE_MATRIX \
MATRIX_ROOT_OUTPUT_DIR "/response_nue_LEE.root" // check this is what theyre actually called.

// ======================================================
// ROOT output + object names
// ======================================================
#define RESP_NAME   "h_resp_LEE"
#define PATH_OUTPUT "forward_folds_output.root"

// ======================================================
// Oscillation MC output prefix
// ======================================================
#define OSC_MC_PREFIX \
BASE_PATH "/output_osc_mc_detail_"

#endif // CONFIG_H
