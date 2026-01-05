# ---------------------------------------------------------
#  # Shared configuration variables for all Makefiles
# ---------------------------------------------------------
#
#  # Directories
MACRO_DIR   = macros
MATRIX_DIR  = matrices
OUTPUT_DIR  = output
DATA_DIR    = data
TEST_DIR    = test
# MiniBooNE tuples live here
MB_DATASET_DIR = /exp/uboone/app/users/jburridg/Geometry/Analysis/MiniBooNEDatasets2023
# Where you want matrices and plots written (project-local is best)
MATRIX_DIR  = matrices
PLOT_DIR    = plots/matrices
# Dataset file pattern
MB_PREFIX   = output_osc_mc_detail_
MB_NFILES   = 10


# Input data
# LARSOFT     = $(DATA_DIR)/larsoft_output.root
#
# # Response matrices
# R_PI0       = $(MATRIX_DIR)/R_pi0.root
# R_NC        = $(MATRIX_DIR)/R_nc.root
# R_NUE       = $(MATRIX_DIR)/R_nue.root
#
# # Output files
# PI0_OUT     = $(OUTPUT_DIR)/category_pi0.root
# NC_OUT      = $(OUTPUT_DIR)/category_nc_delta.root
# NUE_OUT     = $(OUTPUT_DIR)/category_nue.root
# TEMPLATES   = $(OUTPUT_DIR)/templates_reco.root
# MONEYPLOT   = $(OUTPUT_DIR)/moneyplot.pdf
#
