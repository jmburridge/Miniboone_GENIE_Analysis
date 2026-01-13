# -----------------------
# Mode control
# -----------------------
MODE ?= test
ifeq ($(filter $(MODE),run test),)
  $(error MODE must be 'run' or 'test')
endif

# -----------------------
# Base directories
# -----------------------
MACRO_DIR ?= macros

# Put all products under a mode-specific directory
BASE_OUT  ?= output/$(MODE)

# Subdirectories

PLOT_DIR        ?= $(BASE_OUT)/plots/
CACHE_DIR       ?= $(BASE_OUT)/cache/
# Matrices are shared (mode-independent)
MATRIX_ROOT_DIR ?= outputs/matrices/root/
MATRIX_PNG_DIR  ?= outputs/matrices/png/

# -----------------------
# MiniBooNE dataset inputs (for create_all_matrices.C)
# -----------------------
MB_DATASET_DIR ?= /path/to/MiniBooNEDatasets2023
MB_PREFIX      ?= output_osc_mc_detail_
MB_FIRST_FILE  ?= 1
MB_NFILES      ?= 10
MB_TREE_NAME   ?= MiniBooNE_CCQE

# -----------------------
# Truth inputs: GENIE vs NUANCE (selected by MODE)
# -----------------------
GENIE_INPUT_FILE ?= /path/to/genie_or_larsoft_output.root
# optionally: GENIE_TREE_NAME ?= ...
NUANCE_INPUT_FILE ?= /path/to/nuance_output.root
# optionally: NUANCE_TREE_NAME ?= ...

GENIE_PROCESSED_FILE  ?= output/run/cache/truth_genie.root
NUANCE_PROCESSED_FILE ?= output/test/cache/truth_nuance.root

ifeq ($(MODE),test)
  FOLD_INPUT_TRUTH := $(NUANCE_PROCESSED_FILE)
else
  FOLD_INPUT_TRUTH := $(GENIE_PROCESSED_FILE)
endif

# -----------------------
# Matrix inputs (produced in-mode)
# Adjust filenames to match what your matrix macro writes
# -----------------------
PI0_MATRIX_FILE   ?= $(MATRIX_ROOT_DIR)/response_pi0_LEE.root
DELTA_MATRIX_FILE ?= $(MATRIX_ROOT_DIR)/response_ncdelta_LEE.root
NUE_MATRIX_FILE   ?= $(MATRIX_ROOT_DIR)/response_nue_LEE.root
RESP_MATRIX_NAME  ?= h_resp_LEE

# -----------------------
# Downstream products (in-mode)
# -----------------------
TEMPLATES ?= $(BASE_OUT)/cache/templates_reco.root
MONEYPLOT ?= $(PLOT_DIR)/moneyplot.pdf
FOLDED_FILE ?= $(BASE_OUT)/cache/forward_folds_output.root
