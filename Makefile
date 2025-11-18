# ---------------------------------------------------------
# Main GENIE → Forward Folding → Money Plot Pipeline
# ---------------------------------------------------------

# Load shared variables
include config.mk

# Default target: run full pipeline
all: $(MONEYPLOT)

# ---------------------------------------------------------
# Stage 4: Build final "money plot"
# ---------------------------------------------------------
$(MONEYPLOT): $(TEMPLATES) $(MACRO_DIR)/build_moneyplot.C
	@echo "Stage 4: Building final money plot..."
	root -l -b -q '$(MACRO_DIR)/build_moneyplot.C("$(TEMPLATES)", "$(MONEYPLOT)")'

# ---------------------------------------------------------
# Stage 3: Forward folding
# ---------------------------------------------------------
$(TEMPLATES): $(PI0_OUT) $(NC_OUT) $(NUE_OUT) \
              $(R_PI0) $(R_NC) $(R_NUE) \
              $(MACRO_DIR)/forward_fold.C
	@echo "Stage 3: Forward folding categories..."
	root -l -b -q '$(MACRO_DIR)/forward_fold.C("$(OUTPUT_DIR)", "$(MATRIX_DIR)", "$(TEMPLATES)")'

# ---------------------------------------------------------
# Stage 1+2: Process GENIE events → split categories
# ---------------------------------------------------------
$(PI0_OUT) $(NC_OUT) $(NUE_OUT): $(LARSOFT) $(MACRO_DIR)/process_GENIE_events.C
	@echo "Stage 1+2: Processing GENIE/LArSoft outputs..."
	root -l -b -q '$(MACRO_DIR)/process_GENIE_events.C("$(LARSOFT)", "$(OUTPUT_DIR)")'

# ---------------------------------------------------------
# Clean up
# ---------------------------------------------------------
clean:
	rm -f $(OUTPUT_DIR)/*.root $(OUTPUT_DIR)/*.pdf

.PHONY: all clean
