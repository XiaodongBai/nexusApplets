################################################################################
##
## Copyright (c) 2007-2012 Illumina, Inc.
##
## This software is covered by the "Illumina Genome Analyzer Software 
## License Agreement" and the "Illumina Source Code License Agreement", 
## and certain third party copyright/licenses, and any user of this 
## source file is bound by the terms therein (see accompanying files
## Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
## Illumina_Source_Code_License_Agreement.pdf and third party 
## copyright/license notices).
##
## This file is part of the Consensus Assessment of Sequence And VAriation
## (CASAVA) software package.
##
## file generate.mk
##
## brief top-level structure that implements validation --generate mode 
##
## author Roman Petrovski
##
################################################################################

default: all

include $(MAKEFILES_DIR)/Validation/Config.mk

# include validation targets so that the generated chain can be compared 
# to the original if needed. Do it before messing up the ROOT locations below

include $(MAKEFILES_DIR)/Validation/ValidationSummary.mk

include $(MAKEFILES_DIR)/Validation/Demux.mk


Generate.finished: Demux.finished
	$(TOUCH) $(UNPROTECTED_TARGET)

Validation.tar.gz/.finished: tar_params:= \
  $(GENERATE_TAR_PARAMS) $(foreach ign, $(NAME_SHELL_FILE), --exclude '$(ign)')

Validation.tar.gz/.finished: Generate.finished Validation.tar.gz/.sentinel 
	$(TAR) $(tar_params) -cO \
	  -C $(BASE_CALLS_ROOT)/../../../ Data RunInfo.xml \
	  -C $(CURDIR) Unaligned \
	| split $(SPLIT_TAR_PARAMS) - $(dir $@)Validation.tar.gz. \
	$(AND) $(TOUCH) $(UNPROTECTED_TARGET)

all: Validation.tar.gz/.finished
	$(LOG_INFO) Validation dataset generated successfully
