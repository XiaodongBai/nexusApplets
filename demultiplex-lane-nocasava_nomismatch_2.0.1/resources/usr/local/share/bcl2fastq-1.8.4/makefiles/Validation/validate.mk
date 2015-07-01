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
## file validate.mk
##
## brief top-level structure that implements validation --validate mode 
##
## author Roman Petrovski
##
################################################################################

default: all

include $(MAKEFILES_DIR)/Validation/Config.mk

include $(MAKEFILES_DIR)/Validation/Demux.mk

include $(MAKEFILES_DIR)/Validation/ValidationSummary.mk

all: ValidationSummary.txt
	$(CAT) $< >&2
