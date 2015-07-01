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
## file ValidationSummary.mk
##
## brief Rules to generate the top-level validation summary 
##
## author Roman Petrovski
##
################################################################################

include $(MAKEFILES_DIR)/Validation/ValidateDemux.mk

ValidationSummary.txt: Demux.diff
	$(RM) $@.tmp; \
	$(foreach d, $^, (([[ -s $(d) ]] && $(ECHO) -e "$(d) validation    \t[FAIL]") || $(ECHO) -e "$(d) validation    \t[PASS]") \
	  >>$@.tmp &&) $(MV) $@.tmp $(UNPROTECTED_TARGET)
