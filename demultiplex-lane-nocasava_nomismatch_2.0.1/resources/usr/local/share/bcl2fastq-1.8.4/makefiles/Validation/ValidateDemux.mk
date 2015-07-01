################################################################################
##
## Copyright (c) 2007-2011 Illumina, Inc.
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
## file ValidateDemux.mk
##
## brief Rules that perform comparison of Bcl to Fastq conversion data 
##
## author Roman Petrovski
##
################################################################################

ifeq (,$(UNALIGNED_ROOT))
  $(error "UNALIGNED_ROOT is required for Demux validation to run")
endif


Demux.list.ori: UNALIGNED_ROOT:=$(UNALIGNED_ROOT)
Demux.list.ori: Demux.finished $(TEMP_DIR)/.sentinel
	( $(CD) $(UNALIGNED_ROOT) && $(FIND) ./ -type f \
	    $(foreach ign, $(DEMUX_VALIDATE_IGNORE_FILES), -not -name '$(ign)') \
	    $(foreach nos, $(DEMUX_VALIDATE_NOSIZE_FILES), -not -name '$(nos)') -printf '%p %s\n' \
	    $(foreach nos, $(DEMUX_VALIDATE_NOSIZE_FILES), -name '$(nos)') -printf '%p\n') \
	    | $(SORT) | $(INTERSPERSE) >$(SAFEPIPETARGET)

Demux.list.local: Demux.finished $(TEMP_DIR)/.sentinel
	( $(CD) Unaligned && $(FIND) ./ -type f \
	    $(foreach ign, $(DEMUX_VALIDATE_IGNORE_FILES), -not -name '$(ign)') \
	    $(foreach nos, $(DEMUX_VALIDATE_NOSIZE_FILES), -not -name '$(nos)') -printf '%p %s\n' \
	    $(foreach nos, $(DEMUX_VALIDATE_NOSIZE_FILES), -name '$(nos)') -printf '%p\n') \
	    | $(SORT) | $(INTERSPERSE) >$(SAFEPIPETARGET)

Demux.diff: UNALIGNED_ROOT:=$(UNALIGNED_ROOT)
Demux.diff: Demux.list.local Demux.list.ori
	( $(DIFF) Demux.list.local Demux.list.ori || ([[ 1 == $$? ]] && $(EXIT) 0)) \
	  |( $(GREP) -E '^<|^>' \
	    |$(GREP) -vE '^<\s*$$|^>\s*$$' || ([[ 1 == $$? ]] && $(EXIT) 0) ) \
	  |$(SED) 's/^<\s\(.*$$\)/< $(subst /,\/,$(CURDIR))\/Unaligned\/\1/' \
	  |$(SED) 's/^>\s\(.*$$\)/> $(subst /,\/,$(UNALIGNED_ROOT))\/\1/' \
	>$(SAFEPIPETARGET)

