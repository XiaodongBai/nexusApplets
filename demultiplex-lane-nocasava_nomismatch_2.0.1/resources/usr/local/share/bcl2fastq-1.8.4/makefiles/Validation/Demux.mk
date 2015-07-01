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
## file Demux.mk
##
## brief Rules for Bcl to Fastq conversion validation 
##
## author Roman Petrovski
##
################################################################################

ifeq (,$(BASE_CALLS_ROOT))
  $(error "BASE_CALLS_ROOT is required for Demux to run")
endif

Demux.configured: SHELL=$(SHELL_LOG)
Demux.configured:
	$(CMDPREFIX) \
	$(CONFIGURE_BCL_TO_FASTQ) \
	  --input-dir $(BASE_CALLS_ROOT) \
	  $(if $(SAMPLE_SHEET), --sample-sheet $(SAMPLE_SHEET)) \
	  --output-dir ./Unaligned \
	  --skip-variable-metadata \
	  --force \
	  $(CONFIGURE_BCL_TO_FASTQ_PARAMS) \
	$(CP_DEMUX_MK_CMD) \
	$(AND) $(TOUCH) $@ \
	$(AND) ($(LOG_INFO) 'Configure demultiplex [OK]')


ifneq (,$(DEMUX_CONFIG_MK))
include $(DEMUX_CONFIG_MK)
# make sure cat is used instead of copy for ValidationConfig.mk to get the default permission
# bits as opposite to inheriting the installation folder ones which are usually read-only
DEMUX_CONFIG_MK_COPY:=Unaligned/$(notdir $(DEMUX_CONFIG_MK))
$(DEMUX_CONFIG_MK_COPY): Demux.configured Unaligned/.sentinel
	$(CAT) $(DEMUX_CONFIG_MK) >$(SAFEPIPETARGET)
endif

# make sure the submake uses our MAKEFILES_DIR regardless of what it has defined internally.
# this is useful when chaning the top level makefile during debugging
Demux.finished: SHELL=$(SHELL_LOG_OLD)
Demux.finished: Demux.configured $(DEMUX_CONFIG_MK_COPY)
	$(MAKE) -C Unaligned POST_RUN_COMMAND:='$(TOUCH) ../$@' MAKEFILES_DIR:=$(MAKEFILES_DIR) $(EXECUTE_BCL_TO_FASTQ_PARAMS)
