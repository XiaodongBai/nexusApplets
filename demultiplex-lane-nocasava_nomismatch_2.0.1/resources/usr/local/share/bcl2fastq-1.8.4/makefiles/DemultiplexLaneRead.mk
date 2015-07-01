################################################################################
##
## Copyright (c) 2007-2009 Illumina, Inc.
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
## file DemultiplexLaneRead.mk
##
## brief Partial makefile to carry out demultiplexing of qseq files in indexed runs.
##
## author Roman Petrovski
##
################################################################################

ifeq (,$(demuxReadsList))
$(error $(lane): empty list of demux read numbers)
endif

demux_read := $(word 1, $(demuxReadsList))
demuxReadsList := $(wordlist 2, $(words $(demuxReadsList)), $(demuxReadsList))

original_read := $(word 1, $(includedOriginalReadsList))
includedOriginalReadsList := $(wordlist 2, $(words $(includedOriginalReadsList)), $(includedOriginalReadsList))

ANY_BARCODE:=NoIndex
UNKNOWN_BARCODE:=Undetermined
BARCODE_CYCLES_PARAM:=--barcode-cycles "$(BARCODE_CYCLES)"

COMPRESSED_FASTQ_SUFFIX:=$(FASTQ_SUFFIX)$(COMPRESSIONSUFFIX)

l$(lane)_barcode_pos = $(word $(pos), $(l$(lane)_BARCODES))
l$(lane)_sampleid_pos = $(word $(pos), $(l$(lane)_SAMPLEIDS))
l$(lane)_subdir_pos = $(word $(pos), $(l$(lane)_SUBDIRS))

l$(lane)_pos_list := $(shell $(SEQ) 1 $(words $(l$(lane)_SUBDIRS)))

# remove indexes that point to '' subdirs
l$(lane)_pos_list := $(foreach pos, $(l$(lane)_pos_list), $(if $(filter-out '',$(l$(lane)_subdir_pos)), $(pos)) )

l$(lane)_r$(demux_read)_temp_dir := $(TEMP_DIR)/l$(lane)_r$(demux_read)

ifneq (,$(l$(lane)_BARCODES))
l$(lane)_any_barcode_pos = $(l$(lane)_barcode_pos)
else
l$(lane)_any_barcode_pos := $(ANY_BARCODE)
endif
l$(lane)_r$(demux_read)_summary:=$(TEMP_DIR)/L00$(lane)_R$(demux_read)$(DEMUX_SUMMARY_SUFFIX)

# list of targets with pattern to make the demultiplexBcls execute once 
# for each lane rather than once for each target
l$(lane)_r$(demux_read)_pattern_targets:=$(foreach pos, $(l$(lane)_pos_list),\
    $(l$(lane)_subdir_pos)/$(l$(lane)_sampleid_pos)_$(l$(lane)_any_barcode_pos)_L00$(lane)_R$(demux_read)%001$(COMPRESSED_FASTQ_SUFFIX)) \
    $(TEMP_DIR)/L00$(lane)_R$(demux_read)%demux_summary.xml

# list of explicit targets to instantiate the pattern target rule
l$(lane)_r$(demux_read)_results:=$(foreach pos, $(l$(lane)_pos_list),\
    $(l$(lane)_subdir_pos)/$(l$(lane)_sampleid_pos)_$(l$(lane)_any_barcode_pos)_L00$(lane)_R$(demux_read)_001$(COMPRESSED_FASTQ_SUFFIX)) \

#preventing rebuild of everything on make all due to cleanup of Temp
.SECONDARY: $(l$(lane)_r$(demux_read)_summary)

ifeq ($(ANY_BARCODE),$(l$(lane)_BARCODES))
# non-multiplexed lanes are just named differently unknown barcodes.
UNKNOWN_BARCODE:=$(ANY_BARCODE)
# Cleaning the variables that would confuse the demultiplexer otherwise...
BARCODE_CYCLES_PARAM:=
endif

# this is where we list the dependencies on input data (BCL files)
l$(lane)_r$(demux_read)_prereqs:=

# IMPORTANT: this line has to be like this here with spaces the way they are. Make sure
# you've tested all 'if' paths with qmake if you are changing it. Ohterwise it segfaults
filterDirParam := $(if $(FILTER_DIR), --filter-dir $(FILTER_DIR) , )

# move of the result files must happen so that the files *_001$(COMPRESSED_FASTQ_SUFFIX)
# get moved last. Otherwise the termination of a half-way done move will not be recognized
# as an incomplete run when make is restarted (BF-1029).

# Pattern rule for all output files of a lane read:
$(l$(lane)_r$(demux_read)_pattern_targets): lane := $(lane)
$(l$(lane)_r$(demux_read)_pattern_targets): BARCODE_CYCLES_PARAM:= $(BARCODE_CYCLES_PARAM)
$(l$(lane)_r$(demux_read)_pattern_targets): UNKNOWN_BARCODE:=$(UNKNOWN_BARCODE)
$(l$(lane)_r$(demux_read)_pattern_targets): demux_read := $(demux_read)
$(l$(lane)_r$(demux_read)_pattern_targets): original_read := $(original_read) 
$(l$(lane)_r$(demux_read)_pattern_targets): \
  $(foreach lsd, $(sort $(l$(lane)_SUBDIRS)), $(l$(lane)_r$(demux_read)_temp_dir)/$(lsd)/.sentinel)
	$(CMDPREFIX) $(LIBEXEC_DIR)/demultiplexBcls \
        $(DEMUX_OPTIONS) \
        --basecalls-dir $(BASECALLS_DIR) \
        --intensities-dir $(INTENSITIES_DIR) \
        --output-dir $(CURDIR) \
        --output-format fastq \
        $(filterDirParam) \
        --compression $(COMPRESSION) \
        --flow-cell-id $(FLOWCELL) \
        --unknown-barcode '$(UNKNOWN_BARCODE)' \
        --lane $(lane) \
        $(foreach t,$(l$(lane)_TILES), --tile $(t)) \
        $(BARCODE_CYCLES_PARAM) \
        --input-read-number $(original_read) \
        --output-read-number $(demux_read) \
        --read-cycles "$(r$(demux_read)_CYCLES)" \
        $(foreach sa,$(r$(demux_read)_ADAPTER),--adapter $(sa)) \
        --output-summary-path $(l$(lane)_r$(demux_read)_summary) \
        $(foreach sd,$(l$(lane)_SAMPLEIDS),--sample $(sd)) \
        $(foreach sb,$(l$(lane)_BARCODES),--barcode $(sb)) \
        $(foreach sd,$(l$(lane)_SUBDIRS),--sample-dir $(sd)) \
        $(if $(IGNORE_MISSING_BCL),--ignore-missing-bcl) \
        $(if $(IGNORE_MISSING_CTRL),--ignore-missing-control)

# Explicit rule for all files produced by a single lane read as input:
l$(lane)_r$(demux_read).done: $(l$(lane)_r$(demux_read)_results)

ALL_DEMUX_SUMMARIES += $(l$(lane)_r$(demux_read)_summary)

