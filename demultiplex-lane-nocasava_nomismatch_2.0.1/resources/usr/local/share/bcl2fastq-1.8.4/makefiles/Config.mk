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
## file Config.mk
##
## brief Common configuration file for all makefiles.
##
## Defines paths, file names, extensions, etc.
##
## author Come Raczy
##
################################################################################

# The shell pipefail option is required for accurate workflow control
ifeq (,$(MAKESHELL))
SHELL = /bin/bash -o pipefail
else
SHELL = $(MAKESHELL)
endif

# CASAVA installation directories
TOOLS_DIR:=/usr/local/bin
BIN_DIR:=/usr/local/bin
DATA_DIR:=/usr/local/share/bcl2fastq-1.8.4
QCAL_TOOLS_DIR:=/usr/local/bin
ELAND_DIR:=/usr/local/bin
LIBEXEC_DIR:=/usr/local/libexec/bcl2fastq-1.8.4

# Error codes

ERROR_MISSING_QSEQ=1
ERROR_MISSING_QVAL=1
ERROR_MISSING_ELAND_GENOME=1
ERROR_NO_PIPEFAIL=2
ERROR_NO_QTABLE=3
ERROR_NO_PHAGE_ALIGN=4
ERROR_MISSING_PROGRAM=5
ERROR_MISSING_CONFIG=6
ERROR_MISSING_SUMMARY=7
ERROR_MISSING_TILES=8

# System tools
CONFIGURE_DATASET := $(LIBEXEC_DIR)/Alignment/configureDataset.pl
CONFIGURE_BCL_TO_FASTQ := $(BIN_DIR)/configureBclToFastq.pl
COPY_CONFIG := $(LIBEXEC_DIR)/BaseCalls/copyConfig.pl
PRINT_SS_XML := $(LIBEXEC_DIR)/Alignment/printSampleSheetXml.pl

# System templates
BUSTAD_SUMMARY_XSL := $(DATA_DIR)/BustardSummary.xsl

AWK = awk
CAT = cat
CD = cd
CHMOD = chmod
# at the moment all usage of cp is replaced with cat (or tar). The problem with cp is that it
# preserves the permissions of the source file. For example when copying a css from the
# installation folder into the analysis folder, the copy ends up being read-only which is
# not always what is intended. Please make sure you've considered this before uncommenting
# the lines below.
#CP = cp
#CPDIR = cp -R
CUT = cut
DATE = date
DIFF = diff
ECHO = echo
EGREP = egrep
EXIT = exit
FGREP = fgrep
FIND = find
GREP = grep
GUNZIP = gunzip
MKDIR = mkdir
MV = mv
PASTE = paste
PRINTF = printf
PWD = pwd
RM = rm -f
SED = sed
SET=set
SEQ=seq
SORT = sort -T $(SORT_TMPDIR)
TAR = tar
TEE = tee
TOUCH = touch
UNIQ = uniq
XARGS = xargs
XSLTPROC = xsltproc
WC = wc
SLEEP=sleep

AND=&&
OR=||

# Logging

# level 0 - no logging except for error messages. All recipe stderr is suppressed
# level 1 - only text preffixed with known keywords such as (INFO, ERROR, WARNING and similar) is allowed
#           to escape recipes stderr
# level 2 - unfiltered stderr
CASAVA_LOG_LEVEL:=2

LOG_DATE_FORMAT:=%F %T
LOG_DEBUG=[[ $(CASAVA_LOG_LEVEL) != 2 ]] $(OR) 1>&2 $(ECHO)
LOG_WARNING=[[ $(CASAVA_LOG_LEVEL) == 0 ]] $(OR) 1>&2 $(ECHO) -e "WARNING:"
LOG_INFO=[[ $(CASAVA_LOG_LEVEL) == 0 ]] $(OR) 1>&2 $(ECHO) -e "INFO:"

#Condition checking
BOOLEAN_TRUE_WORDS:=Y y yes YES Yes on ON On true TRUE True 1 ok OK Ok
BOOLEAN_FALSE_WORDS:=N n no NO No off OFF Off false FALSE False 0
CHECK_VALID_BOOLEAN=$(if $(filter $(BOOLEAN_TRUE_WORDS) $(BOOLEAN_FALSE_WORDS), $(1)),$(1),\
                                $(error "Incorrect boolean value '$(1)'. Allowed values: $(BOOLEAN_TRUE_WORDS) $(BOOLEAN_FALSE_WORDS)"))
IS_BOOLEAN_TRUE=$(if $(strip $(1)),$(filter $(BOOLEAN_TRUE_WORDS), $(call CHECK_VALID_BOOLEAN,$(1))))

# Global macros
UNPROTECTED_TARGET=$@
SAFEPIPETARGET = $@.tmp && mv $@.tmp $@
SAFEPIPETARGET2=$(subst .,.tmp.,$@) && mv $(subst .,.tmp.,$@) $@
SAFEPIPETARGET3=$(CURDIR)/$@.tmp && mv $(CURDIR)/$@.tmp $(CURDIR)/$@
TOUCH_TARGET=$(TOUCH) $@
CHECK_TARGET=@if [ ! -e $@ ]; then echo "Error: $@ does not exist."; exit 1; fi

# Structure of the analysis folder
TEMP_DIR:=Temp
STATS_DIR:=Stats
PLOT_DIR:=Plots

# Environment variables
SORT_TMPDIR ?= $(TEMP_DIR)

# Analysis files
MISMATCH_THUMB:=Mismatch.htm
PEC_THUMB:=Perfect.htm
RUN_REPORT:=finished.txt

# Analysis file name suffixes

BCL_SUFFIX:=.bcl
FILTER_SUFFIX:=.filter

TILES_SUFFIX:=_tiles.txt
FASTQ_SUFFIX:=.fastq
EXPORT_SUFFIX:=_export.txt.gz
SMT_SUFFIX:=_all.txt
SMT_PNG_SUFFIX:=_all.png
SMT_PNGS_SUFFIX:=_all_pngs.txt
SM_SUFFIX:=_Signal_Means.txt
PEC_IMG_SUFFIX:=_errors.gif
PEC_IMG_THUMB_SUFFIX:=_errors_thumb.gif
DEMUX_SUMMARY_SUFFIX:=_demux_summary.xml
DEMUX_SUMMARY_HTML_SUFFIX:=_demux_summary.htm


#########################################################
# CASAVA Applications


# produce BustardSummary.xml
PIS:=$(LIBEXEC_DIR)/BaseCalls/produceIntensityStats.pl


#Compression format
COMPRESSION:=gzip
COMPRESSIONSUFFIX:=.gz


