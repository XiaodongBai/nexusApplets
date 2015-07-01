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
## file Makefile
##
## brief Master makefile to copy into the validation directory.
##
## author Roman Petrovski
##
################################################################################

# first target needs to be defined in the beginning. Ohterwise includes such as
# Log.mk cause unexpected behavior
firsttarget: default

# Location of the static makefiles
MAKEFILES_DIR:=/usr/local/share/bcl2fastq-1.8.4/makefiles

# Import the global configuration
include $(MAKEFILES_DIR)/Config.mk

include $(MAKEFILES_DIR)/Sentinel.mk

# Import the logging functionalities
include $(MAKEFILES_DIR)/Log.mk

# Import the debug functionalities
include $(MAKEFILES_DIR)/Debug.mk

#Makefiles/config.mk: Makefiles/.sentinel
#	$(LIBEXEC_CONFIGURE_ALIGNMENT) --sample-sheet $< > $(SAFEPIPETARGET)
#
# Import the configuration for the current folder
include config.mk

# force Log.pm to use the level set for make and log to STDERR only.
export CASAVA_LOG_LEVEL

include $(MAKEFILES_DIR)/Validation/$(VALIDATE_ACTION).mk

default: all
