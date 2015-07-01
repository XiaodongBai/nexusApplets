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
## file Debug.mk
##
## brief Partial makefile providing basic debugging support
##
## - Redefines the SHELL
## - provides print finctionalities
##
## author Come Raczy
##
################################################################################

PRINTCMDGOALS:=$(filter print-%, $(MAKECMDGOALS))

ifneq (,$(PRINTCMDGOALS))
# Target to print the value of a variable and exit
print-%: ; @$(error $* is: $($*))

# This will translate the command line "make all print-X" into the dependency
# all: print-X
# as a result, X will get the target-specific value
$(filter-out print-%, $(MAKECMDGOALS)): $(filter print-%, $(MAKECMDGOALS))
endif
