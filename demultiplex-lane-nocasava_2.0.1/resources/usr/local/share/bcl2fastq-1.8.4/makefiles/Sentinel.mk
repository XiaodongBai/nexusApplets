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
## file Log.mk
##
## brief Partial makefile providing basic folder hierarchy creation support
##
## author Roman Petrovski
##
################################################################################

# Pattern rule to create folder hierarchies. If first mkdir fails due to a race
# condition, sleep 5 seconds and retry. Note that patterns dont' work with 
# .SECONDARY on qmake, therefore all sentinels have to be declared as .PRECIOUS.
# As result, they have to be kept when cleaning up Temp or else, everything will
# be rebuilt once the Temp is destroyed.

.PRECIOUS: %/.sentinel
%/.sentinel:
	$(MKDIR) -p $* $(OR) ( $(SLEEP) 5 $(AND) $(MKDIR) -p $* ); $(TOUCH) $@

