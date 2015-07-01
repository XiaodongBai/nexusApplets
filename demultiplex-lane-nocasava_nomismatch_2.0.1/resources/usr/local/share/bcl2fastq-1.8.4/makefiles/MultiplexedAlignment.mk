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
## file MultiplexedAlignment.mk
##
## brief Partial makefile to build each GERALD directory in a demultiplexed output.
##
## author Mauricio Varea
##
################################################################################

ALL_QTABLES += $($(DIRECTORY)_QTABLE_LIST)

.PHONY: $(DIRECTORY)
$(DIRECTORY): $($(DIRECTORY)_QTABLE_LIST)
	$(CD) $@ $(AND) $(MAKE)

