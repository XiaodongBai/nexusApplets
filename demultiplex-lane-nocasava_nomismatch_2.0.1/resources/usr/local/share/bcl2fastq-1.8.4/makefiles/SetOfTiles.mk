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
## file SetOfTiles.mk
##
## brief Partial makefile to partition a lane-based analysis into something that
##       can be analysed on a set of tiles.
##
## This partial makefile encapsulates the mechanisms to split each lane into a
## number of subsets of tiles before processing by ELAND. The size of the sets
## is defined as a number of tiles, controlled by ELAND_FASTQ_FILES_PER_PROCESS
## (customisable on the make command line).
##
## author Mauricio Varea
##
################################################################################

# Create the name for the current set. Sets are sequentially numbered from 1,
# left padded with 0 to a width of 3 (001, 002, 003, etc.)
currentSetNumber:=$(words newSet $(currentSetList))
ifeq (10,$(currentSetNumber))
setPrefix:=0
endif
ifeq (100,$(currentSetNumber))
setPrefix:=
endif
currentSetName:=$(setPrefix)$(currentSetNumber)

# Update the current list of sets
currentSetList:=$(currentSetList) $(currentSetName)

ifneq (,$(ELAND_FASTQ_FILES_PER_PROCESS))
# Get the actual set of tiles
currentSet:=$(wordlist 1, $(ELAND_FASTQ_FILES_PER_PROCESS), $(currentTileList))

# update the current list of tiles (two steps to get rid of the element at position $(ELAND_FASTQ_FILES_PER_PROCESS))
# Note, make version 3.78 incorrectly interpret $(wordlist 2, 1, a) as "a"
# instead of an empty list, which explains the convoluted update
currentTileList:=$(wordlist $(ELAND_FASTQ_FILES_PER_PROCESS), $(words $(currentTileList)), $(currentTileList))
ifeq (1,$(words $(currentTileList)))
currentTileList:=$(wordlist 2, 0, $(currentTileList))
else
currentTileList:=$(wordlist 2, $(words $(currentTileList)), $(currentTileList))
endif

$(lane)_setList := $(currentSetList)
$(lane)_set_$(currentSetName) := $(currentSet)

else # $(ELAND_FASTQ_FILES_PER_PROCESS) == 0
$(error ELAND_FASTQ_FILES_PER_PROCESS has not been specified at configuration time)
endif
