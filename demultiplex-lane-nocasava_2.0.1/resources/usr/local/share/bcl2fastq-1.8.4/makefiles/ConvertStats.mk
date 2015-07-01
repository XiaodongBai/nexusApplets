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
## file ConvertStats.mk
##
## author Roman Petrovski
##
################################################################################
ifeq (,$(paddedTilesList))
$(error $(lane): empty list of demux padded tiles numbers)
endif

paddedTile := $(word 1, $(paddedTilesList))

ifeq (1,$(words $(paddedTilesList)))
paddedTilesList :=
else
paddedTilesList := $(wordlist 2, $(words $(paddedTilesList)), $(paddedTilesList))
endif


NO_CYCLE_CHECK:=--no-cycle-check
repeat:=1

Basecall_Stats_$(FLOWCELL)/SignalMeans/s_$(lane)_$(paddedTile)_all.txt: lane := $(lane)
laneBasePath := $(BASECALLS_DIR)/L00$(lane)

unpaddedTile := $(shell $(PRINTF) "%.0f" $(paddedTile))
Basecall_Stats_$(FLOWCELL)/SignalMeans/s_$(lane)_$(paddedTile)_all.txt: unpaddedTile := $(unpaddedTile)

# IMPORTANT: this line has to be like this here with spaces the way they are. Make sure
# you've tested all 'if' paths with qmake if you are changing it. Ohterwise it segfaults
filterFiles :=  $(if $(FILTER_DIR), \
       $(FILTER_DIR)/s_$(lane)_$(if $(FILTER_PER_READ),$(word 1,$(INCLUDED_ORIGINAL_READS))_)$(paddedTile)$(FILTER_SUFFIX) , \
       $(BASECALLS_DIR)/L00$(lane)/s_$(lane)_$(if $(FILTER_PER_READ),$(word 1,$(INCLUDED_ORIGINAL_READS))_)$(paddedTile)$(FILTER_SUFFIX) )

Basecall_Stats_$(FLOWCELL)/SignalMeans/s_$(lane)_$(paddedTile)_all.txt: \
  $(filterFiles) \
  Basecall_Stats_$(FLOWCELL)/SignalMeans/.sentinel
	$(CMDPREFIX) $(LIBEXEC_DIR)/statsToSignalMeans -l $(lane) -r $(repeat) \
      -i $(BASECALLS_DIR) \
      --filter-file $< \
      --stats-file-name s_$(lane)_$(unpaddedTile).stats \
      $(NO_CYCLE_CHECK) \
      $(if $(IGNORE_MISSING_STATS),--ignore-missing-stats) \
      $(STATS_TO_SIGNAL_MEANS_OPTIONS) \
      $(ALL_ORIGINAL_CYCLES) \
      --signal-means-file $(SAFEPIPETARGET)

include $(MAKEFILES_DIR)/MissingFiles.mk
