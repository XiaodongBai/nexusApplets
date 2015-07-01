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
## file Demultiplexed.mk
##
## brief Partial makefile to carry out demultiplexing of qseq files in indexed runs.
##
## author Roman Petrovski
##
################################################################################
lanes:=$(lanes) l$(lane)

paddedTilesList:=$(l$(lane)_TILES)
include $(foreach paddedTile, $(l$(lane)_TILES), $(MAKEFILES_DIR)/ConvertStats.mk)

Basecall_Stats_$(FLOWCELL)/l$(lane)_Signal_Means.txt: lane:=$(lane)
Basecall_Stats_$(FLOWCELL)/l$(lane)_Signal_Means.txt:  Basecall_Stats_$(FLOWCELL)/SignalMeans/.sentinel \
  $(l$(lane)_TILES:%=Basecall_Stats_$(FLOWCELL)/SignalMeans/s_$(lane)_%_all.txt)
	$(LIBEXEC_DIR)/BaseCalls/Signal_Means.pl \
		$(SM_PARAMS) Basecall_Stats_$(FLOWCELL)/SignalMeans $(l$(lane)_TILES:%=s_$(lane)_%_all.txt) > $(SAFEPIPETARGET)

Basecall_Stats_$(FLOWCELL)/l$(lane)$(SMT_PNGS_SUFFIX): lane:=$(lane)
Basecall_Stats_$(FLOWCELL)/l$(lane)$(SMT_PNGS_SUFFIX): \
  Basecall_Stats_$(FLOWCELL)/tiles.txt Basecall_Stats_$(FLOWCELL)/SignalMeans/.sentinel \
  Basecall_Stats_$(FLOWCELL)/Plots/.sentinel Basecall_Stats_$(FLOWCELL)/Temp/.sentinel \
  $(foreach paddedTile, $(l$(lane)_TILES), Basecall_Stats_$(FLOWCELL)/SignalMeans/s_$(lane)_$(paddedTile)_all.txt) \
  Basecall_Stats_$(FLOWCELL)/config.xml
	$(LIBEXEC_DIR)/BaseCalls/plotIntensity_tiles.pl Basecall_Stats_$(FLOWCELL) \
		tiles.txt s_$(lane) SignalMeans _all.txt _all.png \
	$(AND) $(TOUCH) $(UNPROTECTED_TARGET)

demuxReadsList:= $(DEMUX_READS)
includedOriginalReadsList:=$(INCLUDED_ORIGINAL_READS)
include $(foreach demux_read, $(DEMUX_READS), $(MAKEFILES_DIR)/DemultiplexLaneRead.mk)

# Explicit rule for all files produced by a single lane as input:
l$(lane).done: $(foreach demux_read, $(DEMUX_READS), l$(lane)_r$(demux_read).done)
