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
## file DemultiplexFlowCell.mk
##
## brief Partial makefile to build the structure that turns this into a 
##       BaseCalls dir.
##
## Detailed description comes here.
##
## author Roman Petrovski
##
################################################################################

#include $(MAKEFILES_DIR)/Sentinel.mk

Basecall_Stats_$(FLOWCELL)/Signal_Means.txt: $(foreach l, $(lanes), Basecall_Stats_$(FLOWCELL)/$(l)_Signal_Means.txt)
	$(CAT) $^ > $(SAFEPIPETARGET)

.PRECIOUS: BustardSummary.xsl BustardSummary.xml DemultiplexedBustardSummary.xml DemultiplexedBustardConfig.xml

Basecall_Stats_$(FLOWCELL)/tiles.txt: Basecall_Stats_$(FLOWCELL)/.sentinel
	$(RM) -f $@
	$(ECHO) $(foreach l, $(lanes), $($(l)_TILES:%=$(l:l%=s_%)_%)) > $(SAFEPIPETARGET)


# make sure cat is used instead of copy for target fiel to get the default permission
# bits as opposite to inheriting the source ones which are usually read-only
# BustardSummary.xsl is referenced from BustardSummary.xml. If it is missing, the browsers
# (such as IE) refuse to even display the xml file as xml.
%/BustardSummary.xsl:
	$(CMDPREFIX) $(CAT) $(BUSTAD_SUMMARY_XSL) >$(SAFEPIPETARGET)

# make sure cat is used instead of copy for target fiel to get the default permission
# bits as opposite to inheriting the source ones which are usually read-only
Basecall_Stats_$(FLOWCELL)/Matrix/% : $(BASECALLS_DIR)/Matrix/% \
  Basecall_Stats_$(FLOWCELL)/Matrix/.sentinel
	$(CMDPREFIX) $(CAT) $(wildcard $(BASECALLS_DIR)/Matrix/$*) >$(SAFEPIPETARGET)

# sort below is used purely to make the output of various flavours of make better diffable
Basecall_Stats_$(FLOWCELL)/Matrix: \
  $(foreach f, $(sort $(wildcard $(BASECALLS_DIR)/Matrix/*)), Basecall_Stats_$(FLOWCELL)/Matrix/$(notdir $(f)))
	$(CMDPREFIX) $(foreach wrong_file, \
	             $(filter-out Basecall_Stats_$(FLOWCELL)/Matrix/.sentinel $^, $(wildcard Basecall_Stats_$(FLOWCELL)/Matrix/*)),\
	             $(RM) $(wrong_file) $(AND)) \
	$(TOUCH) $(UNPROTECTED_TARGET)

# make sure cat is used instead of copy for target fiel to get the default permission
# bits as opposite to inheriting the source ones which are usually read-only
Basecall_Stats_$(FLOWCELL)/Phasing/% : $(BASECALLS_DIR)/Phasing/% \
  Basecall_Stats_$(FLOWCELL)/Phasing/.sentinel
	$(CMDPREFIX) $(CAT) $(wildcard $(BASECALLS_DIR)/Phasing/$*) >$(SAFEPIPETARGET)
	
# sort below is used purely to make the output of various flavours of make better diffable
Basecall_Stats_$(FLOWCELL)/Phasing: \
  $(foreach f, $(sort $(wildcard $(BASECALLS_DIR)/Phasing/*)), Basecall_Stats_$(FLOWCELL)/Phasing/$(notdir $(f)))
	$(CMDPREFIX) $(foreach wrong_file, \
	             $(filter-out Basecall_Stats_$(FLOWCELL)/Phasing/.sentinel $^, $(wildcard Basecall_Stats_$(FLOWCELL)/Phasing/*)),\
	             $(RM) $(wrong_file) $(AND)) \
	$(TOUCH) $(UNPROTECTED_TARGET)

Basecall_Stats_$(FLOWCELL)/BustardSummary.xml: \
  Basecall_Stats_$(FLOWCELL)/Matrix Basecall_Stats_$(FLOWCELL)/Phasing \
  Basecall_Stats_$(FLOWCELL)/IVC.htm Basecall_Stats_$(FLOWCELL)/All.htm \
  Basecall_Stats_$(FLOWCELL)/tiles.txt Basecall_Stats_$(FLOWCELL)/config.xml \
  Basecall_Stats_$(FLOWCELL)/BustardSummary.xsl
	$(CMDPREFIX) ( $(CD) Basecall_Stats_$(FLOWCELL) $(AND) $(PIS) . )

DemultiplexedBustardSummary.xml : \
  Basecall_Stats_$(FLOWCELL)/BustardSummary.xml
	$(CMDPREFIX) $(XSLTPROC) \
	  --stringparam INCLUDED_ORIGINAL_READS_PARAM '$(INCLUDED_ORIGINAL_READS)' \
	  --stringparam DEMUX_READS_PARAM '$(DEMUX_READS)' \
	  $(DATA_DIR)/ExcludeReadFromBustardSummary.xsl $< > $(SAFEPIPETARGET)

# Run plotIntensity_for_IVC.pl
Basecall_Stats_$(FLOWCELL)/plotIntensity_for_IVC_finished.txt: \
  Basecall_Stats_$(FLOWCELL)/Signal_Means.txt \
  $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(LIBEXEC_DIR)/BaseCalls/plotIntensity_for_IVC.pl $< Basecall_Stats_$(FLOWCELL) \
		$(AND) $(ECHO) dummy > $@

#Run create_IVC_thumbnail.pl
Basecall_Stats_$(FLOWCELL)/IVC.htm: \
  $(lanes:%=Basecall_Stats_$(FLOWCELL)/%$(SMT_PNGS_SUFFIX)) \
  Basecall_Stats_$(FLOWCELL)/plotIntensity_for_IVC_finished.txt
	$(CMDPREFIX) $(LIBEXEC_DIR)/BaseCalls/create_IVC_thumbnail.pl Basecall_Stats_$(FLOWCELL) > $(SAFEPIPETARGET)

# Build thumbnail of Signal_Means_Tile.pl output
CE:=${LIBEXEC_DIR}/BaseCalls/create_tile_thumbnails.pl
MAX_THUMBNAIL_OPTION:=--maxTiles=20
Basecall_Stats_$(FLOWCELL)/All.htm: \
  Basecall_Stats_$(FLOWCELL)/config.xml \
  $(lanes:%=Basecall_Stats_$(FLOWCELL)/%$(SMT_PNGS_SUFFIX)) \
  Basecall_Stats_$(FLOWCELL)/plotIntensity_for_IVC_finished.txt
	$(CMDPREFIX) ( $(CD) Basecall_Stats_$(FLOWCELL) $(AND) $(CE) all > All.htm \
	    $(AND) $(CE) all $(MAX_THUMBNAIL_OPTION) \
		--link='<a href="All.htm">Full output (Warning: may overload your browser!)</a>') > $(SAFEPIPETARGET)

Basecall_Stats_$(FLOWCELL)/config.xml: $(BASECALLS_DIR)/config.xml Basecall_Stats_$(FLOWCELL)/.sentinel 
	$(COPY_CONFIG) --input-file $< --tiles "$(TILES_FILTER)" > $(SAFEPIPETARGET)

DemultiplexedBustardConfig.xml: Basecall_Stats_$(FLOWCELL)/config.xml
	$(CMDPREFIX) $(XSLTPROC) \
	  --stringparam DEMUX_READS_PARAM '$(DEMUX_READS)' \
	  --stringparam READ_DEMUX_FIRST_CYCLES_PARAM '$(foreach dr, $(DEMUX_READS),$(word 1, $(r$(dr)_DEMUX_CYCLES)))' \
	  --stringparam READ_DEMUX_LAST_CYCLES_PARAM '$(foreach dr, $(DEMUX_READS),$(word $(words $(r$(dr)_DEMUX_CYCLES)), $(r$(dr)_DEMUX_CYCLES)))' \
	  $(DATA_DIR)/ExcludeReadFromBustardConfig.xsl $< > $(SAFEPIPETARGET)


Basecall_Stats_$(FLOWCELL)/Flowcell$(DEMUX_SUMMARY_SUFFIX): ALL_DEMUX_SUMMARIES:=$(ALL_DEMUX_SUMMARIES)
Basecall_Stats_$(FLOWCELL)/Flowcell$(DEMUX_SUMMARY_SUFFIX): \
  Basecall_Stats_$(FLOWCELL)/.sentinel $(ALL_DEMUX_SUMMARIES)
	$(CMDPREFIX) $(ECHO) '<?xml version="1.0"?><Summary/>' $(foreach part_summ, $(ALL_DEMUX_SUMMARIES), \
		| $(XSLTPROC) --stringparam with $(part_summ) $(DATA_DIR)/MergeXmlDocuments.xsl -) \
	> $(SAFEPIPETARGET)

# make sure cat is used instead of copy for target file to get the default permission
# bits as opposite to inheriting the installation folder ones which are usually read-only
Basecall_Stats_$(FLOWCELL)/css/Reports.css: Basecall_Stats_$(FLOWCELL)/css/.sentinel
	$(CMDPREFIX) $(CAT) $(DATA_DIR)/Reports.css >$(SAFEPIPETARGET)

Basecall_Stats_$(FLOWCELL)/Demultiplex_Stats.htm : \
  Basecall_Stats_$(FLOWCELL)/Flowcell$(DEMUX_SUMMARY_SUFFIX) DemultiplexConfig.xml \
  Basecall_Stats_$(FLOWCELL)/css/Reports.css
	$(CMDPREFIX) $(XSLTPROC) \
		--stringparam PROJECTS_ROOT_PATH_PARAM $(CURDIR) \
		--stringparam DEMUX_CONFIG_XML_PATH_PARAM $(CURDIR)/DemultiplexConfig.xml \
		$(DATA_DIR)/DemultiplexSummary.xsl $< >$(SAFEPIPETARGET)

clean_intermediate.done: post_run.done
	$(if $(call IS_BOOLEAN_TRUE, $(KEEP_INTERMEDIARY)), , $(FIND) $(TEMP_DIR) -type f -not -name .sentinel |$(XARGS) $(RM))
