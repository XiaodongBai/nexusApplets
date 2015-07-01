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
## brief Main tweaks for validation framework 
##
## author Roman Petrovski
##
################################################################################

# looks like split command line spec has not settled yet -b300M works on fedora but not on centos
# -b 300m seems to work on both.
SPLIT_TAR_PARAMS:= -d -b20m
GENERATE_TAR_PARAMS:=-z

DEMUX_VALIDATE_IGNORE_FILES:= support.txt moveConfigXml.log $(NAME_SHELL_FILE) .sentinel
DEMUX_VALIDATE_NOSIZE_FILES:= *.gif *.png Demultiplex_Stats.htm Makefile SampleSheet.mk

# introduces empty lines between each line of the file. 
# That way diff does not group lines of the same file together
INTERSPERSE:=sed 's/^\(.*\)$$/\1\n/'
