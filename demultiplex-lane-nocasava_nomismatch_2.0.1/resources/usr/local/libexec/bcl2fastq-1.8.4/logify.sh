#!/usr/bin/env bash
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
## file logifystdin.sh
##
## logging stream formatter
##
## author Roman Petrovski
##
################################################################################

#set -x
set -o pipefail
shopt -s compat31 2>/dev/null

LOG_LEVEL=$1
shift

LOG_DATE_FORMAT=$1
shift

TARGET=$1
shift

case $LOG_LEVEL in
    '' )
        cat ;;
	0 )
		cat >/dev/null;;
	1 )
		grep --line-buffered -i -E '(info|warning|error)' | while read -r l; do
            dt=$(date "+$LOG_DATE_FORMAT") || exit 2
			(echo -en [$dt]"\t[$HOSTNAME]\t[$TARGET]\t" && echo "$l") || exit 2
		done;;
	* )
		while read -r l; do
			dt=$(date "+$LOG_DATE_FORMAT") || exit 2
			(echo -en [$dt]"\t[$HOSTNAME]\t[$TARGET]\t" && echo "$l") || exit 2
		done;;
esac

