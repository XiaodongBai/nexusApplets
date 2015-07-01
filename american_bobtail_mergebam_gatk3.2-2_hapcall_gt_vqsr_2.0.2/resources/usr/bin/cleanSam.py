#!/usr/bin/python

import re
import sys


fo=open(sys.argv[1]+".bwa.unmapped.sam","w")
for line in sys.stdin:
	commentTag = re.compile("^@.*")
	if commentTag.match(line):
		sys.stdout.flush()
		sys.stdout.write(line)
	else:
		splLine = line.split()
		nonMapped = re.compile('[\d+XY]')
		if nonMapped.match(splLine[2]):
			sys.stdout.flush()
			sys.stdout.write(line)
		else:
			fo.write(line)

fo.close()
			
