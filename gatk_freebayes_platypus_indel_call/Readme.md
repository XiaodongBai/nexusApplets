Generate indel genotypes using 3 bwa based variant caller:
	gatk-3.3; Playpus-0.8.1; freebayes-v0.9.21-7-g7dd41db

merge the 3 vcfs together and delete variants had been called in less than 2 caller

extract indels from the above step
