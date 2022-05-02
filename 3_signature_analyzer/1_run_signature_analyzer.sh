#!/bin/bash

# runs the signature-analyzer GPU from getzlab

signatureanalyzer \
	All_to_Prod5.aggregated.ccf.maf \
	-t maf \
	-n 10 \
	--reference pcawg_COMPOSITE \
	--objective poisson \
	--hg_build hg19.2bit \
	-o All_to_Prod5_SA_pcawg_COMPOSITE
