#!/bin/bash

# runs the signature-analyzer GPU from getzlab
# see https://github.com/getzlab/SignatureAnalyzer
# used SHA 19a00a7

signatureanalyzer \
	CTF001_to_058_paper_revisions.maf \
	-t maf \
	-n 100 \
	--reference pcawg_COMPOSITE \
	--objective poisson \
	--hg_build hg19.2bit \
	-o All_to_Prod5_SA_pcawg_COMPOSITE
