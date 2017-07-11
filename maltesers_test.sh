./bin/maltese \
			-T Homo_sapiens \
			-t GRCh38 \
			-a test/annotations \
			--sep "\t"  \
			--format 1,2,3,5,6,4,18,23 \
			-v \
			-o test/MALTESERS \
			test/MATS_output/RI.chr19.txt

cmp test/MALTESERS_expected.csv test/MALTESERS.csv  >/dev/null  && \
	echo The test was successful, maltesers seems to be working || \
	echo Something went wrong, maltesers does not seem to be working 
