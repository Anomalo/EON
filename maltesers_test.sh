./bin/maltese \
			-T Homo_sapiens \
			-t GRCh38 \
			-a test/annotations \
			--sep "\t"  \
			--format 1,2,3,5,6,4,18,23 \
			-v \
			-o test/MALTESERS \
			test/MATS_output/RI.chr19.txt

