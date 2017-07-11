# Motif ALTernative Exons Scanner Enrichment of RNA-Seq

This tool enriches exons that have been spliced as seen in dexseq DEXSeqResults.
It first creates a temporary fasta file with the sequences of the exons spliced, and the rest of the gene sans the spliced exon (background).
It then annotates the sequences extracted with prosite. Finally it calculates the score of each motif by:

score = (motifs in exon / size of exon)/(motifs in background / size of background)

If there is no motifs in the background then the score just gives the count of motifs flagged with an 'N' (number).

## Manual
```
Usage: maltese [options] diff_splicing_file

 Takes as input a dexseq output file and enriches the loci with motifs
overrepresented (compared with the rest of the gene).

Options:
  -h, --help            show this help message and exit
  -V, --version         prints version
  -P PROCESSES, --Processes=PROCESSES
                        number of processors to use

  Required parameters:
    Parameters required to run maltesers

    -T ANNOTATIONS, --taxon=ANNOTATIONS
                        Organisms of study (Mus_musculus, Homo_sapiens...) if
                        no genome files presents,it will try to downloads and
                        generate anotation files of defined taxon
    -t ANNOTATION_VERSION, --taxon_version=ANNOTATION_VERSION
                        defines what genome version to download, default is
                        "GRCm38"
    -a ANNOTATIONSDIR, --annotations=ANNOTATIONSDIR
                        directory with annotations: the gtf file
    -s SEP, --sep=SEP   what separator is present in the input. Also used for
                        output
    -o OUTPUT, --output=OUTPUT
                        output file, Default will output results in current
                        directory
    -F INPUTFORMAT, --format=INPUTFORMAT
                        what format is the input datait takes a string such as
                        "0,8,9,10,12,-" (default dexseq format)each number
                        represents the column where certain information is."en
                        trezID,geneName,chromosome,start,end,strand,pvalue,cha
                        nge"change can be "-" if none present, it is only used
                        for plotting

  Optional parameters:
    -f PLOTFORMAT, --plotFormat=PLOTFORMAT
                        Which format to save the plots, default is pdf
    -v, --verbose       shows you what am I thinking

  Debugging parameters:
    -S, --SkipProsite   Skips the steps leading to analysing the prosite
                        output
    -p, --purge         delete all anotation files
    -m, --Tempfiles     does not erase temprary files
```
## Output

maltesers outputs a few files (with the prefix stipulated in the -o OUTPUT or --output=OUTPU option):
- OUTPUT: The original differential splicing file with six new columns prepended:
  - motif: Motif name.
  - logFold2: The	log2 difference in motif density between the exon and the rest of the gene.
  - ExonCount: The number of amino acids taken by the specific motif in the exon.
  - exonLen: The size of the exon.
  - motifGeneCount: The number of amino acids taken by the specific motif in the gene (minus exon).
  - geneLen: The size of the gene (minus exon).
- OUTPUT_exons.pdf: a boxplot of hits (per logfold2 score) per exon, colored by the change if available.
- OUTPUT_motifs.pdf: a boxplot of hits (per logfold2 score) per exon, colored by the change if available.
- OUTPUT_motifsExon.pdf: clustered heatmap of motifs scores per exon.
- OUTPUT_motifsExon.csv: CSV with the values of the OUTPUT_motifsExon.pdf plot.

temporary files (can be deleted):
- tmp.fasta: a fasta file with the sequences for the exons and the genes for ps_scan.
- tmp.fasta.prosite: prosite output.

As long as the tmp.fasta.prosite the plots can be redone in a different format by reruning maltesers with the same prameters but using the -S (--SkipProsite) and the desired image format -f (--plotFormat) and the plots will be redone without rerunning ps_scan.


## Prequisites

python 2.7

Module | version tested
-------|---------------
biopython| 1.68
matplotlib|1.5.1
scipy|0.17.0

optional for the clustermotif plot:

Module | version tested
-------|---------------
pandas|0.19.2
seaborn|0.7.0

## Installation
The program does not need to be compiled and can run in an UNIX environment.
In order to install simply download repository
```
git clone https://github.com/aLahat/maltese.git

```
or
```
wget https://github.com/aLahat/maltese/archive/dexMotif.zip && unzip dexMotif.zip
```
The program can be run staight from its folder just by calling it 
```
python2.7 [path to maltesers folder]/bin/maltese [options]
```
In order to install maltesers
```
chmod +x [path to maltesers folder]/bin/maltese
```
and then add the bin folder to your path, for example:
```
echo export \${PATH}=[path to maltesers folder]/bin]:\${PATH} >> ~/.bashrc
```
### Testing
To test maltesers, genomes file for human chromosome 19, rMATS output (only chromosome 19), and the expected result of a correct maltesers run are bundled in the test folder.
To test maltesers just run the maltesers_test.sh script.
```
bash maltesers_test.sh
```
The test will run maltesers and then compered if the result is correct.
