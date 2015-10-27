Usage: dexMotif.py [options] [dexseq] 

Takes as input a dexseq output file and enriches the loci with
motifsoverrepresented (compared with the rest of the gene).

Options:
  -h, --help            show this help message and exit
  -A ANNOTATIONS, --annotation=ANNOTATIONS
                        downloads and generates anotation files of defined
                        taxon -A "mus_musculus"
  -t ANNOTATION_VERSION, --taxon_version=ANNOTATION_VERSION
                        defines what genome version to download, default is
                        "GRCm38"
  -p, --purge           delete all anotation files (usefull for changing the
                        taxon of interest)
  -v, --verbose         shows you what am I thinking
  -B, --Bsub            run each dexseq in bsub
  -b BSUB_OPTIONS, --Bsub_options=BSUB_OPTIONS
                        options for bsubs
  -V, --version         prints version



This tool enriches exons that have been spliced as seen in dexseq DEXSeqResults.
It first creates a temporary fasta file with the sequences of the exons spliced, and the rest of the gene sans the spliced exon (background).
It then annotates the sequences extracted with prosite. Finally it calculates the score of each motif by:

score = (motifs in exon / size of exon)/(motifs in background / size of background)

If there is no motifs in the background then the score just gives the count of motifs flagged with an 'N' (number).



