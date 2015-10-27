Usage: dexMotif.py [options]

Takes as input a dexseq output file and enriches the loci with
motifsoverrepresented (compared with the rest of the gene)dexMotif [options]
<dexseq> > <dexseqOutput>

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
