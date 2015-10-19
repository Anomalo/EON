import os
import os.path
from glob import glob

_taxon="mus_musculus"
_version='GRCm38'

def set_taxon(taxon="mus_musculus",version='GRCm38'):
	global _taxon
	global _version
	_taxon = taxon
	_version = version

def chr_url(chromosome):
	url ="ftp://ftp.ensembl.org/pub/release-78/fasta/%(_taxon)s/dna/%(_taxon)s.%(_version)s.dna.chromosome.%(chromosome)s.fa.gz"%locals()
	return url

def get_fa(chromosome, folder="annotations/genome"):
	if not os.path.exists(folder):
		os.makedirs(folder)
	url = chr_url(chromosome, _taxon)
	cmd = " ".join(["wget",
			url,
			"-P",
			folder])
	os.system(cmd)
	filename = url.split("/")[-1]
	cmd = " ".join(["gunzip",
			folder+"/"+filename])
	os.system(cmd)
	'''
	cmd = " ".join(["rm",
			folder+"/"+filename])
	os.system(cmd)
	'''

def chr_filename(chromosome):
	'''
	returns the filename for the chromosome genome
	'''
	files = glob("annotations/genome/*.dna.chromosome.%(chromosome)s.fa"%locals())
	if files == []: return ''
	else : return files[0]

def seq_coords(chromosome, start=0, end=-1, direction='+'):
	'''
	returns the genomic sequence
	'''
	filename = chr_filename(chromosome)
	if filename == '':
		get_fa(chromosome)	
	filename = chr_filename(chromosome)
	f = open(filename, "r")
	seq = f.read().replace("\n","")
	start,end = int(start),int(end)
	if direction == '+':
		return seq[start:end]
	else:
		seqReversed = seq[start:end][::-1]
		convert = {	'A':'T',
				'C':'G',
				'G':'C',
				'T':'A',
				}
		seq = ''
		for i in seqReversed:
			seq+=convert[i]
		return seq


def seqs_coords(coordinates):
	'''
	return a list of sequences based on a list of coords tuples [(chromosome, start,
	end, strand),...]
	'''
	seqs = []
	for chromosome, start, end, strand in coordinates:
		seqs.append(seq_coords(chromosome, start,end, strand))
	return seqs
