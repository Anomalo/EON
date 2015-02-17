from scipy import stats
import glob
import urllib
import os

#Arabidopsis_thaliana_GSEA_GO_sets_all_symbols_September_2013.gmt.zip
#Anaplasma_marginale_str_St_Maries_GSEA_GO_sets_all_symbols_September_2013.gmt.zip
#Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gmt.zip

def getOnlineGMTs():
	'''
	makes a dictionary from http://www.go2msig.org/ with available .gmt (go files) 
	and respective URLS
	'''
	origin = 'http://www.go2msig.org/cgi-bin/prebuilt.cgi?'
	html = urllib.urlopen(origin).read().split('<br>')[3:-2]
	taxids = {}
	for line in html:
		name = line.split('>')[1].split('<')[0]
		url = line.split('"')[1]#.split('"')[1]
		url = 'http://www.go2msig.org/cgi-bin'+url[1:]
		lastURL = urllib.urlopen(url).read().split('a href="')[-1]
		gmtURL = lastURL.split('"')[0]
		gmtName= gmtURL.split('/')[-1].replace('.zip','')
		taxids.update({gmtName:gmtURL})
	return taxids


def updateGMTs(purge = False, prefix = 'annotations/',taxon='Mus musculus'):
	'''
	Checks available .gmt files online and downloads more if available
	purge -- set True to remove all existing .gmt files and redownload them
	prefix -- directory where to save file
	taxon -- taxon of interest
	'''
	gmtsAvail = getOnlineGMTs()
	if purge: os.system('rm -rf '+prefix+'*.gmt')
	
	gmtsDownloaded = glob.glob(prefix+'*.gtm')
	
	for gmt in gmtsDownloaded:
		name = gmt.split('/')[-1]
		if name in gmtsAvail:
			del gmtsAvail[name]
	for name, url in gmtsAvail.items():
		taxon = taxon.upper().replace(' ','_')
		if not taxon in name.upper(): continue
		name = prefix + name
		os.system('clear')
		os.system(' '.join(['wget',
				    url,
				    '-P',
				    prefix]))
		os.system(' '.join(['unzip',
				    name,
				    '-d',
				    prefix]))
		os.system('rm '+name+'.zip')

#updateGMTs(purge = True)
class GO:
	def __init__(self):
		self.readCONFIG()
		self.readGO()
		self.mkGOgenes()
		self.mkGOfreq()
		self.mkGOdef()
		#GOgenes = {} #{gene:[goIDs]...}
		#GOfreq = {} #{GOID: frequency...}

	def mkGOdef(self):
		'''makes a self.GOdef dictionary of {GOID:GO long name}'''
		GOdef = {}
		for line in self.GOTXT:
			items = line.split()
			GOID = items[1].split('=')[-1]
			GOdefinition = items[0]
			GOdef.update({GOID:GOdefinition})
		self.GOdef = GOdef

	def mkGOfreq(self):
		'''makes a self.GOfreq dictionary of {GOID : 0.45 (frequency of that annotation)}'''
		GOfreq = {}
		geneNumber = len(self.GOgenes)
		for line in self.GOTXT:
			items = line.split()
			GOID = items[1].split('=')[-1]
			amount = len(items)-2
			frequency = amount/float(geneNumber)
			GOfreq.update({GOID:frequency})
		self.GOfreq = GOfreq

		
	
	def mkGOgenes(self):
		''' makes a self.GOgenes dictionary of {gene : [GOID, GOID ...]}'''
		GOgenes = {}
		for line in self.GOTXT:
			items = line.split()
			GOlabel = items.pop(0)
			GOID = items.pop(0).split('=')[-1]
			for gene in items:
				gene = gene.upper()
				if not gene in GOgenes:
					GOgenes.update({gene:[]})
				GOgenes[gene].append(GOID)
		self.GOgenes = GOgenes
	def GOgeneNames(self, gene):
		'''returns a list of go names for that gene'''
		goids = self.GOgenes[gene]
		gonames = []
		for id in goids:
			gonames.append(self.GOdef[id])
		return gonames

	def readGO(self):
		''' makes a self.GOTXT wich is a list of all the lines of the .gmt file'''
		f = open(self.GOfile,'r')
		self.GOTXT = f.read().splitlines()
		f.close()


	def readCONFIG(self, fname = 'config.txt'):
		''' makes a self.GOfile, a string of the .gmt file name from the config.txt file'''
		f = open(fname, 'r')
		exec(f.read())
		f.close()
		GO = 'annotations/' + GO
		self.GOfile = GO
	
	def enrich(self, genes, key='GOID', P_cut_off = 0.05):
		'''
		Returns a dictionary of go annotaions present with their P value 
		

		Keyword arguments:
		genes -- a list of genes to enrich
		key --  what key to use for the dictionary to return, default 'GOID' otherwise use GOlabel
		P_cut_off -- the P value to filter with.
		'''
		GOfound = []
		for gene in genes:
			GOfound += self.GOgenes[gene.upper()]
		GOenriched = dict.fromkeys(set(GOfound))
		for GOID in GOenriched:
			x = GOfound.count(GOID)
			n = len(genes)
			frequency = self.GOfreq[GOID]
			allgenes = len(self.GOgenes)
			background = frequency*allgenes			
			P = stats.fisher_exact([[n,x],[allgenes,background]])	
			GOenriched[GOID] = P[-1]
			#print GOID,self.GOdef[GOID],P
		keys = GOenriched.keys()
		for GOID in keys:
			if GOenriched[GOID] > P_cut_off:
				del GOenriched[GOID]
		
		return GOenriched


go = GO()

def enrich(genes,key='GOID', P_cut_off = 0.05):
	return go.enrich(genes, key, P_cut_off)

