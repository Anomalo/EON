import os
def getOBO(dir = 'annotations/'):
	url ='http://geneontology.org/ontology/go.obo'
	os.system(' '.join([	'wget',
				url,
				'-P',
				dir]))

class OBO:
	def __init__(self, fname):
		self.OBO=self.readOBO(fname) 
	
	def findTopGOs(self):
		out = {}
		for key, entry in self.OBO.iteritems():
			if 'is_obsolete' in entry:continue
			if not 'is_a' in entry:
				out[key]=entry
		return out
	
	def stepList(self, GO):
		steps = []
		while True:
			steps.append((GO,self.OBO[GO]['name']))
			if not 'is_a' in self.OBO[GO]:break
			GO = self.OBO[GO]['is_a']
			if type(GO)==list:GO=GO[-0]
			
		return steps
	def readOBO(self, fname):
		f = open(fname,'r')
		obo = f.read()
		f.close()
		obo=obo.split('\n\n')
		OBO={}
		for block in obo:
			lines = block.splitlines()
			if len(lines)<2: continue
			term = lines.pop(0)
			if term != '[Term]':
				continue
			entry={}
			for line in lines:
				key, value = line.split(': ')[:2]
				if key == 'is_a':
					value = value.split()[0]
				if key in entry:
					try: entry[key].append(value)
					except: entry[key]=[entry[key],value]
				else:entry[key]=value
			OBO[entry['id']]=entry
		return OBO
	def levelUp(self, GO,steps=1):
		ladder = self.stepList(GO)
		if len(ladder)<steps:return steps[0]
		else: return ladder[steps]
			
	def levelMax(self, GO,steps=1):
		ladder= self.stepList(GO)
		if len(ladder)<steps:return steps[0]
		else: return ladder[-steps]

