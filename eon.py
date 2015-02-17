#!/usr/bin/python

modules = ['glob','gtf','go','blastgo','fa']
for module in modules:
	print 'importing',module
	exec('import '+module)


def check_files(dir = 'annotations/',taxon='mus musculus'):
	#check gmt
	f = glob.glob(dir+'*.gmt')
	if f == []: i
	#check gtf
	f = glob.glob(dir+'*.gtf')
	if f == []: gtf.getGTF(taxon = 'mus musculus', dir = dir)
	#check gob
	f = glob.glob(dir+'*.gob')
	if f == []: blastgo.genGOB()

	
		

if __name__ == '__main__': 
	check_files()
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e



