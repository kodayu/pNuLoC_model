#!python
import re

pro2uni = {}
uninls = {}
all_pros = []
with open("uniprot-human.tab") as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		pro2uni[data[1]] = data[0]
		if re.search("nuclear localization signal", data[5].lower()):
			motif = data[5].split("MOTIF ")[1:]
			for mot in motif:
				if re.search("nuclear localization signal", mot.lower()):
					tmp1 = mot.split(';')
					tmp2 = tmp1[0].split('..')
					uninls[data[0]+'\t'+data[1]+'\t'+tmp2[0]+'\t'+tmp2[1]] = "Uniprot"
					all_pros.append(data[0]+'\t'+data[1]+'\t'+tmp2[0]+'\t'+tmp2[1])


seqnls = {}
with open("seqnls.txt") as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		orga = data[0].split('_')[1]
		if orga == "HUMAN":
			if data[0] in pro2uni:
				seqnls[pro2uni[data[0]]+'\t'+data[0]+'\t'+data[1]+'\t'+data[2]] = data[3]
				all_pros.append(pro2uni[data[0]]+'\t'+data[0]+'\t'+data[1]+'\t'+data[2])



out = open("seqnls.uniprot.nls.txt", "w")
all_pros = list(set(all_pros))
for i in all_pros:
	evid = []
	if i in uninls:
		evid.append(uninls[i])
	if i in seqnls:
		evid.append(seqnls[i])
	out.write(i+'\t'+';'.join(evid)+'\tNLS'+'\n')

with open("../8.webserver/inuloc.dataset.scores.cutoff.seqnls.uniprot.nlsdb.validnes.nesbase.uniprot.nlsdb") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if data[3] == 'Homo sapiens (Human)':
			if data[-8] != "":
				for i in data[-8].split(';'):
					out.write(data[0]+'\t'+data[1]+'\t'+i.split('...')[0]+'\t'+i.split('...')[1]+'\tValidNES\tNES'+'\n')
			elif data[-6] != "":
				for i in data[-6].split(';'):
					out.write(data[0]+'\t'+data[1]+'\t'+i.split('...')[0]+'\t'+i.split('...')[1]+'\tNESbase\tNES'+'\n')
			elif data[-4] != "":
				for i in data[-4].split(';'):
					out.write(data[0]+'\t'+data[1]+'\t'+i.split('...')[0]+'\t'+i.split('...')[1]+'\tUniProt\tNES'+'\n')


out.close()





