#!python
import re
nlss, ness = [], []
with open("signals.csv") as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split(',')
		if data[2] == '"NLS"':
			nlss.append(data[1].split('"')[1])
		elif data[2] == '"NES"':
			ness.append(data[1].split('"')[1])

hspros = []
with open("../4.Analysis/4.1.hsTruncated/human.dataset.scores") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if data[1] == "Nucleus":
			hspros.append(data[0])

pro2seq = {}
with open("../0.datadeal/eukaryota.whole.dataset") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if data[0] in hspros:
			pro2seq[data[0]] = data[2]

out1 = open("./nlsdb.nls.human.pos.txt", "w")
out2 = open("./nlsdb.nes.human.pos.txt", "w")
for uni in pro2seq:
	tmpres1, tmpres2 = [], []
	for nls in nlss:
		if re.search(nls, pro2seq[uni]):
			for m in re.finditer(nls, pro2seq[uni]):
				tmpres1.append(str(m.start()+1)+'...'+str(m.end()))
		else:
			continue
	for nes in ness:
		if re.search(nes, pro2seq[uni]):
			for m in re.finditer(nes, pro2seq[uni]):
				tmpres2.append(str(m.start()+1)+'...'+str(m.end()))
		else:
			continue
	if tmpres1 != []:
		out1.write(uni+'\t'+';'.join(tmpres1)+'\n')
	if tmpres2 != []:
		out2.write(uni+'\t'+';'.join(tmpres2)+'\n')

out1.close()
out2.close()

def regionseq(pep, region):
	reg = region.split('...')
	outseq = pep[max((int(reg[0])-1)-4, 0):min(int(reg[1])+4, len(pep))]
	return(outseq)

out11 = open("./nlsdb.nls.human.pos.fasta", "w")
with open("./nlsdb.nls.human.pos.txt") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		tmpregion = data[1].split(';')
		peptide = pro2seq[data[0]]
		for i in tmpregion:
			out11.write(">"+data[0]+' '+i+'\n'+regionseq(peptide, i)+'\n')
			
out11.close()

out12 = open("./nlsdb.nes.human.pos.fasta", "w")
with open("./nlsdb.nes.human.pos.txt") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		tmpregion = data[1].split(';')
		peptide = pro2seq[data[0]]
		for i in tmpregion:
			out12.write(">"+data[0]+' '+i+'\n'+regionseq(peptide, i)+'\n')
			
out12.close()
