#!python
import re
def readFa(fa):
	with open(fa,'r') as FA:
		seqName,seq='',''
		while 1:
			line=FA.readline()
			line=line.strip('\n')
			if (line.startswith('>') or not line) and seqName:
				yield((seqName,seq))
			if line.startswith('>'):
				seqName = line[1:]
				seq=''
			else:
				seq+=line
			if not line:break


allfa = readFa("../4.Analysis/4.2.motif/hs.CPT.delta.scores.regions.sigdelta.merged.fasta5")
sub2unis = {}
with open("../9.litingting/merged.cluster.seq.txt") as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		sub2unis[data[0]] = data[1]

AA_list_k=["H","R","K","F","A","L","M","I","W","P","V","C","G","Q","N","S","Y","T","D","E"]
cor = {}
for i in AA_list_k:
	cor[i] = {}

for line in allfa:
	data = list(line)
	data[0] = re.sub("_", ".", data[0])
	if data[0] in sub2unis:
		for i in data[1]:
			if i in cor:
				if sub2unis[data[0]] in cor[i]:
					cor[i][sub2unis[data[0]]] += 1
				else:
					cor[i][sub2unis[data[0]]] = 1
	
out = open("../9.litingting/countAA.cluster.txt", "w")
for i in cor:
	for j in cor[i]:
		out.write(i+'\t'+j+'\t'+str(cor[i][j])+'\n')

out.close()