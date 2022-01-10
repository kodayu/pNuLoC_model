#!python
import random
uni2region = {}

with open("../4.Analysis/4.1.hsTruncated/human.dataset.cpt.sigdelta.sites") as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		tmpuni = data[0]
		if tmpuni in uni2region:
			uni2region[tmpuni].append(int(data[4]))
		else:
			uni2region[tmpuni] = [int(data[4])]

out = open("../4.Analysis/4.2.motif/hs.dataset.highdelta.regions", "w")
for uni in uni2region:
	mins, maxs = [], []
	for i in range(len(sorted(uni2region[uni]))):
		if i == 0:
			mins.append(sorted(uni2region[uni])[i])
		else:
			if sorted(uni2region[uni])[i] - sorted(uni2region[uni])[i-1] == 1:
				tmpmax = sorted(uni2region[uni])[i]
			else:
				maxs.append(sorted(uni2region[uni])[i-1])
				mins.append(sorted(uni2region[uni])[i])
	maxs.append(sorted(uni2region[uni])[-1])
	tmp_res = []
	for i in range(len(mins)):
		tmp_res.append(str(mins[i])+'...'+str(maxs[i]))
	out.write(uni+'\t'+';'.join(tmp_res)+'\n')

out.close()

def regionseq(pep, region):
	reg = region.split('...')
	outseq = pep[max((int(reg[0])-1)-8, 0):min(int(reg[1])+8, len(pep))]
	return(outseq)
	
def randomseq(pep, region):
	reg = region.split('...')
	length = min(int(reg[1])+8, len(pep)) - max((int(reg[0])-1)-8, 0)
	start = random.randint(1, len(pep) - length)
	outseq = pep[start:(start+length)]
	return(outseq)

uni2seq = {}
with open("../0.datadeal/eukaryota.whole.dataset") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		uni2seq[data[0]] = data[2]
		

out1 = open("../4.Analysis/4.2.motif/hs.dataset.highdelta.seq.fasta", "w")
out2 = open("../4.Analysis/4.2.motif/hs.dataset.highdelta.ctr.fasta", "w")
with open("../4.Analysis/4.2.motif/hs.dataset.highdelta.regions") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		tmpregion = data[1].split(';')
		peptide = uni2seq[data[0]]
		for i in tmpregion:
			out1.write(">"+data[0]+'_'+i+'\n'+regionseq(peptide, i)+'\n')
			out2.write(">"+data[0]+'_'+i+'\n'+randomseq(peptide, i)+'\n')
			
out1.close()

















