#!python
#sys.argv[1] = hs.cpt.delta.method1.sigdelta.sites
import random
import sys
uni2region = {}

with open("../4.Analysis/4.1.hsTruncated/"+sys.argv[1]) as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		tmpuni = data[0]
		if tmpuni in uni2region:
			uni2region[tmpuni].append(int(data[4]))
		else:
			uni2region[tmpuni] = [int(data[4])]

out = open("../4.Analysis/4.1.hsTruncated/"+sys.argv[1]+".regions", "w")
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
	outseq = pep[max((int(reg[0])-1), 0):min(int(reg[1]), len(pep))]
	return(outseq)
	
def randomseq(pep, region):
	reg = region.split('...')
	length = min(int(reg[1]), len(pep)) - max((int(reg[0])-1), 0)
	start = random.randint(1, len(pep) - length)
	outseq = pep[start:(start+length)]
	return(outseq)

uni2seq = {}
with open("../0.datadeal/eukaryota.whole.dataset") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		uni2seq[data[0]] = data[2]
		

out11 = open("../4.Analysis/4.2.motif/"+sys.argv[1]+".nls.fasta", "w")
out12 = open("../4.Analysis/4.2.motif/"+sys.argv[1]+".nes.fasta", "w")
out13 = open("../4.Analysis/4.2.motif/"+sys.argv[1]+".unknown.fasta", "w")
out2 = open("../4.Analysis/4.2.motif/"+sys.argv[1]+".ctr.fasta", "w")
out31 = open("../4.Analysis/4.1.hsTruncated/"+sys.argv[1]+".nls.regions", "w")
out32 = open("../4.Analysis/4.1.hsTruncated/"+sys.argv[1]+".nes.regions", "w")
out33 = open("../4.Analysis/4.1.hsTruncated/"+sys.argv[1]+".unknown.regions", "w")
with open("../4.Analysis/4.1.hsTruncated/"+sys.argv[1]+".regions") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		tmpregion = data[1].split(';')
		peptide = uni2seq[data[0]]
		nlss, ness, unknowns = [], [], []
		for i in tmpregion:
			regseq = regionseq(peptide, i)
			hkrsum = regseq.count("H")+regseq.count("K")+regseq.count("R")
			charge = regseq.count("H")+regseq.count("K")+regseq.count("R") - regseq.count("D")-regseq.count("E")
			afilmvsum = regseq.count("A")+regseq.count("F")+regseq.count("I")+regseq.count("L")+regseq.count("M")+regseq.count("V")
			hydrosum = (regseq.count("A")+regseq.count("F")+regseq.count("L")+regseq.count("M")+regseq.count("I")+regseq.count("W")+regseq.count("P")+regseq.count("V"))/len(regseq)
			if (hkrsum > 1) & (charge > 0):
				nlss.append(i)
				out11.write(">"+data[0]+'_'+i+'\n'+regionseq(peptide, i)+'\n')
				out2.write(">"+data[0]+'_'+i+'\n'+randomseq(peptide, i)+'\n')
			elif (afilmvsum > 2) & (hydrosum > 0.3):
				ness.append(i)
				out12.write(">"+data[0]+'_'+i+'\n'+regionseq(peptide, i)+'\n')
				out2.write(">"+data[0]+'_'+i+'\n'+randomseq(peptide, i)+'\n')
			else:
				unknowns.append(i)
				out13.write(">"+data[0]+'_'+i+'\n'+regionseq(peptide, i)+'\n')
				out2.write(">"+data[0]+'_'+i+'\n'+randomseq(peptide, i)+'\n')
		if len(nlss) > 0:
			out31.write(data[0]+'\t'+';'.join(nlss)+'\n')
		if len(ness) > 0:
			out32.write(data[0]+'\t'+';'.join(ness)+'\n')
		if len(unknowns) > 0:
			out33.write(data[0]+'\t'+';'.join(unknowns)+'\n')
		
			
out11.close()
out12.close()
out13.close()
out2.close()
out31.close()
out32.close()
out33.close()
















