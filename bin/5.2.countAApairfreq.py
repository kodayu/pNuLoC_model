#!python
import re
AA_list_k=["H","R","K","F","A","L","M","I","W","P","V","C","G","Q","N","S","Y","T","D","E","-"]
all_seq = []
with open("../4.Analysis/4.2.motif/hs.dataset.highdelta.seq.fasta") as f:
	for line in f:
		line = line.strip('\n')
		if re.search("^>", line):
			continue
		else:
			data = list(line)
			for i in range(len(data)):
				if data[i] not in AA_list_k:
					data[i] == "-"
			all_seq.append("".join(data))
ctr_seq = []
with open("../4.Analysis/4.2.motif/hs.dataset.highdelta.ctr.fasta") as f:
	for line in f:
		line = line.strip('\n')
		if re.search("^>", line):
			continue
		else:
			data = list(line)
			for i in range(len(data)):
				if data[i] not in AA_list_k:
					data[i] == "-"
			ctr_seq.append("".join(data))


all_seq = "_____".join(all_seq)
ctr_seq = "_____".join(ctr_seq)
##one AA
cor_seq, cor_ctr = {}, {}
for i in AA_list_k:
	cor_seq[i] = len(re.findall(i, all_seq))
	cor_ctr[i] = len(re.findall(i, ctr_seq))

out = open("../5.plots/sigdelta.aa.freq.txt", "w")
out.write("AA\tFreq\tFreqCtr\n")
for i in cor_seq:
	out.write(i+'\t'+str(cor_seq[i])+'\t'+str(cor_ctr[i])+'\n')
out.close()
##k-space AA 0
cor_seq, cor_ctr = {}, {}
for i in AA_list_k:
	for j in AA_list_k:
		cor_seq[i+j] = len(re.findall(r''+i+j+'', all_seq))
		cor_ctr[i+j] = len(re.findall(r''+i+j+'', ctr_seq))
out = open("../5.plots/sigdelta.aapair0.freq.txt", "w")
out.write("AA\tFreq\tFreqCtr\n")
for i in cor_seq:
	out.write(i+'\t'+str(cor_seq[i])+'\t'+str(cor_ctr[i])+'\n')
out.close()
##k-space AA 1
cor_seq, cor_ctr = {}, {}
for i in AA_list_k:
	for j in AA_list_k:
		cor_seq[i+j] = len(re.findall(r''+i+"."+j+'', all_seq))
		cor_ctr[i+j] = len(re.findall(r''+i+"."+j+'', ctr_seq))
out = open("../5.plots/sigdelta.aapair1.freq.txt", "w")
out.write("AA\tFreq\tFreqCtr\n")
for i in cor_seq:
	out.write(i+'\t'+str(cor_seq[i])+'\t'+str(cor_ctr[i])+'\n')
out.close()
##k-space AA 2
cor_seq, cor_ctr = {}, {}
for i in AA_list_k:
	for j in AA_list_k:
		cor_seq[i+j] = len(re.findall(r''+i+".."+j+'', all_seq))
		cor_ctr[i+j] = len(re.findall(r''+i+".."+j+'', ctr_seq))

out = open("../5.plots/sigdelta.aapair2.freq.txt", "w")
out.write("AA\tFreq\tFreqCtr\n")
for i in cor_seq:
	out.write(i+'\t'+str(cor_seq[i])+'\t'+str(cor_ctr[i])+'\n')
out.close()







