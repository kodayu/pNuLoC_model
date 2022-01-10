##prepare human dataset, human N/M/C truncated dataset


hscor = {}
with open("../0.datadeal/uniprot.tab") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if data[5] == "9606":
			hscor[data[0]] = data[1]+'\t'+data[3]

out1 = open("../4.Analysis/4.1.hsTruncated/human.dataset", "w")
out2 = open("../4.Analysis/4.1.hsTruncated/human.truncated.dataset", "w")
with open("../0.datadeal/eukaryota.pro.nucleus.class") as f:
	out1.write('UniID\tLocInfo\tSeq\tUniName\tGeneName\n')
	out2.write('UniID\tLocInfo\tSeq\tUniName\tGeneName\tTruncatedTerminal\n')
	for line in f:
		line = line.strip('\n')
		data = line.split('\t')
		if data[0] in hscor:
			out1.write(line+'\t'+hscor[data[0]]+'\n')
			len_seq = len(data[2])
			point1 = int(len_seq/3)
			point2 = int(len_seq*2/3)
			out2.write(data[0]+'\t'+data[1]+'\t'+data[2][point1:]+'\t'+hscor[data[0]]+'\tN\n')
			out2.write(data[0]+'\t'+data[1]+'\t'+data[2][:point1]+data[2][point2:]+'\t'+hscor[data[0]]+'\tM\n')
			out2.write(data[0]+'\t'+data[1]+'\t'+data[2][:point2]+'\t'+hscor[data[0]]+'\tC\n')


out1.close()
out2.close()
		




