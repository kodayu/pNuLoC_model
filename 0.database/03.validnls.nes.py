#!python

out1 = open("./valid.nls.pos.txt", "w")
out2 = open("./valid.nes.pos.txt", "w")
with open("../8.webserver/inuloc.dataset.scores.cutoff.seqnls.uniprot.nlsdb.validnes.nesbase.uniprot.nlsdb") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if data[3] == "Homo sapiens (Human)":
			nls1, nls2 = data[-14].split(";"), data[-12].split(";")
			nes1, nes2, nes3 = data[-8].split(";"), data[-6].split(";"), data[-4].split(";")
			nlss = nls1+nls2
			nlss = [i for i in nlss if i != '']
			ness = nes1+nes2+nes3
			ness = [i for i in ness if i != '']
			if len(nlss) > 0:
				out1.write(data[0]+'\t'+';'.join(nlss)+'\n')
			if len(ness) > 0:
				out2.write(data[0]+'\t'+';'.join(ness)+'\n')

out1.close()
out2.close()		

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

def regionseq(pep, region):
	reg = region.split('...')
	outseq = pep[max((int(reg[0])-1), 0):min(int(reg[1]), len(pep))]
	return(outseq)

out11 = open("./valid.nls.pos.fasta", "w")
with open("./valid.nls.pos.txt") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		tmpregion = data[1].split(';')
		if data[0] not in pro2seq:
			continue
		peptide = pro2seq[data[0]]
		for i in tmpregion:
			out11.write(">"+data[0]+'_'+i+'\n'+regionseq(peptide, i)+'\n')
			
out11.close()

out12 = open("./valid.nes.pos.fasta", "w")
with open("./valid.nes.pos.txt") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		tmpregion = data[1].split(';')
		if data[0] not in pro2seq:
			continue
		peptide = pro2seq[data[0]]
		for i in tmpregion:
			out12.write(">"+data[0]+'_'+i+'\n'+regionseq(peptide, i)+'\n')
			
out12.close()



