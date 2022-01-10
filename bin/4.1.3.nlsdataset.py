#!python

nls2site, nls2source = {}, {}
with open("../0.database/seqnls.uniprot.nls.txt") as f:
	for line in f:
		data = line.strip('\n').split('\t')
		if data[0] in nls2site:
			nls2site[data[0]].append(data[2]+'-'+data[3])
			nls2source[data[0]].append(data[4])
		else:
			nls2site[data[0]] = [data[2]+'-'+data[3]]
			nls2source[data[0]] = [data[4]]


def str2region(str1, regions):
	#example, "fafsfasvfg", [2-3, 5-7]
	str1list = list(str1)
	movesites = []
	for i in regions:
		sites = i.split('-')
		sites = [int(j) for j in sites]
		movesites += range(sites[0], sites[1]+1)
	tmp_list = []
	for i in range(len(str1list)):
		if i+1 in movesites:
			continue
		else:
			tmp_list.append(str1list[i])
	return("".join(tmp_list))

def randomtrunc(str1, length):
	minimal = 0
	maximun = len(str1) - length - 1
	import random
	start = random.randint(minimal, maximun)
	return(str1[:start]+str1[(start+length):])

out = open("../4.Analysis/4.1.hsTruncated/nls.truncated.dataset", "w")
with open("../4.Analysis/4.1.hsTruncated/human.dataset.scores") as f:
	out.write(f.readline().strip('\n')+'\tNLSregions\tSource\tTruncatedSequence\tRandomTruncated\n')
	for line in f:
		line = line.strip('\n')
		data = line.split('\t')
		if data[0] in nls2site:
			seq = str2region(data[2], nls2site[data[0]])
			ranseq = randomtrunc(data[2], 8)
			out.write(line+'\t'+';'.join(nls2site[data[0]])+'\t'+';'.join(nls2source[data[0]])+'\t'+seq+'\t'+ranseq+'\n')


out.close()

