#!python

all_pro_file = "../4.Analysis/4.1.hsTruncated/human.dataset.CodingProbTracy.scores"

pro2pos2score = {}
pro2pos2trunc = {}
pro2pos2delta = {}
with open(all_pro_file) as f:
	f.readline()
	for line in f:
		data = line.strip('\n').split('\t')
		if data[0] not in pro2pos2score:
			pro2pos2score[data[0]] = {}
		pro2pos2score[data[0]][int(data[4])] = float(data[5])


def calc_trunc(p2p2s, uniid, w=16):
	pos = p2p2s[uniid].keys()
	score = p2p2s[uniid].values()
	trunc, delta = {}, {}
	trunc[uniid] = {}
	delta[uniid] = {}
	for i in pos:
		lowlimit = int(max(1,i-(w/2)))
		higlimit = int(min(len(pos)+1,i+(w/2)+1))
		res = [p2p2s[uniid][j] for j in range(lowlimit, higlimit)]
		Strunc = sum(res)/(higlimit-lowlimit)
		trunc[uniid][i] = Strunc
	for i in pos:
		lowlimit = int(max(1,i-(w/2)))
		higlimit = int(min(len(pos)+1,i+(w/2)+1))
		delta[uniid][i] = trunc[uniid][higlimit-1] - trunc[uniid][lowlimit]
	return(trunc, delta)


for pro in pro2pos2score:
	tmp_trunc, tmp_delta = calc_trunc(pro2pos2score, pro, w=16)
	pro2pos2trunc = {**pro2pos2trunc, **tmp_trunc}
	pro2pos2delta = {**pro2pos2delta, **tmp_delta}


out = open("../4.Analysis/4.1.hsTruncated/human.dataset.CodingProbTracy.strunc.delta.scores", "w")
with open(all_pro_file) as f:
	out.write(f.readline().strip('\n')+'\tStrunc\tDelta\n')
	for line in f:
		line = line.strip('\n')
		data = line.split('\t')
		out.write(line+'\t'+str(pro2pos2trunc[data[0]][int(data[4])])+'\t'+str(pro2pos2delta[data[0]][int(data[4])])+'\n')

out.close()



