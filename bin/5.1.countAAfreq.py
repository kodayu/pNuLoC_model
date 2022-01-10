#!python

AA_list_k=["H","R","K","F","A","L","M","I","W","P","V","C","G","Q","N","S","Y","T","D","E","-"]
cor = {}
for i in AA_list_k:
	cor[i] = 0

with open("../0.datadeal/eukaryota.pro.nucleus.class") as f:
	for line in f:
		line = line.strip('\n')
		data = line.split('\t')
		if data[1] != "NonNucle":
			continue
		for i in range(len(data[2])):
			if data[2][i] not in AA_list_k:
				tmp = "-"
			else:
				tmp = data[2][i]
			cor[tmp] += 1

out = open("../0.datadeal/aa.stat.NonNucle.txt", "w")
for i in AA_list_k:
	out.write(i+'\t'+str(cor[i])+'\n')

out.close()






