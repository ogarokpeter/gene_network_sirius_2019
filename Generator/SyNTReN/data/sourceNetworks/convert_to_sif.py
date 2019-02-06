import sys


filename = sys.argv[1]
filename2 = sys.argv[2]

g1 = []
inter = []
g2 = []

all_genes = set()


with open(filename, "r") as f:
	for i, line in enumerate(f):
		ls = line.split('\t')
		if len(ls) == 2:
			all_genes.add(ls[0])
			g1.append(ls[0])
			inter.append('pp')
			g2.append(ls[1])
		else:
			print("fuck at line '", i, "'")

with open(filename2, 'w') as f:
	for i in range(len(g1)):
		f.write(g1[i] + '\t' + inter[i] + '\t' + g2[i])

print(all_genes, len(all_genes))