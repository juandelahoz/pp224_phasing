import sys
import math

gtypes = ['01210','10222','00110']
#gtypes = ['01210']


def possHaps (g):
	pairs = []
	num_het = 0
	
	for snp in g:
		if (snp == '1'):
			num_het = num_het + 1
	ph = ['']*int(math.pow(2, num_het-1))

	i = 0
	j = 0
	het_passed = 1
	while (i < len(g)):
		j = 0
		if (g[i] == '0'):
			while (j < len(ph)):
				ph[j] = ph[j] + '0'
				j = j + 1
		elif (g[i] == '2'):
			while (j < len(ph)):
				ph[j] = ph[j] + '1'
				j = j + 1
		else:
			switch = math.pow(2, num_het-het_passed)
			on = True
			k = 0
			while (j < len(ph)):
				if (on == True):
					ph[j] = ph[j] + '1'
				else:
					ph[j] = ph[j] + '0'
				k = k + 1
				if (k == switch):
					k = 0
					if (on):
						on = False
					else:
						on = True
				j = j + 1
			het_passed = het_passed + 1
		i = i + 1
	for h in ph:
		comp = ''
		i = 0
		while (i < len(g)):
			if (g[i] == '0'):
				comp = comp + '0'
			elif (g[i] == '2'):
				comp = comp + '1'
			else:
				if (h[i] == '1'):
					comp = comp + '0'
				else:
					comp = comp + '1'
			i = i + 1
		pairs.append([h, comp])
	return pairs

total_haps = []
for g in gtypes:
	total_haps.append(possHaps(g))

print(total_haps)

def initialize_probs(listGhap):
	setH = set(listGhap)
	listH = list(setH)
	listpH = [float(1)/float(len(listH))]*len(listH)
	return listH,listpH



def m_step(listH,listpGn,listGhap,listpH,num_ppl):
	for i in range(0,len(listH)):
		haplotype = listH[i]
		matches = [i for i, e in enumerate(listGhap) if e == haplotype]
		freq_sum = 0.0
		for index in matches:
			freq_sum += listpGn[index]
		listpH[i] = freq_sum/(2.0*num_ppl)
	return listpH


num_ppl = max(listG)