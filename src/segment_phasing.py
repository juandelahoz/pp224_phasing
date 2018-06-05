import sys
import math

################################
    ###    FUNCTIONS:    ###

def possHaps (g):
	pairs = []
	num_het = 0
	
	for snp in g:
		if (snp == '1'):
			num_het = num_het + 1
	if (num_het > 0):
		ph = ['']*int(math.pow(2, num_het-1))
	else:
		ph = ['']

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

def threeLists (segs):
	listGhap = []
	listG = []
	listG_h = []
	a = 1 # this is G
	b = 1 # this is G_h
	for g in segs:
		#listGhap.append(possHaps(g))
		b = 1
		for pair in possHaps(g):
			for h in pair:
				listGhap.append(h)
				listG.append(a)
				listG_h.append(b)
			b = b + 1
		a = a + 1

	return listG, listG_h, listGhap

def initialize_probs(listGhap):
	setH = set(listGhap)
	listH = list(setH)
	listpH = [float(1)/float(len(listH))]*len(listH)
	return listH, listpH

def e_step(listG, listGh, listGhap, listH, listpH):

	listpG     = [0] * len(listGhap)    # to store P(h_i |G)
	listpPhase = [0] * len(listGhap)    # to store p(hihj|G)
	for i in range(len(listGhap)):      # look in all possible phases
		h = listGhap[i]
		for j in range(len(listH)):
			if h == listH[j]:           # find each haplotype
				listpG[i] = listpH[j]   # store its probability
				break
		if (i % 2) == 1:
			listpPhase[i-1] = listpPhase[i] = listpG[i-1] * listpG[i]

	listpGn    = [0] * len(listGhap)
	sumpPhase  = 0
	br_bc = 0
	currG = 1

	for i in range(len(listG)):
		if listG[i] == currG:          # while inside each G
			br_fw = i                  # advance forward bracket
			sumpPhase += listpPhase[i] # sum over all probs inside G
		else:                          # normalize for all h in G
			sumpPhase /= 2             # remove double counts (hi + hj)
			for k in range( br_bc, br_fw+1):  # look at h inside each G
				listpGn[k] = listpPhase[k] / sumpPhase

			# update
			currG = listG[i]
			sumpPhase = listpPhase[i]
			br_bc = i

	# normalize for all h in last G
	sumpPhase /= 2
	for k in range( br_bc, br_fw+1):
		listpGn[k] = listpPhase[k] / sumpPhase

	return listpGn

def m_step(listH, listpGn, listGhap, listpH, num_ppl):
	for i in range(0,len(listH)):
		haplotype = listH[i]
		matches = [i for i, e in enumerate(listGhap) if e == haplotype]
		freq_sum = 0.0
		for index in matches:
			freq_sum += listpGn[index]
		listpH[i] = freq_sum/(2.0*num_ppl)
	return listpH

def run_em(listH, listG, listGh, listGhap, listpH, num_iter):
	for i in range(0,num_iter):
		#print("E-step:")
		listpGn = e_step(listG, listGh, listGhap, listH, listpH)
		#print([round(x,2) for x in listpGn])
		#print("M-step:")
		listpH = m_step(listH, listpGn, listGhap, listpH, max(listG))
		#print([round(x,2) for x in listpH])
	listpGn = e_step(listG, listGh, listGhap, listH, listpH)
	return listpGn,listpH

def select_Haps(listG, listGhap, listpGn):

	listMLHGs = [] # most likely haplotypes per genotype
	MLphase  = []  # most likely phase for current genotype
	maxPhase = 0   # likelihood of current phase
	currG = 1

	for i in range(len(listG)):
		if listG[i] == currG:                 # while inside each G
			if (i%2) == 1:                    # every two haplotypes
				if maxPhase < listpGn[i]:     # compare with the past prob
					maxPhase = listpGn[i]     # save highest prob
					MLphase  = [listGhap[i-1],listGhap[i]]  # select hapls
		else:
			listMLHGs.append(MLphase)          # add the best haplotypes
			# update
			currG = listG[i]
			MLphase  = []
			maxPhase = 0

	# add the best haplotypes
	listMLHGs.append(MLphase)

	return listMLHGs

def hamming (s1, s2):
	j = 0
	count = 0
	while (j < len(s1)):
		if (s1[j] != s2[j]):
			count = count + 1
		j = j + 1
	return count

def maxHaplotype (listMLHG, newHapPairs, overlap):
	i = 0 # i is individual we are looking at
	while (i < len(listMLHG)):
		if (len(listMLHG[i][0]) == 0): # case where this is first haplotype added to the list
			listMLHG[i][0] = newHapPairs[i][0]
			listMLHG[i][1] = newHapPairs[i][1]
		else:
			if ((listMLHG[i][0])[-overlap:] == (newHapPairs[i][0])[:overlap]):
				listMLHG[i][0] = listMLHG[i][0] + (newHapPairs[i][0])[overlap:]
				listMLHG[i][1] = listMLHG[i][1] + (newHapPairs[i][1])[overlap:]
			elif ((listMLHG[i][0])[-overlap:] == (newHapPairs[i][1])[:overlap]):
				listMLHG[i][0] = listMLHG[i][0] + (newHapPairs[i][1])[overlap:]
				listMLHG[i][1] = listMLHG[i][1] + (newHapPairs[i][0])[overlap:]
			else:
				if (hamming(listMLHG[i][0][-overlap:], newHapPairs[i][0][:overlap]) + 
					hamming (listMLHG[i][1][-overlap:], newHapPairs[i][1][:overlap])) < (hamming(listMLHG[i][1][-overlap:], newHapPairs[i][0][:overlap]) + 
					hamming (listMLHG[i][0][-overlap:], newHapPairs[i][1][:overlap])):
					listMLHG[i][0] = listMLHG[i][0] + (newHapPairs[i][0])[overlap:]
					listMLHG[i][1] = listMLHG[i][1] + (newHapPairs[i][1])[overlap:]
				else:
					listMLHG[i][0] = listMLHG[i][0] + (newHapPairs[i][1])[overlap:]
					listMLHG[i][1] = listMLHG[i][1] + (newHapPairs[i][0])[overlap:]
		i = i + 1
	return listMLHG

def write_haplotypes(listMLHG, file):
	outfile = open(file, "w")
	for i in range(0,len(listMLHG[0][0])):
		for j in range(0,len(listMLHG)-1):
			outfile.write(listMLHG[j][0][i] + ' ' + listMLHG[j][1][i] + ' ')
		outfile.write(listMLHG[j+1][0][i] + ' ' + listMLHG[j+1][1][i] + '\n')
	outfile.close()

def driver (spot, end):
	segments = []
	if (end == True):
		new_len = len_genome - (spot)
		for gtype in full_genotypes:
			segments.append(gtype[spot:spot+new_len])
	else:
		for gtype in full_genotypes:
			segments.append(gtype[spot:spot+seg_len])
	return segments

################################
    ###    VARIABLES:    ###

input_file = open(sys.argv[1], 'r')
step = int(sys.argv[2])
overlap = int(sys.argv[3])
seg_len = step + overlap
output_file = sys.argv[1][:-4] + "_sol" + ".txt"
individuals = 50
full_genotypes = [''] * individuals
listMLHG = [['',''] for i in range(individuals)]

if int(seg_len) > 14:
	sys.exit("We can't handle segments larger than 14: \
Please, choose a smaller step (" + sys.argv[2] + ") and/or overlap (" + sys.argv[3] + ")")

for line in input_file: # creates full_genotypes
	line = line.strip().split(' ')
	i = 0
	while (i < individuals):
		full_genotypes[i] = full_genotypes[i] + line[i]
		i = i + 1

len_genome = len(full_genotypes[0])
spot = 0   # to walk over the genome

##########################
    ###    RUN!    ###
end = False
while (spot + step) < len_genome:
	if ((spot + seg_len) >= len_genome):
		end = True
	# get segments and advance spot
	segs = driver(spot, end)
	# generate all the needed lists from the genotypes
	listG, listGh, listGhap = threeLists(segs)
	listH, listpH = initialize_probs(listGhap)

	if(end == False):
		print("Phasing..." + "\tpositions: " + str(spot) + "-" + str(spot+seg_len) + 
		"\ttotal_Haplt: " + str(len(listGh)) + "\tunique_Haplt: "  + str(len(listH)))
	else:
		print("Phasing..." + "\tpositions: " + str(spot) + "-" + str(len_genome) + 
		"\ttotal_Haplt: " + str(len(listGh)) + "\tunique_Haplt: "  + str(len(listH)))
		spot = len_genome

	# run the EM algorithm for the current segment
	listpGn, listpH = run_em(listH, listG, listGh, listGhap, listpH, 3)

	# select the most likely haplotypes from the EM probs and the full list
	listMLHGs = select_Haps(listG, listGhap, listpGn)
	# extend the current genotypes 
	listMLHG = maxHaplotype(listMLHG, listMLHGs, overlap)
	spot += step


#print(listMLHG)

write_haplotypes(listMLHG, output_file)

