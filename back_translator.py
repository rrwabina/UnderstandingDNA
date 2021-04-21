import sys,string,mod_translate,mod_seqfiles


#Full IUPAC alphabet for peptide and DNA
alphaDNA = "ACGTRYMKWSBDHVN" 
alphaDNA += alphaDNA.lower()

alphaPEP = "ARNDCQEGHILKMFPSTWYVBZX*"
alphaPEP += alphaPEP.lower()
	 

def degap(s,gapin):
	result = []
	for c in s:
		if not c in gapin: result.append(c)
	return string.join(result,"")

	# +++++++  Start Modifications +++++++

def initialgaps(s,gapin): # count gaps at the beginning of protein sequence
	i = 0
	for c in s:
		if c in gapin: 
			i = i + 1
		else:
			break
	return i
	
	# ++++++++ End Modifications ++++++++

def trim(s,alphabet):
	result = []
	for c in s:
		if c in alphabet: result.append(c)
	return string.join(result,"")
	
def trimseqs(seqs,alphabet):
	for key in seqs.keys():
		s,n = seqs[key]
		#s = s.upper()
		seqs[key] = (trim(s,alphabet) , n)
	
def matchtrans(dnaseqs,pepseqs,gapin,verbose,mtx,allinternal,readthroughstop):
	dnaref = {}
	dnaref_extra = {}
	result = {}
	
	# NOTICE:     We need the handle the situation where more than one DNA
	#             sequence translates to the same peptide sequence.
	#
	# ASSUMPTION: Identical peptide sequences will align exactly the same
	#             way. 
	#
	#
	# EXAMPLE:
	#
	#             dnaSeq17 -> pepSeq17
	#             dnaSeq32 -> pepSeq32   
	#
	#             *) dnaSeq17 and dnaSeq32 differs by a few nucleotides
	# 
	#             *) pepSeq17 and pepSeq32 are exactly the same
	#
	#             Given the assmuption mentioned above, it does NOT
	#             matter if dnaSeq17 gets paired with pepSeq32
	 
	for key in dnaseqs.keys():
		dna,note = dnaseqs[key]
		dna = degap(dna,gapin)
		newpep = mod_translate.translate(dna,mtx,not allinternal,readthroughstop)
		
		# Strip terminal stop-codon
		if newpep.endswith("*"):
			newpep = newpep[:-1]
		
		if verbose > 2:
			warn("DNA sequence "+key+" translated to:\n"+newpep);
		
		if dnaref.has_key(newpep):
			dnaref_extra[key] = newpep
		else:
			dnaref[newpep] = key
	
    for key in pepseqs.keys(): 
        pep, note = pepseqs[key]
        pep = degap(pep,gapin).upper()

		# Strip terminal stop-codon
		if pep.endswith("*"):
			pep = pep[:-1]
		
		if verbose > 2:
			warn("Pep sequence "+key+" degapped: \n"+pep);
		
		if dnaref.has_key(pep):
			result[key] = dnaref.pop(pep)
		else:
			for dnakey in dnaref_extra.keys():
				if pep == dnaref_extra[dnakey]:
					result[key] = dnakey
					dnaref_extra.pop(dnakey)
					break
					
	return result

def matchname(dnaseqs,pepseqs):
	result = {}
	for key in pepseqs.keys():
		if key in dnaseqs.keys():
			result[key] = key
	return result
	
def matchpos(dnaseqs,pepseqs):
	result = {}
	dnakeys = dnaseqs.keys()
	i = 0
	for key in pepseqs.keys():
		if i < len(dnakeys):
			result[key] = dnakeys[i]
			i += 1
	return result
		
def revtrans(dnaseqs,pepseqs,crossref,gapin,gapout,verbose):
	if verbose:
		warn("gapin: '"+gapin+"'")
		warn("gapout: '"+gapout+"'")
	
	newdnaseqs = {}
	error = 0
	for key in pepseqs.keys():
		try:
			# Find the corresponding sequences
			dna, pep, newdna = "","",""  # Just in the case of an exception
			if not key in crossref.keys():
				warn("No cross-reference, skipping peptide sequence "+key)
				continue
#			print key,crossref[key]
			dna,noted = dnaseqs[crossref[key]]
#			dnaName = d_dnames[dna]
			dnaName = crossref[key]
			
			#print dna
			dna = degap(dna,gapin)
			pep,notep = pepseqs[key]
			newdna = ""
			dnap = 0

			# +++++++  Start Modifications +++++++
	
			newpep = ""
			degapped = degap(pep,gapin)
			dnap = -3 * initialgaps(pep,gapin)  # correct start if pep starts with gap characters

			newpep = mod_translate.translate(dna,None,True,False)
			# Strip terminal stop-codon
			if newpep.endswith("*"):
				newpep = newpep[:-1]

			offset = string.find(newpep, degapped) # offset start of rev translation
			if offset < 0:
				warn("Could not match pep:"+key)
			dnap = dnap + 3 * offset

			# +++++++  End Modifications +++++++

				
			# Do the reverse translation for this seq			
			l_dna = []
			for i in range(0,len(pep)):
				c = pep[i]
				if c in gapin:
					l_dna.append(gapout * 3)
				else:
					# Extract codon - keep case from the amino acid
					codon = dna[dnap:dnap+3]
					if c.isupper():
						codon = codon.upper()
					else:
						codon = codon.lower()
						
					l_dna.append(codon)
					dnap = dnap +3
			
			# Everything's cool - add the new seq to the result
			newdna = string.join(l_dna,"")
#			newdnaseqs[key] = (newdna,noted)
			newdnaseqs[dnaName] = (newdna,noted)
		except:
			if verbose:
				warn("Error rev-translating seq:"+key)
				warn("\nLen dna:"+str(len(dna))+" pep:"+str(len(pep))+" newdna:"+str(len(newdna))+"\n")
			error = error +1
	return (newdnaseqs,error)
	
def argerr(arg):
	warn("Error:\nThe parameter "+arg+" must be followed by a value.\n")
	sys.exit(1)
	
def warn(msg):
	sys.stderr.write(msg+"\n")		

def main():
	# Set defaults
	verbose  = 0
	gapin    = "-.~"
	gapout   = "-"
	outfile  = ""
	Idna     = "auto"
	Ipep     = "auto"
	outform  = "fasta"
	matchmet = "trans"
	mtx_file = ""
	mtx      = None
	allinternal = False
	readthroughstop = False
	
	# Quick sanity check
	if len(sys.argv)<3:
		print __doc__
		sys.exit(1)
		
	# Process arguments
	dnafile = ""
	pepfile = ""
	
	argv = sys.argv[1:]
	while (len(argv)>0):
		arg = argv[0]
		
		if arg == "-h" :
			print __doc__
			sys.exit(0)
		
		if arg == "-v"     : verbose = 1
		if arg == "-vv"    : verbose = 2
		if arg == "-vvv"   : verbose = 3
		
		if arg == "-match" :
			if len(argv) == 0: argerr("-match")
			matchmet = argv[1]
			argv = argv[2:]
			continue

		if arg == "-gapin" :
			if len(argv) == 0: argerr("-gapin")
			gapin = argv[1]
			argv = argv[2:]
			continue
			
		if arg == "-gapout" :
			if len(argv) == 0: argerr("-gapout")
			gapout = argv[1][0]			# Use only the first char
			argv = argv[2:]
			continue
	
		if arg == "-Idna" :
			if len(argv) == 0: argerr("-Idna")
			Idna = (argv[1]).lower()
			argv = argv[2:]
			continue
			
		if arg == "-Ipep" :
			if len(argv) == 0: argerr("-Ipep")
			Ipep = (argv[1]).lower()
			argv = argv[2:]
			continue	
				
		if arg == "-O" :
			if len(argv) == 0: argerr("-O")
			outform = (argv[1]).lower()
			argv = argv[2:]
			continue
			
		if arg == "-mtx" :
			if len(argv) == 0: argerr("-mtx")
			mtx_file = argv[1]
			argv = argv[2:]
			continue
			
		if arg == "-allinternal":
			allinternal = True
			
		if arg == "-readthroughstop":
			readthroughstop = True			

		if arg[0] != "-":
			if   dnafile == "" : dnafile = arg
			elif pepfile == "" : pepfile = arg
			else               : outfile = arg						
		
		argv = argv[1:]
	
	# Output extra info if requested
	if verbose:
		warn("verbose level: "+str(verbose))
		warn("dnafile:    "+dnafile+" [format:"+Idna+"]")
		warn("pepfile:    "+pepfile+" [format:"+Ipep+"]")
		if outfile : 
			warn("outfile: "+outfile)
		else       : 
			warn("outfile: None - writing to STDOUT")
		warn("out format:  "+outform)
		warn("Mtxfile:     "+mtx_file)
	
	# Read input files
	try:	
		if Idna == "auto": 
			Idna = mod_seqfiles.autotype(dnafile)
			if verbose: warn("DNA file format appears to be......: "+Idna)
			
		dnaseqs = mod_seqfiles.readfile(dnafile,Idna)
		if verbose: warn("#DNA entries read: "+str(len(dnaseqs)))

		if Ipep == "auto": 
			Ipep = mod_seqfiles.autotype(pepfile)
			if verbose: warn("Peptide file format appears to be..: "+Ipep)

		pepseqs = mod_seqfiles.readfile(pepfile,Ipep)
		if verbose: warn("#pep entries read: "+str(len(pepseqs)))
		
		#if 1:
		if len(dnaseqs) == 0 or len(pepseqs) == 0:
			warn("Error: Bad input.")
			warn("Dna sequences read: "+str(len(dnaseqs)))
			warn("peptide sequence read: "+str(len(pepseqs)))
			sys.exit(1)
			
		if mtx_file:
			try:
				mtx = mod_translate.parseMatrixFile(mtx_file)	
			except:
				warn("Invalid translation matrix: "+mtx_file)
				sys.exit(1)

	except Exception,msg:
		warn("Error reading input files: "+str(msg))
		warn("DNA File type: "+Idna)
		warn("Pep File type: "+Ipep)
		sys.exit(1)
		
	# Remove illegal characters
	trimseqs(dnaseqs,alphaDNA)
	trimseqs(pepseqs,alphaPEP + gapin)

	# Bit of extra info? 
	if verbose > 1:
		warn("DNA names ["+str(len(dnaseqs.keys()))+"] :")
		for key in dnaseqs.keys(): warn(key)
		warn("..")
		
		warn("pep names ["+str(len(pepseqs.keys()))+"] :")
		for key in pepseqs.keys(): warn(key)
		warn("..")

	if verbose:
		warn("Matching DNA and peptide sequences by: "+matchmet)
	
	# Establish cross-references
	if   matchmet == "name":  crossref = matchname(dnaseqs,pepseqs)
	elif matchmet == "pos":   crossref = matchpos(dnaseqs,pepseqs)
	elif matchmet == "trans": crossref = matchtrans(dnaseqs,pepseqs,gapin,verbose,mtx,allinternal,readthroughstop)
	else:
		warn('Match method "'+matchmet+'" not known.')
		sys.exit(1) 
	
	if len(crossref.keys()) <> len(dnaseqs.keys()) <> len(dnaseqs.keys()):
		warn ("Warning: Not all DNA and peptide sequences could be matched.")
		warn (str(len(crossref.keys())+" crossreference(s) could be established."))
	
	if verbose > 1:
		warn("Cross references: (Pep/DNA) sequence names")
		for key in crossref.keys():
			warn(key+" / "+crossref[key]) 

	# Do the reverse translation
	newdnaseqs, error = revtrans(dnaseqs,pepseqs,crossref,gapin,gapout,verbose)
	
	if verbose: warn("#rev-trans DNA entries: "+str(len(newdnaseqs.keys())))
	
	if (error > 0):
		warn ("# errors:"+str(error)+" aborting ...")
		sys.exit(1)
		
	# Output the result
	try:
		if outfile != "": out_stream = open(outfile,"w")
		else:             out_stream = sys.stdout
		
		mod_seqfiles.writestream(out_stream,newdnaseqs,outform,"N")
		
	except Exception, msg:
		warn("Failed to write output."+str(msg))
		if outfile: warn("outfile: "+outfile)
		sys.exit(1)
		
if __name__ == "__main__":
	main()