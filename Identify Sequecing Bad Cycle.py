# Sequencing-Reads
Python code for reading and basic sequencing of DNA reads
###I wrote these algorithms to locate the most frequent position of a poor base quality and  'N'. 


def phred33ToQ(qual):     ##this decodes base qualities from sequencer which are phred33 encoded into a numerical value
    return ord(qual) -33  #ord() converts character to ascii value

readFastQ(filename)        ## already defined this function in my Sequencing Read repository

                          ###Locate and store position of base qualities equal and less than 10 (my abitratry value)###

def badCycle(qual):
	position = []			      ##position of low qual

	for qualities in qual:

		for q in qualities:
			if phred33ToQ(q) <= 10:
				freq.append(phred33ToQ(q))            
				position.append(qualities.index(q))
			
	return position, freq
  
  location = badCycle(quals)
  import statistics
  statistics.mode(location)   

  
  ##stores position of all the 'N' bases and reports most frequent position
  def badCycle2(seqs):
	pos = []                       ##position of N in read
	for read in seqs:
		for i in range(len(read)):
			if read[i]== 'N':
				pos.append(i)             ##stores the index that N occurs
	return pos

position = badCycle2(human_seq)
statistics.mode(position)
