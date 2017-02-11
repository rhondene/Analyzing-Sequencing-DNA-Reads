# Sequencing-Reads
Python code for reading and basic sequencing of DNA reads
##Reading FastQ file 
def readFastq(fileFa):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()                # stores dna sequence
            fh.readline()
            qual = fh.readline().rstrip()               # stores base quality phred33 encoded
            if len(seq) == 0:                           ## to detect when reach end of sequence
                break
            sequences.append(seq)                       # stores each line of dna sequences
            qualities.append(qual)                      #stores base qualities often phred33 encoded

    return sequences, qualities
