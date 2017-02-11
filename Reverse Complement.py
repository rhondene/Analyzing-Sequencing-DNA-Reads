### Outputs Reverse Complement of a Sequence
def revComp(t):
    complement={ 'A': 'T', 'G': 'C', 'T':'A', 'C':'G' , 'N':'N' }    #matches each base with its complement
    tc = ''                           ##stores reverse complement of sequence t
    for base in t:
        tc = complement[base] + tc
    return tc
