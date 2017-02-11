# Sequencing-Reads
Python code for reading and basic sequencing of DNA reads
###Reads sequences from FastA File

def readGenome(fileFa):
  genome = ''
  with open(fileFa, 'r') as f:
      for line in f:           ##loops through each line in f
      if not line[0] == '>':  ##ignores the header symbol
           genome += line.rstrip()
  return genome      
