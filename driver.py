import os
from Bio import AlignIO

#raxml
#os.system('raxmlHPC -p 1 -s "test.fasta" -n "output12.tre" -m GTRGAMMA')

#astral
#os.system('java -jar astral.5.7.8.jar -i test_data/1KP-genetrees.tre -o a.tre')

#branch collapse
#os.system('nw_ed a.tre "i & b <=0" o > o.tre')


#convert to nexus
#AlignIO.convert("test.fasta", "fasta","nn.nex","nexus","DNA")

#svd check inst.txt for the commands
#os.system('paup4c instr.txt')
