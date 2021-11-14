from Bio.Nexus import Nexus
import os
from Bio import AlignIO


x='0000'
file_list=[]
for i in range(1,501):
	file=(x[0:4-len(str(i))]+str(i))
	file='output_'+file+'.nex'
	file_list.append(file)
	#AlignIO.convert(file+'.fas', "fasta",'output_'+file+'.nex',"nexus","DNA")

nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open("COMBINED.nex", "w+"))
