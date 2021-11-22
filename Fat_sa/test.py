from Bio.Nexus import Nexus
import os
from Bio import AlignIO


x='0000'
os. chdir('Fat_sa/')
file_list=[]
for i in range(1,1001):
	file=(x[0:4-len(str(i))]+str(i))
	file=file+'.nex'
	file_list.append(file)

nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open("C://Users//smish//OneDrive//Documents//GitHub//CS581_ASTRAL_SVD//COMBINED.nex", "w+"))
