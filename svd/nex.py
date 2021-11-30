from Bio.Nexus import Nexus
import os
from Bio import AlignIO


location='C://Users//smish//Documents//Astral.5.7.8//Astral//Data_//alignments_miss//25tax-1000gen-0bps-500K-1E-6-rand//'


folder_name=['25tax-1000gen-0bps-500K-1E-6-rand-miss']
sub_folder_name_6=['04','05','07','11','13','16','17','20']


sub_folder_name_7=['03','04','05','08','09','12','15','16','17']

li=['100','1000','500','250','50']
x='0000'



import time

f=open('C:\\Users\\smish\\OneDrive\\Documents\\GitHub\\CS581_ASTRAL_SVD\\timer.txt','w')
start=time.time()
for k in li:
	for j in sub_folder_name_6:
		for i in range(1,int(k)+1):
			file=(x[0:4-len(str(i))]+str(i))
			input_=location+j+'//'+file+'.fas'
			file_list=[]
			file_list.append('output_'+file+'.nex')
			AlignIO.convert(input_, "fasta",'output_'+file+'.nex',"nexus","DNA")


		nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]
		combined = Nexus.combine(nexi)
		combined.write_nexus_data(filename=open("COMBINED_"+k+"_"+j+".nex", "w+"))

end= time.time()

f.write('SVDquartets, '+str(end-start))

f.close()