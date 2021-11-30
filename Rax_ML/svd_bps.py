
import os
from Bio import AlignIO,SeqIO
import pprint
import random
import dendropy
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from dendropy import Tree, TaxonNamespace, Node

from Bio.Nexus import Nexus
import os
from Bio import AlignIO
import copy
import re



x='0000'

location='C://Users//smish//Documents//Astral.5.7.8//Astral//Data_//alignments_miss//25tax-1000gen-0bps-500K-1E-6-rand//'


sub_folder_name_6=['04','05','07','11','13','16','17','20']


sub_folder_name_7=['03','04','05','08','09','12','15','16','17']

li=['100','1000']
#,'500','250','50']
x='0000'

for k in li:
	for j in sub_folder_name_6:
		for i in range(1,int(k)+1):
			file=(x[0:4-len(str(i))]+str(i))
			#file_='output_'+file+'.nex'
			file=location+j+'//'+file
			dictionary_input=SeqIO.to_dict(SeqIO.parse(file+'.fas','fasta'))
			dictionary_output={}


			for i in dictionary_input.keys():
				dictionary_output[i]=dictionary_input[i][:100]
				
			with open(file+'.fasta', "w+") as handle:
				SeqIO.write(dictionary_output.values(),handle,'fasta')