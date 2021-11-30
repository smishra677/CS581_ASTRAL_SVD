import os

location='C://Users//smish//Documents//Astral.5.7.8//Astral//Data_//alignments_miss//25tax-1000gen-0bps-500K-1E-6-rand//'

location1='C://Users//smish//Documents//Astral.5.7.8//Astral//Data_//gene_miss//25tax-1000gen-100bps-500K-1E-6-rand//'


folder_name=['25tax-1000gen-0bps-500K-1E-6-rand-miss']
sub_folder_name_6=['04','05','07','11','13','16','17','20']


sub_folder_name_7=['03','04','05','08','09','12','15','16','17']

li=['500','250','50']
x='0000'

for k in li:
	for j in sub_folder_name_6:
		for i in range(1,int(k)+1):
			file=(x[0:4-len(str(i))]+str(i))
			input_=location+j+'//'+file+'.fasta'
			output_='rax_ml_'+k+'___100'+j+'_'+file+'.tre'
			os.system('raxmlHPC -p 12345 -s "'+input_+'" -n "'+output_+'" -m GTRGAMMA')
