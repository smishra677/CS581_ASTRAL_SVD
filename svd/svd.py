import os


sub_folder_name_6=['04','05','07','11','13','16','17','20']




li=['100','1000']
#,'500','250','50']


x='0000'


for k in li:
	for j in sub_folder_name_6:
		f=open('instrution_'+k+'_'+j+'.txt','w+')
		input_="COMBINED_"+k+"_"+j+".nex"
		print(input_)
		instrution= "exe COMBINED_"+k+"_"+j+".nex;svd showScores=no evalQuartets=all qformat=qmc replace=no;savetrees file='out_svd_"+k+"_"+j+".tre'format='newick' replace='yes';quit;"
		f.write(instrution)
		f.close()







for k in li:
	for j in sub_folder_name_6:
		os.system('paup4c instrution_'+k+'_'+j+'.txt')
