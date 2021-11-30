import shutil


location='C://Users//smish//Documents//Astral.5.7.8//Astral//Data_//species_Alignments_miss//'
folder_name=['25tax-1000gen-100bps-500K-1E-7-rand-miss','25tax-1000gen-100bps-500K-1E-6-rand-miss']
sub_folder_name_6=['//04//','//05//','//07//','//11//','//13//','//16//','//17//','//20//']
sub_folder_name_7=['//03//','//04//','//05//','//08//','//09//','//12//','//15//','//16//','//17//']

file_astral=['astral-keep-miss-50.tre','astral-keep-miss-100.tre','astral-keep-miss-250.tre','astral-keep-miss-500.tre','astral-keep-miss-1000.tre']
file_svd=['svdquartets-keep-miss-50.tre','svdquartets-keep-miss-100.tre','svdquartets-keep-miss-250.tre','svdquartets-keep-miss-500.tre','svdquartets-keep-miss-1000.tre']

'''
for i in folder_name:
	for j in sub_folder_name_7:
		for k in file_astral:
			shutil.copy(location+i+j+k,'astral_7_'+str(j[2:4])+'_'+k.split('-')[3].split('.')[0]+'.tre')


for i in folder_name:
	for j in sub_folder_name_6:
		for k in file_astral:
			shutil.copy(location+i+j+k,'astral_6_'+str(j[2:4])+'_'+k.split('-')[3].split('.')[0]+'.tre')


for i in folder_name:
	for j in sub_folder_name_7:
		for k in file_svd:
			shutil.copy(location+i+j+k,'svd_7_'+str(j[2:4])+'_'+k.split('-')[3].split('.')[0]+'.tre')


for i in folder_name:
	for j in sub_folder_name_6:
		for k in file_svd:
			shutil.copy(location+i+j+k,'svd_6_'+str(j[2:4])+'_'+k.split('-')[3].split('.')[0]+'.tre')


'''
for i in folder_name:
	for j in sub_folder_name_6:
		shutil.copy(location+i+j+'true-species.tre','true-species_6_'+str(j[2:4])+'.tre')


for i in folder_name:
	for j in sub_folder_name_7:
		shutil.copy(location+i+j+'true-species.tre','true-species_7_'+str(j[2:4])+'.tre')
