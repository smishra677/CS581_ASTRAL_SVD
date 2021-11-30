import os

sub_folder_name_6=['04','05','07','11','13','16','17','20']


sub_folder_name_7=['03','04','05','08','09','12','15','16','17']

li=['100','100','500','250','50']
x='0000'




def astral(input_,output):
	os.system('java -jar astral.5.7.8.jar -i'+input_+' -o'+output)


import time

f=open('C:\\Users\\smish\\OneDrive\\Documents\\GitHub\\CS581_ASTRAL_SVD\\timer.txt','w')	
start=time.time()
for k in li:
	for j in sub_folder_name_6:
		astral('RaxML_'+j+'_'+k+'.tre', 'C:\\Users\\smish\\OneDrive\\Documents\\GitHub\\CS581_ASTRAL_SVD\\svd_output\\astral_6_'+j+'_'+k+'.tre')

end= time.time()

f.write('astral, '+str(end-start))

f.close()