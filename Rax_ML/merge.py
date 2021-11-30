


import glob



folder_name=['25tax-1000gen-0bps-500K-1E-6-rand-miss']
sub_folder_name_6=['04','05','07','11','13','16','17','20']


sub_folder_name_7=['03','04','05','08','09','12','15','16','17']

li=['500','250','50']

for k in li:
	for j in sub_folder_name_6:
		read_files = glob.glob("RAxML_bestTree.rax_ml_"+k+"___100"+j+"_*.tre")




		with open("RaxML_"+j+"_"+k+".tre", "wb") as outfile:
		    for f in read_files:
		        with open(f, "rb") as infile:
		            outfile.write(infile.read())
