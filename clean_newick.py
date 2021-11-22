import re


file=open('collapsed_6_07_100.tre')
data=file.read()



NewickTree = [data]

file.close()


print(NewickTree)
pattern = re.compile(r'[,]+([^;:]+)\b')

branch_lengths=[]

for tree in NewickTree:
    branch_lengths = pattern.findall(tree)
    # Do stuff to the list branch_lengths

print(branch_lengths)

new_branch=[]
for i in branch_lengths:
	split=i.split(')')

	if len(split)>=1:
		new_branch.append(''.join(split[:-1])+')')
	else:
		new_branch.append(')')



print(new_branch)

print(','.join(new_branch))
