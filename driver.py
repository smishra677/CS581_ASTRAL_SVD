import os
from Bio import AlignIO,SeqIO
import pprint
import random

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from dendropy import Tree, TaxonNamespace, Node

#raxml
def raxml(input_,output):
	os.system('raxmlHPC -p 1 -s "'+input_+'" -n "'+output+'" -m GTRGAMMA')

#astral
def astral(input_,output):
	os.system('java -jar astral.5.7.8.jar -i'+input_+' -o'+output)

#branch collapse
def branch_collapse(input_,threshold,output):
	os.system('nw_ed '+input_+' "i & b <='+threshold+'" o >'+output)


#convert to nexus
def convert_fatsa_nexus(input_,output):
	AlignIO.convert(input_, "fasta",output,"nexus","DNA")

#svd check inst.txt for the commands
def run_svd():
	os.system('paup4c instr.txt')

#get polytomies
def get_polytomies(input_):
	#detect polytomies got from dendropy
	tax = TaxonNamespace()
	tree = Tree.get(file=open(input_, 'r'),
		            schema="newick",
		            tree_offset=0,
		            taxon_namespace=tax,
		            preserve_underscores=True)
	polytomies = []
	dictionary_={}
	for node in tree.postorder_node_iter():
	    if len(node._child_nodes) > 2:
	    	dictionary_[node]=[]
	    	for children in node._child_nodes:
	    		dictionary_[node].append(children)

	#pprint.pprint(dictionary_)
	return dictionary_



def select_polytomies(dictionary_,Degree):
	selected_polytomies={}

	for node in dictionary_.keys():
		if (len(dictionary_[node]))<=Degree:
			selected_polytomies[node]=dictionary_[node]
		else:
			selected_polytomies[node]=random.sample(dictionary_[node],Degree)

	return selected_polytomies



def get_sequences(input_,dictionary_):
	dictionary_input=SeqIO.to_dict(SeqIO.parse(input_,'fasta'))
	dictionary_output={}

	for j in dictionary_.keys():
		for k in dictionary_[j]:
			if str(k.taxon) !='None':
				tax=str(k.taxon)[1:-1]
				dictionary_output[tax]=dictionary_input[tax]
			
	
	return dictionary_output

	#pprint.pprint(dictionary_input)


def write_fatsa(input_dictionary):

	with open('filtered.fasta', "w+") as handle:
		print(SeqIO.write(input_dictionary.values(),handle,'fasta'))


'''
TODO add the resolved portion of the Tree back to the main Tree)
def regraft(tree,dictionary_):
'''


'''
TODO Compare the tree and generate stats

'''


#Run RAxML
#raxml('test2.fasta','output_raxaml.tre')


#astral('RAxML_result.output_raxaml.tre', 'astral.tre')

branch_collapse('astral.tre','0.1','collapsed.tre')
dict0=get_polytomies('collapsed.tre')
#pprint.pprint(dictn)
dict1=select_polytomies(dict0,5)

for i in dict1.keys():
	send={}
	send[i]=dict1[i]
	print(len(dict1[i]))
	if len(dict1[i])>=4:
		print(send)
		dict2= get_sequences('test.fasta',send)
		print(dict2)
		write_fatsa(dict2)
		convert_fatsa_nexus('filtered.fasta','output.nex')
		run_svd()
		#regraft(tree,dic)


