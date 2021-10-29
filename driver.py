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
	    if len(node._child_nodes) > 3:
	    	dictionary_[node]=[]
	    	for children in node._child_nodes:
	    		dictionary_[node].append(children)

	#pprint.pprint(dictionary_)
	return dictionary_,tree



def select_polytomies(dictionary_,Degree):
	#Selected_polytomies contains all the leaves under the subtree. It is a 2d dictionary[i][j] where i represents the 
	#orginal polytomy and j represents the head_node of the subtree. Input dictionary_ will have multiple polytomies from 
	#the original tree
	selected_polytomies={}
	for node in dictionary_.keys():
		selected_polytomies[node]={}
		for j in dictionary_[node]:
			if j.is_leaf():
				selected_polytomies[node][j]=[j]
			else:
				selected_polytomies[node][j]=[]
				for i in Tree(seed_node=j).postorder_node_iter():
					if i.is_leaf():
						selected_polytomies[node][j].append(i)

	#selected_polymtomie_D contains randomly selected polytomies. If the subtree contains less than D leaves than we select all
	#else we randomize our selection
	selected_polytomies_D={}
	for poly in selected_polytomies.keys():
		selected_polytomies_D[poly]=[]
		for head_ in selected_polytomies[poly].keys():
			if (len(selected_polytomies[poly][head_]))<=Degree:
				selected_polytomies_D[poly]= selected_polytomies_D[poly]+selected_polytomies[poly][head_]
			else:
				selected_polytomies_D[poly]=selected_polytomies_D[poly]+random.sample(selected_polytomies[poly][head_],Degree)

	pprint.pprint(selected_polytomies_D)
	return selected_polytomies_D



def get_sequences(input_,dictionary_):
	dictionary_input=SeqIO.to_dict(SeqIO.parse(input_,'fasta'))
	dictionary_output={}

	for j in dictionary_.keys():
		for k in dictionary_[j]:
			if str(k.taxon) !='None':
				tax=str(k.taxon).strip('"').strip("'")
				dictionary_output[k]=dictionary_input[tax]
			
	
	return dictionary_output

	#pprint.pprint(dictionary_input)


def write_fatsa(input_dictionary):

	with open('filtered.fasta', "w+") as handle:
		SeqIO.write(input_dictionary.values(),handle,'fasta')


#TODO add the resolved portion of the Tree back to the main Tree)
def regraft(tree_,dictionary_):
	tax = TaxonNamespace()
	tree = Tree.get(file=open('out_svd.tre', 'r'),
		            schema="newick",
		            tree_offset=0,
		            taxon_namespace=tax,
		            preserve_underscores=True)


	children_to_prune=[]
	children_to_prune_tree_id=[]
	for i in tree.postorder_node_iter():
		if str(i.taxon) !='None':
			tax=str(i.taxon).strip('"').strip("'")
			children_to_prune.append(tax)


	tree_.prune_taxa_with_labels(children_to_prune)

	node_=list(dictionary_.keys())
	for node in tree_.postorder_node_iter():
		if node in node_:
			node.add_child(tree)

	f=open('out_final.newick','w+')
	f.write(str(tree_))
	f.close()







'''
TODO Compare the tree and generate stats

'''


#Run RAxML
#raxml('test.fasta','output_raxaml.tre')


#astral('RAxML_result.output_raxaml.tre', 'astral.tre')

branch_collapse('out_try.tree','0.1','collapsed.tre')
dict0,tree_=get_polytomies('collapsed.tre')
dict1=select_polytomies(dict0,2)
#pprint.pprint(dict1)
for i in dict1.keys():
	send={}
	send[i]=dict1[i]
	if len(dict1[i])>=4:
		#pprint.pprint(send)
		dict2= get_sequences('test.fasta',send)
		print(dict2)
		write_fatsa(dict2)
		convert_fatsa_nexus('filtered.fasta','output.nex')
		#run_svd()
		regraft(tree_,send)


