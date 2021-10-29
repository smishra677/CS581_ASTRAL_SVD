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
	return dictionary_,tree



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

branch_collapse('out_try.tree','1','collapsed.tre')
dict0,tree_=get_polytomies('collapsed.tre')
#pprint.pprint(dictn)
dict1=select_polytomies(dict0,5)
#pprint.pprint(dict1)
for i in dict1.keys():
	send={}
	send[i]=dict1[i]
	if len(dict1[i])>=4:
		#pprint.pprint(send)
		#dict2= get_sequences('test.fasta',send)
		#print(dict2)
		#write_fatsa(dict2)
		#convert_fatsa_nexus('filtered.fasta','output.nex')
		#run_svd()
		regraft(tree_,send)


