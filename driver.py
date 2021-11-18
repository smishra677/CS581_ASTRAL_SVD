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



	#pprint.pprint(tree)
	return dictionary_,tree



def select_polytomies(tree_,dictionary_,Degree):
	#Selected_polytomies contains all the leaves under the subtree. It is a 2d dictionary[i][j] where i represents the 
	#orginal polytomy and j represents the head_node of the subtree. Input dictionary_ will have multiple polytomies from 
	#the original tree
	for i in tree_.postorder_node_iter():
		print(i)
	print('-----------------------------------------------------')
	tree_x=dendropy.Tree(tree_)
	nds3 = [nd for nd in tree_.postorder_node_iter()]  
	nds4 = [nd for nd in tree_x.postorder_node_iter()] 
	copy_translation={} 
	for i, n in enumerate(nds3):                    
	    copy_translation[n]=nds4[i]

	selected_polytomies={}
	for node in dictionary_.keys():
		selected_polytomies[node]={}
		for j in dictionary_[node]:
			if j.is_leaf():
				selected_polytomies[node][j]=[j]
			else:
				selected_polytomies[node][j]=[]
				new_tree=Tree(seed_node=j)
				for i in new_tree.postorder_node_iter():
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

	#pprint.pprint(selected_polytomies_D)
	for i in tree_x.postorder_node_iter():
		print(i)
	return selected_polytomies_D,tree_x,selected_polytomies,copy_translation



def write_nexus():
	os.system('python ./Fat_sa/test.py')




def get_sequences(file__,dictionary_,m):
	x='0000'
	
	file_list=[]
	for i in range(1,int(m)+1):
		file=(x[0:4-len(str(i))]+str(i))
		file_='output_'+file+'.nex'
		dictionary_input=SeqIO.to_dict(SeqIO.parse(file__+file+'.fas','fasta'))
		dictionary_output={}

	
		for j in dictionary_.keys():
			for k in dictionary_[j]:
				if str(k.taxon) !='None':
					tax=str(k.taxon).strip('"').strip("'")
					if tax not in dictionary_input.keys():
						continue
					else:
						dictionary_output[k]=dictionary_input[tax]

		with open(file__+file+'.fasta', "w+") as handle:
				SeqIO.write(dictionary_output.values(),handle,'fasta')

		AlignIO.convert(file__+file+'.fasta', "fasta",file__+file+'.nex',"nexus","DNA")
		

	write_nexus()


	#pprint.pprint(dictionary_input)


def write_fatsa(input_dictionary):

	with open('filtered.fasta', "w+") as handle:
		SeqIO.write(input_dictionary.values(),handle,'fasta')


#TODO add the resolved portion of the Tree back to the main Tree)
def regraft(tree_,dictionary_,dictionary_1,copy_translation):
	tax = TaxonNamespace()
	tree = Tree.get(file=open('out_svd.tre', 'r'),
		            schema="newick",
		            tree_offset=0,
		            taxon_namespace=tax,
		            preserve_underscores=True)

	taxon_to_remove=set()
	for node in tree.postorder_node_iter():
		tax_svd='581'
		if str(node.taxon)!='None':
			tax_svd=str(node.taxon).strip('"').strip("'")
		for key in dictionary_1.keys():
				for nodes in dictionary_1[key]:
					if str(nodes.taxon) !='None':
						tax=str(nodes.taxon).strip('"').strip("'")
						if tax_svd==tax:
							new_tree=Tree(seed_node=key).clone()
							node.taxon.label=None
							if node.is_leaf():
								if tax not in taxon_to_remove:
									node.add_child(new_tree)
									lis_=set([str(i.taxon).strip('"').strip("'") for i in dictionary_1[key] if str(i.taxon) !='None'])
									taxon_to_remove=taxon_to_remove.union(set(lis_))
								else:
									del node
				



	pprint.pprint(dictionary_1)
	print('----------------')
	print(tree)
	print('----------------')
	for node in dictionary_.keys():
		try:
			for nd in tree_.postorder_node_iter():
				if str(copy_translation[node])==str(nd):
					nd.add_child(tree)
		except:
			continue

	for node in dictionary_1.keys():
		print(node)

	print('-----------------------------------------')
	for node in dictionary_1.keys():
		try:
			for nd in tree_.postorder_node_iter():
				if str(copy_translation[node])==str(nd):
					if len(nd._child_nodes)==0:
						del nd
						break
					else:
						print(nd._child_nodes)
						nd.remove_child(nd._child_nodes)
		except:
			continue

	
	print(tree_)

	
	'''
	node_=list(dictionary_.keys())
	print(node_)
	key = node_[0]

	children_to_prune=[]
	children_to_prune_tree_id=[]
	for i in tree.postorder_node_iter():
		if str(i.taxon) !='None':
			tax=str(i.taxon).strip('"').strip("'")
			children_to_prune.append(tax)

	print(children_to_prune)
	tree_.prune_taxa_with_labels(children_to_prune)
    
	node_=list(dictionary_.keys())
	
	
	for node in tree_.postorder_node_iter():
		if str(node) == str(key):
			print('aaaa')
			node.add_child(tree)
			break
'''
	return tree_

def write_tree(tree_,k,j,l):
	f=open('./svd_output//out_final_'+k+'_'+j+'_'+l+'.tre','w+')
	f.write(str(tree_)+';')
	f.close()



'''
TODO Compare the tree and generate stats


'''

def compareTreesFromPath(treePath1, treePath2):
    print("Comparing {} with {}".format(treePath1, treePath2))
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=treePath1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)
    tr2 = dendropy.Tree.get(path=treePath2,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    return compareDendropyTrees(tr1, tr2)
    #print("RF distance on %d shared leaves: %d" % (nl, fp + fn))
def compareDendropyTrees(tr1, tr2):
    from dendropy.calculate.treecompare \
        import false_positives_and_negatives
    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])
    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)
        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)
        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)
    tr1.update_bipartitions()
    tr2.update_bipartitions()
    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))
    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = float(fp + fn) / (ei1 + ei2)
    return (nl, ei1, ei2, fp, fn, rf)




#Run RAxML
#raxml('0001.fas','output_raxaml.tre')
#astral('raxml-genes.tre', 'astral.tre')
#


def main(astral,file,k,j,l):
	print(astral,file)
	branch_collapse(astral,'0.4','collapsed_'+k+'_'+j+'_'+l+'.tre')
	dict0,tree_=get_polytomies('collapsed'+j+'.tre')


	dict1,tree_,selected_polytomies,copy_translation=select_polytomies(tree_,dict0,1)
	print(dict1)
	for i in dict1.keys():
		send={}
		#print(selected_polytomies)
		send[i]=dict1[i]
		print(len(send[i]))
		if len(send[i])>=4:
			#pprint.pprint(send)
			print('00000')
			get_sequences(file,send,l)
			run_svd()
			tree_=regraft(tree_,send,selected_polytomies[i],copy_translation).clone()

	write_tree(tree_,k,j,l)

	#print(compareTreesFromPath('true-species.tre','astral.tre'))

	#print(compareTreesFromPath('true-species.tre','out_final.tre'))

location='C://Users//smish//Documents//Astral.5.7.8//Astral//Data_//alignments_miss//25tax-1000gen-0bps-500K-1E-6-rand//'


folder_name=['25tax-1000gen-0bps-500K-1E-6-rand-miss']
sub_folder_name_6=['04','05','07','11','13','16','17','20']


sub_folder_name_7=['03','04','05','08','09','12','15','16','17']

li=['100','1000','500','250','50']

for k in li:
	for j in sub_folder_name_6:
		main('./svd_output//astral_6_'+j+'_'+k+'.tre',location+j+'//','6',j,k)


location='C://Users//smish//Documents//Astral.5.7.8//Astral//Data_//alignments_miss//25tax-1000gen-0bps-500K-1E-7-rand//'


for k in li:
	for j in sub_folder_name_7:
		main('./svd_output//astral_7_'+j+'_'+k+'.tre',location+j+'//','7',j,k)


