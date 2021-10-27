import os
from Bio import AlignIO
import pprint


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from dendropy import Tree, TaxonNamespace

#raxml
#os.system('raxmlHPC -p 1 -s "test.fasta" -n "output12.tre" -m GTRGAMMA')

#astral
#os.system('java -jar astral.5.7.8.jar -i test_data/1KP-genetrees.tre -o a.tre')

#branch collapse
#os.system('nw_ed a.tre "i & b <=0.8" o > o.tre')


#convert to nexus
#AlignIO.convert("test.fasta", "fasta","nn.nex","nexus","DNA")

#svd check inst.txt for the commands
#os.system('paup4c instr.txt')
tax = TaxonNamespace()
tree = Tree.get(file=open('o.tre', 'r'),
                schema="newick",
                tree_offset=0,
                taxon_namespace=tax,
                preserve_underscores=True)

#detect polytomies got from dendropy
polytomies = []
dictionary_={}
for node in tree.postorder_node_iter():
    if len(node._child_nodes) > 2:
    	dictionary_[node]=[]
    	for children in node._child_nodes:
    		dictionary_[node].append(children)

pprint.pprint(dictionary_)


'''
dictionary_={}
for node in polytomies:

	for children in node._child_nodes:
		dictionary_[node].append(children)

pprint.pprint(dictionary_)
'''