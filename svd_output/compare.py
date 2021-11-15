import dendropy
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

#f=open('stats_final_1_t.txt','w+')
#f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Method','dataset','nl', 'ei1', 'ei2','fp', 'fn', 'rf'))



print(compareTreesFromPath('astral_6_04_1000.tre','true-species_6_04.tre'))

print(compareTreesFromPath('svd_6_04_1000.tre','true-species_6_04.tre'))


print(compareTreesFromPath('out_final_6_4_1000.tre','true-species_6_04.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_05_1000.tre','true-species_6_05.tre'))

print(compareTreesFromPath('svd_6_05_1000.tre','true-species_6_05.tre'))


print(compareTreesFromPath('out_final_6_5_1000.tre','true-species_6_05.tre'))

print('----------------------------------------------')


print(compareTreesFromPath('astral_6_07_1000.tre','true-species_6_07.tre'))

print(compareTreesFromPath('svd_6_07_1000.tre','true-species_6_07.tre'))


print(compareTreesFromPath('out_final_6_7_1000.tre','true-species_6_07.tre'))

print('----------------------------------------------')


print(compareTreesFromPath('astral_6_11_1000.tre','true-species_6_11.tre'))

print(compareTreesFromPath('svd_6_11_1000.tre','true-species_6_11.tre'))


print(compareTreesFromPath('out_final_6_11_1000.tre','true-species_6_11.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_13_1000.tre','true-species_6_13.tre'))

print(compareTreesFromPath('svd_6_13_1000.tre','true-species_6_13.tre'))


print(compareTreesFromPath('out_final_6_13_1000.tre','true-species_6_13.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_16_1000.tre','true-species_6_16.tre'))

print(compareTreesFromPath('svd_6_16_1000.tre','true-species_6_16.tre'))


print(compareTreesFromPath('out_final_6_16_1000.tre','true-species_6_16.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_17_1000.tre','true-species_6_17.tre'))

print(compareTreesFromPath('svd_6_17_1000.tre','true-species_6_17.tre'))


print(compareTreesFromPath('out_final_6_17_1000.tre','true-species_6_17.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_20_1000.tre','true-species_6_20.tre'))

print(compareTreesFromPath('svd_6_20_1000.tre','true-species_6_20.tre'))


print(compareTreesFromPath('out_final_6_20_1000.tre','true-species_6_20.tre'))

print('----------------------------------------------')
print('----------------------------------------------')
print('----------------------------------------------')

print(compareTreesFromPath('astral_6_04_500.tre','true-species_6_04.tre'))

print(compareTreesFromPath('svd_6_04_500.tre','true-species_6_04.tre'))


print(compareTreesFromPath('out_final_6_4_500.tre','true-species_6_04.tre'))

print('----------------------------------------------')


print(compareTreesFromPath('astral_6_05_500.tre','true-species_6_05.tre'))

print(compareTreesFromPath('svd_6_05_500.tre','true-species_6_05.tre'))


print(compareTreesFromPath('out_final_6_5_500.tre','true-species_6_05.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_07_500.tre','true-species_6_07.tre'))

print(compareTreesFromPath('svd_6_07_500.tre','true-species_6_07.tre'))


print(compareTreesFromPath('out_final_6_7_500.tre','true-species_6_07.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_11_500.tre','true-species_6_11.tre'))

print(compareTreesFromPath('svd_6_11_500.tre','true-species_6_11.tre'))


print(compareTreesFromPath('out_final_6_11_500.tre','true-species_6_11.tre'))

print('----------------------------------------------')



print(compareTreesFromPath('astral_6_13_500.tre','true-species_6_13.tre'))

print(compareTreesFromPath('svd_6_13_500.tre','true-species_6_13.tre'))


print(compareTreesFromPath('out_final_6_13_500.tre','true-species_6_13.tre'))

print('----------------------------------------------')


print(compareTreesFromPath('astral_6_16_500.tre','true-species_6_16.tre'))

print(compareTreesFromPath('svd_6_16_500.tre','true-species_6_16.tre'))


print(compareTreesFromPath('out_final_6_16_500.tre','true-species_6_16.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_17_500.tre','true-species_6_17.tre'))

print(compareTreesFromPath('svd_6_17_500.tre','true-species_6_17.tre'))


print(compareTreesFromPath('out_final_6_17_500.tre','true-species_6_17.tre'))

print('----------------------------------------------')

print(compareTreesFromPath('astral_6_20_500.tre','true-species_6_20.tre'))

print(compareTreesFromPath('svd_6_20_500.tre','true-species_6_20.tre'))


print(compareTreesFromPath('out_final_6_20_500.tre','true-species_6_20.tre'))

print('----------------------------------------------')