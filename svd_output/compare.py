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

f=open('stats_final_1_t.txt','w+')
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Method','dataset','nl', 'ei1', 'ei2','fp', 'fn', 'rf'))





nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_04_1000.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_05_1000.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_07_1000.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_11_1000.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_13_1000.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_16_1000.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_17_1000.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_20_1000.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))





nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_04_1000.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_05_1000.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_07_1000.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_11_1000.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_13_1000.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_16_1000.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_17_1000.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_20_1000.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_04_1000.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_05_1000.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_07_1000.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_11_1000.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_13_1000.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_16_1000.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_17_1000.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_20_1000.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))








nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_03_1000.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_04_1000.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_05_1000.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_08_1000.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_09_1000.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_12_1000.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_15_1000.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_16_1000.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_17_1000.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_03_1000.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_04_1000.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_05_1000.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_08_1000.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_09_1000.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_12_1000.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_15_1000.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_16_1000.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_17_1000.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))










nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_03_1000.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_04_1000.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_05_1000.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_08_1000.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_09_1000.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_12_1000.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_15_1000.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_16_1000.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_17_1000.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','1000',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


'''

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_04_100.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_05_100.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_07_100.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_11_100.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_13_100.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_16_100.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_17_100.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_20_100.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))





nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_04_100.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_05_100.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_07_100.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_11_100.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_13_100.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_16_100.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_17_100.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_20_100.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_04_100.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_05_100.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_07_100.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_11_100.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_13_100.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_16_100.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_17_100.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_20_100.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))








nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_03_100.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_04_100.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_05_100.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_08_100.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_09_100.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_12_100.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_15_100.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_16_100.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_17_100.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_03_100.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_04_100.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_05_100.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_08_100.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_09_100.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_12_100.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_15_100.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_16_100.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_17_100.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))










nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_03_100.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_04_100.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_05_100.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_08_100.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_09_100.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_12_100.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_15_100.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_16_100.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_17_100.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','100',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))



'''

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_04_500.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_05_500.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_07_500.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_11_500.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_13_500.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_16_500.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_17_500.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_20_500.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))





nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_04_500.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_05_500.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_07_500.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_11_500.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_13_500.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_16_500.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_17_500.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_20_500.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_04_500.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_05_500.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_07_500.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_11_500.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_13_500.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_16_500.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_17_500.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_20_500.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))








nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_03_500.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_04_500.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_05_500.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_08_500.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_09_500.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_12_500.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_15_500.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_16_500.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_17_500.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_03_500.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_04_500.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_05_500.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_08_500.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_09_500.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_12_500.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_15_500.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_16_500.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_17_500.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))










nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_03_500.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_04_500.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_05_500.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_08_500.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_09_500.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_12_500.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_15_500.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_16_500.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_17_500.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','500',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))



'''
nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_04_50.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_05_50.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_07_50.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_11_50.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_13_50.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_16_50.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_17_50.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_20_50.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))





nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_04_50.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_05_50.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_07_50.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_11_50.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_13_50.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_16_50.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_17_50.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_20_50.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_04_50.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_05_50.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_07_50.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_11_50.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_13_50.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_16_50.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_17_50.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_20_50.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))








nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_03_50.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_04_50.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_05_50.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_08_50.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_09_50.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_12_50.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_15_50.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_16_50.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_17_50.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_03_50.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_04_50.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_05_50.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_08_50.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_09_50.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_12_50.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_15_50.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_16_50.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_17_50.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))










nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_03_50.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_04_50.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_05_50.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_08_50.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_09_50.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_12_50.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_15_50.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_16_50.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_17_50.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','50',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))




nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_04_250.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_05_250.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_07_250.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_11_250.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_13_250.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_16_250.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_17_250.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_6_20_250.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))





nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_04_250.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_05_250.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_07_250.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_11_250.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_13_250.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_16_250.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_17_250.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_6_20_250.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_04_250.tre','true-species_6_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_05_250.tre','true-species_6_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_07_250.tre','true-species_6_07.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_11_250.tre','true-species_6_11.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_13_250.tre','true-species_6_13.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_16_250.tre','true-species_6_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_17_250.tre','true-species_6_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_6_20_250.tre','true-species_6_20.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))








nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_03_250.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_04_250.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_05_250.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_08_250.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_09_250.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_12_250.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_15_250.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_16_250.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('astral_7_17_250.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('Astral','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))






nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_03_250.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_04_250.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_05_250.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_08_250.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_09_250.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_12_250.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_15_250.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_16_250.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('svd_7_17_250.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('svd','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))










nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_03_250.tre','true-species_7_03.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_04_250.tre','true-species_7_04.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_05_250.tre','true-species_7_05.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_08_250.tre','true-species_7_08.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_09_250.tre','true-species_7_09.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_12_250.tre','true-species_7_12.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_15_250.tre','true-species_7_15.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))


nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_16_250.tre','true-species_7_16.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))

nl, ei1, ei2, fp, fn, rf= list(compareTreesFromPath('out_final_7_17_250.tre','true-species_7_17.tre'))
f.write('{0:20}  {1:14}  {2:14}  {3:10}  {4:10}   {5:10} {6:10}  {7:10}\n'.format('out_final','250',str(nl), str(ei1), str(ei2), str(fp),str(fn),str(rf)))



'''




f.close()