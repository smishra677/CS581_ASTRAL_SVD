# CS581_ASTRAL_SVD



### Introduction

ASTRAL is one of the most cited and most well-known species tree estimation methods in the phylogenetic space. With the increasing optimalities shown in ASTRAL-II and ASTRAL-III, the
race to optimize this tree estimation method has become increasingly focused on optimizing to
its weaknesses rather than strengths \cite{1,2}. Recent findings have shown that ASTRAL tends to show
the most decrease in performance under high ILS and gene tree estimation error (GTEE)
conditions \cite{8, 7}. Thus in solving this problem, it is important to note that another species tree estimation method, SVDquest, is shown to be much more accurate than summary methods in
the high GTEE dataset. SVDquest is another species tree estimation method which runs a
two-phase algorithm that uses SVDquartets implementation in PAUP∗ to compute a set of
quartet trees and uses a quartet amalgamation method on the full set of species \cite{5}. SVDquest
is also shown to be more accurate than ASTRAL in high GTEE conditions but ASTRAL is shown to
have a much higher accuracy in low GTEE datasets \cite{5}.
For this project, we are purposing to create a new pipeline of ASTRAL and SVDquartets+ PAUP*. The main motivation is
to improve ASTRAL by using SVDquartets+ PAUP* on the subsets of ASTRAL that tend to be the least
confidence on, i.e., low support edges, becomes the primary focus for our course project.


### Possible Improvements

In order to improve ASTRAL, we can use ASTRAL as a base method and build off of an initial
tree estimation. We will use ASTRAL-III as the base species tree estimation method and then
create a heuristic threshold marker to signal low support values for branches within this
estimated tree. In doing preliminary research on the iterations of ASTRAL, ASTRAL-III provides
two main improvements compared to ASTRAL-II. The first improvement involves constraining
the bipartition set X in order to grow linearly with the number of species and input genes. The
second, and more consequential to our project proposal, improvement that ASTRAL-III has
compared to ASTRAL-II is the ability to better handle polytomies using techniques to constrain
the optimal tree space and use similarities between gene trees. With the improved run-time,
ASTRAL-III has the feature of removing low-support branches in order to improve accuracy on
estimation.

This feature allows us to move on to the next step of the method and remove low-support
branches thus creating polytomies on the estimated tree output from ASTRAL-III. Using some
threshold to collapse branches (ASTRAL-III paper mentions anything below 10\%) \cite{3}, we will go
through the output tree and remove said branches. The final step involves running SVDquartets+ PAUP* on
the polytomies found from the previous step and resolving those polytomies by running
SVDquartets+ PAUP* on some subset of leaves within those polytomies. Our method will be tested on
datasets with a high gene tree estimation error (GTEE) as ASTRAL-II suffers in accuracy of
datasets with this characteristic.

### Pipeline:

## Data

The actual data can be found at 
https://github.com/smishra677/CS581_ASTRAL_SVD/tree/main/input/25tax-1000gen-0bps-500K-1E-6-rand

The replicates used are 05, 11,13,16, and 20.

## svd_bps.py
There are 1000 alignment files(001-1000) in each replicate folder. These files are in .fas format. We then use the svd_bps.py, which can be found in 
https://github.com/smishra677/CS581_ASTRAL_SVD/tree/main/Rax_ML
svd_bps.py goes over all the 1000 fas files and selects the first 100 sites from the taxon in those files. 
It then creates new alignment files with these 100 sites. The alignment files are named the same 001-1000 but have .fasta formatting.

## rax.py
After this, we run rax.py, which runs raxml for all the alignment files created by svd_bps.py.
rax.py and all the outputs created raxML can be found in  https://github.com/smishra677/CS581_ASTRAL_SVD/tree/main/Rax_ML


## merge.py
After running the raxml we use merge.py to merge all the gene tree files created by rax.py.  The merge happens for a replicate and the model condition.
merge.py can be found in https://github.com/smishra677/CS581_ASTRAL_SVD/tree/main/Rax_ML


## astral.py
Now we run astral using astral.py.  This file can be found in https://github.com/smishra677/CS581_ASTRAL_SVD/tree/main/Rax_ML.


## nex.py
Since we are done with astral we now use nex.py to convert all the fasta files created by nex.py to .nex file and to merge all the nex files concerning a replicate and the model condition. 

## svd.py
We have created a combined nex file; we use svd.py, which first creates all the instructions.txt for the SVDQuartets and runs it according to the instruction. We can find nex.py, svd.py https://github.com/smishra677/CS581_ASTRAL_SVD/tree/main/svd.

## drive.py
Our method, USA, is implemented in driver.py, where it runs on three threshold sizes [0.35,0.4,0.8], and from each subtree, we are selecting one taxon.  Here Newick utility is called to collapse the branches with lower support values than the threshold, and then all the polytomies are extracted, and we run svdQuartets to resolve these polytomies.
driver.py can be found at https://github.com/smishra677/CS581_ASTRAL_SVD.


## display_results.py
The display_results.py script creates the graphical display of all results found in the svd_outputs/FinalData.csv file. Results can be grouped in comparing across methods and across threshold values, as well as comparing runtime results across methods as well. The display_results.py script can be found at https://github.com/smishra677/CS581_ASTRAL_SVD/blob/main/display_results.py.



## comparision.py
There is comparision.py which compares the output of our method with the true species tree. 
comparision.py uses dendropy to compare this and can be comparision.py at https://github.com/smishra677/CS581_ASTRAL_SVD/tree/main/svd_output  








