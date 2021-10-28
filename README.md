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
