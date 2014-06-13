Selaginella_DGE
===============

Differential Gene Expression Analyses on Selaginella

Files
-----

All_by_all:
clust -t 10000 -dt 10 -o All_by_All Selag_counts.txt 

By_leaf:
clust -t 10000 -dt 10 -o By_leaf -g [[0,4,7,11],[1,5,8,12],[2,6,9,13],[3,10,14]] Selag_counts.txt

By_species:
clust -t 10000 -dt 10 -o By_species -g [[0,1,2,3],[4,5,6],[7,8,9,10],[11,12,13,14]] Selag_counts.txt 

GetPvals.py: script to calculate p vlues for all comparisons

Selag_counts.txt: count data for expression analyses


