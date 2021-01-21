This script generates from gbk files the core genome, the persistence and synteny score for each genome.

It must be executed in the directory containing the ".gbk" folder (gbk_coll_125).

Each step of the pipeline can be executed independently of the others by modifying the "operations to do" (true/false booleans) from line 703 to 716.

The 'nb_core' parameter line 697 (by default 50) indicates the number of cores to be used at the maximum during the Blast step.
The 'window_dimension' parameter line 698 (by default 1) indicates the dimension of window to use for synteny score (proporation of the genome, 1 = 1% of the total number of genes).

Mods requiered :

os
sys
pickle
math
re
Bio

Blast requirements :

Blast+
GNU Parallel