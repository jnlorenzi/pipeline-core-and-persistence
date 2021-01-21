This script generates from gbk file the core genome and the persistence score for each genome contains in the repository.

It must be executed in the directory containing the ".gbk" folder (gbk_coll_125).

Each step of the pipeline can be executed independently of the others by modifying the "operations to do" (true/false booleans) from line 589 to 598.

The 'nb_core' parameter line 583 (by default 50) indicates the number of cores to be used at the maximum during the Blast step.