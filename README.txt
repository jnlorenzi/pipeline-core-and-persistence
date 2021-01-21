This script generates from gbk file the core genome and the persistence score for each genome contains in the repository.

It must be executed in the directory containing the ".gbk" folder (gbk_coll_125).

Each step of the pipeline can be executed independently of the others by modifying the "operations to do" (true/false booleans) from line 589 to 598.

The 'nb_core' parameter line 583 (by default 50) indicates the number of cores to be used at the maximum during the Blast step.

The data given in "gbk_coll_125" and "Blast" are a test set. 
The data gbk_col_125 is the input data. The Blast folder contains the Blast outputs on this data set so that you don't have to run it during the test.