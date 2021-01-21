#!/bin/bash

#######################
# Paths to modify before each run, the blast database is not created if already existing.

path_base=$1'/'  # working directory
path_db=$path_base'/Blast_db/'   # database dir
path_data_db=$path_base'/Data_db/'    # db input

#######################


if [ ! -d $path_db ]; then
        mkdir $path_db
fi

for espece in $path_data_db*
do
        temp=${espece/*\//}
        makeblastdb -in $path_data_db$temp -dbtype "prot" -out $path_db${temp/.fa*/}
        echo ""
done


