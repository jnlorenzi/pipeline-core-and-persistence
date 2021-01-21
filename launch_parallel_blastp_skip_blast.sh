#!/bin/bash

# This script creates a blast database for each .fa files contained in the given directory. Then executes all possible blastp (v2.2.31+) and sets the results ( [query file]-vs-[subject file].bl ) in a new directory, './Blast_res/'.
# output format : tabular with the following column organization (see https://www.ncbi.nlm.nih.gov/books/NBK279684/): qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, gaps

# REQUIRE "Data" [= input sequences} and "Data_db" [= database to use]

########################
# Paths to modify before each run, the blast database is not created if already existing.

path_base=$1'/'  # working directory
path_db=$path_base'/Blast_db/'   # database dir
path_data_db=$path_base'/Data_db/'    # db input
path_data=$path_base'/Data/'$2'/'    # input data (.fa files)

########################

path_out=${path_db/db\//res/}
if [ ! -d $path_out ]; then
        mkdir $path_out
fi

format='7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps'
for set1 in $path_data*
do
        if [ ${set1: -3} == ".fa" ]; then
                temp=${set1/*\//}
                for set2 in $path_data_db*
                do
                        if [ ${set2: -3} == ".fa" ]; then
                                temp_set2=${set2/*\//}
                                #
                                if [ $temp == "Streptomyces_ambofaciens_ATCC_23877.fa" ] || [ $temp_set2 == "Streptomyces_ambofaciens_ATCC_23877.fa" ]; then
                                #
                                        out=$path_out${temp/.fa/}-vs-${temp_set2/.fa/.bl}
                                        db=$path_db${temp_set2/.fa/}
                                        if [ ! -f "$out" ]; then
                                                echo "IN: "$set1
                                                echo "DB: "$db
                                                echo "OUT: "$out
                                                echo "RUNNING BLAST..."
                                                echo ""
                                                blastp -query $set1 -db $db -out $out -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps'
                                        fi
                                #        
                                fi
                                #
                        fi
                done
        fi
done
