#!/usr/bin/env python3

#### MODS

import os, sys
import pickle
#import pandas
import math
import re
from Bio import SeqIO
from Bio import Seq

def main_chromosome(path_to_gbk_file):
        """ Returns the main chromosome (the larger) from a genbank file 
        """
        
        main = ''
        for gb_record in SeqIO.parse(open(path_to_gbk_file, 'r'), 'genbank'):
                if not main:
                        main = gb_record
                elif len(main) < len(gb_record):
                        main = gb_record
        return main.name


def int_to_roman(input):
        """
        Convert an integer to Roman numerals.

        Examples:
        >>> int_to_roman(0)
        Traceback (most recent call last):
        ValueError: Argument must be between 1 and 3999

        >>> int_to_roman(-1)
        Traceback (most recent call last):
        ValueError: Argument must be between 1 and 3999

        >>> int_to_roman(1.5)
        Traceback (most recent call last):
        TypeError: expected integer, got <type 'float'>

        >>> for i in range(1, 21): print int_to_roman(i)
        ...
        I
        II
        III
        IV
        V
        VI
        VII
        VIII
        IX
        X
        XI
        XII
        XIII
        XIV
        XV
        XVI
        XVII
        XVIII
        XIX
        XX
        >>> print int_to_roman(2000)
        MM
        >>> print int_to_roman(1999)
        MCMXCIX
        """
   
        if type(input) != type(1):
                raise(TypeError, "expected integer, got %s" % type(input))
        if not 0 < input < 4000:
                raise(ValueError, "Argument must be between 1 and 3999")  
        ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
        nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
        result = ""
        for i in range(len(ints)):
                count = int(input / ints[i])
                result += nums[i] * count
                input -= ints[i] * count
        return result


def position_generator(replicon, feature):
        """ Generates an stantard position identifiant from a replicon and a feature (gbk type)
        """
        strand = str(feature.location)[-2]
        start, end = re.sub('[^0-9:]', '', str(feature.location).split(',')[0]).split(':')
        start = str(int(start) + 1)
        pos = '-'.join([start, end])
        position = '_'.join([ replicon, strand, pos ])
        return position


def position_extractor(CDS_list, annotation):
        '''
        returns the position of each CDS contained in 'CDS_list'. The sequence is researched in all the genome (chromosome + plasmid)
        '''
        pos = {}
        for CDS_id in CDS_list:
                for organism in sorted(annotation):
                        for replicon_type in annotation[organism]:
                                for replicon in annotation[organism][replicon_type]:
                                        if 'CDS' in annotation[organism][replicon_type][replicon] and CDS_id in annotation[organism][replicon_type][replicon]['CDS']:
                                                pos[CDS_id] = annotation[organism][replicon_type][replicon]['CDS'][CDS_id]['position']
                                                continue
        return pos


def short_name_generator(organism_name, dico_short_name):
        """ Generates an unique short name from a an organism name. The dictionary with all the short name is required.
        """
        base = organism_name[0] + str(organism_name.split('_')[1])[0:2]
        short_name = base + str(sum((1 for x in dico_short_name if re.match(base, x))) + 1)
        return short_name


def gbk_to_pickle(path_to_gbk_folder, path_to_short_name = ''):
        """ Generates a dictionary with all the annotation from a gbk file and a dictionary with all the short name.   -> /!\ 2 output <-
        
        !!! Deals only with genome with stricly one chromosome and n plasmids !!!
        
        dictionary structure:
        
                [organism name][replicon type][replicon][sequence type] => specifics fields depending of the sequence type
                
                > organism name: name of organism in collection like "Streptomyces_coelicolor_A3(2)"
                > replicon type: "chromosome" or "plasmid"
                > replicon: "chr" for chromosome or "pls" for plasmid and a number in roman format, like "chrI"
                > sequence type: "source" / "CDS" / "tRNA" / "rRNA". "source" contains dimensions of the current replicon only
                
                [...][source] = replicon_strain_start position-end position. Like "chrI_+_0-8667507"
                
                                    /[position]
                [...][CDS][CDS_name]-[product]
                                    \[translation]
                
                                  /[position]
                [...][tRNA / rRNA]
                                  \[feature]
        Examples:
                >>> annotation['Streptomyces_coelicolor_A3(2)']['chromosome']['chrI']['CDS']['Sco1_1']['position']
                'chrI_+_576-699'
                
                >>> annotation['Streptomyces_coelicolor_A3(2)']['chromosome']['chrI']['CDS']['Sco1_1']['translation']
                'MPASEADDLVAGAQCEVVELDGMAPASRGVLCSRTAGTKA'
                
                >>> annotation['Streptomyces_coelicolor_A3(2)']['chromosome']['chrI']['tRNA']['Sco1_tRNA_1']['product']
                'tRNA-Val-CAC'
                
                >>> annotation['Streptomyces_coelicolor_A3(2)']['plasmid']['plsI']['source']
                'plsI_+_0-31317'

        """
        
        print('\n# .gbk parsing...')
        if path_to_short_name and os.path.exists(path_to_short_name):
                dico_short_name = pickle.load(open(path_to_short_name, 'rb'))
        else:
                dico_short_name = {}
        annotation = {}
        for gbk_file in sorted(os.listdir(path_to_gbk_folder)):
                        if gbk_file.endswith('.gbk'):
                                organism = re.sub('\.gbk$', '', gbk_file)
                                print(organism)
                                short_name = short_name_generator(organism, dico_short_name)
                                dico_short_name[short_name] = organism
                                
                                if os.stat(path_to_gbk_folder + gbk_file).st_size == 0: # If the file is empty, only use the name to construct the short name and skip the annotation extraction
                                        continue
                                main = main_chromosome(path_to_gbk_folder + gbk_file)
                                annotation[organism] = {}
                                annotation[organism]['chromosome'] = {}
                                annotation[organism]['plasmid'] = {}
                                CDS_number, tRNA_number, rRNA_number = 0, 0, 0
                                for gb_record in SeqIO.parse(open(path_to_gbk_folder + gbk_file, 'r'), 'genbank'):
                                        if gb_record.name == main:
                                                replicon = 'chr' + int_to_roman(len(annotation[organism]['chromosome']) + 1)
                                                annotation[organism]['chromosome'][replicon] = {}
                                                
                                                for feature in gb_record.features:
                                                        position = position_generator(replicon, feature)
                                                        handler = position.split('_')
                                                        orientation = handler[1]
                                                        begin, end = handler[-1].split('-')
                                                        begin = int(begin) - 1
                                                        if feature.type not in annotation[organism]['chromosome'][replicon]:
                                                                annotation[organism]['chromosome'][replicon][feature.type] = {}
                                                        if feature.type == 'source':
                                                                annotation[organism]['chromosome'][replicon][feature.type] = position
                                                        elif feature.type == 'CDS':
                                                                CDS_number += 1
                                                                CDS_id = short_name + '_' + str(CDS_number)
                                                                annotation[organism]['chromosome'][replicon][feature.type][CDS_id] = {}
                                                                annotation[organism]['chromosome'][replicon][feature.type][CDS_id]['position'] = position
                                                                annotation[organism]['chromosome'][replicon][feature.type][CDS_id]['translation'] = feature.qualifiers['translation'][0]
                                                                annotation[organism]['chromosome'][replicon][feature.type][CDS_id]['product'] = feature.qualifiers['product'][0]
                                                                if orientation == '+':
                                                                        annotation[organism]['chromosome'][replicon][feature.type][CDS_id]['sequence'] = str(gb_record.seq[int(begin):int(end)])
                                                                else:
                                                                        annotation[organism]['chromosome'][replicon][feature.type][CDS_id]['sequence'] = str(Seq.reverse_complement(gb_record.seq[int(begin):int(end)]))
                                                        elif feature.type == 'tRNA':
                                                                tRNA_number += 1
                                                                tRNA_id = short_name + '_tRNA_' + str(tRNA_number)
                                                                annotation[organism]['chromosome'][replicon][feature.type][tRNA_id] = {}
                                                                annotation[organism]['chromosome'][replicon][feature.type][tRNA_id]['position'] = position
                                                                annotation[organism]['chromosome'][replicon][feature.type][tRNA_id]['product'] = feature.qualifiers['product'][0]
                                                                if orientation == '+':
                                                                        annotation[organism]['chromosome'][replicon][feature.type][tRNA_id]['sequence'] = str(gb_record.seq[int(begin):int(end)])
                                                                else:
                                                                        annotation[organism]['chromosome'][replicon][feature.type][tRNA_id]['sequence'] = str(Seq.reverse_complement(gb_record.seq[int(begin):int(end)]))
                                                        
                                                        elif feature.type == 'rRNA':
                                                                rRNA_number += 1
                                                                rRNA_id = short_name + '_rRNA_' + str(rRNA_number)
                                                                annotation[organism]['chromosome'][replicon][feature.type][rRNA_id] = {}
                                                                annotation[organism]['chromosome'][replicon][feature.type][rRNA_id]['position'] = position
                                                                annotation[organism]['chromosome'][replicon][feature.type][rRNA_id]['product'] = feature.qualifiers['product'][0]
                                                                if orientation == '+':
                                                                        annotation[organism]['chromosome'][replicon][feature.type][rRNA_id]['sequence'] = str(gb_record.seq[int(begin):int(end)])
                                                                else:
                                                                        annotation[organism]['chromosome'][replicon][feature.type][rRNA_id]['sequence'] = str(Seq.reverse_complement(gb_record.seq[int(begin):int(end)]))
                                                        else:
                                                                if feature.type != 'fasta_record':
                                                                        print('\r' + feature.type, end = '')
                                        else:
                                                
                                                replicon = 'pls' + int_to_roman(len(annotation[organism]['plasmid']) + 1)
                                                annotation[organism]['plasmid'][replicon] = {}
                                                
                                                for feature in gb_record.features:
                                                        position = position_generator(replicon, feature)
                                                        handler = position.split('_')
                                                        orientation = handler[1]
                                                        begin, end = handler[-1].split('-')
                                                        begin = int(begin) - 1
                                                        
                                                        if feature.type not in annotation[organism]['plasmid'][replicon]:
                                                                annotation[organism]['plasmid'][replicon][feature.type] = {}
                                                        
                                                        if feature.type == 'source':
                                                                annotation[organism]['plasmid'][replicon][feature.type] = position
                                                        elif feature.type == 'CDS':
                                                                CDS_number += 1
                                                                CDS_id = short_name + '_' + str(CDS_number)
                                                                annotation[organism]['plasmid'][replicon][feature.type][CDS_id] = {}
                                                                annotation[organism]['plasmid'][replicon][feature.type][CDS_id]['position'] = position
                                                                annotation[organism]['plasmid'][replicon][feature.type][CDS_id]['translation'] = feature.qualifiers['translation'][0]
                                                                annotation[organism]['plasmid'][replicon][feature.type][CDS_id]['product'] = feature.qualifiers['product'][0]
                                                                if orientation == '+':
                                                                        annotation[organism]['plasmid'][replicon][feature.type][CDS_id]['sequence'] = str(gb_record.seq[int(begin):int(end)])
                                                                else:
                                                                        annotation[organism]['plasmid'][replicon][feature.type][CDS_id]['sequence'] = str(Seq.reverse_complement(gb_record.seq[int(begin):int(end)]))
                                                                
                                                        
                                                        elif feature.type == 'tRNA':
                                                                tRNA_number += 1
                                                                tRNA_id = short_name + '_tRNA_' + str(tRNA_number)
                                                                annotation[organism]['plasmid'][replicon][feature.type][tRNA_id] = {}
                                                                annotation[organism]['plasmid'][replicon][feature.type][tRNA_id]['position'] = position
                                                                annotation[organism]['plasmid'][replicon][feature.type][tRNA_id]['product'] = feature.qualifiers['product'][0]
                                                                if orientation == '+':
                                                                        annotation[organism]['plasmid'][replicon][feature.type][tRNA_id]['sequence'] = str(gb_record.seq[int(begin):int(end)])
                                                                else:
                                                                        annotation[organism]['plasmid'][replicon][feature.type][tRNA_id]['sequence'] = str(Seq.reverse_complement(gb_record.seq[int(begin):int(end)]))
                                                        
                                                        elif feature.type == 'rRNA':
                                                                rRNA_number += 1
                                                                rRNA_id = short_name + '_rRNA_' + str(rRNA_number)
                                                                annotation[organism]['plasmid'][replicon][feature.type][rRNA_id] = {}
                                                                annotation[organism]['plasmid'][replicon][feature.type][rRNA_id]['position'] = position
                                                                annotation[organism]['plasmid'][replicon][feature.type][rRNA_id]['product'] = feature.qualifiers['product'][0]
                                                                if orientation == '+':
                                                                        annotation[organism]['plasmid'][replicon][feature.type][rRNA_id]['sequence'] = str(gb_record.seq[int(begin):int(end)])
                                                                else:
                                                                        annotation[organism]['plasmid'][replicon][feature.type][rRNA_id]['sequence'] = str(Seq.reverse_complement(gb_record.seq[int(begin):int(end)]))
                                                        
                                                        else:
                                                                if feature.type != 'fasta_record':
                                                                        print('\r' + feature.type, end = '')
        print('\nDone\n')
        return annotation, dico_short_name


def CDS_fasta_file_generator(path_out, annotation, main_chromosome = False):
        """ Generates a fasta file (in path out) whith all the CDS sequence for each organism contains in the annotation dictionary or from simple dictionary structure : annotation[organism][CDS_id]['translation']
        """
        print("# .fa creation")
        if not os.path.exists(path_out):
               os.makedirs(path_out)
        for organism in annotation:
                print(organism)
                with open(path_out + organism + '.fa', 'w') as out_file:
                        if 'chromosome' not in annotation[organism] and 'plasmid' not in annotation[organism]:
                                for CDS_id in sorted(annotation[organism], key = lambda x: int(x.split('_')[1])):
                                                out_file.write('>' + CDS_id + '\n' + annotation[organism][CDS_id]['translation'] + '\n')
                        else:
                                out_list = []
                                for chromosome_id in annotation[organism]['chromosome']:
                                        if 'CDS' in annotation[organism]['chromosome'][chromosome_id].keys():
                                                for CDS_id in sorted(annotation[organism]['chromosome'][chromosome_id]['CDS'], key = lambda x: int(x.split('_')[1])):
                                                        out_list.append('>' + CDS_id + '\n' + annotation[organism]['chromosome'][chromosome_id]['CDS'][CDS_id]['translation'] + '\n')
                                if not main_chromosome:
                                        if 'plasmid' in annotation[organism]:
                                                for contig_id in annotation[organism]['plasmid']:
                                                        if 'CDS' in annotation[organism]['plasmid'][contig_id]:
                                                                for CDS_id in sorted(annotation[organism]['plasmid'][contig_id]['CDS'], key = lambda x: int(x.split('_')[1])):
                                                                        out_list.append('>' + CDS_id + '\n' + annotation[organism]['plasmid'][contig_id]['CDS'][CDS_id]['translation'] + '\n')
                                for l in sorted(out_list, key = lambda x: int((x.split('\n')[0]).split('_')[1])):
                                        out_file.write(l)
        print('\nDone\n')


def CDS_nuc_fasta_file_generator(path_out, annotation, main_chromosome = False):
        """ Generates a fasta file (in path out) whith all the CDS sequence for each organism contains in the annotation dictionary or from simple dictionary structure : annotation[organism][CDS_id]['sequence']
        """
        print("# .fa creation")
        if not os.path.exists(path_out):
               os.makedirs(path_out)
        for organism in annotation:
                print(organism)
                with open(path_out + organism + '.fa', 'w') as out_file:
                        if 'chromosome' not in annotation[organism] and 'plasmid' not in annotation[organism]:
                                for CDS_id in sorted(annotation[organism], key = lambda x: int(x.split('_')[1])):
                                                out_file.write('>' + CDS_id + '\n' + annotation[organism][CDS_id]['sequence'] + '\n')
                        else:
                                out_list = []
                                for chromosome_id in annotation[organism]['chromosome']:
                                        if 'CDS' in annotation[organism]['chromosome'][chromosome_id].keys():
                                                for CDS_id in sorted(annotation[organism]['chromosome'][chromosome_id]['CDS'], key = lambda x: int(x.split('_')[1])):
                                                        out_list.append('>' + CDS_id + '\n' + annotation[organism]['chromosome'][chromosome_id]['CDS'][CDS_id]['sequence'] + '\n')
                                if not main_chromosome:
                                        if 'plasmid' in annotation[organism]:
                                                for contig_id in annotation[organism]['plasmid']:
                                                        if 'CDS' in annotation[organism]['plasmid'][contig_id]:
                                                                for CDS_id in sorted(annotation[organism]['plasmid'][contig_id]['CDS'], key = lambda x: int(x.split('_')[1])):
                                                                        out_list.append('>' + CDS_id + '\n' + annotation[organism]['plasmid'][contig_id]['CDS'][CDS_id]['sequence'] + '\n')
                                for l in sorted(out_list, key = lambda x: int((x.split('\n')[0]).split('_')[1])):
                                        out_file.write(l)
        print('\nDone\n')


def homolog_finder(path_blast_res, identity_threshold = 40, alignement_threshold = 70, evalue_threshold = 1e-10, sequence_length_variation_threshold = 999):
        """ Findq all the homolog genes for each pair of organisms of the collection. BLASTP format 7 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' is required.
        A pair of homologous genes is defined as best hit between two genes which fit thresholds.
        Default thresholds are:
                > % identity >= 40%
                > % alignement >= 70% (with the smallest sequence)
                > E-value <= 1e-10
                > % sequence length variation = 999% (no treshold)
        """
        homolog = {}
        nb_file = sum(1 for f in os.listdir(path_blast_res) if f.endswith('.bl'))
        i = 0
        print('# Homolog dertermination...')
        for blast_res_file in os.listdir(path_blast_res):
                i += 1
                print('\r\t' + str(i) + '/' + str(nb_file) + ' .bl files processed' , end = '')
                if blast_res_file.endswith('.bl'):
                                ref, tar = blast_res_file.split('-vs-')
                                tar = tar[:-3]
                                if ref not in homolog:
                                        homolog[ref] = {}
                                if tar not in homolog[ref]:
                                        homolog[ref][tar] = {}
                                with open(path_blast_res + blast_res_file, 'r') as handler:
                                        #handler = open(path_blast_res + blast_res_file, 'r')
                                        line = handler.readline()
                                        no_hit = False
                                        while line:
                                                if line.startswith('# Query:'):
                                                        query = line[9:].split()[0]
                                                        line = handler.readline()
                                                        no_hit = False
                                                else:
                                                        line = handler.readline()
                                                        continue
                                                while line and line.startswith('#'):    # First hit expecting, if not hit, the next query is analysed
                                                        if line.startswith('# Query:'):
                                                                no_hit = True
                                                                break
                                                        line = handler.readline()
                                                if ref != tar:
                                                        if no_hit or not line:
                                                                homolog[ref][tar][query] = 'NA'
                                                                continue
                                                        line_handler = line.split()
                                                        query, subject, identity, alignement_length, evalue, query_length, subject_length = line_handler[0], line_handler[1], float(line_handler[2]), float(line_handler[3]), float(line_handler[10]), float(line_handler[12]), float(line_handler[13])
                                                        if abs(query_length - subject_length) <= min(query_length, subject_length) * (sequence_length_variation_threshold / 100):
                                                                alignement = (alignement_length / min(query_length, subject_length)) * 100
                                                                if evalue <= evalue_threshold and identity >= identity_threshold and alignement >= alignement_threshold:
                                                                        homolog[ref][tar][query] = subject
                                                                else:
                                                                        homolog[ref][tar][query] = 'NA'
                                                        else:
                                                                homolog[ref][tar][query] = 'NA'
                                                else:
                                                        homolog[ref][tar][query] = query
                                                        continue
        print('\nDone\n')
        return homolog


def core_finder(homolog, condition, reference = '', path = ''):
        """ Returns the core genome with the annotation of the reference from the homolog dictionary as a list of annotation. By default, the core is defined for all organisms.
        
        dictionary structure:
        
                [organism name] = core list
                
                > organism name: name of organism in collection like "Streptomyces_coelicolor_A3(2)"
                > core list = list of ID 
                
        Examples:
                >>> core['Streptomyces_coelicolor_A3(2)']
                ['Sco1_4884', 'Sco1_5243', 'Sco1_3183', 'Sco1_5706', 'Sco1_1644', 'Sco1_1302', 'Sco1_2167', 'Sco1_2592', 'Sco1_1618', 'Sco1_1602', 'Sco1_4231', 'Sco1_4098', 'Sco1_5628', 'Sco1_1377', 'Sco1_5338', ...]
        """
        
        print('\n# Data format conversion...')
        if not os.path.exists(path + '/' + 'ortholog_by_family_redundant' + condition + '.p'):
                ortholog_by_family_redundant = {}
                for ref in homolog:
                        for tar in homolog[ref]:
                                if not ref in ortholog_by_family_redundant:
                                        ortholog_by_family_redundant[ref] = {}
                                if ref != tar:
                                        for CDS_ref in homolog[ref][tar]:
                                                candidat = homolog[ref][tar][CDS_ref]
                                                if candidat in homolog[tar][ref] and homolog[tar][ref][candidat] == CDS_ref:
                                                        if CDS_ref not in ortholog_by_family_redundant[ref]:
                                                                ortholog_by_family_redundant[ref][CDS_ref] = []
                                                                ortholog_by_family_redundant[ref][CDS_ref].append(CDS_ref + ':' + ref)
                                                        ortholog_by_family_redundant[ref][CDS_ref].append(candidat + ':' + tar)
                pickle.dump(ortholog_by_family_redundant, open(path + '/' + 'ortholog_by_family_redundant' + condition + '.p', 'wb'))
        else:
                ortholog_by_family_redundant = pickle.load(open(path + '/' + 'ortholog_by_family_redundant' + condition + '.p', 'rb'))
        limit = len(ortholog_by_family_redundant)
        core = {}
        if not reference:
                list_ref = sorted(ortholog_by_family_redundant)
        else:
                list_ref = [reference]
                if list_ref[0] not in ortholog_by_family_redundant:
                        print('\tInvalid -reference argument: ' + ref + ' - Abort')
        print('\n#Core-genome construction...\n')
        for ref in list_ref:
                core[ref] = []
                i = len(homolog[ref][ref]) - len(ortholog_by_family_redundant[ref])
                nb_CDS = len(homolog[ref][ref])

                for CDS in ortholog_by_family_redundant[ref]:
                        i += 1
                        print('\r\t' + str(i) + '/' + str(nb_CDS) + ' CDS tested - reference: ' + ref, end = '')
                        if len(ortholog_by_family_redundant[ref][CDS]) == limit:
                                add_core = True
                                liste_w = [x.split(':')[0] for x in ortholog_by_family_redundant[ref][CDS]]
                                for w in ortholog_by_family_redundant[ref][CDS]:
                                        prot_w, ref_w = w.split(':')
                                        if ref_w in ortholog_by_family_redundant and prot_w in ortholog_by_family_redundant[ref_w] and len(ortholog_by_family_redundant[ref_w][prot_w]) == limit and (CDS + ':' + ref) in ortholog_by_family_redundant[ref_w][prot_w]:
                                                for e in ortholog_by_family_redundant:
                                                        if prot_w in homolog[ref][e]:
                                                                for j in ortholog_by_family_redundant[e][homolog[ref][e][prot_w]]:
                                                                        if j.split(':')[0] not in liste_w:
                                                                                add_core = False
                                                                                break
                                                        if not add_core:
                                                                break
                                        else:
                                                add_core = False
                                                break
                                if add_core:
                                        core[ref].append(CDS)
                print('\n')
        print('\nDone\n')
        return core


def ortholog_finder(homolog):
        """ Finds all the ortholog genes for each pair of organisms of the collection. A dictionary of homologous genes like dictionary generated by homolog_finder() is required.
        A pair of ortlogous genes is defined as best hit between two genes.
        
        dictionary structure:
        
                [ref name][tar name][ID ref] = ID tar
                
                > ref name: name of organism in collection like "Streptomyces_coelicolor_A3(2)"
                > ref tar: name of organism in collection like "Streptomyces_coelicolor_A3(2)"
                > ID ref: CDS ID of ref organism like: 'chrI_+_576-699'
                > ID tar: CDS ID of tar organism like: 'chrI_+_576-699' if the ortholog exists, else, ID tar = 'NA'
                
        Examples:
                >>> 
        """
        ortholog = {}
        nb_organism = len(homolog.keys())
        i = 0
        print('# Ortholog dertermination...')
        for ref in homolog:
                i += 1 
                print('\r\t' + str(i) + '/' + str(nb_organism), end = '')
                if ref not in ortholog:
                        ortholog[ref] = {}
                for tar in homolog[ref]:
                        if tar not in ortholog[ref]:
                                ortholog[ref][tar] = {}
                        for ref_gene in homolog[ref][tar]:
                                tar_gene = homolog[ref][tar][ref_gene]
                                if tar_gene != 'NA' and homolog[tar][ref][tar_gene] == ref_gene:
                                        ortholog[ref][tar][ref_gene] = tar_gene
                                else:
                                        ortholog[ref][tar][ref_gene] = 'NA'
        print('\nDone\n')
        return ortholog


def persistence_finder(ortholog):
        '''
        Computes gene persitence for each species
        '''
        
        persistence = {}
        nb_species = len(ortholog)
        for species in sorted(ortholog):
                persistence[species] = {}
                for cds in ortholog[species][species]:
                        persistence[species][cds] = 0
                        for target in sorted(ortholog[species]):
                                if ortholog[species][target][cds] != 'NA':
                                        persistence[species][cds] += 1
                        persistence[species][cds] = persistence[species][cds] / nb_species
        return persistence


def seq_prot_extractor(path_out, annotation, dico_to_extract):
        """
        Generates a fasta file (in path out) whith all the CDS sequence in the input list
        """
        print("# .fa creation")
        if not os.path.exists(path_out):
               os.makedirs(path_out)
        for organism in dico_to_extract:
                print(organism)
                with open(path_out + organism + '.fa', 'w') as out_file:
                        if organism in dico_to_extract:
                                out_list = []
                                for chromosome_id in annotation[organism]['chromosome']:
                                        if 'CDS' in annotation[organism]['chromosome'][chromosome_id].keys():
                                                for CDS_id in sorted(annotation[organism]['chromosome'][chromosome_id]['CDS'], key = lambda x: int(x.split('_')[1])):
                                                        if CDS_id in dico_to_extract[organism]:
                                                                # write annotation
                                                                #out_list.append('>' + CDS_id + '\n' + annotation[organism]['chromosome'][chromosome_id]['CDS'][CDS_id]['translation'] + '\n')
                                                                # write position
                                                                start, stop = ((annotation[organism]['chromosome'][chromosome_id]['CDS'][CDS_id]['position']).split('_')[-1]).split('-')
                                                                out_list.append(CDS_id + '\t' + start + '\t' + stop + '\n')
                                if 'plasmid' in  annotation[organism]:
                                        for contig_id in annotation[organism]['plasmid']:
                                                if 'CDS' in annotation[organism]['plasmid'][contig_id]:
                                                        for CDS_id in sorted(annotation[organism]['plasmid'][contig_id]['CDS'], key = lambda x: int(x.split('_')[1])):
                                                                if CDS_id in dico_to_extract[organism]:
                                                                        # write annotation
                                                                        #out_list.append('>' + CDS_id + '\n' + annotation[organism]['plasmid'][contig_id]['CDS'][CDS_id]['translation'] + '\n')
                                                                        # write position
                                                                        start, stop = ((annotation[organism]['plasmid'][contig_id]['CDS'][CDS_id]['position']).split('_')[-1]).split('-')
                                                                        out_list.append(CDS_id + '\t' + start + '\t' + stop + '\n')
                                for l in sorted(out_list, key = lambda x: int((x.split('\t')[0]).split('_')[1])):
                                        out_file.write(l)
        print('\nDone\n')


def conservation_indice_calculator(ortholog, dynamic_window = True, window_dimension = 5, mod = 'all', parallel_mod = False, sub_list = ''):
        """ Returns a dictionary with conservation indice profils for all paires of organism from the ortholog dinctionary. The 'dynamic_window' argument indicates if the value in 'window_length' is a number of CDS or proportion of all the CDS. 3 indices are available ('mod' argument): 'GOC', 'GOC2' and 'Tort'. GOC = Gene Order conservation (Rocha, 2006), GOC2 (Choulet, 2006), Tort = rate of orthologue: number of orthologous gene divided by the number of CDS. 'all' do all the indices.
        By default, a dynamic window with 5% for the 'all' mod is executed.
        """
        dico = {}
        if mod == 'all':
                print('\n# GOC / GOC2 / Tort calculation...\n')
                dico['GOC'] = {}
                dico['GOC2'] = {}
                dico['Tort'] = {}
        elif mod == 'GOC' or mod == 'GOC2' or mod == 'Tort':
                print('\n# ' + mod + ' calculation...\n')
                dico[mod] = {}
        else:
                print('\t\'' + str(mod) + '\' - incorrect \'mod\' argument (\'GOC\'/\'GOC2\'/\'Tort\'/\'all\')')
                return
        if type(window_dimension) != int:
                print('\t\'' + str(window_dimension) + '\' - incorrect \'window_dimension\' argument (positive integer)')
                return
        if window_dimension < 0:
                print('\t\'' + str(window_dimension) + '\' - incorrect \'window_dimension\' argument (positive integer)')
                return
        if type(dynamic_window) != bool:
                print('\t\'' + str(dynamic_window) + '\' - incorrect \'dynamic_window\' argument (boolean: True / False)')
                return
        i = 0
        if not parallel_mod:
                list_ref = ortholog.keys()
        else:
                list_ref = sub_list
        for ref in list_ref:
                if mod == 'all':
                        dico['GOC'][ref] = {}
                        dico['GOC2'][ref] = {}
                        dico['Tort'][ref] = {}
                else:
                        dico[mod][ref] = {}
                for tar in ortholog[ref]:
                        i += 1
                        if not parallel_mod:
                                print('\r\t' + str(i) + '/' + str(len(ortholog) * len(ortholog)), end = '')
                        else:
                                print('\r\t' + str(i) + '/' + str(len(ortholog) * len(list_ref)), end = '')
                        if mod == 'all':
                                dico['GOC'][ref][tar] = []
                                dico['GOC2'][ref][tar] = []
                                dico['Tort'][ref][tar] = []
                        else:
                                dico[mod][ref][tar] = []
                        limit = len(ortholog[ref][tar])
                        if dynamic_window:
                                window_length = math.ceil(limit / 100.0) * window_dimension
                        else:
                                window_length = window_dimension
                        start_window = 0
                        list_CDS_ref = sorted([CDS for CDS in ortholog[ref][tar]], key = lambda x: int(x.split('_')[1]))
                        list_ort_tar = [ortholog[ref][tar][CDS] for CDS in list_CDS_ref]
                        list_CDS_tar = sorted([CDS for CDS in ortholog[tar][tar]], key = lambda x: int(x.split('_')[1]))
                        index_tar_max = len(list_CDS_tar) - 1
                        while start_window < (limit - window_length):
                                window = list_ort_tar[start_window:(start_window + window_length)]
                                if mod == 'all' or mod == 'GOC' or mod == 'GOC2':
                                        number_adjacent = 0
                                        new = True
                                        for index, CDS in enumerate(window):
                                                if CDS != 'NA':
                                                        if index < (len(window) - 1):
                                                                index_tar = list_CDS_tar.index(window[index])
                                                                if index_tar !=  index_tar_max and list_CDS_tar[index_tar + 1] == window[index + 1]:
                                                                        number_adjacent += 1
                                                                        if new:
                                                                                number_adjacent += 1
                                                                                new = False
                                                                        continue
                                                                elif index > 0 and index_tar != 0:
                                                                        if list_CDS_tar[index_tar - 1] == window[index + 1]:
                                                                                number_adjacent += 1
                                                                                if new:
                                                                                        number_adjacent += 1
                                                                                        new = False
                                                                        else:
                                                                                new = True
                                                                                continue
                                                else:
                                                        new = True
                                        if mod == 'all' or mod == 'GOC':
                                                dico['GOC'][ref][tar].append(number_adjacent / len(window))
                                        if mod == 'all' or mod == 'GOC2':
                                                if len([x for x in window if x != 'NA']) != 0:
                                                        dico['GOC2'][ref][tar].append(number_adjacent / len([x for x in window if x != 'NA']))
                                                else:
                                                        dico['GOC2'][ref][tar].append(0)
                                if mod == 'all' or mod == 'Tort':
                                        dico['Tort'][ref][tar].append(len([x for x in window if x != 'NA']) / len(window))
                                start_window += 1
        print('\nDone\n')
        return dico


def profil_dict_to_csv(profil, window_dimension, path_intermediate_data):
        ''' 
        Converts the profil dictionary into a tabulated file. One repository by indice and one file by species.
        '''
        
        for indice in sorted(profil):
                path_out = path_intermediate_data + '/w' + str(window_dimension) + '/profil/' + str(indice) + '/'
                if not os.path.exists(path_out):
                        os.makedirs(path_out)
                print(path_out)
                for species in sorted(profil[indice]):
                        with open(path_out + species + '.tab', 'w') as outfile:
                                for target in sorted(profil[indice][species]):
                                        outfile.write(target + '\t' + '\t'.join([str(x) for x in profil[indice][species][target]]) + '\n')


def main():
    """
        Executes the core and persistence pipeline
    """
    #### CONSTANT

    working_rep = os.getcwd()

    #### PARAMETERS

    condition = '_coll_125'

    nb_core = 50 # Number of core to use for Blast
    window_dimension = 1 # Lengh of the window for the conservation profil construction
    main_chromosome = False # If "True" only the main chromosome is considered



    #### OPERATION TO DO: True / False

    DO_format_annotation =                  True
    DO_blastp =                             False
    DO_homolog =                            True
    DO_ortholog =                           True
    DO_core_genome =                        True
    DO_persistence =                        True
    DO_conservation_profil =                True


    WRITE_core =                            True
    WRITE_profil_conservation =             True
    
    #### PATHS

    path_base = working_rep

    path_gbk = path_base + '/gbk' + condition + '/'
    path_CDS_aa_fa = path_base + '/CDS_aa_fa' + condition + '/'
    path_CDS_nuc_fa = path_base + '/CDS_nuc_fa' + condition + '/'
    path_intermediate_data = path_base + '/intermediate_data/'
    path_blast_rep = path_base + '/Blast/run' + condition + '/'
    path_blast_data = path_blast_rep + '/Data/'
    path_blast_db_data = path_blast_rep + '/Data_db/'
    path_blast_db = path_blast_rep + '/Blast_db/'
    path_blast_output = path_blast_rep + '/Blast_res/'

    #### MAIN


    ### DO

    if DO_format_annotation:
            ''' Annotation extraction: gbk to dict '''
            annotation, dico_short_name_to_organism = gbk_to_pickle(path_gbk, path_base + '/annotation' + condition + '.p')
            pickle.dump(annotation, open(path_base + '/annotation' + condition.replace('_stephanie', '') + '.p', 'wb'))
            pickle.dump(dico_short_name_to_organism, open( path_base + '/short_name_to_organism' + condition + '.p', 'wb'))
            CDS_fasta_file_generator(path_CDS_aa_fa, annotation, main_chromosome = main_chromosome)
            CDS_nuc_fasta_file_generator(path_CDS_nuc_fa, annotation, main_chromosome = main_chromosome)

    if DO_blastp:
            ''' Blastp preparation '''
            if not os.path.exists(path_blast_data):
                    os.makedirs(path_blast_data)
            if not os.path.exists(path_blast_output):
                    os.makedirs(path_blast_output)
            if not os.path.exists(path_blast_db_data):
                    os.makedirs(path_blast_db_data)
            if not os.path.exists(path_blast_db):
                    os.makedirs(path_blast_db)
                    
            os.system('cp ' + path_CDS_aa_fa + '/*' + ' ' + path_blast_data)
            os.system('cp ' + path_CDS_aa_fa + '/*' + ' ' + path_blast_db_data)
            
            ''' Parallel preparation: split into mutltiple sub directory '''
            nb_file = len(os.listdir(path_blast_data))
            nb_by_dir = str(math.ceil(nb_file / nb_core))
            nb_dir = str(math.ceil(nb_file / int(nb_by_dir)))
            
            #Blast DB creation
            os.system('for i in `seq 1 ' + nb_dir + '`; do mkdir -p ' + path_blast_data + '/' + '$i; find ' + path_blast_data + ' -maxdepth 1 -type f | head -n ' + nb_by_dir + ' | xargs -i mv "{}" ' + path_blast_data + '/' + '$i; done')
            
            os.system('./make_blast_db.sh' + ' ' + path_blast_rep)


            ''' Blastp execution '''
            
            #Skip already done BLAST runs
            os.system('parallel ' + './launch_parallel_blastp_skip_blast.sh' + ' ' + path_blast_rep + ' $i ::: `seq ' + nb_dir + '`')

    if DO_homolog:
            ''' Homolog identification '''
            homolog = homolog_finder(path_blast_res = path_blast_output, identity_threshold = 40, alignement_threshold = 70, evalue_threshold = 1e-10)
            pickle.dump(homolog, open(path_base + '/' + 'homolog' + condition + '.p', 'wb'))

    if DO_ortholog:
            ''' Ortholog identification '''
            if os.path.exists(path_base + '/' + 'homolog' + condition + '.p'):
                    homolog = pickle.load(open(path_base + '/' + 'homolog' + condition + '.p', 'rb'))
                    ortholog = ortholog_finder(homolog)
                    pickle.dump(ortholog, open(path_base + '/' + 'ortholog' + condition + '.p', 'wb'))
            else:
                sys.exit("No homolog pickle file found.")

    if DO_core_genome:
            ''' Core genome construction '''
            if os.path.exists(path_base + '/' + 'homolog' + condition + '.p'):
                    homolog = pickle.load(open(path_base + '/' + 'homolog' + condition + '.p', 'rb'))
                    core = core_finder(homolog, condition, reference = '', path = path_base)
                    pickle.dump(core, open(path_base + '/' + 'core' + condition + '.p', 'wb'))
            else:
                sys.exit("No homolog pickle file found.")
                
    if DO_persistence:
            ortholog = pickle.load(open(path_base + '/' + 'ortholog' + condition + '.p', 'rb'))
            annotation = pickle.load(open(path_base + '/annotation' + condition + '.p', 'rb'))
            persistence = persistence_finder(ortholog)        
            if not os.path.exists(path_intermediate_data + '/persistence/'):
                    os.makedirs(path_intermediate_data + '/persistence/')
            for species in persistence:
                    with open(path_intermediate_data + '/persistence/' + species + '.tab', 'w') as outfile:                        
                            pos = position_extractor(persistence[species], annotation)
                            outfile.write('\t'.join([str(pos[x].split('_')[-1].split('-')[0]) for x in sorted(pos, key = lambda x: int(x.split('_')[1]))]) + 
                                        '\n' + 
                                        '\t'.join([str(persistence[species][x]) for x in sorted(persistence[species], key = lambda x: int(x.split('_')[1]))]) + '\n')
    
    if DO_conservation_profil:
        ''' GOC, GOC2 and Tort profil construction '''
        ortholog =  pickle.load(open(path_base + '/' + 'ortholog' + condition + '.p', 'rb'))
        conservation = conservation_indice_calculator(ortholog, dynamic_window = True, window_dimension = window_dimension, mod = 'all')
        pickle.dump(conservation, open(path_base + '/' + 'conservation_indices' + '_w' + str(window_dimension) + condition + '.p', 'wb'))

    ### WRITE

    if not os.path.exists(path_intermediate_data):
            os.makedirs(path_intermediate_data)

    if WRITE_core:
            ''' Creates a file with all the sequence of core '''
            annotation = pickle.load(open(path_base + '/annotation' + condition + '.p', 'rb'))
            core = pickle.load(open(path_base + '/' + 'core' + condition + '.p', 'rb'))
            seq_prot_extractor(path_intermediate_data + '/seq_core/', annotation, core)
    
    if WRITE_profil_conservation:
        ''' Converts the profil dictionary into a tabulated file. One repository by indice and one file by species '''
        profil = pickle.load(open(path_base + '/' + 'conservation_indices' + '_w' + str(window_dimension) + condition + '.p', 'rb'))
        profil_dict_to_csv(profil, window_dimension, path_intermediate_data)


if __name__ == "__main__":
    main()
