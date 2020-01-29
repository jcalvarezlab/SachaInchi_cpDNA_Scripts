#!/usr/bin/env python 
import argparse
import numpy as np
from itertools import product
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable

def get_arguments():
    parser = argparse.ArgumentParser(description='Calculate the codon usage table from a sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('-i', '--input', type=str,
                        help='Path to genbank input')
    inputs.add_argument('-n', '--ncbi', type=str, default="fasta",
                        help='Give an NCBI ID instead of a local file')
    inputs.add_argument('--print', action ="store_true",
                       help='Print genetic code options')
    parser.add_argument('-e','--email', type=str,
                       help='e-mail to use for ncbi searches')
    parser.add_argument('-c', '--code', type=str, default="Standard",
                       help='Genetic code to use, either number or name')
    
    parser.add_argument('-o', '--output', type=str, default="./table.tsv",
                        help='Path to put codon usage table (%(default)s)')
    args = parser.parse_args()
    return args

def extract_fts(record, ft_type):
    return [ft for ft in record.features if ft.type == ft_type]

def init_codon_dict(alphabet):
    return {codon:i for i, codon in enumerate(map("".join, product(alphabet, repeat=3)))}

def codon_usage_table(sequence, features, codon_idx):
    #codon_idx = init_codon_dict()
    mat = np.zeros((64,len(features)), dtype=np.int32) ## there are 64 possible codons
    
    for j,ft in enumerate(features):
        seq_ft = ft.extract(sequence)
        for i in range(0,len(seq_ft),3):
            codon = str(seq_ft[i:i+3])
            if len(codon) != 3:
                name = ft.qualifiers["gene"][0]
                print(f"{name} has {len(seq_ft) % 3} bp not part of a codon, which is abnormal in a coding sequence. Skipping... \n ")
                #continue
            else:
                mat[codon_idx[codon],j] += 1 
    return mat
    
def identify_trnas(features):
    codons = init_codon_dict("ACGU").keys()
    cod2trna = dict(zip(codons,[None]*64)) ## Initialize in None (we have not found any tRNAs yet)
    
    for ft in features:
        name = ft.qualifiers["gene"][0]
        anticodon = name.split("-")[1]
        codon = str(Seq(anticodon, IUPAC.unambiguous_rna).reverse_complement())
        
        if cod2trna[codon] == None:
            cod2trna[codon] = [ft]
        else:
            cod2trna[codon].append(ft)
    
    return cod2trna

def check_gc_by_name(string): ## check if genetic code was given by name or number
        try:
            int(string)
        except ValueError:
            return True
        else:
            return False

def cod2amin(genetic_code, by_name = False, rna = False):
    if rna:
        if by_name:
            cod2aa = CodonTable.unambiguous_rna_by_name[genetic_code].forward_table
            cod2aa.update({cod:"*" for cod in CodonTable.unambiguous_rna_by_name[genetic_code].stop_codons})
        else:    
            cod2aa = CodonTable.unambiguous_rna_by_id[int(genetic_code)].forward_table
            cod2aa.update({cod:"*" for cod in CodonTable.unambiguous_rna_by_id[int(genetic_code)].stop_codons})
    else:
        if by_name:
            cod2aa = CodonTable.unambiguous_dna_by_name[genetic_code].forward_table
            cod2aa.update({cod:"*" for cod in CodonTable.unambiguous_dna_by_name[genetic_code].stop_codons})
        else:    
            cod2aa = CodonTable.unambiguous_dna_by_id[int(genetic_code)].forward_table
            cod2aa.update({cod:"*" for cod in CodonTable.unambiguous_dna_by_id[int(genetic_code)].stop_codons})
    return cod2aa
    
def calc_rscu(matrix, cod_idx, genetic_code):
    cod2aa = cod2amin(genetic_code, by_name = check_gc_by_name(genetic_code))
    aa2cod = {}
    
    for cod, aa in cod2aa.items():
        if aa not in aa2cod:
            aa2cod[aa] = [cod]
        else:
            aa2cod[aa].append(cod)
    
    rscu = {}
    percent = {}
    for cod, aa in cod2aa.items():
        num_aa = len(aa2cod[aa])
        conteo_cod = matrix[cod_idx[cod],:].sum()
        conteo_aa = matrix[list(map(cod_idx.get, aa2cod[aa])),:].sum()
        
        ## num_aa * conteo codon / conteo aminoacidos
        rscu[cod] = num_aa * conteo_cod / conteo_aa
        percent[cod] = rscu[cod] / num_aa
    
    ## In the end Dict: codon --> rscu
    return rscu, percent
    
def save_cod_table(path, matrix, codon_idx, rscu, percent, trnas, genetic_code):
    cod2aa = cod2amin(genetic_code, by_name = check_gc_by_name(genetic_code), rna = True)
    
    with open(path, "w") as file:
        file.write("\t".join(["Codon","Aminoacid","Count","RSCU","%","tRNA"])+ "\n")
        for cod, aa in sorted(cod2aa.items(), key=lambda item: item[1]):
            pre_line = []
            pre_line.append(cod)
            pre_line.append(aa)
            pre_line.append(str(matrix[codon_idx[cod],:].sum()))
            pre_line.append(str(round(rscu[cod.replace("U","T")],3)))
            pre_line.append(str(round(percent[cod.replace("U","T")]*100,2)))
            pre_line.append(" " if trnas[cod] == None else trnas[cod][0].qualifiers["gene"][0])
            line = "\t".join(pre_line) + "\n"
            file.write(line)
    
def main():
    args = get_arguments()
    
    if args.print:
        print("Number\tGenetic Code")
        for name, code in CodonTable.unambiguous_dna_by_name.items():
            print(f"{code.id}\t{name}")
    else:
        ## import necesary libraries and functions
        if args.input:
            rec = list(SeqIO.parse(args.input,"gb")) ## Read sequence from file
        else:
            from Bio import Entrez
            Entrez.email = args.email
            handle = Entrez.efetch(db="nucleotide", id=args.ncbi, rettype="gb", retmode="text")
            rec = list(SeqIO.parse(handle, "gb"))

        aux = {}
        aux["seqs"] = [seq.seq for seq in rec]
        aux["cds"] = list(map(extract_fts, rec, ["CDS"]*len(rec))) ## Extract CDS from sequence
        #cds = extract_fts(rec, "CDS") ## Extract CDS from sequence
        aux["trnas"] = list(map(extract_fts, rec, ["tRNA"]*len(rec))) ## Extract tRNAs from sequence
        #trnas = extract_fts(rec, "tRNA") ## Extract tRNAs from sequence
        
        codon_idx = init_codon_dict("ACGT")
        
        aux["mat"] = [codon_usage_table(seq, cds, codon_idx) for seq, cds in zip(aux["seqs"], aux["cds"])]
        #matrix = codon_usage_table(rec, cds, codon_idx) ## Count codon usage and put in matrix
        matrix = np.concatenate(aux["mat"], axis=1)
        #normal_mat = mat.sum(axis=1)/ mat.sum() ## Sum count and normalize
        
        cod_trnas = identify_trnas((*aux["trnas"])) ## Identify tRNAs in sequence for those codons

        ## Calculate RSCU:
        ## Relative synonymous codon usage (RSCU) is defined as the ratio of the observed 
        ## frequency of codons to the expected frequency given that all the synonymous codons 
        ## for the same amino acids are used equally
        rscu, percent = calc_rscu(matrix, codon_idx, args.code)
        

        ## Create a table and save it to file
        save_cod_table(args.output, matrix, init_codon_dict("ACGU"), rscu, percent, cod_trnas, args.code)
        
if __name__ == "__main__":
    main()
