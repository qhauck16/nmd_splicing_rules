import sys
import gzip
import argparse
import pickle
from statistics import mean, median
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import pyfastx



chrom = 'chr1'
strand = '-'
junc = set((25832158, 25835022), (25837554, 25859280), (25834331, 25837413), (25835846, 25859280), (25835846, 25837413))
start_codons = set((25859355, 25859357))
stop_codons = set((25832088, 25832090), (25834993, 25834995))
gene_name = 'AUNIP'
exonLcutoff=1000
verbose=True


'''
Compute whether there is a possible combination that uses the junction without
inducing a PTC. We start with all annotated stop codon and go backwards.
'''

fa = pyfastx.Fasta('genome.fa')

seed = []

junc.sort()
if strand == "+":
    junc.reverse()
    
"""Quinn Comment: Adds all 'stop codons' to a nested list called seed""" 
##in an individual transcript    
for c in stop_codons:
    if strand == "+":
        seed.append([c[1]])
    else:
        seed.append([c[0]])

# seed starts with just stop codon and then a possible 3'ss-5'ss junction
# without introducing a PTC [stop_codon,3'ss, 5'ss, 3'ss, ..., start_codon]

junc_pass = {}
junc_fail = {}
path_pass = []
proteins = []

dic_terminus = {}

depth = 0

"""Quinn Comment: while our seed length is greater than 0 - which means we have charted all possible paths through 
all junctions ending in a stop codon (or there is an exon longer than 1000 bp and we have no complete paths)"""
while len(seed) > 0:
    new_seed = []
    final_check = []
    depth += 1
    if verbose:
        sys.stdout.write("Depth %s, Seed L = %s\n"%(depth, len(seed)))
    #print(start_codons, [s[-1] for s in seed][-10:], len(junc))
    framepos = {}
                
    for s in seed:
        # first check that the seed paths are good        
        bool_ptc = False
        leftover = ''
        if len(s) > 0:                
            leftover = Seq("")
            allprot = Seq("")

            """Quinn Comment: loop through the exons, calculating lengths"""
            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = exon_coord[1]-exon_coord[0]


                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                startpos = (len(leftover)+exlen+1)%3
                if strand == '+':
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover 
                    prot = seq[startpos:].translate()
                    leftover = seq[:startpos]                                                                                                               
                    allprot = prot+allprot  
                else:
                    seq = leftover+Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    aseq = seq
                    if startpos > 0:
                        leftover = seq[-startpos:]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    prot = seq[startpos:].translate()
                    allprot = prot+allprot

                #found a PTC in this transcript if any element but the last is a stop codon    
                bool_ptc = "*" in allprot[:-1]

        """Quinn Comment: if we found a PTC, add all intron coordinate pairs involved in the transcript to junc_fail"""        
        if bool_ptc:
            #This transcript failed
            for i in range(1, len(s)-1, 2):                                                                                                                  
                j_coord = s[i:i+2]                                                                                                                           
                j_coord.sort()                                                                                                                             
                j_coord = tuple(j_coord)                                                                                                                     
                if j_coord not in junc_fail:                                                                                                                 
                    junc_fail[j_coord] = 0                                                                                                                   
                junc_fail[j_coord] += 1  

            continue
    
        # passed
        """Quinn Comment: if we don't just have a stop codon, create a terminus for this 
        seed at the last 3' splice site or start codon; terminus is last two coordinates and the reading frame, 
        used for dynamic programming later"""
        if len(s) > 2:
            terminus = (s[-2],s[-1],leftover)
            
            if terminus in dic_terminus:
                dic_terminus[terminus].append(tuple(s))
                continue
            else:
                dic_terminus[terminus] = [tuple(s)]
        
        last_pos = s[-1]
        
        """Quinn Comment: check the last position of our seed to see if it is close to a start codon, within a potential exon's length,
            and add your seed plus this start codon to final_check """
        for start in start_codons:                
            #print("start", start, abs(last_pos-start[0]))
            if strand == "+" and last_pos > start[0] and abs(last_pos-start[0]) < exonLcutoff:
                final_check.append(s+[start[0]])
            elif strand == "-" and last_pos < start[1] and abs(last_pos-start[1]) < exonLcutoff:
                final_check.append(s+[start[1]]) 

        """Quinn Comment: add all possible places to go from our last_pos to the seed (nested list)"""
        for j0,j1 in junc:                
            if strand == "+" and last_pos > j1 and abs(last_pos-j1) < exonLcutoff:
                new_seed.append(s+[j1,j0])
            #print("junction", (j0,j1), abs(last_pos-j0))
            if strand == "-" and last_pos < j0 and abs(last_pos-j0) < exonLcutoff: 
                new_seed.append(s+[j0,j1])
                
    """Quinn Comment: Exited from s in seed loop, now we check our final_checks of the full paths, we do not
    eliminate paths based on presence of a PTC, rather we classify full complete paths without PTCs if they exist"""
    # check that the possible final paths are good
    for s in final_check:
        leftover = Seq("")
        allprot = Seq("")
        for i in range(0, len(s)-1, 2):
            exon_coord = s[i:i+2]
            exon_coord.sort()
            exon_coord = tuple(exon_coord)
            exlen = exon_coord[1]-exon_coord[0]
            startpos = (len(leftover)+exlen+1)%3
            if strand == "+":
                seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                leftover = seq[:startpos]  
                prot = seq[startpos:].translate()
                allprot = prot+allprot
            else:
                seq = leftover+Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                if startpos > 0:                                                                                                    
                    leftover = seq[-startpos:]                                    
                else:
                    leftover = Seq("")
                seq = seq.reverse_complement()                                                                                                           
                prot = seq[startpos:].translate()                                                                                                        
                allprot = prot+allprot                    
        bool_ptc = "*" in allprot[:-1]
    
        """Quinn Comment: Classify seed + start codon as a passing path if no PTCs found in previous block of code"""
        if not bool_ptc:
            # all pass
            proteins.append("\t".join([gene_name,chrom,strand, "-".join([str(x) for x in s]), str(allprot)])+'\n')
            #print("ALL PASS %s"%(s))
            path_pass.append(tuple(s))
            for i in range(1, len(s), 2):
                j_coord = s[i:i+2]
                j_coord.sort()
                j_coord = tuple(j_coord)
                if j_coord not in junc_pass:
                    junc_pass[j_coord] = 0
                junc_pass[j_coord] += 1

    seed = new_seed


"""Quinn Comment: OUT OF WHILE LOOP through all possible paths/seeds; 
check all termini to see if they are part of a full path that has been classified as passing"""
while True:
    new_paths = []
    for terminus in dic_terminus:
        terminus_pass = False
        for path_subset in dic_terminus[terminus]:
            for path in path_pass:
                if path[:len(path_subset)] == path_subset:
                    terminus_pass = True
                    break
        #print(terminus, terminus_pass)

        """Quinn Comment: if our terminus is part of a passing path, we want to make sure if is reflected in passing paths and
        add the associate junctions to junc_pass, only if they are not present"""
        if terminus_pass:
            for path_subset in dic_terminus[terminus]:
                if path_subset in path_pass: continue
                new_paths.append(path_subset)
                path_pass.append(path_subset)
                for i in range(1, len(path_subset), 2):
                    j_coord = list(path_subset[i:i+2])
                    j_coord.sort()
                    j_coord = tuple(j_coord)
                    if j_coord not in junc_pass:
                        junc_pass[j_coord] = 0
                        if verbose:
                            sys.stdout.write("junction %s pass\n"%j_coord)
    """Quinn Comment: we could have a new path_pass added, so our while loop checks again to see if there are any new paths 
    that are now going to be passing considering our additions"""
    if len(new_paths) == 0:
        break
        
junc_pass,junc_fail,proteins