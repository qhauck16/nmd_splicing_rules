# SpliceJunctionClassified V0.1 (Updated Jan 2024)
# Written by Yang Li Nov-2023

def tx_by_gene(gtf_annot):
    transcripts_by_gene = {}
    i = 1
    for dic in parse_gtf(gtf_annot):
        if dic['type'] == 'transcript':
            if dic['gene_name'] not in transcripts_by_gene:
                transcripts_by_gene[dic['gene_name']] = {dic['transcript_name']: []}
            else:
                transcripts_by_gene[dic['gene_name']] = transcripts_by_gene[dic['gene_name']] | {dic['transcript_name']: []}
        if dic['type'] == 'CDS':
            insort(transcripts_by_gene[dic['gene_name']][dic['transcript_name']], (dic['start']))
            insort(transcripts_by_gene[dic['gene_name']][dic['transcript_name']], (dic['end']))
    return transcripts_by_gene

def ptc_pos_from_prot(prot, sub):
    to_return = []
    start = 0
    while True:
        start = prot.find(sub, start)
        if start == -1: return to_return
        else:
            to_return.append(start)
        start += 1

def start_proximal(failing_juncs, gene_name, transcripts_by_gene, strand, chrom):
    distances = []
    start_rule = []
    for junc in failing_juncs:
        distance_from_start = []
        possible_transcripts = transcripts_by_gene[gene_name]


        for transcript in possible_transcripts:
            s = list(possible_transcripts[transcript])
            #cast to list of integers
            s = [eval(j) for j in s]
            if s == []:
                continue

            new_junc = False
            failing_junc = False
            for i in range(len(s)-2, 0, -2):
                junction = tuple([s[i-1], s[i]])

                #overlapping 3' or 5' splice site, just replace the junction
                if junc[0] == junction[0] or junc[1] == junction[1]:
                    s[i-1] = junc[0]
                    s[i] = junc[1]
                    new_junc = True
                #if this does not fit in our transcript, move on
                    if not (s[i-1] > s[i-2] and s[i] < s[i+1]):
                        failing_junc = True
                    
                    break  

            #If our overlapping site configuration fails in the transcript, skip this transcript     
            if failing_junc or not new_junc:
                continue

            if strand == '-':
                s.reverse()
            allprot = Seq("")
            leftover = Seq("")

            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = int(exon_coord[1])-int(exon_coord[0])
                if exlen > 1000:
                    break

                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                we must add 1 to length"""
                endpos = (len(leftover)+exlen+1)%3

                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                        leftover = Seq("")
                    else:
                        prot = seq[:-endpos].translate()
                        leftover = seq[-endpos:]                                                                                                               


                    bool_ptc = "*" in prot 

                    allprot = allprot+prot
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()

                    bool_ptc = "*" in prot

                    allprot = allprot+prot
                
                if bool_ptc:
                    ptc_coord = (min(ptc_pos_from_prot(allprot, '*')) + 1)*3
                    distance_from_start.append(ptc_coord)
                    break

        if len(distance_from_start) != 0:        
            distances.append(mode(distance_from_start))
            start_rule.append(junc)

    return start_rule, distances


def nucleotide_rule(failing_juncs, gene_name, transcripts_by_gene, strand, chrom):
    qualifying_prots = {'prot': [], 'gene': [], 'junction': [], 'transcript': []}
    distances = []
    nuc_rule = []
    unique_juncs = []
    unique_juncs_pre_ptc = []
    for junc in failing_juncs:
        junc_break = False
        distances_to_ejc = []
        possible_transcripts = transcripts_by_gene[gene_name]


        for transcript in possible_transcripts:
            s = list(possible_transcripts[transcript])
            #cast to list of integers
            s = [eval(j) for j in s]
            if s == []:
                continue

            new_junc = False
            failing_junc = False
            for i in range(len(s)-2, 0, -2):
                junction = tuple([s[i-1], s[i]])

                #overlapping 3' or 5' splice site, just replace the junction
                if junc[0] == junction[0] or junc[1] == junction[1]:
                    s[i-1] = junc[0]
                    s[i] = junc[1]
                    new_junc = True
                #if this does not fit in our transcript, move on
                    if not (s[i-1] > s[i-2] and s[i] < s[i+1]):
                        failing_junc = True
                    
                    break  

            #If our overlapping site configuration fails in the transcript, skip this transcript     
            if failing_junc or not new_junc:
                continue

            if strand == '-':
                s.reverse()
            allprot = Seq("")
            leftover = Seq("")

            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = int(exon_coord[1])-int(exon_coord[0])
                if exlen > 1000:
                    break

                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                we must add 1 to length"""
                endpos = (len(leftover)+exlen+1)%3
                if junc not in unique_juncs_pre_ptc:
                    unique_juncs_pre_ptc.append(junc)
                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                        leftover = Seq("")
                    else:
                        prot = seq[:-endpos].translate()
                        leftover = seq[-endpos:]                                                                                                               


                    bool_ptc = "*" in prot 
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot) 

                    allprot = allprot+prot
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()

                    bool_ptc = '*' in prot

                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot) 

                    allprot = allprot+prot
                
                if bool_ptc:
                    if junc not in unique_juncs:
                        unique_juncs.append(junc)
                    ptc_coord = [(x + 1)*3 for x in ptc_pos]
                    if i != len(s)-2 and i != len(s) - 4:
                        junc_break = True
                        break
                    else:
                        if i == len(s) - 2:
                            distances_to_ejc.append(-1)
                            qualifying_prots['prot'].append(allprot)
                            qualifying_prots['gene'].append(gene_name)
                            qualifying_prots['junction'].append(junc)
                            qualifying_prots['transcript'].append(possible_transcripts[transcript])
                            break
                        else:
                            distances_to_ejc.append(ptc_prot_len*3-min(ptc_coord))
                            break
                
            if junc_break:
                break
        if not junc_break and len(distances_to_ejc) != 0:
            nuc_rule.append(junc)
            distances.append(mode(distances_to_ejc))
    numbers = [len(failing_juncs), len(unique_juncs_pre_ptc), len(unique_juncs)]

    return nuc_rule, distances, qualifying_prots, numbers

def many_junctions(failing_juncs, gene_name, transcripts_by_gene, strand, chrom):

    before = {}
    after = {}
    for junc in failing_juncs:
        junc_added = False
        possible_transcripts = transcripts_by_gene[gene_name]
        possible_exons_before = []
        possible_exons_after = []
        for transcript in possible_transcripts:
            s = list(possible_transcripts[transcript])
            #cast to list of integers
            s = [eval(j) for j in s]
 

            new_junc = False
            failing_junc = False
            for i in range(len(s)-2, 0, -2):
                junction = tuple([s[i-1], s[i]])

                #overlapping 3' or 5' splice site, just replace the junction
                if junc[0] == junction[0] or junc[1] == junction[1]:
                    s[i-1] = junc[0]
                    s[i] = junc[1]
                    new_junc = True
                #if this does not fit in our transcript, move on
                    if not (s[i-1] > s[i-2] and s[i] < s[i+1]):
                        failing_junc = True
                    
                    break  
            #If our overlapping site configuration fails in the transcript, skip this transcript     
            if failing_junc or not new_junc:
                continue

            if strand == '-':
                s.reverse()
            allprot = Seq("")
            leftover = Seq("")

            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = int(exon_coord[1])-int(exon_coord[0])
                if exlen > 1000:
                    break

                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                we must add 1 to length"""
                endpos = (len(leftover)+exlen+1)%3

                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()
                    leftover = seq[-endpos:]                                                                                                                


                    bool_ptc = "*" in prot 
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot)
                        ptc_coord = (min(ptc_pos) + 1)*3


                    allprot = allprot+prot
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()

                    bool_ptc = "*" in prot
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot) 
                        ptc_coord = (min(ptc_pos) + 1)*3 
                        

                    allprot = allprot+prot

                exons_after = (len(s) - 2)/2 - i/2
                exons_before = i/2
                if bool_ptc and i != len(s)-2 and not (i == len(s)-4 and ptc_prot_len*3 - ptc_coord < 55): #and 55 nt rule not satisfied
                    possible_exons_before.append(exons_before)
                    possible_exons_after.append(exons_after)
                    break
        if len(possible_exons_before) != 0:
            before[junc] = mode(possible_exons_before)
            after[junc] = mode(possible_exons_after)
    return before, after


def long_exon_finder(failing_juncs, gene_name, transcripts_by_gene, strand, chrom):
    ptc_junctions = []
    ptc_exon_lens = []
    ptc_distances = []
    long_exons = []
    for junc in failing_juncs:
        junc_added = False
        possible_transcripts = transcripts_by_gene[gene_name]
        possible_distances = []
        possible_lens = []


        for transcript in possible_transcripts:
            s = list(possible_transcripts[transcript])
            #cast to list of integers
            s = [eval(j) for j in s]
            
            if s == []:
                continue

            new_junc = False
            failing_junc = False
            for i in range(len(s)-2, 0, -2):
                junction = tuple([s[i-1], s[i]])

                #overlapping 3' or 5' splice site, just replace the junction
                if junc[0] == junction[0] or junc[1] == junction[1]:
                    s[i-1] = junc[0]
                    s[i] = junc[1]
                    new_junc = True
                #if this does not fit in our transcript, move on
                    if not (s[i-1] > s[i-2] and s[i] < s[i+1]):
                        failing_junc = True
                    
                    break  

            #If our overlapping site configuration fails in the transcript, skip this transcript     
            if failing_junc or not new_junc:
                continue

            if strand == '-':
                s.reverse()
            allprot = Seq("")
            leftover = Seq("")

            bool_long_exon = False
            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = int(exon_coord[1])-int(exon_coord[0])
                if exlen > 1000:
                    break

                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                we must add 1 to length"""
                endpos = (len(leftover)+exlen+1)%3

                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()
                    leftover = seq[-endpos:]                                                                                                               


                    bool_ptc = "*" in prot 
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot)
                        ptc_coord = (min(ptc_pos) + 1)*3 

                    allprot = allprot+prot
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()

                    bool_ptc = "*" in prot
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot)
                        ptc_coord = (min(ptc_pos) + 1)*3 

                    allprot = allprot+prot
                
                #store a long_exon tag, only add this junction to list if it goes on to cause no PTCs that are not in long exons
            #     if bool_ptc and i != len(s)-2:
            #         ptc_coord = (min(ptc_pos) + 1)*3
            #         #55nt rule
            #         if i == len(s)-4 and ptc_prot_len*3-ptc_coord < 55:
            #             break
            #         elif exlen + 1 > 407:
            #             bool_long_exon = True
            #             break
            #         else:
            #             break
            # if bool_long_exon:
            #     long_exons.append(junc)
            #     ptc_distances.append(ptc_prot_len*3 - ptc_coord)
            #     ptc_exon_lens.append(ptc_prot_len*3)
            #     break
            # else:
                if bool_ptc and i != len(s)-2 and not (i == len(s)-4 and ptc_prot_len*3 - ptc_coord < 55):
            
                    possible_distances.append(ptc_prot_len*3 - ptc_coord)
                    possible_lens.append(ptc_prot_len*3)
                    break

        if len(possible_distances) != 0:
            ptc_distances.append(mode(possible_distances))
            ptc_exon_lens.append(mode(possible_lens))
            ptc_junctions.append(junc)
    return ptc_junctions, ptc_distances, ptc_exon_lens

def check_utrs(junc,utrs):
    '''
    checks if junction is close or within 100bp of UTRs
    '''
    for s1,s2 in list(utrs):
        if abs(junc[0]-s1) < 100 or abs(junc[1]-s2) < 100:
            return True
    return False

def solve_NMD(chrom, strand, junc, start_codons, stop_codons,gene_name, 
              verbose = False, exonLcutoff = 1000):
    '''
    Compute whether there is a possible combination that uses the junction without
    inducing a PTC. We start with all annotated stop codon and go backwards.
    '''

    global fa
    
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
    
    seq_db = {}
    junc_pass = {}
    junc_fail = {}
    path_pass = []
    proteins = []
    
    dic_terminus = {}
    dic_paths = {}

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
                    """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                    we must add 1 to length"""
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
    check all terminus' to see if they are part of a full path that has been classified as passing"""
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
            
    return junc_pass,junc_fail,proteins

def parse_gtf(gtf: str):
    '''Lower level function to parse GTF file
    - gtf: str : path to GTF annotation file
    - returns: dictionary with keys: 
        chrom, source, type, start, end, strand, frame, info, gene_name, transcript_type, transcript_name, gene_type
    '''
    fields = ["chrom", "source", "type", "start", "end", ".","strand","frame", "info"]
    open_gtf = lambda x: gzip.open(x) if ".gz" in x else open(x)
    for ln in open_gtf(gtf):
        ln = ln.decode('ascii') if ".gz" in gtf else ln
        dic = {}
        if ln[0] == "#": continue
        ln = ln.strip().split('\t')
        for i in range(len(fields)):
            dic[fields[i]] = ln[i]

        # add 4 additional fields, parsed from info field
        for ks in ['gene_name', "transcript_type","transcript_name", "gene_type"]:
            info_fields = [{x.split()[0]: x.split()[1].replace('"', '')} 
                          for x in dic['info'].split(';') if len(x.split()) > 1]
            info_fields = {k: v for d in info_fields for k, v in d.items()}
            try: 
                dic[ks] = info_fields[ks]
            except:
                dic[ks] = None # if line is a gene, then wont have transcript info
        yield dic
         

def parse_annotation(gtf_annot: str):
    '''
    Used `parse_gtf` to first parse a gtf file, then extract and return further
    information, including gene coordinates, intron info, and splice site to 
    gene-name dictionary.
    
    - gtf_anno: str : path to GTF annotation file
    - returns: 
        - genes_coords: gene coordinates, grouped by chromosome and strand 
        - introns_info: a dictionary with these keys: junctions (ie all introns),
                        start_codon, stop_codon, utrs, pcjunctions (ie only protein coding introns) 
        - ss2gene: a dictionary with splice site as keys, and gene name as values, eg. ('chr1', 11869): 'DDX11L1'
    '''
    genes_info = {}
    introns_info = {}
    genes_coords = {}
    ss2gene = {}

    for dic in parse_gtf(gtf_annot):
        chrom = dic['chrom']
        gname = dic['gene_name']
        tname = dic['transcript_name'], gname
        anntype  = dic['type']
        if dic['type'] == 'gene': 
            dic['transcript_type'] = "gene"
        if dic['transcript_type'] == "nonsense_mediated_decay" and anntype == "stop_codon": 
            continue 
        if dic['transcript_type'] != "protein_coding" and anntype == "UTR": 
            continue
        """Quinn: CONVERT TO BED FORMAT IS ERROR"""
        start, end = int(dic['start']), int(dic['end']) # convert to BED format
        strand = dic['strand']
    
        if (chrom, strand) not in genes_coords:
            genes_coords[(chrom, strand)] = []

        if tname not in genes_info: # tname is (transcript_name, gene_name)
            genes_info[tname] = {'exons':[],
                                 'start_codon':set(), 
                                 'stop_codon':set(), 
                                 'utrs':set(),
                                 'type':dic['transcript_type']
                                 }
        
        if anntype in ["start_codon", "stop_codon"]:
            genes_info[tname][anntype].add((start,end))
        elif anntype in ["gene"]:
            genes_coords[(chrom,strand)].append(((start,end), gname))
        elif anntype in ['exon']:
            genes_info[tname]['exons'].append((start,end))
            # Store gene info for splice sites
            ss2gene[(chrom, int(dic['start']) - 1)] = dic['gene_name'] # BED
            ss2gene[(chrom, int(dic['end']))] = dic['gene_name'] # BED

        elif anntype in ['UTR']:
            genes_info[tname]['utrs'].add((start,end))

    for gene in genes_info: # this is actually transcript level!!! (transcript_name, gene_name)
        gene_name = gene[1]
        exons = genes_info[gene]['exons']
        exons.sort()

        pc = False # Basically use transcript_type == protein_coding to set flag
        if genes_info[gene]['type'] == "protein_coding": # AGAIN here gene = (transcript_name, gene_name)
            pc = True

        junctions = set() # all junctions/introns
        pcjunctions = set() # only protein coding junctions/introns

        for i in range(len(exons)-1):
            intron = (exons[i][1], exons[i+1][0])
            junctions.add(intron)
            if pc:
                pcjunctions.add(intron)

        if gene_name not in introns_info: # here only gene_name, does not incl. transcript_name
            introns_info[gene_name] = {'junctions':set(),
                                       'start_codon':set(),
                                       'stop_codon':set(),
                                       'utrs':set(),
                                       'pcjunctions':set(),
                                 }

        introns_info[gene_name]['start_codon'] = introns_info[gene_name]['start_codon'].union(genes_info[gene]['start_codon'])
        introns_info[gene_name]['stop_codon'] = introns_info[gene_name]['stop_codon'].union(genes_info[gene]['stop_codon'])
        introns_info[gene_name]['junctions'] = introns_info[gene_name]['junctions'].union(junctions)
        introns_info[gene_name]['utrs'] = introns_info[gene_name]['utrs'].union(genes_info[gene]['utrs'])
        introns_info[gene_name]['pcjunctions'] = introns_info[gene_name]['pcjunctions'].union(pcjunctions)

    return genes_coords, introns_info, ss2gene



# def getmean(lst):
#     return sum(lst)/float(len(lst))

# def getmedian(lst):
#     lst.sort()
#     n = len(lst)
#     if n < 1:
#         return None
#     if n % 2 == 1:
#         return lst[n//2]
#     else:
#         return sum(lst[n//2-1:n//2+1])/2.0

def get_feature(fname, feature = "exon"):
    ss2gene = {}
    if ".gz" in fname:
        F = gzip.open(fname)
    else:
        F = open(fname)
    for ln in F:
        if ln[0] == "#": continue
        ln = ln.split('\t')
        gID = ln[-1].split('gene_name "')[1].split('"')[0]

        if ln[2] != feature:
            continue
        ss2gene[(ln[0], int(ln[3]))] = gID
        ss2gene[(ln[0], int(ln[4]))] = gID
    return ss2gene


def get_overlap_stream(L1, L2, relax = 0):
    '''                                                                                                                                                                
    L1 and L2 are sorted lists of tuples with key, values                                                                                                              
    '''
    i, j = 0, 0
    while i < len(L1) and j < len(L2):
        if L1[i][0][1] < L2[j][0][0]:
            i += 1
            continue
        elif L2[j][0][1] < L1[i][0][0]:
            j += 1
            continue
        else:
            k = 0
            # hits overlapping, check all L2 that may overlap with intron                                                                                              
            while L2[j+k][0][0] <= L1[i][0][1]:
                if overlaps(L1[i][0], L2[j+k][0]):
                    yield L1[i], L2[j+k]
                k += 1
                if j+k == len(L2): break
            i += 1

def overlaps(A: tuple, B: tuple):
    '''                                                                                                                                                                
    Checks if A and B overlaps                                                                                                                                         
    A: tuple : (start, end)
    B: tuple : (start, end)
    '''
    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True


def ClassifySpliceJunction(options):
    '''
        - perind_file: str : path to counts file, e.g. leafcutter_perind.counts.gz
        - gtf_annot: str : Annotation GTF file, for example gencode.v37.annotation.gtf.gz
        - rundir: str : run directory, default is current directory
    '''

    gtf_annot, rundir, outprefix = options.annot, options.rundir, options.outprefix
    verbose = False or options.verbose
    if options.countfile is None:
        perind_file = f"{rundir}/{outprefix}_perind.counts.gz"
    else:
        perind_file = options.countfile

    # read leafcutter perind file and store junctions in dictionary: dic_junc
    # key = (chrom,strand), value = list of junctions [(start,end)]
    dic_junc = {}
    sys.stdout.write(f"Processing junction counts {perind_file}...")
    for ln in gzip.open(perind_file):
        junc_info = ln.decode('ascii').split()[0] # first column
        if junc_info == "chrom": continue # header
        
        chrom, start, end, clu_strand = junc_info.split(":")
        strand = clu_strand.split("_")[-1]
        if (chrom,strand) not in dic_junc: 
            dic_junc[(chrom,strand)] = []
        dic_junc[(chrom,strand)].append((int(start), int(end)))

    sys.stdout.write(" done!\n")
    if verbose:
        sys.stdout.write("Processed: ")
        for chrstrand in dic_junc:
            sys.stdout.write(f"{len(dic_junc[chrstrand])} jxns on {chrstrand[0]} ({chrstrand[1]}).")

    
    # load or parse gtf annotations
    # g_coords: gene coordinates, grouped by chromosome and strand
    # g_info: a dictionary with (transcript_name, gene_name) as keys, and intron info as values
    try: 
        sys.stdout.write("Loading annotations...")
        parsed_gtf = f"{rundir}/{gtf_annot.split('/')[-1].split('.gtf')[0]}_SJC_annotations.pckle"
        with open(parsed_gtf, 'rb') as f:
            g_coords, g_info = pickle.load(f)
        sys.stdout.write(" done!\n")
    except:
        sys.stdout.write("Parsing annotations for the first time...\n")
        g_coords, g_info, ss2gene = parse_annotation(gtf_annot)
        
    
        for chrom,strand in g_coords:
            to_remove_gcoords = set()
            to_remove_ginfo = set()
            for gene in g_coords[(chrom,strand)]:
                if gene[1] in g_info: # check if gene_name is in introns_info keys
                    if len(g_info[gene[1]]['stop_codon']) == 0: # if there are no stop codons
                        to_remove_gcoords.add(gene) # mark gene for remove
                        to_remove_ginfo.add(gene[1]) # mark gene_name for removal

            for g in to_remove_gcoords:
                g_coords[(chrom,strand)].remove(g)
            for g in to_remove_ginfo:
                g_info.pop(g)
        sys.stdout.write("Saving parsed annotations...\n")
        with open(parsed_gtf, 'wb') as f:
            pickle.dump((g_coords, g_info), f)


    txn2gene = f"{rundir}/txn2gene.{gtf_annot.split('/')[-1].split('.gtf')[0]}_SJC_annotations.pckle"
    try:
        sys.stdout.write("Loading txn2gene annotations...")
        with open(txn2gene, 'rb') as f:
            transcripts_by_gene = pickle.load(f)
    except:
        sys.stdout.write("Failed... Making txn2gene annotations...\n")
        transcripts_by_gene = tx_by_gene(gtf_annot)
        with open(txn2gene, 'wb') as f:
            pickle.dump(transcripts_by_gene, f)
    #transcripts_by_gene = tx_by_gene(gtf_annot)
    sys.stdout.write(" done!\n")

    gene_juncs = {}
    for chrom,strand in dic_junc:
        if (chrom,strand) not in g_coords: 
            sys.stderr.write(f"Could not find {chrom} ({strand}) in annotations...\n")
            continue
        juncs = [(x,x) for x in dic_junc[(chrom,strand)]]
        juncs.sort()

        coords = g_coords[(chrom,strand)]
        coords.sort()
        
        # save junctions that overlapping a gene in gene_juncs dictionary: (gene_name, chrom, strand) : [junctions]
        for junc, geneinfo in get_overlap_stream(juncs,coords): 
            info = (geneinfo[1], chrom,strand)
            if info not in gene_juncs:
                gene_juncs[info] = []
            gene_juncs[info].append(junc[0])

    fout = open(f"{rundir}/{outprefix}_junction_classifications.txt",'w')
    fout.write("\t".join(["Gene_name","Intron_coord","Annot","Coding"])+'\n')
    lout = open(f"{rundir}/{outprefix}_long_exon_distances.txt",'w')
    lout.write("\t".join(["Gene_name","Intron_coord","PTC_position","Exon_length"])+'\n')
    eout = open(f"{rundir}/{outprefix}_exon_stats.txt",'w')
    eout.write("\t".join(["Gene_name","Intron_coord","Exons_before","Exons_after"])+'\n')
    nout = open(f"{rundir}/{outprefix}_nuc_rule_distances.txt",'w')
    nout.write("\t".join(["Gene_name","Intron_coord","ejc_distance"])+'\n')
    dout = open(f"{rundir}/{outprefix}_distances_from_start.txt",'w')
    dout.write("\t".join(["Gene_name","Intron_coord","distance_from_start"])+'\n')
    


    nuc_out = open(f"{rundir}/{outprefix}_nuc_rule_proteins.txt",'w')
    junc_numbers = [0,0,0]
    for gene_name, chrom, strand in gene_juncs:

        sys.stdout.write(f"Processing {gene_name} ({chrom}:{strand})\n")
        
        query_juncs = gene_juncs[(gene_name,chrom,strand)] # from LeafCutter perind file
        if gene_name not in g_info: continue
        junctions = g_info[gene_name]['junctions'] # from annotation
        
        # classify all junctions in gene
        junctions = list(junctions.union(query_juncs))

        start_codons = g_info[gene_name]['start_codon'] 
        stop_codons = g_info[gene_name]['stop_codon']

        if verbose:
            sys.stdout.write(f"LeafCutter junctions ({len(query_juncs)}) All junctions ({len(junctions)}) Start codons ({len(start_codons)}) Stop codons ({len(stop_codons)}) \n")

        junc_pass, junc_fail, proteins = solve_NMD(chrom,strand,junctions, 
                                                   start_codons, stop_codons, 
                                                   gene_name)

        junc_fail = set(junc_fail.keys())
        junc_pass = set(junc_pass.keys())
        failing_juncs = junc_fail.difference(junc_pass)

        old_junc_pass = junc_pass
        junc_pass = {}
        junc_pass['normal'] = old_junc_pass
        ptc_junctions, ptc_distances, ptc_exon_lens = long_exon_finder(failing_juncs, gene_name, transcripts_by_gene, strand, chrom)
        exons_before, exons_after = many_junctions(failing_juncs, gene_name, transcripts_by_gene, strand, chrom)
        junc_pass['nuc_rule'], ejc_distances, nuc_prots, numbers = nucleotide_rule(failing_juncs, gene_name, transcripts_by_gene, strand, chrom)
        junc_numbers = [junc_numbers[i]+numbers[i] for i in range(3)]
        start_rule, start_distances = start_proximal(failing_juncs, gene_name, transcripts_by_gene, strand, chrom)
        for j in junctions:

            bool_pass = j in junc_pass['normal'] or j in g_info[gene_name]['pcjunctions']
            bool_fail = j in failing_juncs
            utr = False
            #long_exon = j in junc_pass['long_exon'] and j not in g_info[gene_name]['pcjunctions']
            #nuc_rule = j in junc_pass['nuc_rule'] and j not in g_info[gene_name]['pcjunctions']

            if bool_fail or bool_pass:
                tested = True
            else:
                tested = False
            annotated = j in g_info[gene_name]['junctions']
            #if not bool_pass and annotated:
            #print("%s %s %s junction: %s tested: %s utr: %s coding: %s annotated: %s "%(chrom, strand, gene_name, j, tested,utr, bool_pass, annotated))
            
            fout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]}',
                                  str(annotated), str(bool_pass)])+'\n')
            
        for w in range(len(ptc_junctions)):
            j = ptc_junctions[w]
            if j not in g_info[gene_name]['pcjunctions']:
                lout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]}',
                                    str(ptc_distances[w]), str(ptc_exon_lens[w])])+'\n')
        for j in exons_before:
            if j not in g_info[gene_name]['pcjunctions']:
                eout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]}',
                                    str(exons_before[j]), str(exons_after[j])])+'\n')
        for w in range(len(junc_pass['nuc_rule'])):
            j = junc_pass['nuc_rule'][w]
            if j not in g_info[gene_name]['pcjunctions']:
                nout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]}',
                                    str(ejc_distances[w])])+'\n') 
        for i in range(len(nuc_prots['prot'])):

            nuc_out.write('\t'.join([nuc_prots['gene'][i], str(nuc_prots['transcript'][i]),
                                    str(nuc_prots['junction'][i]), str(nuc_prots['prot'][i])])+'\n') 
        for i in range(len(start_rule)):
            if j not in g_info[gene_name]['pcjunctions']:
                j = start_rule[i]
                dout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]}',
                                    str(start_distances[i])])+'\n')
    print(junc_numbers)        
def boolean_to_bit(bool_vec):
    # Convert boolean vector to string of "1"s and "0"s
    bin_str = ''.join(['1' if b else '0' for b in bool_vec])
    
    # Convert this binary string into an integer
    # bit_num = int(bin_str, 2)
    
    return bin_str
            
def merge_discordant_logics(sjc_file: str):
    '''some junctions have multiple classifications. Use conservative approach
    to merge them.
    '''
    sjc = pd.read_csv(sjc_file, sep = "\t")

    classifier = {
        # each bit represents [ is annotated, is coding, is UTR ]
        '000': 'UP', # UnProductive,
        '001': 'NE', # NEither, not productive, but not considered unprod. due to close to UTR
        '010': 'PR', # PRoductive
        '100': 'UP', # UnProductive
        '101': 'PR', # PRoductive
        '110': 'PR' # PRoductive
        }
    
    # group dt
    sjc = sjc.groupby('Intron_coord').agg('max').reset_index()
    # convert Annotation, Coding, UTR status to binary strings then to SJ categories
    sjc['SJClass'] = sjc.apply(lambda x: boolean_to_bit(x[2:5]), axis=1).map(classifier)
    
    # convert df to dict
    sjc = sjc.set_index('Intron_coord').to_dict(orient='index')
    sjc = {tuple([k.split(':')[0]]) + tuple(k.split(':')[1].split('-')): v for k, v in sjc.items()}
    sjc = {(k[0], int(k[1]), int(k[2])): v for k, v in sjc.items()}


    return sjc
    # sjc is a dcitionary with:
    # - keys: intron coordinates, e.g. ('chr1', 1000, 2000)
    # - values: a dictionary e.g. {'Gene_name': 'DNMBP', 'Annot': False, 'Coding': False, 'UTR': False, 'SJClass': 'UP'})

    
    
    
    


def main(options):

    if options.countfile is None:
        sys.stderr.write("Error: no LeafCutter junction file provided...\npython SpliceJunctionClassifier.py -c Leafcutter_perind.counts.gz\n")
        exit(0)

    if options.genome is None:
        sys.stderr.write("Error: no genome fasta file selected...\npython SpliceJunctionClassifier.py -G genome.fa\n")
        exit(0)

    if options.annot is None:
        sys.stderr.write("Error: no annotation file with gene start and stop codon...\npython SpliceJunctionClassifier.py -A gencode.gtf/ensembl.gtf\n")
        exit(0)
    
    global fa
    sys.stdout.write(f"Loading genome {options.genome} ...")
    fa = pyfastx.Fasta(options.genome)
    sys.stdout.write("done!\n")

    ClassifySpliceJunction(options)

    # for testing
    # sjc_file = f"{options.rundir}/{options.outprefix}_junction_classifications.txt"
    # print(f"Merging discordant logics in {sjc_file}...")
    # # print the first 10 items in the merged dictionary
    # sjc = merge_discordant_logics(sjc_file)
    # print(list(sjc.items())[:10])


if __name__ == "__main__":

    import sys
    import gzip
    import argparse
    import pickle
    from statistics import mean, median
    import numpy as np
    import pandas as pd
    from Bio.Seq import Seq
    import pyfastx
    from bisect import insort
    from statistics import mode
    import time

    parser = argparse.ArgumentParser(description='SpliceJunctionClassifier')

    parser.add_argument("-c", "--countfile", dest="countfile",
                  help="LeafCutter perind counts file, e.g. leafcutter_perind.counts.gz")

    parser.add_argument("-o", "--outprefix", dest="outprefix", default = 'Leaf2',
                      help="output prefix (default: Leaf2)")

    parser.add_argument("-r", "--rundir", dest="rundir", default = '.',
                      help="run directory (default: .)")

    parser.add_argument("-A", "--annotation", dest="annot",
                  help="Annotation GTF file, for example gencode.v37.annotation.gtf.gz")
    
    parser.add_argument("-G", "--genome", dest="genome",
                  help="Reference genome fasta file.")

    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default = False,
                      help="verbose mode")
                
    options = parser.parse_args()
    
    main(options)
