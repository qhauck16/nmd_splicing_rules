# SpliceJunctionClassified V0.1 (Updated Jan 2024)
# Written by Yang Li Nov-2023
# Updates by Quinn Hauck Jun-2024

def ptc_pos_from_prot(prot, sub):
    to_return = []
    start = 0
    while True:
        start = prot.find(sub, start)
        if start == -1: return to_return
        else:
            to_return.append(start)
        start += 1   
   

def check_utrs(junc,utrs):
    '''
    checks if junction is close or within 100bp of UTRs
    '''
    for s1,s2 in list(utrs):
        if abs(junc[0]-s1) < 100 or abs(junc[1]-s2) < 100:
            return True
    return False

def solve_NMD(chrom, strand, junc, start_codons, stop_codons,gene_name, 
              verbose = True, exonLcutoff = 1000):
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

    junc_pass = {'normal':{}, 'long_exon':{}}
    junc_fail = {}
    path_pass = {'normal':[], 'long_exon':[]}
    proteins = []
    short_ptcs = []
    final_check = []

    dic_terminus = {'normal': {}, 'long_exon': {}}

    depth = 0

    """Quinn Comment: while our seed length is greater than 0 - which means we have charted all possible paths through 
    all junctions ending in a stop codon (or there is an exon longer than 1000 bp and we have no complete paths)"""
    while len(seed) > 0:
        new_seed = []
        depth += 1
        if verbose:
            sys.stdout.write("Depth %s, Seed L = %s\n"%(depth, len(seed)))
        #print(start_codons, [s[-1] for s in seed][-10:], len(junc))
        framepos = {}
                    
        for s in seed:
            # first check that the seed paths are good        
            bool_ptc = False
            leftover = ''
            long_ptcs = []
            if len(s) > 0:                
                leftover = Seq("")
                allprot = Seq("")

                """Quinn Comment: loop through the exons, calculating lengths"""
                for i in range(0, len(s)-1, 2):
                    short = []
                    exon_coord = s[i:i+2]
                    exon_coord.sort()
                    exon_coord = tuple(exon_coord)
                    exlen = exon_coord[1]-exon_coord[0]


                    """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                    startpos = (len(leftover)+exlen+1)%3
                    if strand == '+':
                        seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover 
                        #Last exon rule, part of 50-55nt rule
                        if i == 0:
                            prot = seq[startpos:].translate(stop_symbol = '#')
                        #long exon rule
                        elif exlen + 1 > 407:
                            prot = seq[startpos:].translate(stop_symbol = '@')
                            ptc_pos = ptc_pos_from_prot(prot, '@')
                            long_ptcs.append([exon_coord[0] + startpos + 3*x for x in ptc_pos])
                        else:
                            prot = seq[startpos:].translate()

                            #store ptc position for checking with long_exon PTCs later
                            ptc_pos = ptc_pos_from_prot(prot, '*')
                            ptc_coord = [exon_coord[0] + startpos + 3*x for x in ptc_pos]
                            for k in ptc_coord:
                                short.append(k)
                        #50-55nt rule
                        if i == 2:
                            ptc_pos = ptc_pos_from_prot(prot, '@') + ptc_pos_from_prot(prot, '*')
                            close_to_ejc = [(len(prot) - 1)*3 - 3*x - len(leftover) < 56 for x in ptc_pos]
                            if sum(close_to_ejc) == len(close_to_ejc):
                                prot = seq[startpos:].translate(stop_symbol = '#')
                                #don't want to store long or short ptc position if they pass 50-55nt rule
                                long_ptcs = []
                                short = []

                        leftover = seq[:startpos]                                                                                                               
                        allprot = prot+allprot  
                    else:
                        seq = leftover+Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                        if startpos > 0:
                            leftover = seq[-startpos:]
                        else:
                            leftover = Seq("")
                        seq = seq.reverse_complement()
                        
                        if i == 0:
                            prot = seq[startpos:].translate(stop_symbol = '#')
                        elif exlen + 1 > 407:
                            prot = seq[startpos:].translate(stop_symbol = '@')
                            ptc_pos = ptc_pos_from_prot(prot, '@')
                            long_ptcs.append([exon_coord[1] - startpos - 3*x for x in ptc_pos])
                        else:
                            prot = seq[startpos:].translate()
                            ptc_pos = ptc_pos_from_prot(prot, '*')
                            ptc_coord = [exon_coord[1] - startpos - 3*x for x in ptc_pos]
                            for k in ptc_coord:
                                short.append(k)
                        if i == 2:
                            ptc_pos = ptc_pos_from_prot(prot, '@') + ptc_pos_from_prot(prot, '*')
                            close_to_ejc = [(len(prot) - 1)*3 - 3*x - len(leftover) < 56 for x in ptc_pos]
                            if sum(close_to_ejc) == len(close_to_ejc):
                                prot = seq[startpos:].translate(stop_symbol = '#')
                                long_ptcs = []
                                short = []
                        allprot = prot+allprot
                        short_ptcs = short_ptcs + short

                    #found a PTC in this transcript if any element but the last is a stop codon    
                    bool_ptc = "*" in allprot[:-1]
                    bool_long_exon = '@' in allprot[:-1]

                    



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
            seed: terminus is last two coordinates and the reading frame, 
            used for dynamic programming later"""
            if len(s) > 2:
                terminus = (s[-2],s[-1],leftover)
                
                if not bool_long_exon:
                    if terminus in dic_terminus['normal']:
                        dic_terminus['normal'][terminus].append(tuple(s))
                        continue
                    else:
                        dic_terminus['normal'][terminus] = [tuple(s)]
                else:
                    if terminus in dic_terminus['long_exon']:
                        dic_terminus['long_exon'][terminus].append([tuple(s), long_ptcs])
                        continue
                    else:
                        dic_terminus['long_exon'][terminus] = [[tuple(s), long_ptcs]]
            
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
        
        seed = new_seed
                    
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
                #last exon rule, part of 50-55nt rule
                if i == 0:
                    prot = seq[startpos:].translate(stop_symbol = '#')
                elif exlen + 1 > 407:
                    prot = seq[startpos:].translate(stop_symbol = '@')
                    #only allow long_exon tag to persist if PTCs introduced are solely present in a long exon 
                    ptc_pos = ptc_pos_from_prot(prot, '@')
                    ptc_coord = [exon_coord[0] + startpos + 3*x for x in ptc_pos]
                    check = [k in short_ptcs for k in ptc_coord]
                    if sum(check) > 0:
                        prot = seq[startpos:].translate()
                else:
                    prot = seq[startpos:].translate()

                #50-55nt rule, only apply changes if all stop codons in this exon pass this rule
                if i == 2:
                    ptc_pos = ptc_pos_from_prot(prot, '@') + ptc_pos_from_prot(prot, '*')
                    close_to_ejc = [(len(prot) - 1)*3 - 3*x - len(leftover) < 56 for x in ptc_pos]
                    if sum(close_to_ejc) == len(close_to_ejc):
                        prot = seq[startpos:].translate(stop_symbol = '#')
                allprot = prot+allprot

            else:
                seq = leftover+Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                if startpos > 0:                                                                                                    
                    leftover = seq[-startpos:]                                    
                else:
                    leftover = Seq("")
                seq = seq.reverse_complement() 
                if i == 0:
                    prot = seq[startpos:].translate(stop_symbol = '#')                                                                                                        
                elif exlen + 1 > 407:
                    prot = seq[startpos:].translate(stop_symbol = '@')
                    ptc_pos = ptc_pos_from_prot(prot, '@')
                    ptc_coord = [exon_coord[0] - startpos - 3*x for x in ptc_pos]
                    check = [k in short_ptcs for k in ptc_coord]
                    if sum(check) > 0:
                        prot = seq[startpos:].translate()
                else:
                    prot = seq[startpos:].translate()

                if i == 2:
                    ptc_pos = ptc_pos_from_prot(prot, '@') + ptc_pos_from_prot(prot, '*')
                    close_to_ejc = [(len(prot) - 1)*3 - 3*x - len(leftover) < 56 for x in ptc_pos]
                    if sum(close_to_ejc) == len(close_to_ejc):
                        prot = seq[startpos:].translate(stop_symbol = '#')

                allprot = prot+allprot 

        bool_ptc = "*" in allprot[:-1]
        bool_long_exon = '@' in allprot[:-1]
        bool_55nt_rule = '#' in allprot[:-1]
        
        """Quinn Comment: Classify seed + start codon as a passing path if no PTCs found in previous block of code"""
        if not bool_ptc:
            # all pass
            proteins.append("\t".join([gene_name,chrom,strand, "-".join([str(x) for x in s]), str(allprot)])+'\n')
            #print("ALL PASS %s"%(s))
            if bool_long_exon:
                path_pass['long_exon'].append(tuple(s))
            else:
                path_pass['normal'].append(tuple(s))
            for i in range(1, len(s), 2):
                j_coord = s[i:i+2]
                j_coord.sort()
                j_coord = tuple(j_coord)
                if not bool_long_exon: 
                    if j_coord not in junc_pass['normal']:
                        junc_pass['normal'][j_coord] = 0
                    junc_pass['normal'][j_coord] += 1
                else:
                    if j_coord not in junc_pass['long_exon']:
                        junc_pass['long_exon'][j_coord] = 0
                    junc_pass['long_exon'][j_coord] += 1

        

    #remove any paths that rely on a long_exon PTC that is also a short exon PTC
    for k in dic_terminus['long_exon']:
        to_remove = []
        for v in dic_terminus['long_exon'][k]:
            if len(v[0]) < 2: continue
            if sum([x in short_ptcs for x in v[1]]) > 0:
                to_remove.append(v)
        for rem in to_remove:
            dic_terminus['long_exon'][k].remove(rem)

    """Quinn Comment: OUT OF WHILE LOOP through all possible paths/seeds; 
    check all termini to see if they are part of a full path that has been classified as passing"""
    while True:
        new_paths = []
        for terminus in dic_terminus['normal']:
            terminus_pass = False
            for path_subset in dic_terminus['normal'][terminus]:
                for path in path_pass['normal']:
                    if path[:len(path_subset)] == path_subset:
                        terminus_pass = True
                        break
            #print(terminus, terminus_pass)

            """Quinn Comment: if our terminus is part of a passing path, we want to make sure if is reflected in passing paths and
            add the associate junctions to junc_pass, only if they are not present"""
            if terminus_pass:
                subsets_to_check = dic_terminus['normal'][terminus]
                for path_subset in subsets_to_check:
                    if path_subset in path_pass['normal']: continue
                    new_paths.append(path_subset)
                    path_pass['normal'].append(path_subset)
                    for i in range(1, len(path_subset), 2):
                        j_coord = list(path_subset[i:i+2])
                        j_coord.sort()
                        j_coord = tuple(j_coord)
                        if j_coord not in junc_pass['normal']:
                            junc_pass['normal'][j_coord] = 0
                            if verbose:
                                sys.stdout.write("junction pass:" + str(j_coord))

        """Quinn Comment: we could have a new path_pass added, so our while loop checks again to see if there are any new paths 
        that are now going to be passing considering our additions"""
        if len(new_paths) == 0:
            break

    #return long_exon terminus to standard state without PTCs
    for key in dic_terminus['long_exon']:
        to_keep = []
        for path in dic_terminus['long_exon'][key]:
            to_keep.append(path[0])
        dic_terminus['long_exon'][key] = to_keep
    #Now do the same for the long_exon termini and paths
    #We do not care if a terminus is long_exon or not, as long as it ends up leading to a long_exon path
    combined_keys = dic_terminus['long_exon'].keys() | dic_terminus['normal'].keys()
    combined_termini = {key: dic_terminus['long_exon'].get(key, []) + dic_terminus['normal'].get(key, []) for key in combined_keys}  
    while True:
        new_paths = []
        for terminus in combined_termini:
            terminus_pass = False
            for path_subset in combined_termini[terminus]:
                for path in path_pass['long_exon']:
                    if path[:len(path_subset)] == path_subset:
                        terminus_pass = True
                        break
            #print(terminus, terminus_pass)

            """Quinn Comment: if our terminus is part of a passing path, we want to make sure if is reflected in passing paths and
            add the associate junctions to junc_pass, only if they are not present"""
            if terminus_pass:
                subsets_to_check = combined_termini[terminus]
                for path_subset in subsets_to_check:
                    if path_subset in path_pass['long_exon']: continue
                    new_paths.append(path_subset)
                    path_pass['long_exon'].append(path_subset)
                    for i in range(1, len(path_subset), 2):
                        j_coord = list(path_subset[i:i+2])
                        j_coord.sort()
                        j_coord = tuple(j_coord)
                        if j_coord not in junc_pass['long_exon']:
                            junc_pass['long_exon'][j_coord] = 0
                            if verbose:
                                sys.stdout.write("junction long_exon:" + str(j_coord))

        """Quinn Comment: we could have a new path_pass added, so our while loop checks again to see if there are any new paths 
        that are now going to be passing considering our additions"""
        if len(new_paths) == 0:
            break
            
    return junc_pass, junc_fail, proteins

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
            ss2gene[(chrom, int(dic['start']))] = dic['gene_name'] # BED
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

    sys.stdout.write("done!\n")
    if verbose:
        sys.stdout.write("Processed: ")
        for chrstrand in dic_junc:
            sys.stdout.write(f"{len(dic_junc[chrstrand])} jxns on {chrstrand[0]} ({chrstrand[1]}).")

    
    # load or parse gtf annotations
    # g_coords: gene coordinates, grouped by chromosome and strand
    # g_info: a dictionary with (transcript_name, gene_name) as keys, and intron info as values
    try: 
        sys.stdout.write("Loading annotations...\n")
        parsed_gtf = f"{rundir}/{gtf_annot.split('/')[-1].split('.gtf')[0]}_SJC_annotations.pckle"
        with open(parsed_gtf, 'rb') as f:
            g_coords, g_info = pickle.load(f)
        sys.stdout.write("done!\n")
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
    fout.write("\t".join(["Gene_name","Intron_coord","Annot","Coding","UTR","Long_exon"])+'\n')
    
    for gene_name, chrom, strand in gene_juncs:
        if gene_name != 'CCNL2':
            continue
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

        print(chrom)
        print(strand)
        print(junctions)
        print(start_codons)
        print(stop_codons)
        print(gene_name)
        junc_pass, junc_fail, proteins = solve_NMD(chrom,strand,junctions, 
                                                   start_codons, stop_codons, 
                                                   gene_name)
        
        for j in junctions:
            bool_pass = j in junc_pass['normal'] or j in g_info[gene_name]['pcjunctions']
            bool_fail = j in junc_fail or j in junc_pass['long_exon']
            utr = False
            long_exon = j in junc_pass['long_exon'] and not bool_pass
            if not bool_pass:
                # Check that it's not in UTR                
                utr = check_utrs(j,g_info[gene_name]['utrs'])

            if bool_fail or bool_pass:
                tested = True
            else:
                tested = False
            annotated = j in g_info[gene_name]['junctions']
            #if not bool_pass and annotated:
            #print("%s %s %s junction: %s tested: %s utr: %s coding: %s annotated: %s "%(chrom, strand, gene_name, j, tested,utr, bool_pass, annotated))
            
            fout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]}',
                                  str(annotated), str(bool_pass), str(utr), str(long_exon)])+'\n')
        

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
        sys.stderr.write("Error: no LeafCutter junction file provided...\npython SpliceJunctionClassifier.py -j Leafcutter_perind.counts.gz\n")
        exit(0)
    
    global fa
    sys.stdout.write(f"Loading genome {options.genome} ...")
    fa = pyfastx.Fasta(options.genome)
    sys.stdout.write("done!\n")

    ClassifySpliceJunction(options)

    # for testing
    sjc_file = f"{options.rundir}/{options.outprefix}_junction_classifications.txt"
    print(f"Merging discordant logics in {sjc_file}...")
    # print the first 10 items in the merged dictionary
    sjc = merge_discordant_logics(sjc_file)
    print(list(sjc.items())[:10])


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
