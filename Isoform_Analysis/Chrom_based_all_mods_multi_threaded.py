#!/usr/bin/env python3
#Similar to hypermodDistance
import sys
import os
import csv
import re
import itertools
import argparse
import subprocess
#from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import time
import concurrent.futures
#import mpi4py.futures

class Filereader :
    ''' 
    Extracts the users event.align file pertaining to a particular
    gene and stores it as a csv file.
    '''
    def __init__ (self, fname, targets):     #*args):
        '''Contructor: saves attribute fname '''
        self.fname = fname
        self.targets = targets     #('RPL15'), 'SLC25A3')
        self.meta_data = [['Gene', 'Isoform', 'Chromosome', 'BasePosition', 'ReferenceBase', 'MiscalledBase', 'ReadCoverage', '+correct', '-correct', '+error', '-error', 'percent_error']]
        self.slicer = []
        self.chrom = []
        self.pos = [] #could be multiple
        self.mod = []
        self.gene = [] 
        self.gene_ID = []
        self.direct_mod = [] #could be multiple
        self.IVT_mod = [] #could be multiple
        self.kmer = [] #could be multiple
        self.nsites = [] #could be multiple
        
        #self.experiments = [('Direct_1st', 'hela_D1.fq')] #older version 
        self.experiments = [('HeLa_merged_D1_D2', 'HeLa_D_1_2_merged.fq')] 
        #print(self.targets)
    def openf (self):
        '''Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname) 
    
    def prepareMeta (self):
        '''Parses gene information from csv file'''
        with self.openf() as File:
            file_reader = csv.reader(File, delimiter=',')
            n = range(0,616) #150, 395
            interestingrows=[row for idx, row in enumerate(file_reader) if idx in n]
            for line in interestingrows:
                #if line[0] == 'Annotation':
            #for line in file_reader:
                if line[0].lower() == 'gene':
                    print("First line out")
                    continue
                elif 'All' in self.targets:
                    T = re.split('/', line[0])
                    P = re.split(',', line[1])
                    if len(T) > 1:
                        Gene = T[1]
                    elif len(T) == 1:
                        Gene = T[0]
                    else:
                        Gene = T[0]
                    pos = P
                    print(pos)
                    chrom = line[2]#chrom = line[7]
                    start = line[11]#start = line[8]
                    end = line[12]#end = line[9]
                    D_mod = re.split(',', line[3])
                    IVT_mod = re.split(',', line[4])
                    kmer = re.split(',', line[9])
                    n_sites = re.split(',', line[10])
                    #b = f'"{a}"'
                    #Gene = f'"{Gene}"'
                    GENE = f'"{Gene}"' #added to make sure correct gene ID is pulled
                    print("GENE: " + str(Gene))
                   # with open('output.bam', 'w') as fh_bam:
                   #     param_1 = ["samtools", "view", "-h", "-Sb", "../IVT_1st_hg38.sorted.bam", Slicer]
                   #     output1 = Popen(param_1, stdin=PIPE, stdout=fh_bam)
                   #     output1.wait()
                   #Popen("mycmd" + " myarg", shell=True).wait()
                   #c = "grep {} ../gencode.v36.annotation.gtf > EID.txt".format(Gene)
                    subprocess.Popen("grep '{}' ../gencode.v36.annotation.gtf > EID.txt".format(GENE), shell=True).wait()
                    grab_ID = 'EID.txt'
                    with open(grab_ID, 'r') as f:
                        ###############
                        #n = range(0,1)
                        #interestingrows=[row for idx, row in enumerate(f) if idx in n]
                        #for line in interestingrows:
                        ###############
                        for line in f:
                            EID = re.split(';|\t| ', line)
                            EID = EID[9]
                            EID = re.split('"', EID)
                            EID = EID[1]
                    print("Gene ID that was parsed: " + str(EID))
                    #EID = line[1]
                    #Gene = line[0]
                    #pos1 = line[3]
                    #pos2 = line[4]
                    self.slicer.append(chrom + ':' + start + '-' + end)
                    self.gene.append(Gene)
                    self.chrom.append(chrom)
                    self.gene_ID.append(EID)
                    self.pos.append(pos) #could be multiple
                    self.direct_mod.append(D_mod) #could be multiple
                    self.IVT_mod.append(IVT_mod) #could be multiple
                    self.kmer.append(kmer) #could be multiple
                    self.nsites.append(n_sites)
                    subprocess.Popen("rm EID.txt".format(Gene), shell=True).wait()
                elif line[0] in self.targets:
                    T = re.split('/', line[0])
                    P = re.split(',', line[1])
                    if len(T) > 1:
                        Gene = T[1]
                    elif len(T) == 1:
                        Gene = T[0]
                    else:
                        Gene = T[0]
                    pos = P
                    print(pos)
                    chrom = line[2]#chrom = line[7]
                    start = line[11]#start = line[8]
                    end = line[12]#end = line[9]
                    D_mod = re.split(',', line[3])
                    IVT_mod = re.split(',', line[4])
                    kmer = re.split(',', line[9])
                    n_sites = re.split(',', line[10])
                    Gene = "'" + '"' + Gene + '"' + "'" #added to make sure correct gene ID is pulled                   
                   # with open('output.bam', 'w') as fh_bam:
                   #     param_1 = ["samtools", "view", "-h", "-Sb", "../IVT_1st_hg38.sorted.bam", Slicer]
                   #     output1 = Popen(param_1, stdin=PIPE, stdout=fh_bam)
                   #     output1.wait()
                   #Popen("mycmd" + " myarg", shell=True).wait()
                   #c = "grep {} ../gencode.v36.annotation.gtf > EID.txt".format(Gene)
                    subprocess.Popen("grep {} ../gencode.v36.annotation.gtf > EID.txt".format(Gene), shell=True).wait()
                    grab_ID = 'EID.txt'
                    with open(grab_ID, 'r') as f:
                        for line in f:
                            EID = re.split(';|\t| ', line)
                            EID = EID[9]
                            EID = re.split('"', EID)
                            EID = EID[1]
                    #EID = line[1]
                    #Gene = line[0]
                    #pos1 = line[3]
                    #pos2 = line[4]
                    self.slicer.append(chrom + ':' + start + '-' + end)
                    self.gene.append(Gene)
                    self.chrom.append(chrom)
                    self.gene_ID.append(EID)
                    self.pos.append(pos) #could be multiple
                    self.direct_mod.append(D_mod) #could be multiple
                    self.IVT_mod.append(IVT_mod) #could be multiple
                    self.kmer.append(kmer) #could be multiple
                    self.nsites.append(n_sites)
                    subprocess.Popen("rm EID.txt".format(Gene), shell=True).wait()
            print(self.gene) #order is perserved in gene
            print(self.targets)
        for i in self.experiments:
            print(i[0])
    
    def highestExpression (self, idx1, idx2):
        #os.system("mkdir Merged_Isoform_Data_Updated")
        dominant_isoforms = []
        #for frag, site, gen in zip(self.slicer, self.mod, self.gene):
        batchdir=f"/scratch/makhamreh.a/Batches/{idx1}_{idx2}_Batch/"
        os.system( f"mkdir /scratch/makhamreh.a/Batches/{idx1}_{idx2}_Batch" )
        for frag, gID, gen, P, D_mod, IVT_mod, kmer, nsites in zip(self.chrom[idx1:idx2], self.gene_ID[idx1:idx2], self.gene[idx1:idx2], self.pos[idx1:idx2], self.direct_mod[idx1:idx2], self.IVT_mod[idx1:idx2], self.kmer[idx1:idx2], self.nsites[idx1:idx2]):
            for experiment in self.experiments:
                all_isoforms = []
                subprocess.Popen( "samtools view -h -Sb ../{}.bam {} > {}/output.bam".format(experiment[0], frag, batchdir), shell=True).wait()
                subprocess.Popen( f"samtools index {batchdir}/output.bam", shell=True).wait()
                subprocess.Popen( f"samtools view -h -F 256 -q 10 {batchdir}/output.bam > {batchdir}/output.sorted.bam", shell=True).wait()
                subprocess.Popen( f"python ../bin/bam2Bed12.py --input_bam {batchdir}/output.sorted.bam > {batchdir}/output.bed", shell=True).wait()
            
                os.chdir(batchdir)
                #subprocess.call('ls', shell=True, cwd='path/to/wanted/dir/')
                subprocess.Popen( f"python /home/makhamreh.a/rna_MODS/flair/flair/flair.py correct -q output.bed -g /home/makhamreh.a/rna_MODS/flair/flair/GRCh38.p10.genome.fa -f /home/makhamreh.a/rna_MODS/flair/flair/gencode.v36.annotation.gtf", shell=True).wait()
                subprocess.Popen( f"python /home/makhamreh.a/rna_MODS/flair/flair/flair.py collapse -g /home/makhamreh.a/rna_MODS/flair/flair/GRCh38.p10.genome.fa -r /home/makhamreh.a/rna_MODS/flair/flair/{experiment[1]} -q flair_all_corrected.bed -f /home/makhamreh.a/rna_MODS/flair/flair/gencode.v36.annotation.gtf --generate_map", shell=True).wait()
                #os.system( "samtools view -h -o output_{}.sam output.bam".format(experiment[0]))
                
                #ids = 'flair.collapse.isoform.read.map.txt'
                #all_isoforms.append(self.isoformReads(ids))
                
                subprocess.Popen(f"grep {gID} flair.collapse.isoform.read.map.txt > updated.read.map.txt", shell=True).wait()
                ids = 'updated.read.map.txt'
                all_isoforms = self.isoformReads(ids)
                high = self.domFinder(all_isoforms)
                #################
                P_Site = [] 
                splice_Distance = []
                location_Length = []
                Bias = []
                perc_mod = []
                n_sites = []
                
                for d_mod, ivt_mod in zip(D_mod, IVT_mod):
                    #p_diff = float(d_mod) - float(ivt_mod)
                    #perc_Diff.append(round(p_diff,3))
                    d_mod = float(d_mod)
                    perc_mod.append(d_mod)
                for n in nsites:
                    n_sites.append(int(n))
                    
                
                if "ENST" in str(high):
                    subprocess.Popen("grep {} /home/makhamreh.a/rna_MODS/flair/flair/gencode.v36.annotation.gtf > GTF_Info.txt".format(high), shell=True).wait()
                    meta_transcript = 'GTF_Info.txt'
                    
                    with open(meta_transcript, 'r') as f:
                        lines = f.readlines()
                        for position in P:
                            floating_pos = 0
                            for idx, line in enumerate(lines):
                                x = re.split(';|\t', line)
                                location = x[2]
                                if location == 'transcript':
                                    full_transcript_start = int(x[3])
                                    full_transcript_end = int(x[4])
                                print("Gene: " + str(gen))
                                print("location: " + str(location))
                                print("Full start: " + str(full_transcript_start))
                                print("Full end: " + str(full_transcript_end))
                                loc_start = int(x[3])
                                loc_end = int(x[4])
                                strand_direction = x[6]
                                print("loc_start: " + str(loc_start) )
                                print("loc_end: " + str(loc_end))
                            
                            ###Here is where we break the multiple positions down
                            
                                if int(position) < full_transcript_start or int(position) > full_transcript_end:
                                    print("HIT")
                                    P_Site.append(0)
                                    splice_Distance.append(0)
                                    location_Length.append(0)
                                    Bias.append(0)
                                    break 
                                #Fix code to deal with positions being considered twice if on start or stop codon, they are also considered UTR or CDS!!!!   
                                if location == 'CDS' or location == 'start_codon' or location == 'stop_codon':
                                    if int(position) >= loc_start and int(position) <= loc_end:
                                        floating_pos = 1
                                        location_length = loc_end - loc_start
                                        P_Site.append(location)
                                        print("Hit: " + str(location))
                                        location_Length.append(location_length)
                                        x1 = abs(int(position)-loc_start)
                                        x2 = abs(int(position)-loc_end)
                                        if x1 < x2 and strand_direction == '+':
                                            splice_Distance.append(x1)
                                            bias = '5 prime end'
                                            Bias.append(bias)
                                        elif x2 < x1 and strand_direction == '+':
                                            splice_Distance.append(x2)
                                            bias = '3 prime end'
                                            Bias.append(bias)
                                        elif x1 < x2 and strand_direction == '-':
                                            splice_Distance.append(x1)
                                            bias = '3 prime end'
                                            Bias.append(bias)
                                        elif x2 < x1 and strand_direction == '-':
                                            splice_Distance.append(x2)
                                            bias = '5 prime end'
                                            Bias.append(bias)
                                        else:
                                            splice_Distance.append(x1)
                                            bias = '5 prime end'
                                            Bias.append(bias)
                                elif location == 'UTR': 
                                    if int(position) >= loc_start and int(position) <= loc_end:
                                        floating_pos = 1
                                        location_length = loc_end - loc_start
                                        location_Length.append(location_length)
                                        x1 = abs(int(position)-loc_start)
                                        x2 = abs(int(position)-loc_end)
                                        if strand_direction == '+' and x1 < x2:
                                            splice_Distance.append(x1)
                                            if abs(loc_start - full_transcript_start) < abs(loc_start - full_transcript_end):
                                                location = '5 prime ' + location
                                                bias = '5 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                            else:
                                                location = '3 prime ' + location
                                                bias = '5 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                        elif strand_direction == '+' and x2 < x1:
                                            splice_Distance.append(x2)
                                            if abs(loc_start - full_transcript_start) < abs(loc_start - full_transcript_end):
                                                location = '5 prime ' + location
                                                bias = '3 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                            else:
                                                location = '3 prime ' + location
                                                bias = '3 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                        elif strand_direction == '-' and x1 < x2:
                                            splice_Distance.append(x1)
                                            if abs(loc_start - full_transcript_start) < abs(loc_start - full_transcript_end):
                                                location = '3 prime ' + location
                                                bias = '3 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                            else:
                                                location = '5 prime ' + location
                                                bias = '3 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                        elif strand_direction == '-' and x2 < x1:
                                            splice_Distance.append(x2)
                                            if abs(loc_start - full_transcript_start) < abs(loc_start - full_transcript_end):
                                                location = '3 prime ' + location
                                                bias = '5 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                            else:
                                                location = '5 prime ' + location
                                                bias = '5 prime end'
                                                P_Site.append(location)
                                                Bias.append(bias)
                                        else:
                                            P_Site.append(0)
                                            splice_Distance.append(0)
                                            location_Length.append(0)
                                            Bias.append(0)
                                else: 
                                    if line == lines[-1] and floating_pos == 0 and (int(position) < full_transcript_start or int(position) > full_transcript_end):
                                        P_Site.append(0)
                                        splice_Distance.append(0)
                                        location_Length.append(0)
                                        Bias.append(0)
                                    
# =============================================================================
#                                     P_Site.append(0)
#                                     splice_Distance.append(0)
#                                     location_Length.append(0)
#                                     Bias.append(0)
# =============================================================================
                                
                else:
                    P_Site.append(0)
                    splice_Distance.append(0)
                    location_Length.append(0)
                    Bias.append(0)
            

            dominant_isoforms.append([gen, high, P_Site, splice_Distance, location_Length, Bias, perc_mod, kmer, n_sites])
                
                #os.system("rm GTF_Info.txt")                 
                #################
                #dominant_isoforms.append(high)
            os.chdir("/home/makhamreh.a/rna_MODS/flair/flair/isoform_bam_prep/")
        return dominant_isoforms
            
            
                
    def isoformReads (self, ID):
        isoforms = []
        with open(ID, 'r') as f:
                for line in f:
                    x = re.split(',|\t|_', line)
                    transcript = x[0]
                    gene = x[1]
                    reads = x[2:]
                    isoforms.append((transcript, gene, reads))
        return isoforms
    
    def domFinder (self, isoforms_Dir1):
        previous_sz = 0
        dom = 0
        if len(isoforms_Dir1) == 1:
            dom = isoforms_Dir1[0][0]
        else:
            for j in isoforms_Dir1:
                print("Current Isoform: " + str(j[0]))
                current_transcript = j[0]
                current_sz = len(j[2])#number of reads in isoform
                print(" A "+ str(current_sz))
                print(" B " + str(previous_sz))
                if current_sz > previous_sz:
                    dom = current_transcript
                    print(type(dom), dom)
                    previous_sz = current_sz
                    print("Switch " + str(dom))
                else:
                    continue
        return dom
         
                    
        
# =============================================================================
#     def compareIsoforms (self, isoforms_Dir1, isoforms_Dir2, isoforms_IVT):
#         ISO = []
#         for i in range(0, len(isoforms_Dir1), 1):
#             trans_check = isoforms_Dir1[i][0]
#             for t in range(0, len(isoforms_Dir2), 1):
#                 if isoforms_Dir2[t][0] == trans_check:
#                     ISO.append([trans_check, isoforms_Dir1[i][2], isoforms_Dir2[t][2]])
# 
#         for i in range(0, len(isoforms_IVT), 1):
#             trans_check = isoforms_IVT[i][0]
#             for t in range(0, len(ISO), 1):
#                 if ISO[t][0] == trans_check:
#                     ISO[t].insert(3, isoforms_IVT[i][2])
#         for i in ISO:
#             if len(i) < 4:
#                 i.insert(3, None)
#         return ISO
# =============================================================================
    
        
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fname', type = str, default = False,
                        help = 'What file are you parsing?')
    parser.add_argument('--targets', nargs = '*', type = str, default = False,
                        help = 'What genes do you want to analyze?')
    parser.add_argument('--batches', nargs = '*', type = str, default = False,
                        help = 'How do you wish to batch the list of gene targets?')
    #parser.add_argument('--kmer', type = str, default = False,
    #                    help = 'What Kmer are you parsing out from your file?')
    args = parser.parse_args()
    main_prep = Filereader(args.fname, args.targets)
    prep_meta = main_prep.prepareMeta()
    field = ['Gene','Transcript ID', 'P_location', 'P_Splice_Distance', 'Length of CDS/UTR', 'Splice Bias', 'Percent_Diff', 'Kmer', 'nsites']
    with open('dominant_isoform_combined.csv', 'w') as f:
        write = csv.writer(f)
        write.writerow(field)
    print(args.fname)
    print(args.targets)
    print(args.batches)
    
    if args.batches: #If the user wishes to run their job in batches
        batches = args.batches[0].split(",")
        i1 = int(batches[0])
        i2 = int(batches[1])
        batch_sz = int(batches[2])
        batch_List = []
        x = range(i1, i2, batch_sz)
        for idx, n in enumerate(x):
            if n != x[-1]:
                Slice = (n,n+batch_sz)
                batch_List.append(Slice) #batch list contains tuples of each batch
            else:
                Slice = (x[-2]+batch_sz, i2)
                batch_List.append(Slice)
        print(f"Prepared batch list: {batch_List}")
        start_process = time.time()
        with concurrent.futures.ProcessPoolExecutor() as executer:
        #with mpi4py.futures.MPIPoolExecutor() as executer:
            results = executer.map(main_prep.highestExpression, [i[0] for i in batch_List], [i[1] for i in batch_List])
            for result in results:
                with open('dominant_isoform_combined.csv', 'a') as f:
                    write = csv.writer(f)
                    for i in result:
                        write.writerow([i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8]])
        end_process = time.time()
        print(f"Batching time: {end_process-start_process}")
                
        #batching = map(Filereader.highestExpression, [i[0] for i in batch_List], [i[1] for i in batch_List])
        #batching = Filereader.highestExpression(idx1, idx2)
        #result = Filereader(args.fname, args.targets)
    #print('Targets: {}'.format(args.targets))
    #p = result.prepareMeta()
    #x = result.highestExpression()
    #s = result.filteredAlign()
    #field = ['Gene','Transcript ID', 'P1_location', 'P1_Splice_Distance', 'Length of CDS/UTR', 'Splice Bias']
    #with open('dominant_isoform_combined.csv', 'w') as f:
    #    write = csv.writer(f)
    #    write.writerow(field)
        #write.writerow(x)
    #    for i in x:
    #        write.writerow([i[0], i[1], i[2], i[3], i[4], i[5]])
        
    #print(p)
    #p, q = result.Print()
    #print("You have {} instances of this kmer in your sample set\
    #      and {} instances of this kmer in your control".format(p, q))
if __name__ == '__main__':
    main()






