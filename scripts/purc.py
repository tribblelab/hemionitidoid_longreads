#!/usr/bin/env python

logo = """
-------------------------------------------------------------
|                            PURC                           |
|        Pipeline for Untangling Reticulate Complexes       |
|                        version 2                          |
|        https://bitbucket.org/peter_schafran/purc          |
|                                                            |
|      Fay-Wei Li & Carl J Rothfels & Peter W Schafran      |
-------------------------------------------------------------
"""

usage = """
You need to provide a configuration file.

Usage: ./purc.py configuration_file
Example: ./purc.py ppp_configuration.txt
For more info, try: ./purc.py -help

"""

citation = """
This script relies heavily on VSEARCH, DADA2, MAFFT, LIMA, and BLAST.
If this script assisted with a publication, please cite the following papers
(or updated citations, depending on the versions of VSEARCH, etc., used).

PURC v1
-Rothfels, C.J., K.M. Pryer, and F-W. Li. 2017. Next-generation polyploid
phylogenetics: rapid resolution of hybrid polyploid complexes using PacBio
single-molecule sequencing. New Phytologist 213: 413-429.

PURC v2
-Schafran, P., F-W. Li, and C.J. Rothfels. PURC v2 provides improved sequence
inference for polyploid phylogenetics and other manifestations of the
multiple-copy problem. Submitted.

BLAST
-Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, et al. 2009.
BLAST+: Architecture and applications. BMC Bioinformatics 10: 421.

Cutadapt
-Martin, M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads.
EMBnet.journal 17:10-12.

DADA2
-Callahan, B.J., P.J. McMurdie, M.J. Rosen, A.W. Han, A.J.A. Johnson, and Susan P. Holmes.
2016. DADA2: High-resolution sample inference from Illumina amplicon data.
Nature Methods 13: 581-583.

MAFFT
-Katoh, K., and D.M. Standley. 2013. MAFFT multiple sequence alignment software version 7:
improvements in performace and usability. Molecular Biology and Evolution 30: 772-780.

VSEARCH
-Rognes, T., T. Flouri, B. Nichols, C. Quince, and F. Mahé. 2016. VSEARCH: a versatile
open source tool for metagenomics. PeerJ 4:e2584.


Deprecated PURC v1 dependencies:
MUSCLE
-Edgar, R.C. 2004. MUSCLE: Multiple sequence alignment with high accuracy and high throughput.
Nucleic Acids Research 32:1792-1797.

USEARCH/UCLUST
-Edgar, R.C. 2010. Search and clustering orders of magnitude faster than BLAST.
Bioinformatics 26(19), 2460-2461.

UCHIME:
-Edgar, R.C., B.J. Haas, J.C. Clemente, C. Quince, R. Knight. 2011.
UCHIME improves sensitivity and speed of chimera detection, Bioinformatics 27(16), 2194-2200.
"""

import re
import sys
import os
import glob
import subprocess
import shutil
import time
import datetime
import io
import fileinput
import platform
import pickle
from os import listdir
from os.path import isfile, join

try:
    from Bio import SeqIO
    from Bio import AlignIO
    from Bio.Align import AlignInfo
except:
    sys.exit('ERROR: could not import BioPython; please install BioPython first')

def convertTime(time): # returns a string with time.time converted to minutes, hours, or days
    if time > 60:
        time = time/60
        if time > 60:
            time = time/60
            if time > 24:
                time = time/24
                convertedTime = "%.2f days" % time
            else:
                convertedTime = "%.2f hours" % time
        else:
            convertedTime = "%.2f minutes" % time
    else:
        convertedTime = "%.2f seconds" % time
    return convertedTime

def parse_fasta(infile):
    """Reads in a fasta, returns a list of biopython seq objects"""
    AllSeq = SeqIO.parse(infile, 'fasta')
    return [i for i in AllSeq]

def parse_fastq(infile):
    """Reads in a fastq, returns a list of biopython seq objects"""
    AllSeq = SeqIO.parse(infile, 'fastq')
    return [i for i in AllSeq]

def rename_fasta(infile):
    """Makes sure there are no ';' or '='' or '/'' characters in the fasta file, so that they won't confuse blast"""
    prefix = infile.replace('.fasta', '').replace('.fas', '').replace('.fa', '').replace('.txt', '')
    path = "/".join(prefix.split("/")[:-1])
    filename = prefix.split("/")[-1]
    outfile = "%s/tmp/%s_renamed.fasta" %(Output_folder, filename)
    #outfile = "%s/tmp/%s_renamed.fasta" %(Output_folder, prefix)
    sed_cmd = "sed 's/;/_/g;s/=/_/g;s/\//_/g' %s > %s" % (infile, outfile)
    process = subprocess.Popen(sed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.wait()
    (out, err) = process.communicate()
    return outfile

def rename_fastq(infile):
    """Makes sure there are no ';' or '='' or '/'' characters in the fastq file"""
    replace = [';', '=', '/']
    prefix = infile.replace('.fastq', '').replace('.fq', '')
    path = "/".join(prefix.split("/")[:-1])
    filename = prefix.split("/")[-1]
    outfile = "%s/tmp/%s_renamed.fastq" %(Output_folder, filename)
    with open(infile, "r") as open_infile:
        with open(outfile, "w") as open_outfile:
            linecount = 0
            for line in open_infile:
                if linecount % 4 == 0:
                    for char in replace:
                        line = line.replace(char, '_')
                open_outfile.write(line)
                linecount += 1
    return outfile

def check_fasta(infile):
    """Check if this is a good fasta file, if so return True, if not return False"""
    for line in open(infile, 'r'):
        if line.startswith('>'):
            return True
        else:
            return False

def count_seq_from_fasta(infile):
    """Yep. Just returns the number of sequence contain in a fasta file"""
    cmdLine = "grep '>' %s | wc -w" % infile
    process = subprocess.Popen(cmdLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    process.wait()
    (out, err) = process.communicate()
    out = out.strip('\n').replace(' ', '')
    return out

def count_seq_from_fastq(infile):
    cmdLine = "awk 'END {print NR/4}'"
    process = subprocess.Popen(cmdLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    process.wait()
    (out, err) = process.communicate()
    out = out.strip('\n').replace(' ', '')
    return out

def convert_fastq_to_fasta(infile):
    prefix = infile.replace('.fastq', '').replace('.fq', '')
    outfile = prefix + ".fasta"
    SeqIO.convert(infile, 'fastq', outfile, 'fasta')
    return outfile

def subset_fasta_seqs_from_fastq(fasta_file, fastq_file):
    prefix = ".".join(fasta_file.split(".")[:-1])
    outfile = "%s.fastq" % prefix
    AllSeq = SeqIO.parse(fastq_file, "fastq")
    seqNames = []
    with open(fasta_file, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                seqNames.append(line.strip(">\n").split("|")[-1])
    outseqs = [i for i in AllSeq if i.id in seqNames]
    SeqIO.write(outseqs, outfile, "fastq")

def ReverseComplement(seq):
    """Returns reverse complement sequence, ignores gaps"""
    seq = seq.replace(' ','')
    seq = seq[::-1] # Reverse the sequence
    basecomplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y', 'Y':'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W', 'H': 'D', 'D': 'H', 'B': 'V', 'V': 'B'} # Make a dictionary for complement
    letters = list(seq)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)

def reorientSeqs(sequence_file, refseq_blastDB, num_threads): # function needed after lima demultiplexing to orient sequences so all are facing the same way. BLAST-based demultiplexing does this within its function
    print("Reorienting sequences based on references...")
    filePrefix = ".".join(sequence_file.split(".")[:-1])
    basename = filePrefix.split("/")[-1]
    fileExt = sequence_file.split(".")[-1]
    if fileExt in ["fastq", "fq"]:
        fastq_sequence_file = sequence_file
        fasta_sequence_file = convert_fastq_to_fasta(sequence_file)
    else:
        fasta_sequence_file = sequence_file
    # use BLAST against ref sequences to determine orientation of reads
    blastoutfile = "%s/tmp/blast_reorient_out.txt" % Output_folder
    #print("Reorienting: %s" % fasta_sequence_file)
    #print("BlastDB: %s" % refseq_blastDB)
    blastCMD = "blastn -query %s -db %s -num_threads %s -max_target_seqs 1 -out %s -outfmt 6 -max_hsps 1" % (fasta_sequence_file, refseq_blastDB, num_threads, blastoutfile)
    process = subprocess.Popen(blastCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    stdout,stderr = process.communicate()
    reorientList = []
    with open(blastoutfile, "r") as open_blast_file:
        for line in open_blast_file:
            seqname = line.strip("\n").split("\t")[0]
            sstart = int(line.strip("\n").split("\t")[8])
            send = int(line.strip("\n").split("\t")[9])
            evalue = line.strip("\n").split("\t")[10]
            if send < sstart:
                reorientList.append(seqname)
    if len(reorientList) >= 1:
        with open("%s/tmp/reorient_list.txt" % Output_folder, "w") as open_orient_list:
            for name in reorientList:
                open_orient_list.write("%s\n" % name)
        print("\tReorienting %s sequences..." % len(reorientList))
        fastaSeqs = SeqIO.parse(fasta_sequence_file, 'fasta')
        with open("%s/tmp/%s.reoriented.fa" %(Output_folder, basename), "w") as open_out_fasta:
            print("\tWriting reoriented sequences to %s/tmp/%s.reoriented.fa" %(Output_folder, basename))
            for read in fastaSeqs:
                if read.id in reorientList:
                    open_out_fasta.write(">%s\n" % read.id)
                    open_out_fasta.write("%s\n" % ReverseComplement(read.seq))
                else:
                    open_out_fasta.write(">%s\n" % read.id)
                    open_out_fasta.write("%s\n" % read.seq)
        reorientedSequences = "%s/tmp/%s.reoriented.fa" %(Output_folder, basename)
        if fileExt in ["fastq", "fq"]:
            with open(fastq_sequence_file, "r") as open_fastq_file:
                with open("%s/tmp/%s.reoriented.fq" %(Output_folder,basename), "w") as open_out_fastq:
                    print("\tWriting reoriented sequences to %s/tmp/%s.reoriented.fq" % (Output_folder, basename))
                    lineCount = 0
                    reverseLine = 0
                    for line in open_fastq_file:
                        if lineCount == 0 or (lineCount % 4) == 0:
                            if line.strip("@\n") in reorientList:
                                reverseLine = 1
                                open_out_fastq.write(line)
                            else:
                                reverseLine = 0
                                open_out_fastq.write(line)
                        elif lineCount == 1 or ((lineCount-1)%4) == 0:
                            if reverseLine == 1:
                                reversed_seq = ReverseComplement(line.strip("\n"))
                                open_out_fastq.write("%s\n" % reversed_seq)
                            else:
                                open_out_fastq.write(line)
                        elif lineCount == 2 or ((lineCount-2)%4) ==  0:
                            open_out_fastq.write(line)
                        elif lineCount == 3 or ((lineCount-3)%4) == 0:
                            if reverseLine == 1:
                                reversed_qscore = line.strip("\n")[::-1]
                                open_out_fastq.write("%s\n" % reversed_qscore)
                            else:
                                open_out_fastq.write(line)
            reorientedSequences = "%s/tmp/%s.reoriented.fq" %(Output_folder, basename)
    else:
        reorientedSequences = sequence_file
    return reorientedSequences

def checkDuplicateBC(barcode_seq_filename):
    """Checks for duplicate seqs (including reverse complements in barcode file, which will cause lima to fail)"""
    BCdict = SeqIO.parse(barcode_seq_filename, 'fasta')
    BCseqdict = {} # Key: BC sequence. Entry: list of BCs with sequence (including reverse complements)
    duplicateBClist = []
    deduplicatedBClist = []
    BCpairdict = {}
    dupesFound = 0
    for i in BCdict:
        rc = ReverseComplement(i.seq)
        if i.seq not in BCseqdict and rc not in BCseqdict:
            BCseqdict.update({i.seq : [i.id]})
        else:
            dupesFound = 1
            duplicateBClist.append(i.id)
            if i.seq in BCseqdict:
                print("WARNING: Duplicate barcodes found! %s and %s " %(i.id, BCseqdict[i.seq][0]))
                BCseqdict[i.seq].append(i.id)
                duplicateBClist.append(BCseqdict[i.seq][0])
                BCpairdict.update({i.id : BCseqdict[i.seq][0]})
                BCpairdict.update({BCseqdict[i.seq][0] : i.id})
            elif rc in BCseqdict:
                print("WARNING: Duplicate reverse complement barcodes found! %s and %s " %(i.id, BCseqdict[rc][0]))
                BCseqdict[rc].append(i.id)
                duplicateBClist.append(BCseqdict[rc][0])
                BCpairdict.update({i.id : BCseqdict[rc][0]})
                BCpairdict.update({BCseqdict[rc][0] : i.id})
    duplicateBClist = set(duplicateBClist)
    if dupesFound == 1:
        print("Remove duplicates to run lima, or change config file to use BLAST demultiplexing (Lima_override = 1)")
        print("Will attempt to run by splitting same/different dual barcodes then rejoining before annotation")
    for seq in BCseqdict:
        if len(BCseqdict[seq]) > 1:
            deduplicatedBClist.append(BCseqdict[seq][0])
    BCdict = SeqIO.parse(barcode_seq_filename, 'fasta') # easiest way to get back to top of file
    nonduplicatedBC = [i for i in BCdict if i.id not in deduplicatedBClist]
    BCdict = SeqIO.parse(barcode_seq_filename, 'fasta') # easiest way to get back to top of file
    deduplicatedBC = [i for i in BCdict if i.id in deduplicatedBClist]
    SeqIO.write(nonduplicatedBC, "%s/tmp/nonduplicated_BCs.fasta" % Output_folder, "fasta")
    SeqIO.write(deduplicatedBC, "%s/tmp/deduplicated_BCs.fasta" % Output_folder, "fasta")
    return dupesFound,BCpairdict

def makeBlastDB(inFileName, outDBname):
    """Makes a blast database from the input file"""

    # remove any '-' in the sequence, so that BLAST won't freak out
    path = "/".join(inFileName.split("/")[:-1])
    basename = inFileName.split("/")[-1]
    seq_no_hyphen = open("%s/tmp/%s.nohyphen.fasta" % (Output_folder, basename), 'w')
    makeDBfail = "FALSE"
    for i in parse_fasta(inFileName):
        new_seq = str(i.seq).replace('-','')
        seq_no_hyphen.write('>' + str(i.id) + '\n' + new_seq + '\n')
        if len(i.id) > 50:
            print("ERROR: %s has sequence name too long for BLAST. Reduce to < 50 characters" % inFileName)
            print("Offending line: %s" % i.id)
            makeDBfail = "TRUE"
    if makeDBfail == "TRUE":
        sys.exit(1)
    seq_no_hyphen.close()
    makeblastdb_cmd = "makeblastdb -in %s/tmp/%s.nohyphen.fasta -dbtype nucl -parse_seqids -out %s" % (Output_folder, basename, outDBname)
    process = subprocess.Popen(makeblastdb_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text = True)
    (out, err) = process.communicate()
    #print err, out
    if verbose_level in [1, 2]:
        log.write(str(err))
        log.write(str(out))
    return

def BlastSeq(inputfile, outputfile, databasefile, num_threads=1, evalue=0.0000001, max_target=1, outfmt='6 qacc sacc nident mismatch length pident bitscore'):
    """Calls blastn. The output format can be changed by outfmt. Requires the blast database to be made already"""
    blastn_cLine = "blastn -query %s -task blastn -num_threads %s -db %s -out %s -evalue %s -max_target_seqs %d -outfmt '%s'" % (inputfile, num_threads, databasefile, outputfile, evalue, max_target, outfmt)
    process = subprocess.Popen(blastn_cLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text = True)
    process.wait()
    check_cmdLine = "wc -l %s | awk '{print $1}'" % outputfile
    proc = subprocess.Popen(check_cmdLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text = True)
    proc.wait()
    (out, err) =proc.communicate()
    if int(out.strip("\n")) == 0:
        print("ERROR: BLAST output is empty.\n\tInput file: %s\n\tSearch database: %s\n\tOutput file: %s\n" %(inputfile, databasefile, outputfile))
        print("If input file is correct, check that BLAST database was created correctly.")
        sys.exit(1)
    return

def CheckChimericLoci(inputfile_raw_sequences, outputfile_blast, outputfile_goodSeqs, outputfile_chimeras, databasefile, SeqDict, SplitChimera=False):
    """Blasts each input sequence to the reference database to detect concatemers (i.e. a sequence that matches two different loci)
    Return "chimera_dict", in which the sequence name is the key and [locus_name1, locus_name2] is the value
    If SplitChimera = True, then will split the concatemer into two loci, based on the coordinates returned from BLAST.
        NOTE: to ensure the BC are split together with each locus, the orientation/strand also matters.
    """
    BlastSeq(inputfile_raw_sequences, outputfile_blast, databasefile, num_threads=num_threads, evalue=1e-100, max_target=100, outfmt='6 qacc sacc length pident evalue qstart qend qlen sstrand')

    chimera_blast = open(outputfile_blast, 'r') # Read the blast result
    chimeras = open(outputfile_chimeras, 'w') # File to save chimera sequences (concatemers)
    non_chimeras = open(outputfile_goodSeqs, 'w') # File to save non-chimeric, good sequences

    loci_info_dict = {}
    chimera_dict = {}

    # Go through the blast result, and see if any sequence matches to two different loci
    for each_rec in chimera_blast:
        each_rec = each_rec.strip('\n')
        seq_name = each_rec.split('\t')[0]
        First = True
        if seq_name not in list(chimera_dict.keys()):
            try:
                locus_name = re.search('LOCUS=(\w+)/', each_rec.split('\t')[1], re.IGNORECASE).group(1) # The names are in the format "locus=X/group=XY/ref_taxon=XYZ"
            except:
                sys.exit('ERROR in parsing locus annotations in the reference sequences; should be in the format of >locus=X/group=XY/ref_taxon=XYZ')
            locus_name = locus_name.upper()
            locus_start = each_rec.split('\t')[5]
            locus_end = each_rec.split('\t')[6]
            locus_direction = each_rec.split('\t')[-1]
            if seq_name not in list(loci_info_dict.keys()):
                loci_info_dict[seq_name] = [(locus_name, int(locus_start), int(locus_end), locus_direction)] # e.g. ('IBR', 41, 906, 'plus')

            else:
                if locus_name != loci_info_dict[seq_name][0][0]:
                    chimera_dict[seq_name] = [locus_name, loci_info_dict[seq_name][0][0]]
                    if First:
                        loci_info_dict[seq_name].append((locus_name, int(locus_start), int(locus_end), locus_direction))
                        First = False

    # Get the non-chimeric seq, and save them
    non_chimeras_list = list(set(list(SeqDict.keys())) - set(list(chimera_dict.keys())))
    for each_non_chimera in non_chimeras_list:
        non_chimeras.write('>' + each_non_chimera + '\n' + str(SeqDict[each_non_chimera].seq) + '\n')

    if SplitChimera:
        for each_chimera in chimera_dict:
            loci_info_list = sorted(loci_info_dict[each_chimera], key=lambda LIST: LIST[1]) # loci_info_list as [('IBR', 41, 906, 'plus'), ('APP', 956, 1920, 'minus')]

            new_seq_name1 = each_chimera + 'ERRchimera_' + loci_info_list[0][0] + 'of' + loci_info_list[0][0] + '+' + loci_info_list[1][0]
            new_seq_name2 = each_chimera + 'ERRchimera_' + loci_info_list[1][0] + 'of' + loci_info_list[0][0] + '+' + loci_info_list[1][0]

            # The chimeric sequences can be in any orientation, and each requires different ways to split the sequence
            if loci_info_list[0][3] == 'plus' and loci_info_list[1][3] == 'plus':
                locus_seq1 = str(SeqDict[each_chimera].seq)[:loci_info_list[0][2]]
                locus_seq2 = str(SeqDict[each_chimera].seq)[loci_info_list[0][2]:]

            elif loci_info_list[0][3] == 'plus' and loci_info_list[1][3] == 'minus':
                locus_seq1 = str(SeqDict[each_chimera].seq)[:loci_info_list[0][2]]
                locus_seq2 = str(SeqDict[each_chimera].seq)[loci_info_list[0][2]:]

            elif loci_info_list[0][3] == 'minus' and loci_info_list[1][3] == 'plus':
                midpoint = int((loci_info_list[0][2] + loci_info_list[1][1]) / 2)
                locus_seq1 = str(SeqDict[each_chimera].seq)[:midpoint]
                locus_seq2 = str(SeqDict[each_chimera].seq)[midpoint:]

            elif loci_info_list[0][3] == 'minus' and loci_info_list[1][3] == 'minus':
                locus_seq1 = str(SeqDict[each_chimera].seq)[:loci_info_list[1][1]]
                locus_seq2 = str(SeqDict[each_chimera].seq)[loci_info_list[1][1]:]

            non_chimeras.write('>' + new_seq_name1 + '\n' + locus_seq1 + '\n')
            non_chimeras.write('>' + new_seq_name2 + '\n' + locus_seq2 + '\n')

    for each_chimera in chimera_dict:
        new_seq_name = each_chimera + 'ERRchimera_' + chimera_dict[each_chimera][0] + '+' + chimera_dict[each_chimera][1]
        chimeras.write('>' + new_seq_name + '\n' + str(SeqDict[each_chimera].seq) + '\n')

    return chimera_dict # as {seq1: [locus1, locus2], seq2: [locus1, locus3]}

def lima_output_rename(Output_folder, Output_prefix, fileExt, barcode_seq_file): # renaming sequences to match PURC _1_bc_trimmed.fa format
    with open(barcode_seq_file, "r") as open_barcode_file:
        barcodeList = []
        for line in open_barcode_file:
            if line.startswith(">"):
                barcodeList.append(line.strip(">\n"))
    if fileExt in ["fasta", "fa", "fas", "fna", "faa"]:
        with open("%s/%s.demux.%s" %(Output_folder, Output_prefix, fileExt), "r") as open_seq_file:
            with open("%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix), "w") as open_outfile:
                for line in open_seq_file:
                    if line.startswith(">"):
                        splitline = line.strip(">\n").split(" ")
                        readName = splitline[0]
                        for item in splitline:
                            if item.split("=")[0] == "bc":
                                barcodes = item.split("=")[1]
                                barcode1 = int(barcodes.split(",")[0])
                                barcode1Name = barcodeList[barcode1]
                                if len(barcodes.split(",")) == 2:
                                    barcode2 = int(barcodes.split(",")[1])
                                    barcode2Name = barcodeList[barcode2]
                        if Lima_barcode_type == "single-side":
                            newReadName = ">%s|%s\n" %(barcode1Name, readName)
                        else:
                            newReadName = ">%s^%s|%s\n" %(barcode1Name, barcode2Name, readName)
                        for char in [';', '=', '/']:
                            newReadName = newReadName.replace(char, "_")
                        open_outfile.write(newReadName)
                    else:
                        open_outfile.write(line)
    elif fileExt in ["fastq", "fq"]:
        with open("%s/%s.demux.%s" %(Output_folder, Output_prefix, fileExt), "r") as open_seq_file:
            with open("%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix), "w") as open_outfile:
                lineCount = 0
                for line in open_seq_file:
                    if lineCount == 0 or (lineCount % 4) == 0:
                        splitline = line.strip("@\n").split(" ")
                        readName = splitline[0]
                        for item in splitline:
                            if item.split("=")[0] == "bc":
                                barcodes = item.split("=")[1]
                                barcode1 = int(barcodes.split(",")[0])
                                barcode1Name = barcodeList[barcode1]
                                if len(barcodes.split(",")) == 2:
                                    barcode2 = int(barcodes.split(",")[1])
                                    barcode2Name = barcodeList[barcode2]
                        if Lima_barcode_type == "single-side":
                            newReadName = ">%s|%s\n" %(barcode1Name, readName)
                        else:
                            newReadName = ">%s^%s|%s\n" %(barcode1Name, barcode2Name, readName)
                        for char in [';', '=', '/']:
                            newReadName = newReadName.replace(char, "_")
                        open_outfile.write(newReadName)
                    elif lineCount == 1 or ((lineCount-1)%4) == 0:
                        open_outfile.write(line)
                    lineCount += 1

def lima_dual(inputfile_raw_sequences, Lima_barcode_type, barcode_seq_file, Output_folder, Output_prefix):
# commands for creating biosample.csv file for lima -- not currently working in lima
#    for mapFile in mapping_file_list:
#        biosamplecsv = ".".join(mapFile.split(".")[:-1]) + ".biosample.csv"
#        with open(mapFile, "r") as open_infile:
#            with open(biosamplecsv, "w") as open_outfile:
#                open_outfile.write("Barcodes,Bio Sample\n")
#                for line in open_infile:
#                    splitline = line.strip("\n").split("\t")
#                    open_outfile.write("%s--%s,%s\n" % (splitline[0], splitline[1], splitline[2]))
    prefix = ".".join(inputfile_raw_sequences.split(".")[:-1])
    fileExt = inputfile_raw_sequences.split(".")[-1]
    cmdLine = "lima --different --ccs --min-score 80 %s %s %s/%s.demux.%s" % (inputfile_raw_sequences, barcode_seq_file, Output_folder, Output_prefix, fileExt)
    process = subprocess.Popen(cmdLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text = True)
    stdout,stderr = process.communicate()
    print(stdout)
    print(stderr)
    lima_output_rename(Output_folder, Output_prefix, fileExt, barcode_seq_file)

def lima_symmetric(inputfile_raw_sequences, Lima_barcode_type, barcode_seq_file, Output_folder, Output_prefix):
    prefix = ".".join(inputfile_raw_sequences.split(".")[:-1])
    fileExt = inputfile_raw_sequences.split(".")[-1]
    cmdLine = "lima --same --ccs --min-score 80 %s %s %s/%s.demux.%s" % (inputfile_raw_sequences, barcode_seq_file, Output_folder, Output_prefix, fileExt)
    process = subprocess.Popen(cmdLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text = True)
    stdout,stderr = process.communicate()
    print(stdout)
    print(stderr)
    lima_output_rename(Output_folder, Output_prefix, fileExt, barcode_seq_file)

def lima_singleend(inputfile_raw_sequences, Lima_barcode_type, barcode_seq_file, Output_folder, Output_prefix):
    prefix = ".".join(inputfile_raw_sequences.split(".")[:-1])
    fileExt = inputfile_raw_sequences.split(".")[-1]
    cmdLine = "lima --single-side --ccs --min-score 80 %s %s %s/%s.demux.%s" % (inputfile_raw_sequences, barcode_seq_file, Output_folder, Output_prefix, fileExt)
    process = subprocess.Popen(cmdLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text = True)
    stdout,stderr = process.communicate()
    print(stdout)
    print(stderr)
    lima_output_rename(Output_folder, Output_prefix, fileExt, barcode_seq_file)

def DeBarcoder(inputfile_raw_sequences, databasefile, SeqDict, Output_folder, Output_prefix):
    """Blasts the raw sequences against the barcode blast database, identifies the barcode, adds the barcode ID to the
    sequence name, removes the barcode from sequence; returns a list of barcode names.
    Note that the range of acceptable bc starting point is hardcoded here, e.g. "if barcode_info_dict[each_seq][1] < 15:", "elif barcode_info_dict[each_seq][1] > len(str(SeqDict[each_seq].seq))-40:"
    """

    BlastSeq(inputfile_raw_sequences, Output_folder + '/blast_barcode_out.txt', databasefile, num_threads=num_threads, evalue=1, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')

    bc_blast = open(Output_folder + '/blast_barcode_out.txt', 'r') # Read the blast result
    bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'w') # For writing the de-barcoded sequences
    bc_leftover = open(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'w') # For saving those without barcodes
    bc_toomany = open(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa', 'w') # For saving those more than one barcode

    seq_withbc_list = [] # A list containing all the seq names that have barcodes
    seq_withbc_morethanone_list = [] # A list containing all the seq names that have more than one barcode
    seq_withoutbc_list = [] # A list containing all the seq names that do not have barcode identified by BLAST

    barcode_info_dict = {} # {seq_name1: [BC01, 0, 12], seq_name2: [BC08, 0, 12]}; barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi]

    # Go through the blast output file, and complete the barcode_info_dict, seq_withbc_list, and seq_withbc_morethanone_list
    for each_rec in bc_blast:
        each_rec = each_rec.strip('\n')
        seq_name = each_rec.split('\t')[0]
        barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
        barcode_start_posi = int(each_rec.split('\t')[5])
        barcode_end_posi = int(each_rec.split('\t')[6])

        if seq_name not in list(barcode_info_dict.keys()) and seq_name not in seq_withbc_morethanone_list:
            barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi]
            seq_withbc_list.append(seq_name)
        elif seq_name in seq_withbc_morethanone_list:
            continue
        else: # means that this seq has more than one barcode, then take out this seq record from seq_withbc_list, but append it to seq_withbc_morethanone_list
            del barcode_info_dict[seq_name]
            seq_withbc_list.remove(seq_name)
            seq_withbc_morethanone_list.append(seq_name)

    # De-barcode and write sequences
    for each_seq in seq_withbc_list:
        new_seq_name = str(barcode_info_dict[each_seq][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name

        # Check the orientation of the sequence; if the barcode is in the 3' end, reverse complement the seq
        if barcode_info_dict[each_seq][1] < 15: # The start position _should_ be 1, but this allows for some slop
            new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][2]:]) # "barcode_end_posi" to the end of the sequence
        elif barcode_info_dict[each_seq][1] > len(str(SeqDict[each_seq].seq))-40: # Those barcodes that are at the end of the sequences, so need to be reversecomplemented
            new_seq_trimmed = ReverseComplement(str(SeqDict[each_seq].seq[:barcode_info_dict[each_seq][1]-1])) # the beginning of the sequence to "barcode_start_posi" - 1
        else: # Those barcodes that are at the middle of the sequences
            new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][2]:])
            new_seq_name = new_seq_name + "ERRmidBC"
        bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')

    # Write the sequences with multiple barcodes
    for each_seq in seq_withbc_morethanone_list:
        bc_toomany.write('>' + str(each_seq) + '\n' + str(SeqDict[each_seq].seq) + '\n')

    # Write the sequences without identified barcode to bc_leftover
    seq_withoutbc_list = list(set(list(SeqDict.keys())) - set(seq_withbc_list) - set(seq_withbc_morethanone_list))
    for seq_withoutbc in seq_withoutbc_list:
        bc_leftover.write('>' + str(seq_withoutbc) + '\n' + str(SeqDict[seq_withoutbc].seq) + '\n')


    bc_blast.close()
    bc_toomany.close()
    bc_leftover.close()
    bc_trimmed.close() #this is the file that now has all the sequences, labelled with the barcode, and the barcodes themselves removed

    return

def DeBarcoder_ends(SeqDict, databasefile, Output_folder, Output_prefix, search_range=25):
    """This function looks for barcodes at the ends of the sequence, only. So it avoids internal barcodes
    (helpful in part because some of our loci had regions that were similar to some primer sequences)
    """
    bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'w') # For writing the de-barcoded sequences
    bc_leftover = open(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'w') # For saving those without barcodes
    bc_toomany = open(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa', 'w') # For saving those more than one barcode

    # Get 5' and 3' end sequences, the length of seq is determined by search_range
    F_ends = open('tempF', 'w')
    R_ends = open('tempR', 'w')
    for each_rec in sorted(SeqDict):
        seq_to_search_F = str(SeqDict[each_rec].seq)[:search_range]
        F_ends.write('>' + str(each_rec) + '\n' + seq_to_search_F + '\n')

        seq_to_search_R = ReverseComplement(str(SeqDict[each_rec].seq))[:search_range]
        R_ends.write('>' + str(each_rec) + '\n' + seq_to_search_R + '\n')

    F_ends.close()
    R_ends.close()
    BlastSeq('tempF', Output_folder + '/blast_barcodeF_out.txt', databasefile, num_threads=num_threads, evalue=1, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')
    BlastSeq('tempR', Output_folder + '/blast_barcodeR_out.txt', databasefile, num_threads=num_threads, evalue=1, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')

    seq_withbc_list = [] # A list containing all the seq names that have barcodes
    seq_withbc_morethanone_list = [] # A list containing all the seq names that have more than one barcode
    seq_withoutbc_list = [] # A list containing all the seq names that do not have barcode identified by BLAST
    barcode_info_dict = {} # {seq_name1: [BC01, 0, 12], seq_name2: [BC08, 0, 12]}; barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi]

    bc_blast_F = open(Output_folder + '/blast_barcodeF_out.txt', 'r')
    for each_rec in bc_blast_F:
        each_rec = each_rec.strip('\n')
        seq_name = each_rec.split('\t')[0]
        barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
        barcode_start_posi = int(each_rec.split('\t')[5])
        barcode_end_posi = int(each_rec.split('\t')[6])
        if seq_name not in list(barcode_info_dict.keys()):
            barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi, '+']
            seq_withbc_list.append(seq_name)
        else: # means that this seq has more than one barcode, then take out this seq record from seq_withbc_list, but append it to seq_withbc_morethanone_list
            del barcode_info_dict[seq_name]
            seq_withbc_list.remove(seq_name)
            seq_withbc_morethanone_list.append(seq_name)

    bc_blast_R = open(Output_folder + '/blast_barcodeR_out.txt', 'r')
    for each_rec in bc_blast_R:
        each_rec = each_rec.strip('\n')
        seq_name = each_rec.split('\t')[0]
        barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
        barcode_start_posi = int(each_rec.split('\t')[5])
        barcode_end_posi = int(each_rec.split('\t')[6])
        if seq_name not in list(barcode_info_dict.keys()) and seq_name not in seq_withbc_morethanone_list:
            barcode_info_dict[seq_name] = [barcode_name, barcode_start_posi, barcode_end_posi, '-']
            seq_withbc_list.append(seq_name)
        elif seq_name in seq_withbc_morethanone_list:
            continue
        else: # means that this seq has more than one barcode, then take out this seq record from seq_withbc_list, but append it to seq_withbc_morethanone_list
            del barcode_info_dict[seq_name]
            seq_withbc_list.remove(seq_name)
            seq_withbc_morethanone_list.append(seq_name)

    # De-barcode and write sequences
    for each_seq in seq_withbc_list:
        new_seq_name = str(barcode_info_dict[each_seq][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name

        #check the orientation of the sequence; if the barcode is in the 3' end, reverse complement the seq
        if barcode_info_dict[each_seq][-1] == '+': # bc on the 5' end
            new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][2]:])
        elif barcode_info_dict[each_seq][-1] == '-': # bc on the 3' end
            new_seq_trimmed = ReverseComplement(str(SeqDict[each_seq].seq))[barcode_info_dict[each_seq][2]:]

        bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')

    # Write the sequences with multiple barcodes
    for each_seq in seq_withbc_morethanone_list:
        bc_toomany.write('>' + str(each_seq) + '\n' + str(SeqDict[each_seq].seq) + '\n')

    # Write the sequences without identified barcode to bc_leftover
    seq_withoutbc_list = list(set(list(SeqDict.keys())) - set(seq_withbc_list) - set(seq_withbc_morethanone_list))
    for seq_withoutbc in seq_withoutbc_list:
        bc_leftover.write('>' + str(seq_withoutbc) + '\n' + str(SeqDict[seq_withoutbc].seq) + '\n')

    os.remove('tempF')
    os.remove('tempR')
    bc_blast_F.close()
    bc_blast_R.close()
    bc_toomany.close()
    bc_leftover.close()
    bc_trimmed.close() #this is the file that now has all the sequences, labelled with the barcode, and the barcodes themselves removed

def DeBarcoder_dual(inputfile_raw_sequences, databasefile, SeqDict):
    """Blasts the raw sequences against the barcode blast database, identifies the barcode, adds the barcode ID to the
    sequence name, removes the barcode from sequence; deal with barcodes at both primers.
    Note that the range of acceptable bc starting point is hardcoded here, e.g. "if barcode_info_dict[each_seq][0][1] < 5" and "elif barcode_info_dict[each_seq][0][1] > len(str(SeqDict[each_seq].seq))-30".
    """

    BlastSeq(inputfile_raw_sequences, Output_folder + '/blast_barcode_out.txt', databasefile, num_threads=num_threads, evalue=1, max_target=2, outfmt='6 qacc sacc length pident bitscore qstart qend')

    bc_blast = open(Output_folder + '/blast_barcode_out.txt', 'r') # Read the blast result
    bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'w') # For writing the de-barcoded sequences
    bc_leftover = open(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'w') # For saving those without barcodes
    bc_onlyF = open(Output_folder + '/' + Output_prefix + '_1_trashBin_onlyF_bc.fa', 'w')
    bc_onlyR = open(Output_folder + '/' + Output_prefix + '_1_trashBin_onlyR_bc.fa', 'w')
    bc_toomany = open(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa', 'w') # For saving those more than one barcode
    bc_invalid = open(Output_folder + '/' + Output_prefix + '_1_trashBin_invalid_bc.fa', 'w') # For saving those having FF or RR barcodes

    seq_withbc_list = [] # A list containing all the seq names that have barcodes
    seq_withoutbc_list = [] # A list containing all the seq names that do not have barcode identified by BLAST
    seq_onlyF_list = []
    seq_onlyR_list = []
    barcode_info_dict = {} # {seq_name1: [(BCF01, 0, 12), (BCR02, 459, 511)], seq_name2: [(BC08, 0, 12)]}; barcode_info_dict[seq_name] = [(barcodeF_name, barcode_start_posi, barcode_end_posi), (barcodeR_name, barcode_start_posi, barcode_end_posi)]
                            # each seq has a list of tuples, each tuple contains the (barcode name, start, end)
    # Go through the blast output file, and complete the barcode_info_dict, and seq_withbc_list
    for each_rec in bc_blast:
        each_rec = each_rec.strip('\n')
        seq_name = each_rec.split('\t')[0]
        barcode_name = each_rec.split('\t')[1] # E.g. BC01, BC24...
        barcode_start_posi = int(each_rec.split('\t')[5])
        barcode_end_posi = int(each_rec.split('\t')[6])
        seq_withbc_list.append(seq_name)
        # {seq_name1: [(BCF01, 0, 12), (BCR02, 459, 511)], seq_name2: [(BC08, 0, 12)]}
        try:
            barcode_info_dict[seq_name].append((barcode_name, barcode_start_posi, barcode_end_posi))
        except:
            barcode_info_dict[seq_name] = [(barcode_name, barcode_start_posi, barcode_end_posi)]

    # De-barcode and write sequences
    for each_seq in set(seq_withbc_list):
        # When more than three barcodes are present in seq
        if len(barcode_info_dict[each_seq]) >= 3:
            bc_toomany.write('>' + str(each_seq) + '\n' + str(SeqDict[each_seq].seq) + '\n')
            continue
        # When only one barcode is present in seq
        elif len(barcode_info_dict[each_seq]) == 1:
            new_seq_name = str(barcode_info_dict[each_seq][0][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name
            # Check the orientation of the sequence; if the barcode is in the 3' end, reverse complement the seq
            if barcode_info_dict[each_seq][0][1] < 5: # The start position _should_ be 1, but this allows for some slop
                new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][0][2]:]) # "barcode_end_posi" to the end of the sequence
            elif barcode_info_dict[each_seq][0][1] > len(str(SeqDict[each_seq].seq))-30: # Those barcodes that are at the end of the sequences, so need to be reversecomplemented
                new_seq_trimmed = ReverseComplement(str(SeqDict[each_seq].seq[:barcode_info_dict[each_seq][0][1]-1])) # the beginning of the sequence to "barcode_start_posi" - 1
            else: # Those barcodes that are at the middle of the sequences
                new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][0][2]:])
                new_seq_name = new_seq_name + "ERRmidBC"
            # Check where barcode is located, on F or R primer
            if barcode_info_dict[each_seq][0][0][2] == 'F':
                bc_onlyF.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
            elif barcode_info_dict[each_seq][0][0][2] == 'R':
                bc_onlyR.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
        # When both barcodes are present in seq
        elif len(barcode_info_dict[each_seq]) == 2:
            # Check if barcodes are in the same category (F, R); don't want these
            if barcode_info_dict[each_seq][0][0][2] == barcode_info_dict[each_seq][1][0][2]:
                new_seq_name = str(barcode_info_dict[each_seq][0][0]) + '^' + str(barcode_info_dict[each_seq][1][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BC01|sequence_name
                bc_invalid.write('>' + new_seq_name + '\n' + str(SeqDict[each_seq].seq) + '\n')
                continue

            # Re-order the list, so that the F barcode tuple is at the first position in the list
            if barcode_info_dict[each_seq][0][0][2] == 'R':
                barcode_info_dict[each_seq] = barcode_info_dict[each_seq][::-1]

            new_seq_name = str(barcode_info_dict[each_seq][0][0]) + '^' + str(barcode_info_dict[each_seq][1][0]) + '|' + str(each_seq) # Add the barcode ID to the sequence name: BCF01|BCR02|sequence_name

            # Trim the barcodes
            # BCF - F primer - Seq - R primer - BCR
            if barcode_info_dict[each_seq][0][1] < 5:
                new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][0][2]:barcode_info_dict[each_seq][1][1]-1])
                if len(new_seq_trimmed) > 0: # PWS 2019 Sept 9 added if statement to avoid printing empty sequences
                    bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
            # BCR - R primer - Seq - F primer - BCF
            elif barcode_info_dict[each_seq][0][1] > len(str(SeqDict[each_seq].seq))-30: # When the seq is reverse complemented
                new_seq_trimmed = str(SeqDict[each_seq].seq[barcode_info_dict[each_seq][1][2]:barcode_info_dict[each_seq][0][1]-1])
                new_seq_trimmed = ReverseComplement(new_seq_trimmed)
                if len(new_seq_trimmed) > 0: # PWS 2019 Sept 9 added if statement to avoid printing empty sequences
                    bc_trimmed.write('>' + new_seq_name + '\n' + new_seq_trimmed + '\n')
            else:
                new_seq_name = new_seq_name + "ERRmidBC"
                # Do something here #

    # Save those without barcode
    seq_withoutbc_list = list(set(list(SeqDict.keys())) - set(seq_withbc_list))
    for seq_withoutbc in seq_withoutbc_list:
        bc_leftover.write('>' + str(seq_withoutbc) + '\n' + str(SeqDict[seq_withoutbc].seq) + '\n')

    bc_blast.close()
    bc_trimmed.close()
    bc_leftover.close()
    bc_onlyF.close()
    bc_onlyR.close()
    bc_toomany.close()
    bc_invalid.close()

def DeBarcoder_SWalign(SeqDict, barcode_seq_filename, Output_folder, Output_prefix, search_range=25):
    """Use Smith-Waterman local alignment to find barcodes. This approach is slow; recommended as a last resort"""
    #sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1))
    sw = LocalAlignment(NucleotideScoringMatrix(2, -1))

    barcode_seq = parse_fasta(barcode_seq_filename)
    out_stat = open(Output_folder + '/SWalign_barcode_out.txt', 'w')
    bc_trimmed = open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', 'a') # For writing the de-barcoded sequences

    count_matches = 0
    for each_rec in sorted(SeqDict):
        seq_to_search_F = str(SeqDict[each_rec].seq)[:search_range]
        seq_to_search_R = ReverseComplement(str(SeqDict[each_rec].seq))[:search_range]
        Match = None
        # Go through each barcode, do alignment, and see if there is a match or not
        for each_bc in barcode_seq:
            # Forward
            aln = sw.align(each_bc.seq, seq_to_search_F)
            aln_len = aln.r_end - aln.r_pos
            mismatches = len(each_bc.seq) - aln_len + aln.mismatches
            if mismatches <= 3:
                if not Match or aln.score > Match[1].score:
                    Match = (each_bc.id, aln, mismatches, '+', each_bc.id)
            # Reverse
            aln = sw.align(each_bc.seq, seq_to_search_R)
            aln_len = aln.r_end - aln.r_pos
            mismatches = len(each_bc.seq) - aln_len + aln.mismatches
            if mismatches <= 3:
                if not Match or aln.score > Match[1].score:
                    Match = (each_bc.id, aln, mismatches, '-', each_bc.id)

        if Match:
            count_matches = count_matches + 1
            out_stat.write(str(each_rec) + '\t' + Match[0] + '\t' + str(Match[1].q_end) + '\t' + str(Match[2]) + '\t' + Match[3] + '\n')
            new_seq_name = Match[4] + '|' + str(each_rec) + '_recycled'
            if Match[3] == '+':
                bc_trimmed.write('>' + new_seq_name + '\n' + str(SeqDict[each_rec].seq)[Match[1].q_end:] + '\n')
            else:
                bc_trimmed.write('>' + new_seq_name + '\n' + ReverseComplement(str(SeqDict[each_rec].seq)[Match[1].q_end:]) + '\n')
    return count_matches

def doCutAdapt(Fprims, Rprims, InFile, OutFile, minimum_len=50):
    """Removes the primers using the Cutadapt program."""

    # Build the forward and reverse primer portions of the cutadapt command lines by adding each primer in turn
    F_cutadapt_cline = ''
    for each_primer in Fprims:
        each_primer = each_primer.upper() #in case the primers aren't uppercase (lower case nucs will confuse ReverseComplement)
        F_cutadapt_cline = F_cutadapt_cline + '-g ^' + each_primer + ' '
        # e.g.,  F_cutadapt_cline as '-g ^GGACCTGGSCTYGCTGARGAGTG -g ^TCTGCMCATGCMATTGAAAGAGAG -g ^GAGYGTTTGGAATGTYTCWTTCCTYGG'
        """ The version with the "^" included anchors the fprimer so that it has to occur at the beginning of the sequence.
        This has the advantage of not picking up primers in the middle (in which case the first part--presumably including the
        original barcode--would be erased). But it may leave behind too many primers (primers have to be within the error rate
        of full-length) which may mess up the clustering. And R primers can still be found in the middle of the sequence. """

    R_cutadapt_cline = ''
    for each_primer in Rprims:
        each_primer = each_primer.upper() #in case the primers aren't uppercase (lower case nucs will confuse ReverseComplement)
        R_cutadapt_cline = R_cutadapt_cline + '-a ' + ReverseComplement(each_primer) + ' '

    # Build the complete cutadapt command line
    cutadapt_cline = '%s -O 15 -e 0.05 -n 2 -m %s %s %s %s > %s' % (Cutadapt, minimum_len, F_cutadapt_cline, R_cutadapt_cline, InFile, OutFile)
    '''Here's where you can control some of the primer-detection settings. Pay particular attention to the error rate: if that
    is too large, primers may match internal regions of the sequence reads, resulting in that region and everything downstream
    being trimmed away. If the sequences at the end of the pipeline are too short, that's probably what happened '''
    # -O: Minimum overlap length. If the overlap between the read and the primer is shorter than -O, the read is not modified.
    # -e: Maximum allowed error rate (no. of errors divided by the length of the matching region)
    # -n: Try to remove primers at most -n times.

    process = subprocess.Popen(cutadapt_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    process.wait()
    (out, err) = process.communicate()
    if verbose_level in [1, 2]:
        log.write('**Primer-trimming results**\n')
        log.write(str(err))

def SplitBy(annotd_seqs_file, split_by = "locus-taxon", Multiplex_perBC_flag=True):
    """Uses the annotated sequences to split sequences into different files based on splits_list
    (could be by barcode or by locus, etc); returns a dictionary of seq counts for each subgroup"""

    unsplit_seq = parse_fasta(annotd_seqs_file)
    splits_file_dic = {}
    splits_count_dic = {}
    splits_list = []

    for each_seq in unsplit_seq:
        # Find the annotation for the split of interest.
        # e.g., BC01, BC02, ... or gapCp, PGIC, ...

        if split_by == "taxon":
            split = str(each_seq.id).split('|')[0]
        elif split_by == "locus":
            split = str(each_seq.id).split('|')[1]
        elif split_by == "taxon-locus":
            split = str(each_seq.id).split('|')[0] + "_" + str(each_seq.id).split('|')[1]
        elif split_by == "locus-taxon": # same as above, but folders/files labeled with locus name before taxon name
            split = str(each_seq.id).split('|')[1] + "_" + str(each_seq.id).split('|')[0]
        elif split_by == "barcode" and Multiplex_perBC_flag:
            split = str(each_seq.id).split('|')[3]
        elif split_by == "group" and Multiplex_perBC_flag: # where group is the let set to identify the taxonomic groups, e.g, A, B, C, ...
            split = str(each_seq.id).split('|')[2]
        elif split_by == "barcode" and not Multiplex_perBC_flag:
            split = str(each_seq.id).split('|')[2]
        else:
            sys.exit('Error: Attempting to split sequences by an inappropriate (absent?) category.\n')

        if split not in splits_list: # add split identifier to the list of splits
            splits_list.append(split)
            os.makedirs(split, exist_ok=True)
            os.chdir(split)
            seq_file = split + '.fa'
            file_handle = open(seq_file, 'w') # Only have to do this once, for the first time that split occurs
            splits_file_dic[split] = file_handle # e.g., {BC01:BC01.fa, BC02:BC02.fa, ...}
        else:
            os.chdir(split)

        file_handle = splits_file_dic[split] #use that split as the key to find the corresponding output file
        file_handle.write('>' + str(each_seq.id) + '\n' + str(each_seq.seq) + '\n') #write the seq

        try:
            splits_count_dic[split] += 1
        except:
            splits_count_dic[split] = 1

        os.chdir('..')
    return splits_count_dic #as {'BC01': 150, 'BC02': 156} for example

def makeMapDict(mapping_file, locus, Multiplex_perBC_flag=True, DualBC_flag=False): #danger? Locus used somewhere else?
    """Constructs a dictionary, MapDict, that maps each barcode/locus combination to a specific accession
    This function is used in annotateIt"""

    map = open(mapping_file, 'r') #open the mapping file
    log.write("Working on %s..." % mapping_file)
    if Multiplex_perBC_flag and not DualBC_flag: # when there are multiple individuals sharing the same barcodes
        MapDict = {}
        for each_rec in map:
            each_rec = each_rec.strip('\n')
            map_barcode_name = each_rec.split('\t')[0]
            map_group_name = each_rec.split('\t')[1].upper()
            map_taxon_name = each_rec.split('\t')[2]
            key = map_barcode_name + '_' + map_group_name #BC01_A, BC01_B, BC010_C...
            MapDict[key] = map_taxon_name
    elif not Multiplex_perBC_flag and DualBC_flag:
        MapDict = {}
        for each_rec in map:
            each_rec = each_rec.strip('\n')
            map_barcodeF_name = each_rec.split('\t')[0]
            map_barcodeR_name = each_rec.split('\t')[1]
            map_taxon_name = each_rec.split('\t')[2]
            key = map_barcodeF_name + '^' + map_barcodeR_name #BCF01^BCR01, BCF02^BCR02, ...
            MapDict[key] = map_taxon_name
    elif not Multiplex_perBC_flag and not DualBC_flag: # when a barcode points to one individual; no multiplex per barcode; ignore "group" info and reads a different mapping file format
        MapDict = {}
        for each_rec in map:
            each_rec = each_rec.strip('\n')
            map_barcode_name = each_rec.split('\t')[0]
            map_taxon_name = each_rec.split('\t')[1]
            key = map_barcode_name #BC01_A, BC01_B, BC010_C...
            MapDict[key] = map_taxon_name

    map.close()
    return MapDict

def annotateIt(filetoannotate, outFile, failsFile, Multiplex_perBC_flag=True, DualBC_flag=False, verbose_level=0):
    """Uses the blast results (against the reference sequence database) to assign locus and taxon, and write sequences
    for a particular locus as specified by map_locus; returns a dictionary containing taxon-locus seq counts"""

    # Blasts each sequence in the input file (e.g., BC01.fa) against the reference sequences
    BlastSeq(filetoannotate, Output_folder + '/blast_refseq_out.txt', BLAST_DBs_folder + '/' + refseq_databasefile, num_threads=num_threads, evalue=0.0000001, max_target=1, outfmt='6 qacc sacc length pident evalue qstart qend qlen')

    # Reads the  sequences as a dict
    SeqDict = SeqIO.index(filetoannotate, 'fasta')

    # Using blast matches to the reference sequences, and barcode <-> taxon mapping files, to assign
    # each seq to a particular locus and taxon
    dictOfMapDicts = {} # A dictionary to store all of the map dictionaries
    for each_file, each_locus in zip(mapping_file_list, locus_list):
        dictOfMapDicts[each_locus.upper()] = makeMapDict(each_file, each_locus, Multiplex_perBC_flag, DualBC_flag) # Note the Multiplex_perBC and DualBC flags (as True/False)

    refseq_blast = open(Output_folder + '/blast_refseq_out.txt', 'r')
    annotated_seqs = open(outFile, "w")
    no_matches = open(failsFile, "w")
    groupsList = []
    locusList = []
    LocusTaxonCountDict = {}
    seq_processed_list = []
    for each_rec in refseq_blast:
        each_rec = each_rec.strip('\n')
        seq_name = each_rec.split('\t')[0] # The un-annotated sequence name, e.g., "BC02|m131213_174801_42153_c100618932550000001823119607181400_s1_p0/282/ccs;ee=7.2;"
        refseq_name = each_rec.split('\t')[1].replace(' ','').upper() # The best-hit reference sequence name, e.g., "locus=PGI/group=C/ref_taxon=C_diapA_BC17" ##Need to change this format

        # Get the key for retrieving taxon_name in dictOfMapDicts[locus_name]
        if Multiplex_perBC_flag:
            try:
                group_name = re.search('GROUP=(\w)/', refseq_name, re.IGNORECASE).group(1)
            except:
                sys.exit('ERROR in parsing group annotations in the reference sequences; should be in the format of >locus=X/group=XY/ref_taxon=XYZ')
            try:
                locus_name = re.search('LOCUS=(\w+)/', refseq_name, re.IGNORECASE).group(1) # The names are in the format "locus=X/group=XY/ref_taxon=XYZ"
            except:
                sys.exit('ERROR in parsing locus annotations in the reference sequences; should be in the format of >locus=X/group=XY/ref_taxon=XYZ')
            key = seq_name.split('|')[0] + '_' + group_name # Grabbing the barcode from the source seq, and the group from the matching ref seq.
            #i.e., gets the unique identifier that can link to a specific sample; i.e. BC01_A, BC01_B, BC01_C...
            if not group_name in groupsList: #keeping track of which groups are found, as a way of potentially diagnosing errors
                groupsList.append(group_name)
            if not locus_name in locusList: #keeping track of which loci are found, as a way of potentially diagnosing errors
                locusList.append(locus_name)
        else:
            try:
                locus_name = re.search('LOCUS=(\w+)/', refseq_name, re.IGNORECASE).group(1)
                if not locus_name in locusList:
                    locusList.append(locus_name)
                #i.e., gets the unique barcode that can link to a specific sample; i.e. BC01, BC02, BC03...
            except:
                sys.exit('ERROR in parsing locus annotations in the reference sequences; should be in the format of >locus=X/group=XY/ref_taxon=XYZ')
            key = seq_name.split('|')[0]
        try: #use try/except to avoid the error when the key is not present in MapDict
            taxon_name = dictOfMapDicts[locus_name][key]
            #getting to the dict corresponding to this locus, and then finding that taxon that matches the barcode+group (the key)

            if Multiplex_perBC_flag:
                new_seq_name = taxon_name + '|' + locus_name + '|' + group_name + '|' + seq_name.replace(seq_name_toErase, '')
            else:
                new_seq_name = taxon_name + '|' + locus_name + '|' + seq_name.replace(seq_name_toErase, '')

            if seq_name not in seq_processed_list:
                annotated_seqs.write('>' + new_seq_name + '\n' + str(SeqDict[seq_name].seq) + '\n')
                try:
                    LocusTaxonCountDict[taxon_name, locus_name] += 1 #as {('C_mem_6732', 'PGI'): 2, ('C_mem_6732', 'IBR'): 4} for example
                except:
                    LocusTaxonCountDict[taxon_name, locus_name] = 1 #initiate the key and give count = 1
                seq_processed_list.append(seq_name)
        except:
            log.write("\tThe combo '" + str(key) + "' wasn't found in " + str(locus_name) + '\n')
            if Multiplex_perBC_flag:
                new_seq_name = locus_name + '|' + group_name + '|' + seq_name.replace(seq_name_toErase, '')
            else:
                new_seq_name = locus_name + '|' + seq_name.replace(seq_name_toErase, '')
            if seq_name not in seq_processed_list:
                no_matches.write('>' + new_seq_name + '\n' + str(SeqDict[seq_name].seq) + '\n')
                seq_processed_list.append(seq_name)
            continue

    seq_no_hit = list(set(SeqDict.keys()) - set(seq_processed_list))
    log.write("\tThere are " + str(len(seq_no_hit)) + " sequences that failed to match any of the reference sequences -- these are likely contaminants and added to the 'unclassifiable' output fasta file\n")
    for each_rec in seq_no_hit:
        no_matches.write('>' + each_rec + '\n' + str(SeqDict[each_rec].seq) + '\n')

    refseq_blast.close()
    annotated_seqs.close()
    no_matches.close()

    if verbose_level in [1, 2]:
        log.write("The groups found are " + ', '.join(groupsList) + "\nAnd the loci found are " + ', '.join(locusList) + "\n")
    return LocusTaxonCountDict #as {('C_mem_6732', 'PGI'): 2, ('C_mem_6732', 'IBR'): 4} for example

def sortIt_length(file, verbose_level=0):
    log.write("sortIt_length")
    """Sorts clusters by seq length"""
    outFile = re.sub(r"(.*)\..*", r"\1_Sl.fa", file) # Make the outfile name by cutting off the extension of the infile name, and adding "_S1.fa"
    vsearch_cline = "%s --sortbylength %s --output %s --threads %s" %(Vsearch, file, outFile, num_threads)
    process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=log, shell=True, text=True)
    process.wait()
    (out, err) = process.communicate() #the stdout and stderr
    savestdout = sys.stdout
    if verbose_level == 2:
        log.write('\n**vsearch-sorting output on' + str(file) + '**\n')
        log.write(str(err))
    return outFile # having it spit out the outfile name, if necessary, so that it can be used to call downstream stuff and avoid complicated glob.globbing

def sortIt_size(file, thresh, round, verbose_level=0):
    log.write("sortIt_size\n")
    """Sorts clusters by size, and removes those that are smaller than a particular size
    (sent as thresh -- ie, sizeThreshold).
    "round" is used to annotate the outfile name with S1, S2, etc. depending on which sort this is"""
    outFile = re.sub(r"(.*)\.fa", r"\1Ss%s.fa" %(round), file)
    logFile = outFile + ".sortIt_size.log"
    vsearch_cline = "%s --sortbysize %s --output %s --minsize %d --log %s --threads %s" %(Vsearch, file, outFile, thresh, logFile, num_threads)
    process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=log, shell=True, text=True)
    process.wait()
    (out, err) = process.communicate() #the stdout and stderr
    savestdout = sys.stdout
    if verbose_level == 2:
        log.write('\n**vsearch-sorting output on' + str(file) + '**\n')
        log.write(str(err))
    return outFile

def align_and_consensus(inputfile, output_prefix):
    """Align a fasta and output the consensus sequence; Note that consensus threshold can be modified"""
    output_alignment = output_prefix.split(';')[0] + '_aligned.fa'
    output_consensus = output_prefix.split(';')[0] + '_consensus.fa'
    #muscle_cline = '%s -in %s -out %s' % (Muscle, inputfile, output_alignment)
    mafft_cline = '%s --auto --thread %s %s > %s' %(Mafft, num_threads, inputfile, output_alignment)
    #process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=log, shell=True, text=True)
    process = subprocess.Popen(mafft_cline, stdout=subprocess.PIPE, stderr=log, shell=True, text=True)
    process.wait()
    (out, err) = process.communicate() #the stdout and stderr
    savestdout = sys.stdout

    alignment = AlignIO.read(output_alignment, 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.gap_consensus(ambiguous='N',threshold=0.51)
    output = open(output_consensus, 'w')
    output.write('>' + output_prefix + '\n' + str(consensus).replace('-','') + '\n')
    return

def clusterIt(file, clustID, round, previousClusterToCentroid_dict, verbose_level=0):
    log.write("clusterIt\n")
    """The clustering step, using the clustID value"""
    outFile = re.sub(r"(.*)\.fa", r"\1C%s_%s.fa" %(round, clustID), file) # The rs indicate "raw" and thus python's escaping gets turned off
    outClustFile = re.sub(r"(.*).fa", r"\1clusts%s.uc" %(round), file)
    logFile = outFile + ".clusterIt.log"
    if round == 1:
        vsearch_cline = "%s --cluster_fast %s --id %f --gapopen 3I/1E --consout %s --uc %s --sizeout --threads %s --log %s" % (Vsearch, file, clustID, outFile, outClustFile, num_threads, logFile)
    elif round > 1:
        vsearch_cline = "%s --cluster_fast %s --id %f --gapopen 3I/1E --consout %s --uc %s --sizein --sizeout --threads %s --log %s" % (Vsearch, file, clustID, outFile, outClustFile, num_threads, logFile)
    process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=log, shell=True, text=True)
    (out, err) = process.communicate() #the stdout and stderr
    savestdout = sys.stdout
    if verbose_level == 2:
        log.write('\n**vsearch-clustering output on' + str(file) + '**\n')
        log.write(str(err))

    # remove "centroid=" instances from output files
    with fileinput.input(files = outFile, inplace = True) as fa:
        for line in fa:
            line = line.strip("\n")
            if line.startswith(">"):
                line = line.replace("centroid=", "")
                print(line)
            else:
                print(line)

    uc = open(outClustFile, 'r')
    ClusterToCentroid_dict = {} # {Cluster0:}
    for line in uc:
        if line.startswith('C'):
            cluster_name = 'Cluster' + str(line.split('\t')[1])
            centroid_seq_name = line.split('\t')[-2]

            if round > 1:
                #ClusterToCentroid_dict[cluster_name] = previousClusterToCentroid_dict[centroid_seq_name.split(';')[0]]
                ClusterToCentroid_dict[cluster_name] = centroid_seq_name
            elif round == 1:
                ClusterToCentroid_dict[cluster_name] = centroid_seq_name

    if round == 4:
        try:
            clustered_seqs = parse_fasta(outFile)
            renamed_clustered_seqs = open('temp', 'a')
            for seq in clustered_seqs: # seq as 'Cluster15;size=1;'
                new_seq_id = str(seq.id).replace(str(seq.id).split(';')[0], ClusterToCentroid_dict[str(seq.id).split(';')[0]])
                renamed_clustered_seqs.write('>' + new_seq_id + '\n' + str(seq.seq) + '\n')
            os.remove(outFile)
            os.rename('temp', outFile)
        except:
            log.write(outFile + ' is empty; perhaps sizeThreshold too high\n')

    return outFile, ClusterToCentroid_dict, outClustFile

def deChimeIt(file, round, abskew=1.9, verbose_level=0):
    log.write("deChimeIt\n")
    """Chimera remover"""
    outFile = re.sub(r"(.*)\.fa", r"\1dCh%s.fa" %(round), file) # The rs indicate "raw" and thus python's escaping gets turned off
    outFile_uchime = re.sub(r"(.*)\.fa", r"\1dCh%s.uchime" %(round), file) # The rs indicate "raw" and thus python's escaping gets turned off
    logFile = outFile + ".deChimeIt.log"
    vsearch_cline = "%s --uchime_denovo %s --abskew %s --nonchimeras %s --uchimeout %s --log %s" % (Vsearch, file, abskew, outFile, outFile_uchime, logFile)

    ## To make the chimera-killing more stringent, change the vsearch command line (above) to something like:
    # vsearch_cline = "%s -uchime_denovo %s -abskew 1.1 -minh 0.2 -xn 3 -dn 0.5 -nonchimeras %s -uchimeout %s" % (Vsearch, file, outFile, outFile_uchime)

    process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=log, shell=True, text=True)
    (out, err) = process.communicate() #the stdout and stderr
    savestdout = sys.stdout
    if verbose_level in [1, 2]:
        log.write('\n**Uchime output on ' + str(file) + ' written to ' + str(outFile) + ' and ' + str(outFile_uchime) + '**\n')
        log.write(str(err))

    # count number of chimera seq found by parsing 'outFile_uchime'
    chimera_count = 0
    uchime_out = open(outFile_uchime, 'r')
    for line in uchime_out:
        line = line.strip('\n')
        if line.split('\t')[-1] == 'Y':
            chimera_count = chimera_count + 1
    return outFile, chimera_count

#def muscleIt(file, verbose_level=0):
#    """Aligns the sequences using MUSCLE"""
#    outFile = ".".join(file.split(".")[:-1]) + ".aligned.fa"
#    muscle_cline = '%s -in %s -out %s' % (Muscle, file, outFile)
#    process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
#    (out, err) = process.communicate() #the stdout and stderr
#    savestdout = sys.stdout
#    if verbose_level == 2:
#        log.write('\n**Muscle output on ' + str(file) + '**\n')
#        log.write(str(err))
#    return outFile

def mafftIt(file, verbose_level=0):
    """Aligns the sequences using MAFFT"""
    outFile = ".".join(file.split(".")[:-1]) + ".aligned.fa"
    mafft_cline = '%s --auto --thread %s %s > %s' % (Mafft, num_threads, file, outFile)
    process = subprocess.Popen(mafft_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    (out, err) = process.communicate() #the stdout and stderr
    savestdout = sys.stdout
    if verbose_level == 2:
        log.write('\n**MAFFT output on ' + str(file) + '**\n')
        log.write(str(err))
    return outFile

def IterativeClusterDechimera(annotd_seqs_file, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2, verbose_level = 1):
    log.write("IterativeClusterDechimera\n")
    ## Split sequences into separate files/folders for each locus ##
    sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
    locusCounts = SplitBy(annotd_seqs_file = annoFileName, split_by = "locus", Multiplex_perBC_flag = Multiplex_per_barcode)

    sys.stderr.write('Clustering/dechimera-izing seqs...\n')
    log.write('#Sequence clustering/dechimera-izing#\n')
    all_folders_loci = list(locusCounts.keys()) # SplitBy makes a dictionary where the keys are the subcategories (and thus also the
        # folders) and they correspond to the counts for each.
    LocusTaxonCountDict_clustd = {} # {('C_dia_5316', 'ApP'): 28} for example, to store clusted seq count
    LocusTaxonCountDict_chimera = {} # {('C_dia_5316', 'ApP'): [1,0,0,0,0]} for example, to store the chimerc seq count for each chimera-killing step

    ## Go through each locus ##
    for locus_folder in all_folders_loci: # locus_folder = locus name
        os.chdir(locus_folder)
        sys.stderr.write('\nWorking on: ' + locus_folder + '...\n')
        if verbose_level in [1,2]:
            log.write('\nWorking on ' + str(locus_folder) + ' ...\n')

        if not os.stat(locus_folder + ".fa").st_size == 0: # ie, the file is not empty
            ## Split sequences into separate taxon folders ##
            taxonCounts = SplitBy(annotd_seqs_file = locus_folder + ".fa", split_by = "taxon", Multiplex_perBC_flag = Multiplex_per_barcode)
            all_folders_taxon = list(taxonCounts.keys())

            for taxon_folder in all_folders_taxon:
                if verbose_level in [1,2]:
                    log.write("Working on " + taxon_folder + '\n')
                os.chdir(taxon_folder)

                if verbose_level in [1,2]:
                    log.write("\tFirst clustering\n")
                    log.write("\nAttempting to sort: " + taxon_folder + ".fa\n")

                sorted_length = sortIt_length(file = taxon_folder + ".fa", verbose_level = verbose_level)
                clustered1, previousClusterToCentroid_dict, outClustFile1 = clusterIt(file = sorted_length, previousClusterToCentroid_dict = '', clustID = clustID, round = 1, verbose_level = verbose_level)

                if verbose_level in [1,2]:
                    log.write("\tFirst chimera slaying expedition\n")
                deChimered1, chimera_count1 = deChimeIt(file = clustered1, round = 1, abskew = abskew, verbose_level = verbose_level)

                if verbose_level in [1,2]:
                    log.write("\tSecond clustering\n")
                sorted_size1 = sortIt_size(file = deChimered1, thresh = sizeThreshold, round = 1, verbose_level = verbose_level)
                clustered2, previousClusterToCentroid_dict, outClustFile2 = clusterIt(file = sorted_size1, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID2, round = 2, verbose_level = verbose_level)

                if verbose_level in [1,2]:
                    log.write("\tSecond chimera slaying expedition\n")
                deChimered2, chimera_count2 = deChimeIt(file = clustered2, round = 2, abskew = abskew, verbose_level = verbose_level)

                if verbose_level in [1,2]:
                    log.write("\tThird clustering\n")
                sorted_size2 = sortIt_size(file = deChimered2, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
                clustered3, previousClusterToCentroid_dict, outClustFile3 = clusterIt(file = sorted_size2, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID3, round = 3, verbose_level = verbose_level)

                if verbose_level in [1,2]:
                    log.write("\tThird chimera slaying expedition\n")
                deChimered3, chimera_count3 = deChimeIt(file = clustered3, round = 3, abskew = abskew, verbose_level = verbose_level)

                if verbose_level in [1,2]:
                    log.write("\Fourth clustering\n")
                sorted_size3 = sortIt_size(file = deChimered3, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
                clustered4, previousClusterToCentroid_dict, outClustFile4 = clusterIt(file = sorted_size3, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID3, round = 4, verbose_level = verbose_level)

                if verbose_level in [1,2]:
                    log.write("\tThird chimera slaying expedition\n")
                deChimered4, chimera_count4 = deChimeIt(file = clustered4, round = 4, abskew = abskew, verbose_level = verbose_level)

                sorted_size4 = sortIt_size(file = deChimered4, thresh = sizeThreshold2, round = 4, verbose_level = verbose_level)

                ## Collect all sequences from each cluster and re-consensus ##
                ClusterToSeq_dict4 = {}
                with open(outClustFile4, 'r') as file_to_read:
                    for line in file_to_read:
                        if line.startswith("H") or line.startswith("C"):
                            key = 'Cluster' + line.split('\t')[1]
                            seq = line.split("\t")[8].split(";")[0]
                            try:
                                ClusterToSeq_dict4[key].append(seq)
                            except:
                                ClusterToSeq_dict4[key] = [seq]

                with open(outClustFile3, 'r') as file_to_read:
                    for line in file_to_read:
                        if line.startswith("H"):
                            seq = line.split("\t")[8].split(";")[0]
                            centroid = line.strip("\n").split("\t")[9].split(";")[0]
                            for key in ClusterToSeq_dict4.keys():
                                if centroid in ClusterToSeq_dict4[key]:
                                    ClusterToSeq_dict4[key].append(seq)

                with open(outClustFile2, 'r') as file_to_read:
                    for line in file_to_read:
                        if line.startswith("H"):
                            seq = line.split("\t")[8].split(";")[0]
                            centroid = line.strip("\n").split("\t")[9].split(";")[0]
                            for key in ClusterToSeq_dict4.keys():
                                if centroid in ClusterToSeq_dict4[key]:
                                    ClusterToSeq_dict4[key].append(seq)

                with open(outClustFile1, 'r') as file_to_read:
                    for line in file_to_read:
                        if line.startswith("H"):
                            seq = line.split("\t")[8].split(";")[0]
                            centroid = line.strip("\n").split("\t")[9].split(";")[0]
                            for key in ClusterToSeq_dict4.keys():
                                if centroid in ClusterToSeq_dict4[key]:
                                    ClusterToSeq_dict4[key].append(seq)

                # Go through the first clustering uc file
                #ClusterToSeq_dict1 = {}
                #file_to_read = open(outClustFile1, 'r')
                #for line in file_to_read:
                #    line = line.strip('\n')
                #    if line.startswith('H') or line.startswith('C'):
                #        key = 'Cluster' + line.split('\t')[1]
                #        seq = line.split('\t')[8].split(';')[0]
                #        try:
                #            ClusterToSeq_dict1[key].append(seq)
                #        except:
                #            ClusterToSeq_dict1[key] = [seq]
                #file_to_read.close()

                # Go through the second clustering uc file
                #ClusterToSeq_dict2 = {}
                #file_to_read = open(outClustFile2, 'r')
                #for line in file_to_read:
                #    line = line.strip('\n')
                #    if line.startswith('H') or line.startswith('C'):
                #        key = 'Cluster' + line.split('\t')[1]
                #        #seqs = ClusterToSeq_dict1[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
                #        seq = line.split('\t')[8].split(';')[0]
                #        #for seq in seqs:
                #        try:
                #            ClusterToSeq_dict2[key].append(seq)
                #        except:
                #            ClusterToSeq_dict2[key] = [seq]
                #file_to_read.close()

                # Go through the third clustering uc file
                #ClusterToSeq_dict3 = {}
                #file_to_read = open(outClustFile3, 'r')
                #for line in file_to_read:
                #    line = line.strip('\n')
                #    if line.startswith('H') or line.startswith('C'):
                #        key = 'Cluster' + line.split('\t')[1]
                #        #seqs = ClusterToSeq_dict2[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
                #        seq = line.split('\t')[8].split(';')[0]
                #        #for seq in seqs:
                #        try:
                #            ClusterToSeq_dict3[key].append(seq)
                #        except:
                #            ClusterToSeq_dict3[key] = [seq]
                #file_to_read.close()

                # Go through the forth clustering uc file
                #ClusterToSeq_dict4 = {}
                #file_to_read = open(outClustFile4, 'r')
                #for line in file_to_read:
                #    line = line.strip('\n')
                #    if line.startswith('H') or line.startswith('C'):
                #        key = 'Cluster' + line.split('\t')[1]
                #        #seqs = ClusterToSeq_dict3[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
                #        seq = line.split('\t')[8].split(';')[0]
                #        #for seq in seqs:
                #        try:
                #            ClusterToSeq_dict4[key].append(seq)
                #        except:
                #            ClusterToSeq_dict4[key] = [seq]
                #file_to_read.close()

                ## Align and re-consensus of all constituent seq for each cluster ##
                SeqDict = SeqIO.index(taxon_folder + ".fa", 'fasta')
                for each_cluster in ClusterToSeq_dict4:
                    if len(ClusterToSeq_dict4[each_cluster]) >= int(sizeThreshold2):
                        cluster_seq_file = open(each_cluster, 'w')
                        for seq in ClusterToSeq_dict4[each_cluster]:
                            cluster_seq_file.write('>' + seq + '\n' + str(SeqDict[seq].seq) + '\n')
                        cluster_seq_file.close()
                        new_seq_name = taxon_folder + '_' + each_cluster + ';size=' + str(len(ClusterToSeq_dict4[each_cluster])) + ';'
                        align_and_consensus(each_cluster, new_seq_name)

                ## Put all consensus seq into one file *_Cluster_Finalconsensus.fa ##
                all_consensus_seq = open(taxon_folder + '_Cluster_Finalconsensus.fa', 'w')
                files_to_add_reconsensus = glob.glob("*_consensus.fa")
                for file in files_to_add_reconsensus:
                    shutil.copyfileobj(open(file,'r'), all_consensus_seq) #Add each file to the final output
                all_consensus_seq.close()

                ## Do final clustering and chimera-killing ##
                vsearch_cline = "%s --sortbysize %s --output %s" %(Vsearch, taxon_folder + '_Cluster_Finalconsensus.fa', taxon_folder + '_Cluster_FinalconsensusSs.fa')
                process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
                process.wait()
                (out, err) = process.communicate() #the stdout and stderr

                vsearch_cline = "%s --cluster_fast %s --id %f --gapopen 3I/1E --consout %s --uc %s --sizein --sizeout" % (Vsearch, taxon_folder + '_Cluster_FinalconsensusSs.fa', clustID4, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.uc')
                process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
                process.wait()
                (out, err) = process.communicate() #the stdout and stderr

                vsearch_cline = "%s --uchime_denovo %s --abskew %s --nonchimeras %s --uchimeout %s" % (Vsearch, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.fa', abskew, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.uchime')
                process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
                process.wait()
                (out, err) = process.communicate() #the stdout and stderr

                ## Count number of chimera seq found by parsing 'outFile_uchime' ##
                chimera_count5 = 0
                uchime_out = open(taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.uchime', 'r')
                for line in uchime_out:
                    line = line.strip('\n')
                    if line.split('\t')[-1] == 'Y':
                        chimera_count5 = chimera_count5 + 1
                uchime_out.close()
                LocusTaxonCountDict_chimera[taxon_folder, locus_folder] = [chimera_count1, chimera_count2, chimera_count3, chimera_count4, chimera_count5] # {('C_dia_5316', 'ApP'): [1,0,0,0,0]} for example

                ## Count clustered seq and store in LocusTaxonCountDict_clustd as {('C_dia_5316', 'ApP'): 28} for example ##
                try:
                    clustered_seq_file = parse_fasta(taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa')
                    for each_seq in clustered_seq_file:
                        try:
                            LocusTaxonCountDict_clustd[taxon_folder, locus_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
                        except:
                            LocusTaxonCountDict_clustd[taxon_folder, locus_folder] = 1
                except:
                    if verbose_level in [1,2]:
                        log.write(str(all_consensus_seq) + 'is an empty file\n')

                ## Rename sequences in the final fasta: add taxon name ##
                with open(taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', "r") as renamingFile:
                    with open(taxon_folder + '_OTUs.fa', "w") as renamedFile:
                        for line in renamingFile:
                            if line.startswith(">"):
                                line = re.sub(r"centroid=","",line)
                                line = re.sub(r";seqs=\d+;", ";", line)
                            renamedFile.write(line)
                #sed_cmd = "sed 's/>/>%s_/g' %s > %s" % (taxon_folder, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh_renamed.fa')
                #process = subprocess.Popen(sed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                #process.wait()
                #(out, err) = process.communicate() #the stdout and stderr

                ## Remove intermediate files if requsted ##
                if remove_intermediates == 1:
                    filt_to_remove = list(set(glob.glob('*')) - set(glob.glob('*_OTUs.fa')) - set(glob.glob(taxon_folder+'.fa')) - set(glob.glob(taxon_folder+'.fastq')) - set(glob.glob(taxon_folder+'.R')) - set(glob.glob('*ASVs.fa')) - set(glob.glob('*.pdf')) - set(glob.glob('*.log')))
                    for file in filt_to_remove:
                        os.remove(file)

                os.chdir("..") # To get out of the current barcode folder and ready for the next one
        os.chdir("..") # To get out of the current locus folder and ready for the next one
    log.write('\t...done\n\n')

    ## Put all the sequences together ##
    sys.stderr.write('\n\nPutting all the sequences together...\n\n')
    for locus_folder in all_folders_loci: # Looping through each of the locus folders
        outputfile_name_reconsensus = Output_prefix + '_4_' + str(locus_folder) + '_OTUs.fa'
        outputfile_reconsensus = open(outputfile_name_reconsensus, 'w')

        os.chdir(locus_folder)
        taxonForThisLocus = glob.glob("*")
        for taxon_folder in taxonForThisLocus: # have to go into each barcode folder in each locus folder
            if os.path.isdir(taxon_folder): # the glob might have found some files as well as folders
                os.chdir(taxon_folder)
                files_to_add_reconsensus = glob.glob("*_OTUs.fa")
                for file in files_to_add_reconsensus:
                    shutil.copyfileobj(open(file,'r'), outputfile_reconsensus) #Add each file to the final output

                os.chdir('..')
        os.chdir('..')
        outputfile_reconsensus.close()

    return LocusTaxonCountDict_clustd, LocusTaxonCountDict_chimera

def writeASV(sample, Fprimer, Rprimer, minLen, maxLen, maxEE): # Writes a custom R script for each sample to run DADA2
    with open("%s_DADA2.R" % sample, "w") as rscript:
        rscript.write('''print(format(Sys.time()))

# load libraries
library(dada2);packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(gridExtra); packageVersion("gridExtra")

path <- "."

fns <- list.files(path, pattern="fastq", full.names=TRUE)
Fprimer <- "%s"
Rprimer <- "%s"
rc <- dada2:::rc
print("Forward Primer:")
print(Fprimer)
print("Reverse Primer:")
print(Rprimer)

nops <- file.path(path, "noprimers", basename(fns))
prim <- removePrimers(fns, nops, primer.fwd=Fprimer, primer.rev=dada2:::rc(Rprimer), orient=TRUE)
if (prim[2] == 0){
print("WARNING: No primers found. Continuing with raw sequences.")
nops <- fns
}
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)

minLen <- %s
maxLen <- %s
maxEE <- %s

# Calculate outlier boundaries if min/max length not specified in config file
if (minLen == 0){
iqr <- IQR(lens)
q25 <- quantile(lens, 0.25)
minLen <- round(q25-(1.5*iqr),0)
}
if (maxLen == 0){
iqr <- IQR(lens)
q75 <- quantile(lens, 0.75)
maxLen <- round(q75+(1.5*iqr),0)
}

print("Minimum length:")
print(minLen)
print("Maximum length:")
print(maxLen)
print("Maximum expected errors:")
print(maxEE)

lowerBound <- min(c(minLen, min(lens)))
upperBound <- max(c(maxLen, max(lens)))

# round down to nearest 10 base number
while ((lowerBound %s 10) != 0){
    lowerBound <- lowerBound-1
}
# round up to nearest 10 base number
while ((upperBound %s 10) != 0){
    upperBound <- upperBound+1
}

pdf("%s_read_lengths.pdf")
hist(lens, 100, xlim = c(lowerBound, upperBound), xlab = "Length (bp)", main = "%s", xaxt="n")
axis(1, at=seq(lowerBound , upperBound, by=10))
abline(v= c(minLen, maxLen), lty=c(2,2))
dev.off()

filts <- file.path(path, "noprimers", "filtered", paste(basename(fns), ".gz", sep=""))
track <- filterAndTrim(nops, filts, minQ=3, minLen=minLen, maxLen=maxLen, maxN=0, rm.phix=FALSE, maxEE=maxEE)
print("Reads filtered:")
print(track)

drp <- derepFastq(filts, verbose=TRUE)
err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
pdf("%s_error_profile.pdf")
plotErrors(err)
dev.off()
''' % (Fprimer, Rprimer, minLen, maxLen, maxEE, "%%", "%%", sample, sample, sample))
        if useOTUpriors == "TRUE":
            otuPriors = [i.seq for i in SeqIO.parse("%s_OTUs.fa" % sample, 'fasta')]
            priors = ', '.join(['"{}"'.format(value) for value in otuPriors])
            rscript.write("dd2 <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE, priors = c(%s))" % priors)
        else:
            rscript.write("dd2 <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)")

        rscript.write('''
st <- makeSequenceTable(dd2)

# Identify chimeras
st.chim <- isBimeraDenovo(st, multithread=TRUE, verbose=TRUE)
st.chim.seq <- which(st.chim == TRUE)
st.nochim <- removeBimeraDenovo(st, method="consensus", multithread=TRUE, verbose=TRUE)
print("Chimeric Reads:")
print(sum(st.chim, na.rm = TRUE))
print("Non-chimeric reads:")
print(dim(st.nochim))
print("Percent non-chimeric reads:")
print(sum(st.nochim)/sum(st)*100)

print("SUMMARY")
print(cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sum(dd2$denoised), ASVs=dim(st.nochim)[2]))

# Output final ASV seqs
if (length(st.nochim) > 0){
    outputNames <- vector()
    for (i in 1:dim(st.nochim)[2]){
        newName <- paste("%s_ASV", i, sep = "")
        size <- paste("size=", st.nochim[i], sep = "")
        newName <- paste(newName, size, sep = "_")
        outputNames <- c(outputNames, newName)
    }
    uniquesToFasta(st.nochim, "%s_ASVs.fa", ids=outputNames)
} else { print("WARNING: No ASVs found!")}

# Output chimeras if any
if (length(st.chim.seq) > 0){
    chimeraNames <- vector()
    for (i in 1:dim(st.chim.seq)[2]){
        newName <- paste("%s_chimera", i, sep = "")
        chimeraNames <- c(chimeraNames, newName)
    }
    uniquesToFasta(st.chim.seq, "%s_chimeras.fa", ids=chimeraNames)
} else { print("No chimeras detected")}

print(format(Sys.time()))
quit()
''' % (sample, sample, sample, sample))

def dada(annotd_seqs_file, raw_fastq_sequences, Forward_primer, Reverse_primer, locus_list, minLen, maxLen, maxEE, RscriptPath, verbose_level = 1):
    log.write("DADA2\n")
    ## Split sequences into separate files/folders for each locus ##
    sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
    locusCounts = SplitBy(annotd_seqs_file = annoFileName, split_by = "locus", Multiplex_perBC_flag = Multiplex_per_barcode)

    sys.stderr.write('Sorting sequences into ASVs...\n')
    log.write('#Sorting sequences into ASVs#\n')
    all_folders_loci = list(locusCounts.keys()) # SplitBy makes a dictionary where the keys are the subcategories (and thus also the
        # folders) and they correspond to the counts for each.
    LocusTaxonCountDict_clustd = {} # {('C_dia_5316', 'ApP'): 28} for example, to store clusted seq count
    LocusTaxonCountDict_chimera = {} # {('C_dia_5316', 'ApP'): [1,0,0,0,0]} for example, to store the chimerc seq count for each chimera-killing step

    ## Go through each locus ##
    for locus_folder in locus_list: # locus_folder = locus name
        locusIndex = locus_list.index(locus_folder)
        try:
            os.chdir(locus_folder)
        except:
            log.write("WARNING: %s directory not found" % locus_folder)
            continue
        sys.stderr.write('\nWorking on: ' + locus_folder + '...\n')
        #sys.stderr.write("Forward primer: %s\n" % Forward_primer[locusIndex])
        #sys.stderr.write("Reverse primer: %s\n" % Reverse_primer[locusIndex])
        if verbose_level in [1,2]:
            log.write('\nWorking on ' + str(locus_folder) + ' ...\n')
            log.write("Forward primer: %s\n" % Forward_primer[locusIndex])
            log.write("Reverse primer: %s\n" % Reverse_primer[locusIndex])
        if not os.stat(locus_folder + ".fa").st_size == 0: # ie, the file is not empty
            ## Split sequences into separate taxon folders ##
            taxonCounts = SplitBy(annotd_seqs_file = locus_folder + ".fa", split_by = "taxon", Multiplex_perBC_flag = Multiplex_per_barcode)
            all_folders_taxon = list(taxonCounts.keys())

            for taxon_folder in all_folders_taxon:
                if verbose_level in [1,2]:
                    log.write("Working on " + taxon_folder + '\n')
                os.chdir(taxon_folder)
                subset_fasta_seqs_from_fastq("%s.fa" % taxon_folder, raw_fastq_sequences)
                if mode != 3:
                    writeASV(taxon_folder, Forward_primer[locusIndex], Reverse_primer[locusIndex], minLen, maxLen, maxEE)
                    with open("%s_DADA2.log" % taxon_folder, "w") as logfile:
                        dadaCMD = "%s %s_DADA2.R" %(RscriptPath, taxon_folder)
                        process = subprocess.Popen(dadaCMD, stdout=logfile, stderr=logfile, shell=True, text=True)
                        process.communicate()
                    ## Count ASVs and store in LocusTaxonCountDict_clustd as {('C_dia_5316', 'ApP'): 28} for example ##
                    try:
                        clustered_seq_file = parse_fasta(taxon_folder + '_ASVs.fa')
                        for each_seq in clustered_seq_file:
                            try:
                                LocusTaxonCountDict_clustd[taxon_folder, locus_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
                            except:
                                LocusTaxonCountDict_clustd[taxon_folder, locus_folder] = 1
                    except:
                        print("WARNING: No ASVs found for %s" % taxon_folder)
                        log.write("WARNING: No ASVs found for %s" % taxon_folder)

                ## Count chimeras
                try:
                    chimera_seq_file = parse_fasta(taxon_folder + '_chimeras.fa')
                    for each_seq in chimera_seq_file:
                        try:
                            LocusTaxonCountDict_chimera[taxon_folder, locus_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
                        except:
                            LocusTaxonCountDict_chimera[taxon_folder, locus_folder] = 1
                except:
                    if verbose_level in [2]:
                        log.write("No chimeras detected in %s\n" % taxon_folder)
                    LocusTaxonCountDict_chimera[taxon_folder, locus_folder] = 0
                os.chdir("..")
        os.chdir("..")
    log.write('\t...done\n\n')

    if mode == 3:
        sys.stderr.write("PURC stopped after annotation\n")
        exit(0)
    ## Put all the sequences together ##
    sys.stderr.write('\n\nPutting all the sequences together...\n\n')
    for locus_folder in all_folders_loci: # Looping through each of the locus folders
        outputfile_name_reconsensus = Output_prefix + '_4_' + str(locus_folder) + '_ASVs.fa'
        outputfile_reconsensus = open(outputfile_name_reconsensus, 'w')

        os.chdir(locus_folder)
        taxonForThisLocus = glob.glob("*")
        for taxon_folder in taxonForThisLocus: # have to go into each barcode folder in each locus folder
            if os.path.isdir(taxon_folder): # the glob might have found some files as well as folders
                os.chdir(taxon_folder)
                files_to_add_reconsensus = glob.glob("*_ASVs.fa")
                for file in files_to_add_reconsensus:
                    shutil.copyfileobj(open(file,'r'), outputfile_reconsensus) #Add each file to the final output

                os.chdir('..')
        os.chdir('..')
        outputfile_reconsensus.close()

    return LocusTaxonCountDict_clustd, LocusTaxonCountDict_chimera

'''
Smith-Waterman aligner, modified from github.com/mbreese/swalign
'''
class ScoringMatrix(object):
    '''
    Read scoring matrix from a file or string

    Matrix should be space-delimited in a format like:

      A C G T
    A 1 0 0 0
    C 0 1 0 0
    G 0 0 1 0
    T 0 0 0 1

    Rows and Columns must be in the same order

    '''
    def __init__(self, filename=None, text=None, wildcard_score=0):
        assert filename or text

        if filename:
            fs = open(filename)
        else:
            fs = io.StringIO(text)
        self.scores = []
        self.bases = None
        self.wildcard_score = wildcard_score
        for line in fs:
            if line[0] == '#':
                continue

            if not self.bases:
                self.bases = line.split()
                self.base_count = len(self.bases)
            else:
                cols = line.split()
                self.scores.extend([int(x) for x in cols[1:]])
        fs.close()

    def score(self, one, two, wildcard=None):
        if self.wildcard_score and wildcard and (one in wildcard or two in wildcard):
            return self.wildcard_score

        one_idx = 0
        two_idx = 0
        for i, b in enumerate(self.bases):
            if b == one:
                one_idx = i
            if b == two:
                two_idx = i

        return self.scores[(one_idx * self.base_count) + two_idx]

class IdentityScoringMatrix(object):
    def __init__(self, match=1, mismatch=-1):
        self.match = match
        self.mismatch = mismatch

    def score(self, one, two, wildcard=None):
        if wildcard and (one in wildcard or two in wildcard):
            return self.match

        if one == two:
            return self.match
        return self.mismatch

NucleotideScoringMatrix = IdentityScoringMatrix

class Matrix(object):
    def __init__(self, rows, cols, init=None):
        self.rows = rows
        self.cols = cols
        self.values = [init, ] * rows * cols

    def get(self, row, col):
        return self.values[(row * self.cols) + col]

    def set(self, row, col, val):
        self.values[(row * self.cols) + col] = val

class LocalAlignment(object):
    def __init__(self, scoring_matrix, gap_penalty=-1, gap_extension_penalty=-1, gap_extension_decay=0.0, prefer_gap_runs=True, verbose=False, globalalign=False, wildcard=None, full_query=False):
        self.scoring_matrix = scoring_matrix
        self.gap_penalty = gap_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.gap_extension_decay = gap_extension_decay
        self.verbose = verbose
        self.prefer_gap_runs = prefer_gap_runs
        self.globalalign = globalalign
        self.wildcard = wildcard
        self.full_query = full_query

    def align(self, ref, query, ref_name='', query_name='', rc=False):
        orig_ref = ref
        orig_query = query
        ref = ref.upper()
        query = query.upper()
        matrix = Matrix(len(query) + 1, len(ref) + 1, (0, ' ', 0))
        for row in range(1, matrix.rows):
            matrix.set(row, 0, (0, 'i', 0))
        for col in range(1, matrix.cols):
            matrix.set(0, col, (0, 'd', 0))
        max_val = 0
        max_row = 0
        max_col = 0
        # calculate matrix
        for row in range(1, matrix.rows):
            for col in range(1, matrix.cols):
                mm_val = matrix.get(row - 1, col - 1)[0] + self.scoring_matrix.score(query[row - 1], ref[col - 1], self.wildcard)
                ins_run = 0
                del_run = 0
                if matrix.get(row - 1, col)[1] == 'i':
                    ins_run = matrix.get(row - 1, col)[2]
                    if matrix.get(row - 1, col)[0] == 0:
                        # no penalty to start the alignment
                        ins_val = 0
                    else:
                        if not self.gap_extension_decay:
                            ins_val = matrix.get(row - 1, col)[0] + self.gap_extension_penalty
                        else:
                            ins_val = matrix.get(row - 1, col)[0] + min(0, self.gap_extension_penalty + ins_run * self.gap_extension_decay)
                else:
                    ins_val = matrix.get(row - 1, col)[0] + self.gap_penalty
                if matrix.get(row, col - 1)[1] == 'd':
                    del_run = matrix.get(row, col - 1)[2]
                    if matrix.get(row, col - 1)[0] == 0:
                        # no penalty to start the alignment
                        del_val = 0
                    else:
                        if not self.gap_extension_decay:
                            del_val = matrix.get(row, col - 1)[0] + self.gap_extension_penalty
                        else:
                            del_val = matrix.get(row, col - 1)[0] + min(0, self.gap_extension_penalty + del_run * self.gap_extension_decay)

                else:
                    del_val = matrix.get(row, col - 1)[0] + self.gap_penalty

                if self.globalalign or self.full_query:
                    cell_val = max(mm_val, del_val, ins_val)
                else:
                    cell_val = max(mm_val, del_val, ins_val, 0)
                if not self.prefer_gap_runs:
                    ins_run = 0
                    del_run = 0
                if del_run and cell_val == del_val:
                    val = (cell_val, 'd', del_run + 1)
                elif ins_run and cell_val == ins_val:
                    val = (cell_val, 'i', ins_run + 1)
                elif cell_val == mm_val:
                    val = (cell_val, 'm', 0)
                elif cell_val == del_val:
                    val = (cell_val, 'd', 1)
                elif cell_val == ins_val:
                    val = (cell_val, 'i', 1)
                else:
                    val = (0, 'x', 0)
                if val[0] >= max_val:
                    max_val = val[0]
                    max_row = row
                    max_col = col

                matrix.set(row, col, val)
        # backtrack
        if self.globalalign:
            # backtrack from last cell
            row = matrix.rows - 1
            col = matrix.cols - 1
            val = matrix.get(row, col)[0]
        elif self.full_query:
            # backtrack from max in last row
            row = matrix.rows - 1
            max_val = 0
            col = 0
            for c in range(1, matrix.cols):
                if matrix.get(row, c)[0] > max_val:
                    col = c
                    max_val = matrix.get(row, c)[0]
            col = matrix.cols - 1
            val = matrix.get(row, col)[0]
        else:
            # backtrack from max
            row = max_row
            col = max_col
            val = max_val
        op = ''
        aln = []
        path = []
        while True:
            val, op, runlen = matrix.get(row, col)

            if self.globalalign:
                if row == 0 and col == 0:
                    break
            elif self.full_query:
                if row == 0:
                    break
            else:
                if val <= 0:
                    break
            path.append((row, col))
            aln.append(op)
            if op == 'm':
                row -= 1
                col -= 1
            elif op == 'i':
                row -= 1
            elif op == 'd':
                col -= 1
            else:
                break
        aln.reverse()
        if self.verbose:
            self.dump_matrix(ref, query, matrix, path)
            print(aln)
            print((max_row, max_col), max_val)
        cigar = _reduce_cigar(aln)
        return Alignment(orig_query, orig_ref, row, col, cigar, max_val, ref_name, query_name, rc, self.globalalign, self.wildcard)
    def dump_matrix(self, ref, query, matrix, path, show_row=-1, show_col=-1):
        sys.stdout.write('      -      ')
        sys.stdout.write('       '.join(ref))
        sys.stdout.write('\n')
        for row in range(matrix.rows):
            if row == 0:
                sys.stdout.write('-')
            else:
                sys.stdout.write(query[row - 1])

            for col in range(matrix.cols):
                if show_row == row and show_col == col:
                    sys.stdout.write('       *')
                else:
                    sys.stdout.write(' %5s%s%s' % (matrix.get(row, col)[0], matrix.get(row, col)[1], '$' if (row, col) in path else ' '))
            sys.stdout.write('\n')
def _reduce_cigar(operations):
    count = 1
    last = None
    ret = []
    for op in operations:
        if last and op == last:
            count += 1
        elif last:
            ret.append((count, last.upper()))
            count = 1
        last = op

    if last:
        ret.append((count, last.upper()))
    return ret
def _cigar_str(cigar):
    out = ''
    for num, op in cigar:
        out += '%s%s' % (num, op)
    return out
class Alignment(object):
    def __init__(self, query, ref, q_pos, r_pos, cigar, score, ref_name='', query_name='', rc=False, globalalign=False, wildcard=None):
        self.query = query
        self.ref = ref
        self.q_pos = q_pos
        self.r_pos = r_pos
        self.cigar = cigar
        self.score = score
        self.r_name = ref_name
        self.q_name = query_name
        self.rc = rc
        self.globalalign = globalalign
        self.wildcard = wildcard
        self.r_offset = 0
        self.r_region = None
        self.orig_query = query
        self.query = query.upper()
        self.orig_ref = ref
        self.ref = ref.upper()
        q_len = 0
        r_len = 0
        self.matches = 0
        self.mismatches = 0
        i = self.r_pos
        j = self.q_pos
        for count, op in self.cigar:
            if op == 'M':
                q_len += count
                r_len += count
                for k in range(count):
                    if self.query[j] == self.ref[i]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    i += 1
                    j += 1
            elif op == 'I':
                q_len += count
                j += count
                self.mismatches += count
            elif op == 'D':
                r_len += count
                i += count
                self.mismatches += count
        self.q_end = q_pos + q_len
        self.r_end = r_pos + r_len
        if self.mismatches + self.matches > 0:
            self.identity = float(self.matches) / (self.mismatches + self.matches)
        else:
            self.identity = 0
    def set_ref_offset(self, ref, offset, region):
        self.r_name = ref
        self.r_offset = offset
        self.r_region = region

    @property
    def extended_cigar_str(self):
        qpos = 0
        rpos = 0
        ext_cigar_str = ''
        working = []
        for count, op in self.cigar:
            if op == 'M':
                for k in range(count):
                    if self.query[self.q_pos + qpos + k] == self.ref[self.r_pos + rpos + k]:
                        ext_cigar_str += 'M'
                    else:
                        ext_cigar_str += 'X'
                qpos += count
                rpos += count

            elif op == 'I':
                qpos += count
                ext_cigar_str += 'I' * count
            elif op == 'D':
                rpos += count
                ext_cigar_str += 'D' * count
            working = _reduce_cigar(ext_cigar_str)
        out = ''
        for num, op in working:
            out += '%s%s' % (num, op)
        return out

    @property
    def cigar_str(self):
        return _cigar_str(self.cigar)

    def dump(self, wrap=None, out=sys.stdout):
        i = self.r_pos
        j = self.q_pos
        q = ''
        m = ''
        r = ''
        qlen = 0
        rlen = 0
        for count, op in self.cigar:
            if op == 'M':
                qlen += count
                rlen += count
                for k in range(count):
                    q += self.orig_query[j]
                    r += self.orig_ref[i]
                    if self.query[j] == self.ref[i] or (self.wildcard and (self.query[j] in self.wildcard or self.ref[i] in self.wildcard)):
                        m += '|'
                    else:
                        m += '.'
                    i += 1
                    j += 1
            elif op == 'D':
                rlen += count
                for k in range(count):
                    q += '-'
                    r += self.orig_ref[i]
                    m += ' '
                    i += 1
            elif op == 'I':
                qlen += count
                for k in range(count):
                    q += self.orig_query[j]
                    r += '-'
                    m += ' '
                    j += 1
            elif op == 'N':
                q += '-//-'
                r += '-//-'
                m += '    '
        if self.q_name:
            out.write('Query: %s%s (%s nt)\n' % (self.q_name, ' (reverse-compliment)' if self.rc else '', len(self.query)))
        if self.r_name:
            if self.r_region:
                out.write('Ref  : %s (%s)\n\n' % (self.r_name, self.r_region))
            else:
                out.write('Ref  : %s (%s nt)\n\n' % (self.r_name, len(self.ref)))
        poslens = [self.q_pos + 1, self.q_end + 1, self.r_pos + self.r_offset + 1, self.r_end + self.r_offset + 1]
        maxlen = max([len(str(x)) for x in poslens])
        q_pre = 'Query: %%%ss ' % maxlen
        r_pre = 'Ref  : %%%ss ' % maxlen
        m_pre = ' ' * (8 + maxlen)
        rpos = self.r_pos
        if not self.rc:
            qpos = self.q_pos
        else:
            qpos = self.q_end

        while q and r and m:
            if not self.rc:
                out.write(q_pre % (qpos + 1))  # pos is displayed as 1-based
            else:
                out.write(q_pre % (qpos))  # revcomp is 1-based on the 3' end

            if wrap:
                qfragment = q[:wrap]
                mfragment = m[:wrap]
                rfragment = r[:wrap]

                q = q[wrap:]
                m = m[wrap:]
                r = r[wrap:]
            else:
                qfragment = q
                mfragment = m
                rfragment = r
                q = ''
                m = ''
                r = ''
            out.write(qfragment)
            if not self.rc:
                for base in qfragment:
                    if base != '-':
                        qpos += 1
            else:
                for base in qfragment:
                    if base != '-':
                        qpos -= 1

            if not self.rc:
                out.write(' %s\n' % qpos)
            else:
                out.write(' %s\n' % (qpos + 1))
            out.write(m_pre)
            out.write(mfragment)
            out.write('\n')
            out.write(r_pre % (rpos + self.r_offset + 1))
            out.write(rfragment)
            for base in rfragment:
                if base != '-':
                    rpos += 1
            out.write(' %s\n\n' % (rpos + self.r_offset))
        out.write("Score: %s\n" % self.score)
        out.write("Matches: %s (%.1f%%)\n" % (self.matches, self.identity * 100))
        out.write("Mismatches: %s\n" % (self.mismatches,))
        out.write("CIGAR: %s\n" % self.cigar_str)


################################################ Setup ################################################
ts = time.time()
time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
print("Start Time: %s" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

if len(sys.argv) < 2:
    sys.exit(usage)

elif sys.argv[1] in ['-help', '-h', '-citation']:
    sys.exit(usage + citation)

else:
    try:
        configuration = open(sys.argv[1], 'r')
    except:
        sys.stderr.write('Error: Cannot open the configuration file\n')
        sys.exit(usage)
    ppp_location = os.getcwd()

    ## Setting up defaults ##
    parameterDict = {"mode" : 1, "Check_chimeras" : False, "Multiplex_per_barcode" : False, "Dual_barcode" : False, "Search_ends_only" : True, "Recycle_bc" : False, "Align" : 0, "minLen" : 0, "maxLen" : 0, "maxEE" : 10, "clustID" : 0.997, "clustID2" : 0.995, "clustID3" : 0.990, "clustID4" : 0.997, "sizeThreshold" : 1, "sizeThreshold2" : 4, "abskew" : 1.9, "verbose_level" : 0, "num_threads" : 1, "barcode_databasefile" : 'barcode_blastdb', "refseq_databasefile" : 'refseq_blastdb',"log_file" : 'purc_log_' + time_stamp + '.txt', "useOTUpriors" : 'False', "Clustering_method" : "OTU" }
    mode = 1
    Check_chimeras = False
    Multiplex_per_barcode = False
    Dual_barcode = False
    Search_ends_only = True
    Recycle_bc = False
    Align = 0
    minLen = 0
    maxLen = 0
    maxEE = 10
    clustID = 0.997
    clustID2 = 0.995
    clustID3 = 0.990
    clustID4 = 0.997
    sizeThreshold = 1
    sizeThreshold2 = 4
    abskew = 1.9
    verbose_level = 0
    num_threads = 1
    barcode_databasefile = 'barcode_blastdb'
    refseq_databasefile = 'refseq_blastdb'
    seq_name_toErase = ''
    Vsearch = 'vsearch'
    Cutadapt = 'cutadapt'
    #Muscle = 'muscle'
    Mafft = 'mafft'
    RscriptPath = 'Rscript'
    log_file = 'purc_log_' + time_stamp + '.txt'
    useOTUpriors = "FALSE"

    ## Read-in the parameters and settings ##
    for line in configuration:
        line = line.strip(' \t\n\r').replace(' ', '').replace('\t', '').split('#')[0]
        setting_line = line.find('=')
        if setting_line != -1:
            setting = line.split('=')
            setting_name = setting[0]
            setting_argument = setting[1]
            if setting_name == 'Mode':
                mode = int(setting_argument)
                parameterDict['mode'] = setting_argument
            elif setting_name == 'Input_sequence_file':
                raw_sequences = setting_argument
                parameterDict['raw_sequences'] = setting_argument
            elif setting_name == 'Input_sequence_dir':
                demux_input_dir = setting_argument + "/"
                parameterDict['demux_input_dir'] = demux_input_dir
            elif setting_name == "Align":
                Align = int(setting_argument)
                parameterDict['Align'] = setting_argument
            elif setting_name == 'Output_prefix':
                Output_prefix = setting_argument
                parameterDict['Output_prefix'] = setting_argument
            elif setting_name == 'Output_folder':
                Output_folder = os.path.abspath(setting_argument)
                parameterDict['Output_folder'] = setting_argument
                print("Output folder: %s" % Output_folder)
                BLAST_DBs_folder = Output_folder + '_BlastDBs'
            elif setting_name == 'Barcode_blastDB':
                barcode_databasefile = setting_argument
                parameterDict['barcode_databasefile'] = setting_argument
            elif setting_name == 'RefSeq_blastDB':
                refseq_databasefile = setting_argument
                parameterDict['refseq_databasefile'] = setting_argument
            elif setting_name == 'Log_file':
                if setting_argument == '':
                    log_file = 'purc_log_' + time_stamp + '.txt'
                else:
                    log_file = setting_argument
                    parameterDict['log_file'] = setting_argument
            elif setting_name == 'Locus_name':
                locus_list = setting_argument.upper().replace(' ', '').replace('\t', '').split(',') #needs the upper() now that LocusTaxonCountDict_unclustd has the loci in uppercase
                parameterDict['locus_list'] = locus_list
            elif setting_name == 'Locus-barcode-taxon_map':
                mapping_file_list = setting_argument.replace(' ', '').replace('\t', '').split(',')
                parameterDict['mapping_file_list'] = mapping_file_list
				#mapping_file_list = []
                #for mapfile in mapping_file_list_tmp:
                #    mapfile = os.path.abspath(mapfile)
                #    if not os.path.isfile(mapfile):
                #        sys.exit("Error: could not find " + mapfile)
                #    else:
                #        mapping_file_list.append(mapfile)
                #parameterDict['mapping_file_list'] = mapping_file_list
            elif setting_name == 'Vsearch':
                if setting_argument.startswith('Dependencies/'):
                    Vsearch = ppp_location + '/' + setting_argument
                else:
                    Vsearch = setting_argument
            elif setting_name == 'Cutadapt':
                if setting_argument.startswith('Dependencies/'):
                    Cutadapt = ppp_location + '/' + setting_argument
                else:
                    Cutadapt = setting_argument
            #elif setting_name == 'Muscle':
            #    if setting_argument.startswith('Dependencies/'):
            #        Muscle = ppp_location + '/' + setting_argument
            #    else:
            #        Muscle = setting_argument
            elif setting_name == 'MAFFT':
                if setting_argument.startswith('Dependencies/'):
                    Mafft = ppp_location + '/' + setting_argument
                else:
                    Mafft = setting_argument
            elif setting_name == 'Rscript':
                RscriptPath = setting_argument
            elif setting_name == "Lima_override":
                Lima_override = setting_argument
                parameterDict['Lima_override'] = setting_argument
            elif setting_name == "minLen":
                minLen = setting_argument
                parameterDict['minLen'] = setting_argument
            elif setting_name == "maxLen":
                maxLen = setting_argument
                parameterDict['maxLen'] = setting_argument
            elif setting_name == "maxEE":
                maxEE = setting_argument
                parameterDict['maxEE'] = setting_argument
            elif setting_name == "Use_OTU_priors":
                useOTUpriors = str(setting_argument).upper()
                parameterDict['useOTUpriors'] = useOTUpriors
            elif setting_name == 'clustID1':
                clustID = float(setting_argument)
                parameterDict['clustID'] = clustID
            elif setting_name == 'clustID2':
                clustID2 = float(setting_argument)
                parameterDict['clustID2'] = clustID2
            elif setting_name == 'clustID3':
                clustID3 = float(setting_argument)
                parameterDict['clustID3'] = clustID3
            elif setting_name == 'clustID4':
                clustID4 = float(setting_argument)
                parameterDict['clustID4'] = clustID4
            elif setting_name == 'sizeThreshold1':
                sizeThreshold = float(setting_argument)
                parameterDict['sizeThreshold'] = sizeThreshold
            elif setting_name == 'sizeThreshold2':
                sizeThreshold2 = float(setting_argument)
                parameterDict['sizeThreshold2'] = sizeThreshold2
            elif setting_name == 'abundance_skew':
                abskew = str(setting_argument)
                parameterDict['abskew'] = abskew
            elif setting_name == 'Forward_primer':
                Forward_primer = setting_argument.replace(' ', '').replace('\t', '').split(',')
                parameterDict['Forward_primer'] = Forward_primer
            elif setting_name == 'Reverse_primer':
                Reverse_primer = setting_argument.replace(' ', '').replace('\t', '').split(',')
                parameterDict['Reverse_primer'] = Reverse_primer
            elif setting_name == 'seq_name_toErase':
                seq_name_toErase = setting_argument
                parameterDict['seq_name_toErase'] = seq_name_toErase
            elif setting_name == 'Verbose_level':
                verbose_level = int(setting_argument)
                parameterDict['verbose_level'] = verbose_level
            elif setting_name == 'Threads':
                num_threads = int(setting_argument)
                parameterDict['num_threads'] = num_threads
            elif setting_name == 'Remove_intermediates':
                remove_intermediates = int(setting_argument)
                parameterDict['remove_intermediates'] = remove_intermediates
            elif setting_name == 'in_Barcode_seq_file':
                barcode_seq_filename = os.path.abspath(setting_argument)
                parameterDict['barcode_seq_filename'] = setting_argument
                #if not os.path.isfile(barcode_seq_filename):
                #    sys.exit("Error: could not find " + barcode_seq_filename)
            elif setting_name == 'in_RefSeq_seq_file':
                refseq_filename = os.path.abspath(setting_argument)
                parameterDict['refseq_filename'] = setting_argument
                if not os.path.isfile(refseq_filename):
                    sys.exit("Error: could not find " + refseq_filename)
            elif setting_name == 'Dual_barcode':
                if setting_argument == '0':
                    Dual_barcode = False
                    Lima_barcode_type = "single-side"
                    parameterDict['Dual_barcode'] = Dual_barcode
                    parameterDict['Lima_barcode_type'] = Lima_barcode_type
                elif setting_argument == '1':
                    Dual_barcode = True
                    Lima_barcode_type = "different"
                    parameterDict['Dual_barcode'] = Dual_barcode
                    parameterDict['Lima_barcode_type'] = Lima_barcode_type
                elif setting_argument == '2':
                    Dual_barcode = True
                    Lima_barcode_type = "same"
                    parameterDict['Dual_barcode'] = Dual_barcode
                    parameterDict['Lima_barcode_type'] = Lima_barcode_type
                else:
                    sys.exit('Error: incorrect setting of Dual_barcode')
            elif setting_name == 'Multiplex_per_barcode':
                if setting_argument == '0':
                    Multiplex_per_barcode = False
                    parameterDict['Multiplex_per_barcode'] = Multiplex_per_barcode
                elif setting_argument == '1':
                    Multiplex_per_barcode = True
                    parameterDict['Multiplex_per_barcode'] = Multiplex_per_barcode
                else:
                    sys.exit('Error: incorrect setting of Multiplex_per_barcode')
            elif setting_name == 'Barcode_detection':
                if setting_argument == '0':
                    Search_ends_only = False
                    parameterDict['Search_ends_only'] = Search_ends_only
                elif setting_argument == '1':
                    Search_ends_only = True
                    parameterDict['Search_ends_only'] = Search_ends_only
                else:
                    sys.exit('Error: incorrect setting of Barcode_detection')
            elif setting_name == 'Recycle_no_barcoded_seq':
                if setting_argument == '0':
                    Recycle_bc = False
                    parameterDict['Recycle_bc'] = Recycle_bc
                elif setting_argument == '1':
                    Recycle_bc = True
                    parameterDict['Recycle_bc'] = Recycle_bc
                else:
                    sys.exit('Error: incorrect setting of Recycle_no_barcoded_seq')
            elif setting_name == 'Recycle_chimeric_seq':
                if setting_argument == '0':
                    Recycle_chimera = False
                    parameterDict['Recycle_chimera'] = Recycle_chimera
                elif setting_argument == '1':
                    Recycle_chimera = True
                    parameterDict['Recycle_chimera'] = Recycle_chimera
                else:
                    sys.exit('Error: incorrect setting of Recycle_chimeric_seq')
            elif setting_name == 'Clustering_method':
                if setting_argument == '0':
                    Clustering_method = "ASV"
                    parameterDict['Clustering_method'] = Clustering_method
                elif setting_argument == '1':
                    Clustering_method = "OTU"
                    parameterDict['Clustering_method'] = Clustering_method
                elif setting_argument == '2':
                    Clustering_method = "BOTH"
                    parameterDict['Clustering_method'] = Clustering_method
                else:
                    sys.exit('Error: incorrect setting of Clustering_method')
    # Check is sequence file type is compatible with clustering method
    fileExt = raw_sequences.split(".")[-1]
    if fileExt == "gz":
        fileExt = raw_sequences.split(".")[-2]
    if fileExt in ["fasta", "fa", "fas", "fna", "faa"] and Clustering_method in ["ASV", "BOTH"]:
        print("Error: ASV inference requires a FASTQ file. Your input appears to be FASTA format. Change to 'Clustering_method = 1' in config file.")
        sys.exit(1)

    # Check if settings are compatible
    if Clustering_method != "BOTH" and useOTUpriors == "TRUE":
        print("Error: Clustering_method = 2 must be set to use OTU priors for ASV inference.")
        sys.exit(1)

    if mode == 3:
        sys.stderr.write("Running annotation-only mode...\n")

    # Check if dependencies are in place
    sys.stderr.write('Checking dependencies...\n')
    # Check if muscle can be executed
    #muscle_cline = '%s -version' % (Muscle)
    #print(muscle_cline)
    #process = subprocess.Popen(muscle_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    #(out, err) = process.communicate() #the stdout and stderr
    #if not str(out).startswith("MUSCLE"):
    #    sys.exit("Error: could not execute Muscle")

    # Check if MAFFT can be executed
    if not shutil.which(Mafft):
        sys.exit("Error: could not execute MAFFT")
    else:
        mafftPath = shutil.which(Mafft)

    if not os.access(mafftPath, os.X_OK):
        if not shutil.which("mafft"):
            sys.exit("Error: could not execute MAFFT")
        elif os.access(shutil.which("mafft"), os.X_OK):
            print("WARNING: MAFFT executable found at %s may be different than one specified in config file" % (shutil.which("mafft")))
        else:
            sys.exit("Error: could not execute MAFFT")

    # Check if Vsearch can be executed
    vsearch_cline = '%s --version' % (Vsearch)
    #print(vsearch_cline)
    process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    (out, err) = process.communicate() #the stdout and stderr
    if not str(err).startswith("vsearch"):
        sys.exit("Error: could not execute vsearch")

    # Check if cutadapt can be executed
    cutadapt_cline = '%s --help' % (Cutadapt)
    #print(cutadapt_cline)
    process = subprocess.Popen(cutadapt_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    (out, err) = process.communicate() #the stdout and stderr
    #print out
    if not str(out).startswith("cutadapt"):
        if not str(out).startswith("Usage"): #for older cutadapt version
            sys.exit("Error: could not execute Cutadapt")

    # Check if blast can be execuated
    blast_cline = 'blastn -version'
    #print(blast_cline)
    process = subprocess.Popen(blast_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    (out, err) = process.communicate() #the stdout and stderr
    if not str(out).replace(' ','').startswith("blastn"):
        sys.exit("Error: could not execute BLAST")

    # Check if Rscript can be executed
    Rscript_cline = 'Rscript --version'
    #print(Rscript_cline)
    process = subprocess.Popen(Rscript_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    (out, err) = process.communicate() #the stdout and stderr
    if not str(err).replace(' ','').startswith("R"):
        sys.exit("Error: could not execute Rscript")

    # Check for input files if using pre-demultiplexed files
    if mode == 2:
        print("Checking demultiplexed input files...")
        try:
            demux_input_files = [file for file in listdir(parameterDict['demux_input_dir']) if isfile(join(parameterDict['demux_input_dir'], file))]
            if len(demux_input_files) == 0:
                print("\tERROR: No files found in %s" % parameterDict['demux_input_dir'])
                exit(1)
            else:
                parameterDict['demux_input_files'] = demux_input_files
                print("\tFound files:")
                for i in demux_input_files:
                    print("\t%s" % i)
        except:
            print("\tERROR: Could not parse %s" % parameterDict['demux_input_dir'])
            exit(1)

	# Check for other input files
    if mode != 2:
		# Raw sequences
        if not os.path.isfile(parameterDict['raw_sequences']):
            sys.exit("Error: could not find sequence file " + parameterDict['raw_sequences'])
        else:
            raw_sequences = os.path.abspath(parameterDict['raw_sequences'])
            print("Input sequence file: %s" % parameterDict['raw_sequences'])

		# Map files
        for file in parameterDict['mapping_file_list']:
            mapfile = os.path.abspath(file)
            if not os.path.isfile(mapfile):
                sys.exit("Error: could not find map file " + mapfile)

		# Barcode file
        if not os.path.isfile(parameterDict['barcode_seq_filename']):
            sys.exit("Error: could not find barcode file " + barcode_seq_filename)


    ## Make output folder or read in previously completed steps ##
    if os.path.exists(Output_folder): # overwrite existing folder
        try:
            with open("%s/tmp/run_settings.pkl" % Output_folder, "rb") as previous_settings_file:
                previousSettingsDict = pickle.load(previous_settings_file)
                if parameterDict == previousSettingsDict:
                    with open("%s/tmp/checkpoint.txt" % Output_folder, "r") as checkpoint_file:
                        checkpoints_complete = []
                        for line in checkpoint_file:
                            checkpoints_complete.append(line.strip("\n"))
                else:
                    print("Error: Config file settings don't match those used before. Change config file back to continue previous run.")
                    sys.exit(1)
        except:
            shutil.rmtree(Output_folder)
            os.makedirs(Output_folder)
            os.makedirs("%s/tmp" % Output_folder)
            checkpoints_complete = []
            with open("%s/tmp/run_settings.pkl" % Output_folder, "wb") as settings_file:
                pickle.dump(parameterDict, settings_file)
    else:
        os.makedirs(Output_folder)
        os.makedirs("%s/tmp" % Output_folder)
        checkpoints_complete = []
        with open("%s/tmp/run_settings.pkl" % Output_folder, "wb") as settings_file:
            pickle.dump(parameterDict, settings_file)

    log = open("%s/%s" %(Output_folder,log_file), 'w')
    log.write("Start Time: %s\n" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    log.write(logo + '\n')
    log.write("PURC called with: \n\t" + str(sys.argv) + "\n")
    log.write('Vsearch location: ' + str(Vsearch) + '\n')
    log.write('Cutadapt location: ' + str(Cutadapt) + '\n')
    #log.write('Muscle location: ' + str(Muscle) + '\n')
    log.write("Settings for this run:\n" + "\tSequence file:\t" + str(raw_sequences) + "\n\tLoci:\t" + '\t'.join(locus_list) + '\n')
    log.write("\tMapping files: " + ', '.join(mapping_file_list) + '\n')
    log.write("\tForward primers: " + ', '.join(Forward_primer) + '\n')
    log.write("\tReverse primers: " + ', '.join(Reverse_primer) + '\n')

    if Dual_barcode:
        log.write("\tExpecting barcodes at each end of the sequence\n")
    else:
        log.write("\tExpecting barcodes on forward primers only\n")
    if Multiplex_per_barcode:
        log.write("\tExpecting barcodes to be shared across multiple taxa (genera, etc)\n")
    else:
        log.write("\tExpecting each barcode to be used for only a single taxon\n")
    if Search_ends_only:
        log.write("\tBarcodes will be looked for in the terminal 25 bases of each sequence, only\n")
    else:
        log.write("\tThe full sequence will be searched for primers; internal primers may be found\n")
    if Recycle_bc:
        log.write('''\tIf the BLAST-based approach doesn't find a barcode in a particular sequence
        the sequence will be re-searched using a Smith-Waterman pairwise alignment approach\n''')
        #try:
        #    import swalign
        #except:
        #    sys.exit("Error: could not import swalign; maybe turn off the recycle barcode option")
    else:
        log.write("\tPrimer-searching will be done using the BLAST-based approach, only\n")
    if Recycle_chimera:
        log.write("\tInter-locus chimeras will be split into their component pieces and fed into the pipeline\n")
    else:
        log.write("\tInter-locus chimeras will be removed from the analysis\n")

    log.write("\tSimilarity cut-off for clustering:\t" + str(clustID) + "\t" + str(clustID2) + "\t" + str(clustID3) + '\n')
    log.write("\tCluster size for retention:\t" + str(sizeThreshold) + '\n')




################################################ RUN YO!!! ########################################
sys.stderr.write('Renaming sequences...\n')
fileExt = raw_sequences.split(".")[-1]
filePrefix = ".".join(raw_sequences.split(".")[:-1])
if fileExt == "gz":
    fileType = filePrefix.split(".")[-1]
    if os.path.isfile(filePrefix):
        print("WARNING: Unzipped file with same name already exists. Will use what I assume is the uncompressed version of the read file...")
        #gunzipCmd = "gunzip -c -k %s > tmp_sequences.%s" %(raw_sequences, fileType)
        #process = subprocess.Popen(gunzipCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        #process.wait()
        #raw_sequences = os.path.abspath("tmp_sequences.%s" % fileType)
        raw_sequences = filePrefix
    else:
        gunzipCmd = "gunzip -k %s" % raw_sequences
        process = subprocess.Popen(gunzipCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        process.wait()
        raw_sequences = filePrefix

if mode != 2:
	fileType = raw_sequences.split(".")[-1]
	if fileType in ["fasta", "fas", "fa", "fna", "faa"]:
	    try:
	        fasta_sequences = rename_fasta(raw_sequences)
	    except:
	        sys.exit('ERROR: failed to rename ' + raw_sequences)
	elif fileType in ["fastq", "fq"]:
	    try:
	        fastq_sequences = rename_fastq(raw_sequences)
	    except:
	        sys.exit("ERROR: failed to rename %s" % raw_sequences)
	    try:
	        fasta_sequences = convert_fastq_to_fasta(fastq_sequences)
	    except:
	        sys.exit("ERROR: failed to convert %s" % raw_sequences)
	else:
	    sys.exit("ERROR: Sequence file type not recognized. Expects standard file extensions (.fasta, .fa, .fastq, .fq)")


if os.path.exists(BLAST_DBs_folder): # overwrite existing folder
    shutil.rmtree(BLAST_DBs_folder)
os.makedirs(BLAST_DBs_folder)
os.chdir(BLAST_DBs_folder)
try:
    makeBlastDB(refseq_filename, refseq_databasefile) # and one of the reference sequences
    if mode != 2:
        makeBlastDB(barcode_seq_filename, barcode_databasefile) # one of the barcodes
except:
    sys.exit('ERROR: failed to make blast database')
os.chdir('..')

## Read sequences ##
if mode != 2:
    sys.stderr.write('Reading sequences...\n')
    try:
        SeqDict = SeqIO.index(fasta_sequences, 'fasta') # Read in the raw sequences as dictionary, using biopython's function
    except:
        sys.exit('ERROR: failed to index ' + fasta_sequences)
    count_total_input_sequences = len(SeqDict)
    sys.stderr.write('\t' + str(count_total_input_sequences) + ' sequences read\n')


if mode == 0 and "concatemerCheck" not in checkpoints_complete: # QC mode.
    ## Check chimeras ##
    log.write('\n#Concatemers Removal#\n')
    sys.stderr.write('Checking for inter-locus chimeric sequences (concatemers)...\n')
    if not Recycle_chimera:
        chimeras_file = Output_prefix + '_0_chimeras.fa'
        non_chimeras_file = Output_prefix + '_0_nonchimeras.fa'
        chimera_dict = CheckChimericLoci(fasta_sequences, Output_folder + '/' + 'blast_full_refseq_out.txt', Output_folder + '/' + non_chimeras_file, Output_folder + '/' + chimeras_file, BLAST_DBs_folder + '/' + refseq_databasefile, SeqDict, SplitChimera=False)
    else:
        chimeras_file = Output_prefix + '_0_chimeras.fa'
        non_chimeras_file = Output_prefix + '_0_nonchimeras+split.fa'
        chimera_dict = CheckChimericLoci(fasta_sequences, Output_folder + '/' + 'blast_full_refseq_out.txt', Output_folder + '/' + non_chimeras_file, Output_folder + '/' + chimeras_file, BLAST_DBs_folder + '/' + refseq_databasefile, SeqDict, SplitChimera=True)

    count_chimeric_sequences = len(chimera_dict)
    fasta_sequences = Output_folder + '/' + non_chimeras_file
    SeqDict = SeqIO.index(fasta_sequences, 'fasta')
    sys.stderr.write('\t' + str(count_chimeric_sequences) + ' inter-locus chimeric sequences found\n')
    log.write('\t' + str(count_chimeric_sequences) + ' inter-locus chimeric sequences found\n')
    if not check_fasta(fasta_sequences):
        sys.exit('Error: concatemers-removal returned no sequence')
    else:
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("concatemerCheck\n")
elif mode == 0 and "concatemerCheck" in checkpoints_complete:
    log.write('\n#Concatemers Removal#\n')
    log.write('Reusing previous concatemer results...\n')
    print('Reusing previous concatemer results...\n')


## Remove barcodes ##
if mode != 2:
	log.write('\n#Barcode Removal# %s \n' % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
	if "barcodeRemoval" in checkpoints_complete:
	    log.write('Reusing previous barcode trimmed files...\n')
	    print('Reusing previous barcode trimmed files...\n')
	else:
	    if platform.system() == 'Linux' and Lima_override == "0":
	        print("Demultiplexing with Lima")
	        dupesFound, BCpairdict = checkDuplicateBC(barcode_seq_filename) # Can't have duplicate barcodes (including reverse complements) in lima
	        if dupesFound == 0:
	            limaOutputPrefix = "%s.demux.lima" % Output_prefix
	            if Lima_barcode_type == "same":
	                lima_symmetric(raw_sequences, Lima_barcode_type, barcode_seq_filename, Output_folder, Output_prefix)
	            elif Lima_barcode_type == "different":
	                lima_dual(raw_sequences, Lima_barcode_type, barcode_seq_filename, Output_folder, Output_prefix)
	            elif Lima_barcode_type == "single-side":
	                lima_singleend(raw_sequences, Lima_barcode_type, barcode_seq_filename, Output_folder, Output_prefix)
	        elif dupesFound == 1: # work around if duplicated barcodes are found. Split duplicated and nonduplicated barcodes and check sequences separately.
	            lima_dual(raw_sequences, Lima_barcode_type, "%s/tmp/nonduplicated_BCs.fasta" % Output_folder, Output_folder, "nonduplicated_BCs")
	            lima_symmetric(raw_sequences, Lima_barcode_type, "%s/tmp/deduplicated_BCs.fasta" % Output_folder, Output_folder, "deduplicated_BCs")
	            nondupeseqs = SeqIO.parse("%s/nonduplicated_BCs_1_bc_trimmed.fa" % Output_folder, "fasta")
	            nondupeseqIDs = []
	            for seq in nondupeseqs:
	                nondupeseqIDs.append(seq.id.split("|")[1])
	            dedupeseqs = SeqIO.parse("%s/deduplicated_BCs_1_bc_trimmed.fa" % Output_folder, "fasta")
	            dedupeseqIDs = []
	            for seq in dedupeseqs:
	                dedupeseqIDs.append(seq.id.split("|")[1])

	            shared_seq_list = list(set(nondupeseqIDs) & set(dedupeseqIDs))
	            if len(shared_seq_list) > 0:
	                print("WARNING: %s sequences identified with both symmetric and asymmetric barcodes" % len(shared_seq_list))

	            # Retrieve sequences from asymmetric demultiplexing that weren't also found in symmetric demultiplexing
	            nondupeseqs = SeqIO.parse("%s/nonduplicated_BCs_1_bc_trimmed.fa" % Output_folder, "fasta")
	            write_nondupeseqs = [i for i in nondupeseqs if i.id not in nondupeseqIDs]
	            SeqIO.write(write_nondupeseqs, "%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix), "fasta")
	            # Retrieve sequences from symmetric demultiplexing that weren't also found in asymmetric demultiplexing. Rename second barcode in pair to original names so can be handled downstream.
	            with open("%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix), "a") as open_bc_trimmed_file:
	                #dedupeseqs = SeqIO.parse("%s/deduplicated_BCs_1_bc_trimmed.fa" % Output_folder, "fasta")
	                with open("%s/deduplicated_BCs_1_bc_trimmed.fa" % Output_folder, "r") as open_dedup_bcs:
	                    writeOut = 0
	                    for line in open_dedup_bcs:
	                        if line.startswith(">") and line.strip(">\n").split("|")[1] not in shared_seq_list:
	                            writeOut = 1
	                            barcodes = line.strip(">\n").split("|")[0]
	                            seqid = line.strip(">\n").split("|")[1]
	                            barcode1 = barcodes.split("^")[0]
	                            barcode2 = BCpairdict[barcode1]
	                            newReadName = ">%s^%s|%s\n" %(barcode1, barcode2, seqid)
	                            open_bc_trimmed_file.write(newReadName)
	                        elif line.startswith(">") and line.strip(">\n").split("|")[1] in shared_seq_list:
	                            writeOut = 0
	                        elif writeOut == 1:
	                            open_bc_trimmed_file.write(line)
	        reorientedSequences = reorientSeqs("%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix), "%s/%s" %(BLAST_DBs_folder,refseq_databasefile), num_threads) # Seqs need to have same orientation for OTU clustering (otherwise get 1 OTUs for each orientation of same sequence). DADA2 includes a reorientation step so not necessary here
	        mvCMD = "cp %s %s/%s_1_bc_trimmed.fa " %(reorientedSequences, Output_folder, Output_prefix)
	        process = subprocess.Popen(mvCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	        process.communicate()
	    # BLAST-based demultiplexing
	    else:
	        print("Demultiplexing with BLAST")
	        if Dual_barcode:
	            sys.stderr.write('Removing dual barcodes...\n')
	            DeBarcoder_dual(fasta_sequences, BLAST_DBs_folder + '/' + barcode_databasefile, SeqDict)
	        else:
	            sys.stderr.write('Removing barcodes...\n')
	            if not Search_ends_only:
	                DeBarcoder(fasta_sequences, BLAST_DBs_folder + '/' + barcode_databasefile, SeqDict, Output_folder, Output_prefix)
	            else:
	                DeBarcoder_ends(SeqDict, BLAST_DBs_folder + '/' + barcode_databasefile, Output_folder, Output_prefix, search_range=25)
	    count_seq_w_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa')
	    count_seq_wo_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa')
	    count_seq_w_toomany_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa')
	    sys.stderr.write('\t' + str(count_seq_w_bc) + ' sequences have barcode\n')
	    sys.stderr.write('\t' + str(count_seq_wo_bc) + ' sequences have no barcode\n')
	    sys.stderr.write('\t' + str(count_seq_w_toomany_bc) + ' sequences have too many barcodes\n')
	    if count_seq_w_bc == str(0):
	        sys.exit('Error: barcode-removal returned no sequence')
	    else:
	        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
	            open_checkpoint_file.write("barcodeRemoval\n")

	    if Recycle_bc:
	        sys.stderr.write('Looking for barcodes in the no-barcode sequences, using Smith-Waterman pairwise alignment...\n')
	        SeqDict_no_bc = SeqIO.index(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa', 'fasta') # Read in the raw sequences as dictionary, using biopython's function
	        count_seq_recycled = DeBarcoder_SWalign(SeqDict_no_bc, barcode_seq_filename, Output_folder, Output_prefix, search_range=25)
	        sys.stderr.write('\t' + str(count_seq_recycled) + ' sequences recycled from ' + str(count_seq_wo_bc) + ' sequences\n')
	log.write('\t...done\n\n')

## Use pre-demultiplexed files
if mode == 2:
	log.write('#Reading Demultiplexed Files# %s \n' % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
	# Convert any files from FASTQ to FASTA
	with open(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', "w") as outfile:
		fastaFiles = []
		demux_sample_map = {}
		for file in parameterDict['demux_input_files']:
			if not check_fasta(demux_input_dir + file):
				convertedFile = convert_fastq_to_fasta(demux_input_dir + file)
				fastaFiles.append(convertedFile)
			else:
				fastaFiles.append(demux_input_dir + file)
		for path in fastaFiles:
			with open(path, "r") as f:
				file = path.split("/")[-1]
				sample = ".".join(file.split(".")[:-1])
				for line in f:
					if line.startswith(">"):
						newline = ">" + sample + "|" + line.split(">")[1]
						outfile.write(newline)
					else:
						outfile.write(line)
	reorientedSequences = reorientSeqs("%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix), "%s/%s" %(BLAST_DBs_folder,refseq_databasefile), num_threads) # Seqs need to have same orientation for OTU clustering (otherwise get 1 OTUs for each orientation of same sequence). DADA2 includes a reorientation step so not necessary here
	mvCMD = "cp %s %s/%s_1_bc_trimmed.fa " %(reorientedSequences, Output_folder, Output_prefix)
	process = subprocess.Popen(mvCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	process.communicate()
	# remove characters from sequence ids that will break later steps
	try:
		fasta_sequences = rename_fasta("%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix))
	except:
		sys.exit("ERROR: failed to rename %s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix))
	mvCMD = "cp %s/tmp/%s_1_bc_trimmed_renamed.fasta %s/%s_1_bc_trimmed.fa " %(Output_folder, Output_prefix, Output_folder, Output_prefix)
	process = subprocess.Popen(mvCMD, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
	process.communicate()
	count_total_input_sequences = count_seq_from_fasta("%s/%s_1_bc_trimmed.fa" %(Output_folder, Output_prefix))

## Remove primers ##
log.write('#Primer Removal# %s \n' % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
if "primerRemoval" in checkpoints_complete:
    log.write('Reusing previous primer trimmed files...\n')
    print('Reusing previous primer trimmed files...\n')
else:
    if mode != 3:
        if Clustering_method == "OTU" or Clustering_method == "BOTH":
            sys.stderr.write('Removing primers...\n')
            primer_trimmed_file = Output_prefix + '_2_pr_trimmed.fa'
            doCutAdapt(Fprims = Forward_primer, Rprims = Reverse_primer, InFile = Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa', OutFile = Output_folder + '/' + primer_trimmed_file)
            count_seq_pr_trimmed = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_2_pr_trimmed.fa')
            sys.stderr.write('\t' + str(count_seq_pr_trimmed) + ' sequences survived after primer-trimming\n')
            log.write('\t...done\n\n')
            if count_seq_pr_trimmed == str(0):
                sys.exit('Error: primer-removal returned no sequence')
            else:
                with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
                    open_checkpoint_file.write("primerRemoval\n")
        elif Clustering_method == "ASV":
            primer_trimmed_file = Output_prefix + '_1_bc_trimmed.fa'
            with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
                open_checkpoint_file.write("primerRemoval\n")
    else:
        primer_trimmed_file = Output_prefix + '_1_bc_trimmed.fa'
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("primerRemoval\n")

## Annotate the sequences with the taxon and locus names, based on the reference sequences ##
log.write('#Sequence Annotation# %s \n' % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
if "seqAnnotating" in checkpoints_complete:
    log.write('Reusing previous annotated files...\n')
    print('Reusing previous annotated files...\n')
    annoFileName = Output_prefix + '_3_annotated.fa'
    with open("%s/tmp/LocusTaxonCountDict_unclustd.pkl" % Output_folder, "rb") as picklefile:
        LocusTaxonCountDict_unclustd = pickle.load(picklefile)
else:
    if "primerRemoval" in checkpoints_complete:
        if Clustering_method == "OTU" or Clustering_method == "BOTH":
            primer_trimmed_file = Output_prefix + '_2_pr_trimmed.fa'
        elif Clustering_method == "ASV":
            primer_trimmed_file = Output_prefix + '_1_bc_trimmed.fa'
    sys.stderr.write('Annotating seqs...\n')
    toAnnotate = primer_trimmed_file
    annoFileName = Output_prefix + '_3_annotated.fa'
    if mode != 2:
        LocusTaxonCountDict_unclustd = annotateIt(filetoannotate = Output_folder + '/' + toAnnotate, outFile = Output_folder + '/' + annoFileName, Multiplex_perBC_flag = Multiplex_per_barcode, DualBC_flag = Dual_barcode, failsFile = Output_folder + '/' + Output_prefix + '_3_unclassifiable.fa')
    elif mode == 2:
        with open("%s/%s_2_pr_trimmed.fa" %(Output_folder, Output_prefix), "r") as infile, open("%s/%s" %(Output_folder, annoFileName), "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    newline = line.split("|")[0] + "|" + parameterDict['locus_list'][0] + "|" + line.split("|")[1]
                    outfile.write(newline)
                else:
                    outfile.write(line)
    count_seq_annotated = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_3_annotated.fa')
    count_seq_unclassifiable = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_3_unclassifiable.fa')
    sys.stderr.write('\t' + str(count_seq_annotated) + ' sequences annotated\n')
    sys.stderr.write('\t' + str(count_seq_unclassifiable) + ' sequences cannot be classified\n')
    log.write('\t...done\n\n')
    if count_seq_annotated == str(0):
        sys.exit('Error: annotation step returned no sequence')
    else:
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("seqAnnotating\n")
        if mode != 2:
            with open("%s/tmp/LocusTaxonCountDict_unclustd.pkl" % Output_folder, "wb") as picklefile:
                pickle.dump(LocusTaxonCountDict_unclustd, picklefile)

	#Does this block do anything? Do not see where mode would be set to 2 before I added the new mode to skip demultiplexing 2025-03-07 PWS
    #if mode == 2:
    #    log.write("PURC completed!\n")
    #    log.write("%s" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    #    print("PURC completed!")
    #    print("%s" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    #    sys.exit(0)

## Iterative clustering and chimera-killing ##
os.chdir(Output_folder) # move into the designated output folder

# force split annotating to FASTQ files if exiting early
if mode == 3:
	Clustering_method = "ASV"

if Clustering_method == "OTU":
    if "otuClustering" in checkpoints_complete:
        log.write('Reusing previous OTU clustering results...\n')
        print('Reusing previous OTU clustering results...\n')
        with open("%s/tmp/OTU_LocusTaxonCountDict_clustd.pkl" % Output_folder, "rb") as picklefile:
            LocusTaxonCountDict_clustd = pickle.load(picklefile)
        with open("%s/tmp/OTU_LocusTaxonCountDict_chimera.pkl" % Output_folder, "rb") as picklefile:
            LocusTaxonCountDict_chimera = pickle.load(picklefile)
    else:
        otuStartTime = time.time()
        LocusTaxonCountDict_clustd, LocusTaxonCountDict_chimera = IterativeClusterDechimera(annoFileName, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2)
        otuStopTime = time.time()
        otuRunTime = otuStopTime - otuStartTime
        print("OTU Runtime: %s" % convertTime(otuRunTime))
        log.write("OTU Runtime: %s\n" % convertTime(otuRunTime))
        log.write("OTU stop time: %s\n" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("otuClustering\n")
        with open("%s/tmp/OTU_LocusTaxonCountDict_clustd.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(LocusTaxonCountDict_clustd, picklefile)
        with open("%s/tmp/OTU_LocusTaxonCountDict_chimera.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(LocusTaxonCountDict_chimera, picklefile)
elif Clustering_method == "ASV":
    if "asvInference" in checkpoints_complete:
        log.write('Reusing previous ASV inference results...\n')
        print('Reusing previous ASV inference results...\n')
        with open("%s/tmp/ASV_LocusTaxonCountDict_clustd.pkl" % Output_folder, "rb") as picklefile:
            LocusTaxonCountDict_clustd = pickle.load(picklefile)
        with open("%s/tmp/ASV_LocusTaxonCountDict_chimera.pkl" % Output_folder, "rb") as picklefile:
            LocusTaxonCountDict_chimera = pickle.load(picklefile)
    else:
        asvStartTime = time.time()
        LocusTaxonCountDict_clustd, LocusTaxonCountDict_chimera = dada(annoFileName, fastq_sequences, Forward_primer, Reverse_primer, locus_list, minLen, maxLen, maxEE, RscriptPath)
        asvStopTime = time.time()
        asvRunTime = asvStopTime - asvStartTime
        print("ASV Runtime: %s" % convertTime(asvRunTime))
        log.write("ASV Runtime: %s\n" % convertTime(asvRunTime))
        log.write("ASV stop time: %s\n" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("asvInference\n")
        with open("%s/tmp/ASV_LocusTaxonCountDict_clustd.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(LocusTaxonCountDict_clustd, picklefile)
        with open("%s/tmp/ASV_LocusTaxonCountDict_chimera.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(LocusTaxonCountDict_chimera, picklefile)
elif Clustering_method == "BOTH" and useOTUpriors == "FALSE":
    # ASV
    if "asvInference" in checkpoints_complete:
        log.write('Reusing previous ASV inference results...\n')
        print('Reusing previous ASV inference results...\n')
        with open("%s/tmp/ASV_LocusTaxonCountDict_clustd.pkl" % Output_folder, "rb") as picklefile:
            ASV_LocusTaxonCountDict_clustd = pickle.load(picklefile)
        with open("%s/tmp/ASV_LocusTaxonCountDict_chimera.pkl" % Output_folder, "rb") as picklefile:
            ASV_LocusTaxonCountDict_chimera = pickle.load(picklefile)
    else:
        asvStartTime = time.time()
        ASV_LocusTaxonCountDict_clustd, ASV_LocusTaxonCountDict_chimera = dada(annoFileName, fastq_sequences, Forward_primer, Reverse_primer, locus_list, minLen, maxLen, maxEE, RscriptPath)
        asvStopTime = time.time()
        asvRunTime = asvStopTime - asvStartTime
        print("ASV Runtime: %s" % convertTime(asvRunTime))
        log.write("ASV Runtime: %s\n" % convertTime(asvRunTime))
        log.write("ASV stop time: %s\n" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("asvInference\n")
        with open("%s/tmp/ASV_LocusTaxonCountDict_clustd.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(ASV_LocusTaxonCountDict_clustd, picklefile)
        with open("%s/tmp/ASV_LocusTaxonCountDict_chimera.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(ASV_LocusTaxonCountDict_chimera, picklefile)
    # OTU
    if "otuClustering" in checkpoints_complete:
        log.write('Reusing previous OTU clustering results...\n')
        print('Reusing previous OTU clustering results...\n')
        with open("%s/tmp/OTU_LocusTaxonCountDict_clustd.pkl" % Output_folder, "rb") as picklefile:
            OTU_LocusTaxonCountDict_clustd = pickle.load(picklefile)
        with open("%s/tmp/OTU_LocusTaxonCountDict_chimera.pkl" % Output_folder, "rb") as picklefile:
            OTU_LocusTaxonCountDict_chimera = pickle.load(picklefile)
    else:
        otuStartTime = time.time()
        OTU_LocusTaxonCountDict_clustd, OTU_LocusTaxonCountDict_chimera = IterativeClusterDechimera(annoFileName, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2)
        otuStopTime = time.time()
        otuRunTime = otuStopTime - otuStartTime
        print("OTU Runtime: %s" % convertTime(otuRunTime))
        log.write("OTU Runtime: %s\n" % convertTime(otuRunTime))
        log.write("OTU stop time: %s\n" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("otuClustering\n")
        with open("%s/tmp/OTU_LocusTaxonCountDict_clustd.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(OTU_LocusTaxonCountDict_clustd, picklefile)
        with open("%s/tmp/OTU_LocusTaxonCountDict_chimera.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(OTU_LocusTaxonCountDict_chimera, picklefile)
    # Combine outputs
    LocusTaxonCountDict_clustd = {}
    for locus_taxon in OTU_LocusTaxonCountDict_clustd:
        LocusTaxonCountDict_clustd[locus_taxon] = {"OTU" : OTU_LocusTaxonCountDict_clustd[locus_taxon]}
    for locus_taxon in ASV_LocusTaxonCountDict_clustd:
        try:
            LocusTaxonCountDict_clustd[locus_taxon].update({"ASV" : ASV_LocusTaxonCountDict_clustd[locus_taxon]})
        except:
            LocusTaxonCountDict_clustd[locus_taxon] = {"ASV" : ASV_LocusTaxonCountDict_clustd[locus_taxon]}

elif Clustering_method == "BOTH" and useOTUpriors == "TRUE":
    # OTU
    if "otuClustering" in checkpoints_complete:
        log.write('Reusing previous OTU clustering results...\n')
        print('Reusing previous OTU clustering results...\n')
        with open("%s/tmp/OTU_LocusTaxonCountDict_clustd.pkl" % Output_folder, "rb") as picklefile:
            OTU_LocusTaxonCountDict_clustd = pickle.load(picklefile)
        with open("%s/tmp/OTU_LocusTaxonCountDict_chimera.pkl" % Output_folder, "rb") as picklefile:
            OTU_LocusTaxonCountDict_chimera = pickle.load(picklefile)
    else:
        otuStartTime = time.time()
        OTU_LocusTaxonCountDict_clustd, OTU_LocusTaxonCountDict_chimera = IterativeClusterDechimera(annoFileName, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2)
        otuStopTime = time.time()
        otuRunTime = otuStopTime - otuStartTime
        print("OTU Runtime: %s" % convertTime(otuRunTime))
        log.write("OTU Runtime: %s\n" % convertTime(otuRunTime))
        log.write("OTU stop time: %s\n" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("otuClustering\n")
        with open("%s/tmp/OTU_LocusTaxonCountDict_clustd.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(OTU_LocusTaxonCountDict_clustd, picklefile)
        with open("%s/tmp/OTU_LocusTaxonCountDict_chimera.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(OTU_LocusTaxonCountDict_chimera, picklefile)
    # ASV
    if "asvInference" in checkpoints_complete:
        log.write('Reusing previous ASV inference results...\n')
        print('Reusing previous ASV inference results...\n')
        with open("%s/tmp/ASV_LocusTaxonCountDict_clustd.pkl" % Output_folder, "rb") as picklefile:
            ASV_LocusTaxonCountDict_clustd = pickle.load(picklefile)
        with open("%s/tmp/ASV_LocusTaxonCountDict_chimera.pkl" % Output_folder, "rb") as picklefile:
            ASV_LocusTaxonCountDict_chimera = pickle.load(picklefile)
    else:
        asvStartTime = time.time()
        ASV_LocusTaxonCountDict_clustd, ASV_LocusTaxonCountDict_chimera = dada(annoFileName, fastq_sequences, Forward_primer, Reverse_primer, locus_list, minLen, maxLen, maxEE, RscriptPath)
        asvStopTime = time.time()
        asvRunTime = asvStopTime - asvStartTime
        print("ASV Runtime: %s" % convertTime(asvRunTime))
        log.write("ASV Runtime: %s\n" % convertTime(asvRunTime))
        log.write("ASV stop time: %s\n" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        with open("%s/tmp/checkpoint.txt" % Output_folder, "a") as open_checkpoint_file:
            open_checkpoint_file.write("asvInference\n")
        with open("%s/tmp/ASV_LocusTaxonCountDict_clustd.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(ASV_LocusTaxonCountDict_clustd, picklefile)
        with open("%s/tmp/ASV_LocusTaxonCountDict_chimera.pkl" % Output_folder, "wb") as picklefile:
            pickle.dump(ASV_LocusTaxonCountDict_chimera, picklefile)
    # Combine outputs
    LocusTaxonCountDict_clustd = {}
    for locus_taxon in OTU_LocusTaxonCountDict_clustd:
        LocusTaxonCountDict_clustd[locus_taxon] = {"OTU" : OTU_LocusTaxonCountDict_clustd[locus_taxon]}
    for locus_taxon in ASV_LocusTaxonCountDict_clustd:
        try:
            LocusTaxonCountDict_clustd[locus_taxon].update({"ASV" : ASV_LocusTaxonCountDict_clustd[locus_taxon]})
        except:
            LocusTaxonCountDict_clustd[locus_taxon] = {"ASV" : ASV_LocusTaxonCountDict_clustd[locus_taxon]}
if mode == 3:
    log.write("PURC stopped early after annotation")
    print("PURC completed!")
    print("Stop Time: %s" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    exit(0)
## Producing a summary ##
log.write('#Run Summary# %s \n\n' % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
count_output = open(Output_prefix + '_5_counts.xls', 'w')
count_output.write('Total input sequences:\t' + str(count_total_input_sequences) + '\n')
log.write('Total input sequences:\t' + str(count_total_input_sequences) + '\n')
if mode == 0:
    count_chimeric_sequences = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_0_chimeras.fa')
    count_output.write('Concatemers (multi-locus seqs):\t' + str(count_chimeric_sequences) + '\n')
    log.write('Concatemers (multi-locus seqs):\t' + str(count_chimeric_sequences) + '\n')

count_seq_w_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_bc_trimmed.fa')
count_seq_wo_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_trashBin_no_bc.fa')
count_seq_w_toomany_bc = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_1_trashBin_tooMany_bc.fa')
count_seq_annotated = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_3_annotated.fa')
count_seq_unclassifiable = count_seq_from_fasta(Output_folder + '/' + Output_prefix + '_3_unclassifiable.fa')

count_output.write('Sequences with barcodes:\t' + str(count_seq_w_bc) + '\n')
if not platform.system() == 'Linux' and Lima_override == "0":
	count_output.write('Sequences without barcodes:\t' + str(count_seq_wo_bc) + '\n')
	count_output.write('Sequences with too many barcodes:\t' + str(count_seq_w_toomany_bc) + '\n')
count_output.write('Sequences annotated:\t' + str(count_seq_annotated) + '\n')
count_output.write('Sequences that cannot be classified:\t' + str(count_seq_unclassifiable) + '\n')

log.write('Sequences with barcodes:\t' + str(count_seq_w_bc) + '\n')
if not platform.system() == 'Linux' and Lima_override == "0":
	log.write('Sequences without barcodes:\t' + str(count_seq_wo_bc) + '\n')
	log.write('Sequences with too many barcodes:\t' + str(count_seq_w_toomany_bc) + '\n')
log.write('Sequences annotated:\t' + str(count_seq_annotated) + '\n')
log.write('Sequences that cannot be classified:\t' + str(count_seq_unclassifiable) + '\n')

# Get lists of taxa and loci
taxon_list = []
locus_list = []
for taxon_locus in list(LocusTaxonCountDict_clustd.keys()):
    taxon_list.append(taxon_locus[0])
    locus_list.append(taxon_locus[1])
taxon_list = sorted(set(taxon_list))
locus_list = sorted(set(locus_list))

# Output read count per taxon per locus
log.write("\n**Raw reads per accession per locus**\n")
count_output.write('\n**Raw reads per accession per locus**\n')
count_output.write('\t' + '\t'.join(locus_list) + '\n')
log.write('\t' + '\t'.join(locus_list) + '\n')
for each_taxon in sorted(set(taxon_list)):
    count_output.write(each_taxon + '\t')
    log.write(each_taxon + '\t')
    for each_locus in sorted(locus_list):
        try:
            count_output.write(str(LocusTaxonCountDict_unclustd[each_taxon, each_locus]) + '\t')
            log.write(str(LocusTaxonCountDict_unclustd[each_taxon, each_locus]) + '\t')
        except:
            count_output.write('0\t')
            log.write('0\t')
    count_output.write('\n')
    log.write('\n')

# Output clustered seq count
if Clustering_method != "BOTH":
    count_output.write('\n**Final clustered sequences per accession per locus**\n')
    log.write('\n**Final clustered sequences per accession per locus**\n')
    count_output.write('\t' + '\t'.join(locus_list) + '\n')
    log.write('\t' + '\t'.join(locus_list) + '\n')
    copynumber_totals = {}
    for each_taxon in sorted(set(taxon_list)):
        count_output.write(each_taxon + '\t')
        log.write(each_taxon + '\t')
        for each_locus in sorted(locus_list):
            try:
                count = LocusTaxonCountDict_clustd[each_taxon, each_locus]
                count_output.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
                log.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
            except:
                count = 0
                count_output.write('0\t')
                log.write('0\t')
            try:
                copynumber_totals[each_locus] += count
            except:
                copynumber_totals[each_locus] = count
        count_output.write('\n')
        log.write('\n')
elif Clustering_method == "BOTH":
    count_output.write('\n**Final clustered sequences per accession per locus comparing OTUs and ASVs**\n')
    count_output.write('\t' + '\t\t'.join(locus_list) + '\n')
    count_output.write('\t' + "OTU\tASV\t"*len(locus_list) + '\n')
    log.write('\n**Final clustered sequences per accession per locus comparing OTUs and ASVs**\n')
    log.write('\t' + '\t\t'.join(locus_list) + '\n')
    log.write('\t' + "OTU\tASV\t"*len(locus_list) + '\n')
    copynumber_totals = {}
    for each_taxon in sorted(set(taxon_list)):
        count_output.write(each_taxon + '\t')
        log.write(each_taxon + '\t')
        for each_locus in sorted(locus_list):
            try:
                otu_count = LocusTaxonCountDict_clustd[each_taxon, each_locus]["OTU"]
                count_output.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]["OTU"]) + '\t')
                log.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]["OTU"]) + '\t')
            except:
                otu_count = 0
                count_output.write('0\t')
                log.write('0\t')
            try:
                asv_count = LocusTaxonCountDict_clustd[each_taxon, each_locus]["ASV"]
                count_output.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]["ASV"]) + '\t')
                log.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]["ASV"]) + '\t')
            except:
                asv_count = 0
                count_output.write('0\t')
                log.write('0\t')
            try:
                copynumber_totals[each_locus]["ASV"] += asv_count
            except:
                try:
                    copynumber_totals[each_locus].update({"ASV" : asv_count})
                except:
                    copynumber_totals[each_locus] = {"ASV" : asv_count}
            try:
                copynumber_totals[each_locus]["OTU"] += otu_count
            except:
                try:
                    copynumber_totals[each_locus].update({"OTU" : otu_count})
                except:
                    copynumber_totals[each_locus] = {"OTU" : otu_count}
        count_output.write('\n')
        log.write('\n')

# Output cluster count per locus
if Clustering_method != "BOTH":
    count_output.write('\n**Allele/copy/cluster/whatever count by locus**\n')
    log.write('\n**Allele/copy/cluster/whatever count by locus**\n')
    for each_locus in sorted(locus_list):
        if Clustering_method == "OTU":
            file_name = Output_prefix + '_4_' + str(each_locus) + '_OTUs.fa'
        elif Clustering_method == "ASV":
            file_name = Output_prefix + '_4_' + str(each_locus) + '_ASVs.fa'
        try: #I'm hoping that this will stop the program from crashing if a locus has no sequences
            seq_no = len(parse_fasta(file_name))
            count_output.write(str(each_locus) + '\t' + str(seq_no) + '\n')
            log.write(str(each_locus) + '\t' + str(seq_no) + '\n')
        except:
            count_output.write(str(each_locus) + '\t0\n')
            log.write(str(each_locus) + '\t0\n')
elif Clustering_method == "BOTH":
    count_output.write('\n**Allele/copy/cluster/whatever count by locus comparing OTUs and ASVs**\n')
    log.write('\n**Allele/copy/cluster/whatever count by locus comparing OTUs and ASVs**\n')
    count_output.write("\tOTU\tASV\n")
    log.write("\tOTU\tASV\n")
    for each_locus in sorted(locus_list):
        count_output.write("%s\t%s\t%s\n" % (each_locus, copynumber_totals[each_locus]["OTU"], copynumber_totals[each_locus]["ASV"]))

# Output chimeric seq count
count_output.write('\n**Chimeric clusters/sequence count by locus**\n')
for each_locus in sorted(locus_list):
    count_output.write(each_locus + '\n')
    #log.write(each_locus + '\t')
    for each_taxon in sorted(set(taxon_list)):
        try:
            if Clustering_method == "OTU":
                count_output.write('\t' + each_taxon + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][0]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][1]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][2]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][3]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][4]))
            elif Clustering_method == "ASV":
                count_output.write('\t' + each_taxon + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus]))
        except:
            count_output.write('\t' + each_taxon)

        count_output.write('\n')
count_output.close()


# Create Proportion List

otuFastas = glob.glob("*_OTUs.fa")
asvFastas = glob.glob("*_ASVs.fa")
otuLocusList = []
asvLocusList = []

sizeDict = {}
if len(otuFastas) >= 1:
    for file in otuFastas:
        locus = file.removeprefix("%s_4_" % Output_prefix).removesuffix("_OTUs.fa")
        otuLocusList.append(locus)
        with open(file, "r") as infile:
            for line in infile:
                if line.startswith(">"):
                    sample = line.strip(">\n").split("_Cluster")[0]
                    sizeTmp = line.strip(">\n").split(";")[1]
                    size = int(sizeTmp.split("=")[1])
                    try:
                        sizeDict[sample]["OTU"][locus].append(size)
                    except:
                        try:
                            sizeDict[sample]["OTU"].update({locus : [size]})
                        except:
                            try:
                                sizeDict[sample].update({"OTU" : {locus : [size]}})
                            except:
                                sizeDict.update({sample : {"OTU" : {locus : [size]}}})
if len(asvFastas) >= 1:
    for file in asvFastas:
        locus = file.removeprefix("%s_4_" % Output_prefix).removesuffix("_ASVs.fa")
        asvLocusList.append(locus)
        with open(file, "r") as infile:
            for line in infile:
                if line.startswith(">"):
                    sample = line.strip(">\n").split("_ASV")[0]
                    size = int(line.strip(">\n").split("_size=")[1])
                    try:
                        sizeDict[sample]["ASV"][locus].append(size)
                    except:
                        try:
                            sizeDict[sample]["ASV"].update({locus : [size]})
                        except:
                            try:
                                sizeDict[sample].update({"ASV" : {locus : [size]}})
                            except:
                                sizeDict.update({sample : {"ASV" : {locus : [size]}}})
    for locus in asvLocusList:
        for sample in sizeDict:
            if locus not in sizeDict[sample]["ASV"].keys():
                sizeDict[sample]["ASV"].update({locus : [0]})
with open("%s_5_proportions.tsv" % Output_prefix, "w") as outfile:
    otuLocusString = ""
    if len(otuLocusList) >= 1:
        otuLocusCounter = 1
        for otuLocus in sorted(otuLocusList):
            if otuLocusCounter < len(otuLocusList):
                otuLocusString = "%s%s-OTUs\t" %(otuLocusString, otuLocus)
            else:
                otuLocusString = "%s%s-OTUs" %(otuLocusString, otuLocus)
            otuLocusCounter += 1
    asvLocusString = ""
    if len(asvLocusList) >= 1:
        asvLocusCounter = 1
        for asvLocus in sorted(asvLocusList):
            if asvLocusCounter < len(asvLocusList):
                asvLocusString ="%s%s-ASVs\t" %(asvLocusString, asvLocus)
            else:
                asvLocusString ="%s%s-ASVs" %(asvLocusString, asvLocus)
            asvLocusCounter+=1
    if len(otuLocusList) >= 1 and len(asvLocusList) >= 1:
        outfile.write("Sample\t%s\t%s\n" %(otuLocusString, asvLocusString))
    elif len(otuLocusList) >= 1 and len(asvLocusList) == 0:
        outfile.write("Sample\t%s\n" %(otuLocusString))
    elif len(otuLocusList) == 0 and len(asvLocusList) >= 1:
        outfile.write("Sample\t%s\n" %(asvLocusString))
    else:
        print("WARNING: No OTUs or ASVs found in results")
    for sample in sorted(sizeDict.keys()):
        outfile.write("%s\t" % sample)
        try:
            for locus in sorted(otuLocusList):
                try:
                    counter = 1
                    for i in sizeDict[sample]['OTU'][locus]:
                        if counter > 1:
                            outfile.write(",")
                        outfile.write("%s" %i)
                        counter += 1
                    outfile.write("\t") #### TODO: Fix this section so trailing tabs are created when only OTUs or ASVs are present
                except:
                    outfile.write("0\t")
        except:
            outfile.write("\t" * len(otuLocusList))
        try:
            locusCounter = 0
            for locus in sorted(asvLocusList):
                locusCounter += 1
                try:
                    counter = 1
                    for i in sizeDict[sample]["ASV"][locus]:
                        if counter > 1:
                            outfile.write(",")
                        outfile.write("%s" %i)
                        counter += 1
                    if locusCounter < len(sizeDict[sample]["ASV"].keys()):
                        outfile.write("\t")
                except:
                    if locusCounter < len(sizeDict[sample]["ASV"].keys()):
                        outfile.write("0\t")
                    elif locusCounter == len(sizeDict[sample]["ASV"].keys()):
                        outfile.write("0")
        except:
            outfile.write("\n")
        outfile.write("\n")
# Make plot of proportions
with open("tmp/proportions.R", "w") as propRscript:
    propRscript.write('''
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
props <- read.delim("%s_5_proportions.tsv", sep = "\t", header = TRUE)
len <- dim(props)[1]
width <- dim(props)[2]
plotList <- list()
counter = 0
for (row in 1:len){
  counter = counter + 1
  sample <- props[row,1]
  titlePlot <- ggplot() +
    annotate("text", x = -1000,y = 1,size = 5,label = sample)+
    theme_void()
  plotList[[counter]] <- titlePlot
  for (col in c(2,3)){
    counter = counter + 1
    if (col == 2){
      clusterMethod <- "OTUs"
    } else if (col == 3) { clusterMethod <- "ASVs"}
    else if (col == 4) {clusterMethod <- "ASVs No Priors"}
    numbers <- props[row,col]
    for (j in strsplit(numbers,",")){
      newList <- c()
      for (k in j){
        newList <- c(newList,as.numeric(k))
        data <- data.frame(
          group=LETTERS[1:length(newList)],
          value=sort(newList, decreasing = TRUE)
          )
      }
    }
  # Get the positions
  df2 <- data %s>%s
      mutate(csum = rev(cumsum(rev(value))),
             pos = value/2 + lead(csum, 1),
             pos = if_else(is.na(pos), value/2, pos))
  if (counter < 5){
  pie <- ggplot(df2, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0, direction = -1) +
    theme_void() +
    theme(legend.position="none") +
    ggtitle(clusterMethod) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_label_repel(data = df2, aes(y = pos, label = value), size = 4.5, nudge_x = 1, show.legend = FALSE) +
    scale_fill_brewer(palette = "Pastel1")
  } else {
    pie <- ggplot(data, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=1) +
      coord_polar("y", start=0, direction = -1) +
      theme_void() +
      theme(legend.position="none") +
      ggtitle("") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_label_repel(data = df2, aes(y = pos, label = value), size = 4.5, nudge_x = 1, show.legend = FALSE) +
      scale_fill_brewer(palette = "Pastel1")
  }
  plotList[[counter]] <- pie
  }
}
relwidths = c(3)
for (i in 1:width){
  relwidths = c(relwidths,1)
}
plotsPerPage = width*6
pdf("%s_5_proportions.pdf", 8.5, 11)
for (i in seq(1, length(plotList), plotsPerPage)) {
  print(plot_grid(plotlist = plotList[i:(i+(plotsPerPage-1))], ncol = width, rel_widths = relwidths))
}
dev.off()
    ''' % (Output_prefix, "%", "%", Output_prefix))
propPlotCMD = "%s tmp/proportions.R" %(RscriptPath)
process = subprocess.Popen(propPlotCMD, stdout=log, stderr=log, shell=True, text=True)
process.communicate()


## Aligning the sequences ##
if Align == 1: # Aligning can be turned on/off in the configuration file
    fastas = glob.glob("*_OTUs.fa")
    for file in fastas:
        sys.stderr.write("Aligning " + file + "\n")
        log.write("Aligning " + file + "\n")
        outFile = mafftIt(file, verbose_level)
    fastas = glob.glob("*_ASVs.fa")
    for file in fastas:
        sys.stderr.write("Aligning " + file + "\n")
        log.write("Aligning " + file + "\n")
        outFile = mafftIt(file, verbose_level)
log.write("PURC completed!\nStop Time: %s" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
print("PURC completed!")
print("Stop Time: %s" % datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
