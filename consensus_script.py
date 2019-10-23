#! /usr/bin/env python3

import argparse
import os
import subprocess as sb
import tempfile
import uuid
import xml.etree.ElementTree as ET
from io import StringIO

import requests
from Bio import SeqIO
from Bio import pairwise2

mlst_database = "https://pubmlst.org/data/dbases.xml"

MINIMAP = "minimap2 -a -x ava-ont -t %s %s %s "
MINIMAP_S = "minimap2 -a -x map-ont --secondary=no -t %s %s %s "
MINIMAP_SY = "minimap2 -a -x map-ont -t %s %s %s "
SAMSORT = "samtools sort -o %s -O BAM"
SAMINDEX = "samtools index %s"
RACON = "racon -t %s -f %s %s %s"
JELLYFISH_COUNT= "jellyfish count -m 100 -s 100M -t 10 -o /dev/stdout -C %s"
JELLYFISH_DUMP= "jellyfish dump -o %s /dev/stdin"
SPADES = "spades.py -s %s -o %s --only-assembler"
PORECHOP = "porechop -i %s -t %s -o %s"
NANOPOLISHV = "nanopolish variants --consensus -o %s -w %s -r %s -b %s -g %s -t %s --min-candidate-frequency 0.1 -p 1 " \
              "--fix-homopolymers"
NANOPOLISHI = "nanopolish index -d %s %s "
NANOPOLISHVA = "nanopolish vcf2fasta -g %s %s"
BWAI = "bwa index %s"
BWA = "bwa mem -t %s %s %s"
BCFTOOLS_MP = "bcftools mpileup -Ou -f %s %s"
BCFTOOLS_CALL = "bcftools call -mv -Oz -o %s --ploidy 1"
BCFTOOLS_IN = "bcftools index %s"
BCFTOOLS_CO = "bcftools consensus -f %s -o %s %s"
FREEBAYES = "freebayes -f %s -p 1 %s"
KROCUS_D = "krocus_database_downloader --species %s"


def setting():
    parser = argparse.ArgumentParser(prog='quantify and detect pathogens', usage='%(prog)s [options]')
    parser.add_argument("-o","--output", nargs="?", default="output")
    parser.add_argument("-f","--fastq", nargs="?", default="", required=True)
    parser.add_argument("-t", "--threads", nargs="?", default="1")
    parser.add_argument("-w", "--workdir", nargs="?", default="./")
    parser.add_argument("-s", "--species", nargs="?", default="")
    parser.add_argument("-f5", "--fast5", nargs="?", default="", required=True)
    args = parser.parse_args()
    return args


def racon():
    args = setting()
    cwd = os.getcwd()
    if args.species != "":
        req = requests.get(mlst_database)
        root = ET.fromstring(req.text)
        mlst_all = []
        for child in root:
            if args.species in child.text:
                mlst_comb = child.findall("mlst/database/profiles")
                mlst_location = child.findall("mlst/database/loci")
                for item in mlst_comb:
                    mlst_pattern_ftp =item.find('url').text
                for item in mlst_location:
                    mlst_ftps = item.findall('locus/url')
                    for mlst_ftp in mlst_ftps:
                        mlst_seq = requests.get(mlst_ftp.text).text
                        for record in SeqIO.parse(StringIO(mlst_seq), "fasta"):
                            mlst_all.append(record)

        mlst_complete = "/tmp/" + uuid.uuid4().hex
        with open (mlst_complete, "w") as fh:
            SeqIO.write(mlst_all, fh, "fasta")
        patternsST=[]
        mlst_pattern = requests.get(mlst_pattern_ftp).text.split("\n")
        for combination in mlst_pattern:
            if combination.startswith("ST"):
                orderST = combination.split("\t")[1:-1]
            else:
                patternsST.append([combination.split("\t")[0],"".join(combination.split("\t")[1:-1])])


        print("RUNNING MINIMAP")
        m = MINIMAP_SY % (args.threads, mlst_complete, args.fastq)
        print(m)
        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        out = minimap.communicate()[0].decode().split("\n")
        id_to_keep = []
        for read in out:
            fields = read.split("\t")
            if not read.startswith("@") and len(fields) > 18 and not fields[2].startswith("*") and not fields[9].startswith("*"):
                id_to_keep.append(fields[0])
        id_to_search = list(set(id_to_keep))

        records = [record for record in SeqIO.parse(args.fastq, "fastq") for name in id_to_search if name in record.id]

        read_mapping = "/tmp/" + uuid.uuid4().hex
        with open (read_mapping, "w") as fh:
            SeqIO.write(records, fh, "fastq")
    else:
        read_mapping=args.fastq

    print("RUNNING PORECHOP")
    reads_porechop = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
    porechop_cmd = PORECHOP % (read_mapping, args.threads, reads_porechop.name)
    porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
    porechop.communicate()
    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
    reads = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP % (args.threads, reads_porechop.name, reads_porechop.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    r = RACON % (args.threads, reads_porechop.name, sam.name, reads_porechop.name)
    print(r)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
    racon_cmd.communicate()
    sam = tempfile.NamedTemporaryFile(suffix=".sam")
    print("RUNNING MINIMAP")
    m = MINIMAP % (args.threads, reads.name, reads.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    r = RACON % (args.threads, reads.name, sam.name, reads.name)
    print(r)
    output = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    jfc_out = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
    racon_cmd.communicate()
    print("RUNNING JELLYFISH")
    jfc = JELLYFISH_COUNT % output.name
    jfd = JELLYFISH_DUMP % jfc_out.name
    print(jfc)
    print(jfd)
    jellyfishc = sb.Popen(jfc, shell=True, cwd=cwd, stdout=sb.PIPE)
    jellyfishd = sb.Popen(jfd, shell=True, cwd=cwd, stdin=jellyfishc.stdout)
    jellyfishd.communicate()
    count = 0
    kmer = "/tmp/" + uuid.uuid4().hex + ".fasta"
    with open(kmer, "w") as fh:
        for record in SeqIO.parse(jfc_out.name, "fasta"):
            if int(record.id) > 10 and len(record.seq) == 100:
                repetitive = 0
                while 10 >= repetitive:
                    repetitive += 1
                    count += 1
                    record.id = "kmer_" + str(count)
                    #print (record.id)
                    SeqIO.write(record, fh, "fasta")
    print(kmer)
    tmp_dir = tempfile.mkdtemp(dir="/tmp")
    print("RUNNING SPADES")
    spa = SPADES % (kmer, tmp_dir)
    print(spa)
    spade = sb.Popen(spa, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
    spade.communicate()
    assembled_contigs = os.path.join(tmp_dir, "contigs.fasta")
    bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP_S % (args.threads, assembled_contigs, read_mapping)
    ss = SAMSORT % bam.name
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
    samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=minimap.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
    samtools_sort.communicate()
    si = SAMINDEX % bam.name
    samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
    samtoos_index.communicate()
    ni = NANOPOLISHI % (args.fast5, read_mapping)
    print(ni)
    nanopolish_index = sb.Popen(ni, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
    nanopolish_index.communicate()
    regions = []
    print("RUNNING NANOPOLISH")
    for record in SeqIO.parse(assembled_contigs, "fasta"):
        region = record.id + ":0-" + str(len(record.seq))
        vcf = tempfile.NamedTemporaryFile(prefix="polished.", suffix=".vcf", delete=False)
        regions.append(vcf.name)
        nv = NANOPOLISHV % (vcf.name, region, read_mapping, bam.name, assembled_contigs, args.threads)
        print(nv)
        nanopolish_var = sb.Popen(nv, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        nanopolish_var.communicate()
    nva = NANOPOLISHVA % (assembled_contigs, " ".join(regions))
    print(nva)
    mlst = os.path.join(assembled_contigs + ".MLST.fasta")
    with open(mlst, "w") as fh:
        nanopolish_vcf = sb.Popen(nva, shell=True, cwd=cwd, stdout=fh, stderr=sb.PIPE)
    nanopolish_vcf.communicate()
    bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
    bi = BWAI % assembled_contigs
    bwa_index = sb.Popen(bi, shell=True, cwd="/tmp")
    bwa_index.communicate()
    bm = BWA % (args.threads, mlst, output.name)
    ss = SAMSORT % bam.name
    print(bm)
    bwa_mem = sb.Popen(bm, shell=True, cwd="/tmp", stdout=sb.PIPE, stderr=sb.PIPE)
    samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=bwa_mem.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
    samtools_sort.communicate()
    si = SAMINDEX % bam.name
    samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
    samtoos_index.communicate()
    vcf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
    bcfm = BCFTOOLS_MP % (mlst, bam.name) #FREEBAYES
    bcfc = BCFTOOLS_CALL % vcf.name
    print(bcfm)
    print(bcfc)
    bcftools_mpile = sb.Popen(bcfm, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
    bcftools_call = sb.Popen(bcfc, shell=True, cwd=cwd, stdin=bcftools_mpile.stdout)
    bcftools_call.communicate()
    bcfi = BCFTOOLS_IN % vcf.name
    print(bcfi)
    bcftools_index = sb.Popen(bcfi, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
    bcftools_index.communicate()
    mlst_done = os.path.join(cwd, args.fastq + ".MLST.done.fasta")
    bcfco = BCFTOOLS_CO % (mlst, mlst_done, vcf.name)
    print(bcfi)
    bcftools_call = sb.Popen(bcfco, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
    bcftools_call.communicate()
    if args.species != "":
        combinations = [(mlst, sample) for mlst in SeqIO.parse(mlst_complete, "fasta") for sample in SeqIO.parse(mlst_done, "fasta")]
        best_match = {}

        for pair in combinations:
            aln = pairwise2.align.globalms(pair[0].seq, pair[1].seq, 5, 0, 0, 0)[0]
            aln_rev = pairwise2.align.globalms(pair[0].seq, pair[1].reverse_complement().seq, 5, 0, 0, 0)[0]
            if aln_rev[2] > aln[2]:
                aln = aln_rev
            max_identity = len(pair[0].seq)*5
            identity = aln[2]/max_identity*100
            id_name = pair[0].id.split("_")[0]
            if id_name in best_match:
                if best_match[id_name][0][2] < aln[2]:
                    best_match[id_name] = [aln, pair[0].id, identity]
            else:
                best_match[id_name] = [aln, pair[0].id, identity]
        sample_st = []
        for st in orderST:
            sample_st.append(best_match[st][1].split("_")[1])
        for st in patternsST:
            if st[1] in "".join(sample_st):
                pattern = st[1]
        if not pattern == "":
            print(pattern)
        else:
            print("NO ST FOUND")

    print("DONE")




if __name__ == '__main__':
    racon()
