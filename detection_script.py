#! /usr/bin/env python3

import argparse
import gzip
import os
import subprocess as sb
from pathlib import Path

import matplotlib as mpl
import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa

if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import shutil
import sys
import re
import time
import datetime


refseq_file = "assembly_summary_refseq"
assembly_file = "assembly_summary_genbank.txt"


WGET = 'wget ftp://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/%s'
WGET_GENOME = "wget %s"
MINIMAP = "minimap2 -a --secondary=no --MD -x map-ont -t %s %s %s > %s"
MINIMAP_INDEX = "minimap2 -t %s -d %s %s"


types = [("kingdom", "k__"),
         ("phylum", "p__"),
         ("class", "c__"),
         ("order", "o__"),
         ("family", "f__"),
         ("genus", "g__"),
         ("species", "s__"),
         ("tribu", "t__")]

def setting():
    parser = argparse.ArgumentParser(prog='quantify and detect pathogens', usage='%(prog)s [options]')
    parser.add_argument("-hs","--host_specie", nargs="?", default="")
    parser.add_argument("-ps","--pathogens_species", nargs="?", default="")
    parser.add_argument("-o","--output", nargs="?", default="output")
    parser.add_argument("-f","--fastq", nargs="?", default="", required=True)
    parser.add_argument("-g", "--genomes", nargs="?", default="")
    parser.add_argument("-t", "--threads", nargs="?", default="1")
    parser.add_argument("-d", "--NCBIdatabase", nargs="?", default="refseq",
                        help="The parameter indicate from with database to downlaod the genomes. "
                             "It can have assembly or refseq. Default is refseq")
    parser.add_argument("-w", "--workdir", nargs="?", default="./")
    parser.add_argument("-u", "--update", action='store_true')
    args = parser.parse_args()
    return args


def analysis():
    args = setting()
    cwd = args.workdir #os.getcwd()
    ncbi = NCBITaxa()
    home = str(Path.home())
    pathogens = args.pathogens_species.split(",")
    file_combined_fastq = os.path.join(os.getcwd(), args.fastq)
    if not os.path.isfile(file_combined_fastq):
        fastq_files = [os.path.join(file_combined_fastq, f) for f in listdir(file_combined_fastq) if isfile(join(file_combined_fastq, f)) and f.endswith("fastq")]
        k = file_combined_fastq.rfind("/")
        file_combined_fastq = file_combined_fastq[:k] + ".fastq" + file_combined_fastq[k + 1:]
        with open(file_combined_fastq, 'wb') as wfd:
            for file in fastq_files:
                with open(file, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024 * 1024 * 10)

    reads_fastq = []
    if file_combined_fastq.endswith("fastq") or file_combined_fastq.endswith("fq"):
        for record in SeqIO.parse(file_combined_fastq, "fastq"):
            reads_fastq.append(str(record.id))
    elif file_combined_fastq.endswith("fasta") or file_combined_fastq.endswith("fa"):
        for record in SeqIO.parse(file_combined_fastq, "fasta"):
            reads_fastq.append(str(record.id))
    else:
        print("Not known reads file format")

    number_reads = len(reads_fastq)

    if args.host_specie == "" and args.pathogens_species == "":
        species = ""
    elif args.host_specie == "" and not args.pathogens_species == "":
        species = pathogens
    elif not args.host_specie == "" and args.pathogens_species == "":
        species = [args.host_specie]
    else:
        species = [args.host_specie] + pathogens

    species.sort()
    name_database = "_".join(species).replace(" ", "_")
    genome_db = os.path.join(cwd, name_database + ".fasta")
    genome_db_id = os.path.join(cwd, name_database + ".txt")
    all_genomes = False
    if "refseq" in args.NCBIdatabase:
        table_file = "assembly_summary_refseq.txt"
    if "assembly" in args.NCBIdatabase:
        all_genomes = True
        table_file = "assembly_summary_genbank.txt"
    if os.path.exists(os.path.join(cwd,table_file)):
        os.remove(os.path.join(cwd,table_file))
    cmd = WGET % table_file
    wget = sb.Popen(cmd, shell=True, stdout=sb.PIPE, stderr=sb.PIPE, cwd=cwd)
    wget.communicate()



    sys.stdout.write("### UPDATING THE DATABASE\n")
    # This part checks for a new version of the taxdump.tar.gz; the code looks for a new version every day
    ete = os.path.expanduser("~/.etetoolkit/taxa.sqlite.traverse.pkl")
    modified = os.path.getmtime(ete)
    modificationTime = time.strftime('%m', time.localtime(modified))
    today = datetime.date.today()
    month = today.strftime("%m")
    if modificationTime != month:
        ncbi.update_taxonomy_database()
    dict_species = {}

    # here we set if is an not know pathogen or we have and idea of which pathogen to investigate
    with open(os.path.join(cwd, table_file), "r") as fh:
        descendants_all = []
        for specie in species:
            name2taxid = ncbi.get_name_translator([specie])
            if args.host_specie in specie:
                plant = name2taxid[specie]
            for key in name2taxid[specie]:
                descendants = ncbi.get_descendant_taxa(key, collapse_subspecies=True)
                for sstaxa in descendants:
                    descendants_all.append(str(sstaxa))
        for line in fh:
            if not line.startswith("#"):
                if line.split("\t")[6] in descendants_all:# and "subsp" in line:
                    ssname = " ".join([line.split("\t")[7].split(" ")[0], line.split("\t")[7].split(" ")[1]])
                    tax = line.split("\t")[6]
                    ftp = line.split("\t")[19]
                    genome = ftp.split("/")[-1] + "_genomic.fna.gz"
                    ftp_genome = os.path.join(ftp, genome)
                    path_genome = os.path.join(cwd, genome)
                    #species_assembly = " ".join([line.split("\t")[7]].split(" ")[0], [line.split("\t")[7]].split(" ")[1])
                    if ssname in dict_species:
                        dict_species[ssname] = dict_species[ssname] + [(ftp_genome, path_genome, tax, genome, ssname)]
                    else:
                        dict_species[ssname] = [(ftp_genome, path_genome, tax, genome, ssname)]
    db_file = os.path.join(home, ".db_monica." + name_database)
    if all_genomes:
        print("DOWNLOADING MULTIPLE GENOMES FOR THE SAME SPECIES")
        genomes_select = [name for specie in dict_species for name in dict_species[specie]]
    else:
        print("DOWNLOADING ONE GENOME FOR SPECIES")
        genomes_select = [dict_species[specie][-1] for specie in dict_species]
    print("I WILL DOWNLOAD %s GENOMES" % str(len(genomes_select)))
    if not os.path.exists(db_file) or not os.path.exists(genome_db):
        with open(genome_db, "w") as output_handle, open(genome_db_id, "w") as output_handle_id:
            with open(db_file, "w") as fh:
                for names in genomes_select:
                    ftp_genome, path_genome, tax, genome, ssname = names
                    if genome.startswith("GC"):
                        genome_used = cwd + genome + "\n"
                        fh.write(genome_used)
                        if not os.path.exists(path_genome):
                            cmd = WGET_GENOME % ftp_genome
                            wget_gen = sb.Popen(cmd, shell=True, stdout=sb.PIPE, stderr=sb.PIPE, cwd=cwd)
                            wget_gen.communicate()
                        with gzip.open(path_genome, "rt") as handle:
                            print("PARSING " + genome + " GENOME")
                            for record in SeqIO.parse(handle, "fasta"):
                                record.id = tax + "_" + str(record.id)
                                record.description = genome.split(".")[0]
                                SeqIO.write(record, output_handle, "fasta")
                                output_handle_id.write(str(record.name) + "%" + str(record.description) + "\n")

    sys.stdout.write("### PREPARING FOR MAPPING\n")
    genome_to_contig = {}
    with open(genome_db_id, "r") as fhtxt:
        for record in fhtxt: #txt SeqIO.parse(genome_db, "fasta"):
            line = record.split("%")
            genome_to_contig[line[0]] = line[1].rsplit()
    genome_to_species= {}
    with open(os.path.join(cwd, table_file), "r") as fh:
        for line in fh:
            line = line.rstrip().split("\t")
            genome = line[0].split(".")[0]
            if len(line) > 9 and not line[0].startswith("#"):
                subspecies = line[7].split(" ")[:2]
                subspecie = "_".join(subspecies) #+ " " + line[8].split("=")[1:]
                tribu = "_".join(line[8].split("=")[1:])
                genome_to_species[genome] = subspecie + "-" + tribu
    sam_output = file_combined_fastq + ".sam"
    cmd = MINIMAP % (str(args.threads), genome_db, file_combined_fastq, sam_output)
    sys.stdout.write("RUNNING MINIMAP2\n")
    minimap = sb.Popen(cmd, shell=True, cwd=cwd)
    minimap.communicate()
    reads_dict = {}
    count = 0
    with open(sam_output) as fh:
        for sam in fh:
            if sam != "" and not sam.startswith("@"):
                fields = sam.split("\t")
                if not fields[2] == "*":
                    for entry in fields:
                        if entry.startswith("MD"):
                            md = entry.split(":")[-1]
                            mismatch = len(re.findall("[A-Z]", md))
                            match = sum([int(number) for number in re.sub('[A-Z]|\^', ',', md).split(",") if number != "" and number.isdigit()])
                            if match > 0:
                                if mismatch > 0:
                                    iden = (match - mismatch) / match * 100
                                    if fields[0] in reads_dict:
                                        if iden == reads_dict[fields[0]][0]:
                                            if reads_dict[fields[0]][1].startswith(fields[2].split("_")[0]):
                                                continue
                                            else:
                                                count += 1
                                                reads_dict.pop(fields[0], None)
                                        elif iden > reads_dict[fields[0]][0]:
                                            reads_dict[fields[0]] = (iden, fields[2], fields[0])
                                    else:
                                        reads_dict[fields[0]] = (iden, fields[2], fields[0])
                                else:
                                    iden = 100
                                    if fields[0] in reads_dict:
                                        if iden == reads_dict[fields[0]][0]:
                                            if reads_dict[fields[0]][1].startswith(fields[2].split("_")[0]):
                                                continue
                                            else:
                                                count += 1
                                                reads_dict.pop(fields[0], None)
                                        elif iden > reads_dict[fields[0]][0]:
                                            reads_dict[fields[0]] = (iden, fields[2], fields[0])
                                    else:
                                        reads_dict[fields[0]] = (iden, fields[2], fields[0])

    out_file = file_combined_fastq + ".reads.txt"
    with open(out_file, "w") as csv:
        for key in reads_dict:
            csv.write("\t".join([reads_dict[key][1], reads_dict[key][2]]) + " \n")
    print(count)
    count = {}
    number_reads_mapped = 0
    for read in reads_dict:
        match = reads_dict[read][1].split("_")
        if len(match) > 1:
            number_reads_mapped += 1
            if all_genomes:
                contig = match[1] #+ "_" + match[2]
            else:
                contig = match[1] + "_" + match[2]
            genome_map = genome_to_contig[contig]
            species_ss = genome_to_species[genome_map[0]]
            uniq_name = match[0] + "_" + species_ss
            if not uniq_name in count:
                count[uniq_name] = 1
            else:
                count[uniq_name] = count[uniq_name] + 1
    print("Name sample: " + file_combined_fastq)
    print("Number reads:" + str(number_reads))
    print("Number reads mapped:" + str(number_reads_mapped) + "\nPercentage of reads mapped:" + str(
        number_reads_mapped/number_reads * 100) + " %\n")
    header = []
    reads_mapped = []
    partial_tree = []
    for clade in types:
        header.append(clade[0])
        reads_mapped.append("")
    header.append("A")
    reads_mapped.append(str(number_reads-number_reads_mapped))
    total = [header] + [reads_mapped]
    tribu_dict = {}
    sorted_list = []
    for value in count:
        key = value.split("_")[0]
        if not str(key).startswith(str(plant[0])):
            sorted_list.append((value[1],(count[value]/number_reads_mapped*100)))
            lineage = ncbi.get_lineage(int(key))
            a = ncbi.get_rank(lineage)
            tribu = value.split("-")[1]
            tribu_dict["tribu"] = tribu
            tree = []
            for match in types:
                combination = [match[1]]
                if match[0] in tribu_dict:
                    combination.append("".join([tribu_dict[match[0]]]))
                else:
                    for tax in a:
                        if match[0].startswith(a[tax]) and match[0].endswith(a[tax]):
                            combination.append(ncbi.get_taxid_translator([int(tax)])[tax].replace(" ","_"))
                tree.append("".join(combination))
            tree.append(str(count[value]))
            partial_tree = partial_tree + [tree]
    partial_tree.sort()
    total = total + partial_tree
    out_file = file_combined_fastq + ".txt"
    with open(out_file, "w") as csv:
        for line in total:
            csv.write(",".join(line) + " \n")
    plot_circ(out_file, file_combined_fastq)
    print("done")


def plot_circ(out_file, file_combined_fastq):

    fig, ax = plt.subplots()
    size = 0.3

    aa = pd.read_csv(out_file)
        #'/nanopore/data/xylella_1/in/20190115_1528_MN24941_FAK46406_0417e368/demultiplexed_fast5/barcode08_guppy.fastq.txt')
    # aa = pd.read_csv('/nanopore/monica/monica/libs/species.csv')
    df = aa.dropna()  # Drop the "no hits" line

    # Do the summing to get the values for each layer
    def nested_pie(df):
        cols = df.columns.tolist()
        outd = {}
        gb = df.groupby(cols[0], sort=False).sum()
        names = gb.index.values
        values = gb.values
        outd[0] = {'names': names, 'values': values}
        for lev in range(1, 8):
            gb = df.groupby(cols[:(lev + 1)], sort=False).sum()
            names = gb.index.levels[lev][gb.index.labels[lev]].values
            values = gb.values
            outd[lev] = {'names': names, 'values': values}
        return outd

    a = ['r', 'g', 'b', 'k', 'y', 'm', 'c', 'g']  # colors=plt.style.library["bmh"]["axes.prop_cycle"]
    outd = nested_pie(df)
    ax.pie(outd[6]['values'].sum(axis=1), labels=outd[6]['names'], radius=1, colors=a,
           wedgeprops=dict(width=size, edgecolor='w'))

    ax.pie(outd[7]['values'].flatten(), labels=outd[7]['names'], radius=1 + size, colors=a,
           wedgeprops=dict(width=size, edgecolor='w'))

    plt.savefig(file_combined_fastq + '.png')
    return ()
    #
    # aa = pd.read_csv(out_file)
    # df = aa.dropna()  # Drop the "no hits" line
    # df['A'] = np.random.rand(len(df)) * 100 + 1
    #
    # # Do the summing to get the values for each layer
    # def nested_pie(df):
    #
    #     cols = df.columns.tolist()
    #     outd = {}
    #     gb = df.groupby(cols[0], sort=False).sum()
    #     name = gb.index.values
    #     value = gb.values
    #     outd[0] = {'names': name, 'values': value}
    #     for lev in range(1, 7):
    #         gb = df.groupby(cols[:(lev + 1)], sort=False).sum()
    #         outd[lev] = {'names': gb.index.levels[lev][gb.index.labels[lev]].tolist(), 'values': gb.values}
    #     return outd
    # #a=mp.Set1(np.arange(10)/10.)
    # a = ['r', 'g', 'b', 'k', 'y', 'm', 'c']  # colors=plt.style.library["bmh"]["axes.prop_cycle"]
    # outd = nested_pie(df)
    # #diff = 1 / 7.0
    #
    # # This first pie chart fill the plot, it's the lowest level
    # plt.pie(outd[6]['values'], labels=outd[6]['names'], labeldistance=0.9, colors=a)
    # ax = plt.gca()
    # # For each successive plot, change the max radius so that they overlay
    # for i in np.arange(5, -1, -1):
    #     ax.pie(outd[i]['values'], labels=outd[i]['names'], radius=np.float(i + 1) / 7.0,
    #            labeldistance=((2 * (i + 1) - 1) / 14.0) / ((i + 1) / 7.0),
    #            colors=a)
    # ax.set_aspect('equal')
    # return plt
    #



if __name__ == '__main__':
    analysis()

