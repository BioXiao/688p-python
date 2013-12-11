#!/bin/env python2
# -*- coding: utf-8 -*-
"""
CBMG 688P case study: 3'UTR motif analysis

This script brings together some of the concepts we have discussed in class
up to this point to complete a fairly complex and useful task: scanning a
genome to look for for interesting 3'UTR motifs 

Basic process:

    1) Define our goals for the script -- what do we want it to do? What
       should the inputs and outputs be?
    2) Break up task into modular tasks: e.g. "load genome", "compute 3'UTR
       coordinates", etc.
    3) Decide how to store things internally: using lists? dicts? classes? etc.
    4) Code each part of the process, testing along the way
    5) Validate output.

References
----------
- Najafabadi, H. S., Lu, Z., Macpherson, C., Mehta, V., Adoue, V., 
  Pastinen, T., & Salavati, R. (2013). Global identification of conserved 
  post-transcriptional regulatory programs in trypanosomatids. 
  Nucleic acids research, 1–10. doi:10.1093/nar/gkt647
"""
import os
import re
import sys
import csv
import StringIO
import urllib2
from Bio import SeqIO
from Bio import SeqUtils
from Bio.Alphabet import IUPAC

def main():
    """Main application body"""
    # Genome sequence and annotations
    genome = load_file('http://tritrypdb.org/common/downloads/release-6.0/TcruziCLBrenerEsmeraldo-like/fasta/data/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta')
    annotations = load_file('http://tritrypdb.org/common/downloads/release-6.0/TcruziCLBrenerEsmeraldo-like/gff/data/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like.gff')

    # 3'UTR motifs from supplementary table 2 in Najafabadi et al. (2013)
    motifs = load_motifs('najafabadi_table_s1_2013.csv')

    # Load genome sequence
    chromosomes = load_fasta(genome)

    # Prase annotations and return 3'UTR coordinates
    genes = get_utr_coords(annotations, utr_length=500)

    # For each gene, return a list of the motifs that are present in its 3'UTR
    for gene in genes[0:2]:
        utr_seq = get_3utr_seq(chromosomes, gene)

        # check each motif to see if it is present
        utr3_motifs = []

        for motif in motifs:
            matches = SeqUtils.nt_search(utr_seq, motif)[1:]

            # save matched motif
            if len(matches) > 0:
                utr3_motifs.append(motif)

        # output results
        print("%s: %s" % (gene['id'], ", ".join(utr3_motifs)))

def load_motifs(motif_file):
    """
    Returns a list of 3'UTR motifs.

    Arguments
    ---------
    motif_file: string
        Filepath to the supplementary table from Najafabadi et. al.

    Return
    ------
    out: list
        List of motifs from the table.
    """
    with open(motif_file, 'r') as fp:
        reader = csv.DictReader(fp)

        # Convert sequence to DNA and add to list
        motifs = [row['Sequence'].replace('U', 'T') for row in reader]

    return motifs

def get_utr_coords(filepath, utr_length=500):
    """
    Parses a GFF file and returns the estimated 3'UTR coordinates for each gene
    in the file.
    """
    # parse GFF file and get a list of the results
    gene_rows = load_gff(filepath)

    # For each gene from the above list, create a dictionary containing
    # the information needed to find the gene UTRs
    genes = []

    for row in gene_rows:
        # chromosome number
        match = re.search('\d+', row['seqid'])
        chromosome = int(match.group())

        # gene id
        match = re.search('ID=[^;]+;', row['attributes'])
        gene_id = match.group()[3:-1]

        # 3'UTR location
        # since the gene start and end ranges correspond to the CDS in this
        # case, we can simply find the N bases immediately following the CDS.
        if (row['strand'] == '+'):
            utr3_start = int(row['end'])
            utr3_end = int(row['end']) + utr_length
        else:
            utr3_start = int(row['start']) - 500
            utr3_end = int(row['start'])

        # create a dictionary representation of the gene
        genes.append({
            "id": gene_id,
            "chromosome": chromosome,
            "strand": row['strand'],
            "start": int(row['start']),
            "end": int(row['end']),
            "utr3_start": utr3_start,
            "utr3_end": utr3_end
        })

    return genes

def load_gff(filepath):
    """
    Loads a GFF file and returns a csv.DictReader instance corresponding
    to the gene rows in the file.
    """
    # GFF fields
    colnames = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand',
                'phase', 'attributes']

    # get lines from file
    with open(filepath, 'r') as fp:
        lines = fp.readlines()

    # filter out non-gene entries
    gene_rows = [ x for x in lines if 'gene\t' in x]

    # Next, let's create a StringIO buffer -- this is similar to the file
    # and url handles we have seen so far. We can then pass this to a csv
    # reader instance in the same way we have seen for actual files

    # First though, let's collapse the rows back into a single string
    csv_str = "".join(gene_rows)
    str_buffer = StringIO.StringIO(csv_str)

    return csv.DictReader(str_buffer, fieldnames=colnames, delimiter='\t')

def get_3utr_seq(chromosomes, gene):
    """
    Returns the sequence for the 3'UTR of a given gene
    """
    # get chromosome SeqRecord
    chromosome = chromosomes[gene['chromosome']]

    # get Seq for 3'UTR range
    seq = chromosome[gene['utr3_start'] - 1:gene['utr3_end']].seq

    # positive strand
    if gene['strand'] == '+':
        return str(seq)
    else:
        return str(seq.reverse_complement())

def load_fasta(filepath):
    """
    Loads a genome FASTA file and returns dictionary of chromosome sequences
    indexed by chromosome number.
    """
    chromosomes = {}

    seqs = SeqIO.parse(filepath, format='fasta', 
                       alphabet=IUPAC.ambiguous_dna)

    # iterate over seqs and add to chromosome dictionary
    for seq in seqs:
        # determine chromosome number
        match = re.search('\d+', seq.name)
        chromosome_number = int(match.group())

        chromosomes[chromosome_number] = seq

    return chromosomes

def load_file(uri):
    """
    Checks to see if a file exists in the current working directory and
    attempts to download otherwise.

    Arguments
    ---------
    uri : string
        A filepath or URI from which the file can be downloaded
    """
    # get filename
    filename = os.path.basename(uri)

    # If filepath specified and file exists, return filepath as-is
    if os.path.isfile(uri):
        return uri
    # Check to see if file exists in current working directory
    elif os.path.isfile(filename):
        return filename
    # Otherwise, check to see if URI is a URL
    elif uri.startswith(('http', 'ftp')):
        # retrieve remote file contents
        print("Downloading %s" % uri)
        handle = urllib2.urlopen(uri)
        contents = handle.read()

        with open(filename) as fp:
            fp.write(contents)
        return filename
    # if it's note a URL or a valid filepath, raise and exception
    else:
        raise Exception("Invalid URI specified: %s" % uri)

if __name__ == "__main__":
    sys.exit(main())

