#!/bin/env python2
# -*- coding: utf-8 -*-
"""
CBMG 688P case study: 3'UTR motif analysis
Keith Hughitt
2013/12/11

Overview
--------
This script brings together some of the concepts we have discussed in class
up to this point to complete a fairly complex and useful task: scanning a
genome to look for for interesting 3'UTR motifs 

The script expects a genome sequence (as a FASTA file), gene annotations (GFF),
and a simple CSV file containing a list of short DNA or RNA motifs that should
be searched for. The lines in the GFF file corresponding to genes are extracted
and an approximate 3'UTR location is determine by grabbing the next N (where N
is something like 500) bases after the coding sequence of the gene. The 3'UTR
sequence for each gene is retrieved from the FASTA file, and then each scanned
for the presence of each motif in the input list. The matches for each gene
are then printed to the screen.

Currently, the script has been hard-coded to look for Genome and annotation
data related to an organism that we work on in my lab (http://en.wikipedia.org/wiki/Trypanosoma_cruzi).
This was done primarily for convenience since I am more familiar with this
data. The organism has been sequenced, and the ~65 Megabase haploid sequence 
along with some annotation data are available at http://tritrypdb.org/tritrypdb/.

The script assumes that the 3'UTR boundaries have not been well-defined, and
attempts to make a guess by taking a region of a specified length following
the each CDS. This, of course, is not perfect and is only meant for the
demonstrative purposes in this course. A better approach would be to attempt
to use existing 3'UTR information or make a more sophisticated guess at the
3'UTR location based on other properties of the surrounding sequence, etc.

Answers
-------
This script contains example implmentations for two of the homework tasks 
(command-line argument parsing and CSV output.) There are many different ways
these can be approached, however, so this is by no means the only or the best
solution.

Example Usage
-------------
python utr_motif_analysis_answers.py \
    -g http://tritrypdb.org/common/downloads/release-6.0/TcruziCLBrenerEsmeraldo-like/fasta/data/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \
    -a http://tritrypdb.org/common/downloads/release-6.0/TcruziCLBrenerEsmeraldo-like/gff/data/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like.gff \
    --utr-length=700 results2.csv

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
import argparse
from Bio import SeqIO
from Bio import SeqUtils
from Bio.Alphabet import IUPAC

def main():
    """Main application body"""
    # Parse command-line arguments
    args = parse_args()

    # Genome sequence and annotations
    genome = load_file(args.input_genome)
    annotations = load_file(args.input_annotations)

    # 3'UTR motifs from supplementary table 2 in Najafabadi et al. (2013)
    motifs = load_motifs('najafabadi_table_s1_2013.csv')

    # Load genome sequence
    chromosomes = load_fasta(genome)

    # Parse annotations and return 3'UTR coordinates
    genes = get_utr_coords(annotations, utr_length=args.utr_length)

    # Create a list to store output rows
    output = []

    # For each gene, return a list of the motifs that are present in its 3'UTR
    num_genes = len(genes)

    for i, gene in enumerate(genes):
        utr_seq = get_3utr_seq(chromosomes, gene)

        print('Processing gene %d/%d' % (i + 1, num_genes))

        # check each motif to see if it is present
        utr3_motifs = []

        for motif in motifs:
            matches = SeqUtils.nt_search(utr_seq, motif)[1:]

            # save matched motif
            if len(matches) > 0:
                utr3_motifs.append(motif)

        output.append([gene['id']] + utr3_motifs)

    # output results
    with open(args.output, 'w') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(output)

def parse_args():
    """Parses command-line arguments and returns a Namespace instance
    containing the user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--input-genome', required=True,
                        help='Input FASTA file to use')
    parser.add_argument('-a', '--input-annotations', required=True,
                        help='Input GFF annotation file to use')
    parser.add_argument('-u', '--utr-length', type=int, default=500,
                        help='Estimated size of the UTR (default: 500)')
    parser.add_argument('output', default='output.csv',
                        help='Filepath to write results to.')
    return parser.parse_args()

def get_utr_coords(filepath, utr_length=500):
    """
    Parses a GFF file and returns the estimated 3'UTR coordinates for each gene
    in the file.

    Parameters
    ----------
    filepath: string
        Location of a GFF file to process.
    utf_length: int
        Length of 3'UTR sequence to assume (Default: 500)

    Returns
    -------
    out: A list containing dictionary representations of each gene in the GFF
         file.
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

def get_3utr_seq(chromosomes, gene):
    """
    Returns the sequence for the 3'UTR of a given gene.

    Parameters
    ----------
    chromosomes: dict
        A dictionary containing SeqRecord entries for each chromosome in the
        input genome as outputted from the `load_fasta` function.
    gene: dict
        A dictionary representation of a single gene including the basic
        information required to determine the 3'UTR sequence.

    Returns
    -------
    out: string
        3'UTR sequence string.

    See Also
    --------
    * gene_utr_coords
    * load_fasta
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

def load_gff(filepath):
    """
    Loads a GFF file and returns a csv.DictReader instance corresponding
    to the gene rows in the file.

    Parameters
    ----------
    filepath: string
        Filepath to a GFF file to be processed.

    Returns
    -------
    out: csv.DictReader
        Returns a `DictReader` instance representing the gene-related fields
        of the input GFF file.

    See Also
    --------
    * http://www.sanger.ac.uk/resources/software/gff/spec.html
    * http://gmod.org/wiki/GFF
    * http://docs.python.org/2/library/stringio.html
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

def load_fasta(filepath):
    """
    Loads a genome FASTA file and returns dictionary of chromosome sequences
    indexed by chromosome number.

    Parameters
    ----------
    filepath: string
        Location of a FASTA genome file to be loaded. The FASTA file should
        contain one entry for each chromosome in the genome.

    Returns
    -------
    out: dict
        A dictionary of SeqRecord objects indexed by chromosome number.
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

def load_motifs(motif_file):
    """
    Returns a list of 3'UTR motifs.

    Parameters
    ----------
    motif_file: string
        Filepath to the supplementary table from Najafabadi et. al.

    Returns
    -------
    out: list
        List of motifs from the table.
    """
    with open(motif_file, 'r') as fp:
        reader = csv.DictReader(fp)

        # Convert sequence to DNA and add to list
        motifs = [row['Sequence'].replace('U', 'T') for row in reader]

    return motifs

def load_file(uri):
    """
    Checks to see if a file exists either at the specified location or in the
    the current working directory and attempts to download otherwise. The
    filepath to the matched file is then returned.

    Parameters
    ----------
    uri : string
        A filepath or URI from which the file can be downloaded

    Returns
    -------
    out : filepath
        Location of the requested file.
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

        with open(filename, 'w') as fp:
            fp.write(contents)
        return filename
    # if it's note a URL or a valid filepath, raise and exception
    else:
        raise Exception("Invalid URI specified: %s" % uri)

if __name__ == "__main__":
    sys.exit(main())

