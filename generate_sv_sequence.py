import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
from warnings import warn
import re
import argparse
from pathlib import Path
import csv

# Functions to extract the sequence for structural variants of interest 
def test_key(dict_test, key_test):
    """
    Function to test if a key exists in a dict

    Parameters:
        - dict_test: dict. Dictionary to test
        - key_test: str, int, or float. Key to test
    Returns:
        - Empty string if the key does not exist or the input dictionary
    Raises:

    """
    try:
        tested_dict=dict_test[key_test]
        return tested_dict
    except KeyError:
        return ''

# Function to extract chromosme mate from the ALT field of a VcfRecord object
def extract_chr_mate(record):
    """
    Function to extract the position end of an SV based on the breakpoint mate

    Parameters:
        - record: VariantRecord object from Pysam. Contains information for a single variant from a structural variant VCF
    
    Returns:
        - str, int. The chromosome and position of the breakend mate
    Raises:
    """
    mate_pos=record.alts[0]
    chr_end=re.sub(r'(.*[\[\]])([A-z0-9]+)(:.+)', r'\2', mate_pos)
    pos_end=int(re.sub(r'(.*[\[\]][A-z0-9]+:)([0-9]+)(.*)', r'\2', mate_pos))

    return chr_end,pos_end

# Function to extract the start-end position of a variant
def extract_sv_startend(record, caller):
    """
    Function to extract the start and end position, reference and alternative alleles and SV type of a structural variant VCF from either Sniffles, Savana or Svim calls

    Parameters:
        - record: VariantRecord object from Pysam. Contains information for a single variant from a structural variant VCF
        - caller: Str. One of 'sniffles', 'svim' or 'savana'   
    Returns:
        - dict. In the following order, chromosome start, position start, position end, alt allele, reference allele, SV type and chromosome end.
    Raises:
        - Exception if the caller input is neither of 'sniffles', 'svim' or 'savana'

    """
    if caller == 'sniffles':
        return {'chr_start':record.contig,
                'startpos':record.pos,
                'endpos':record.stop,
                'alt':record.alts,
                'ref':record.ref,
                'svtype':record.info['SVTYPE'],
                'chr_end':test_key(record.info, 'CHR2')}
    if caller == 'svim':
        svtype=record.info['SVTYPE']
        # Simplify duplication types
        if bool(re.search('DUP', svtype)):
            svtype='DUP'
        # Extract destination chromosome and position for BNDs
        if svtype == 'BND':
            chr_end,endpos=extract_chr_mate(record)
            altallele=record.alts
        elif svtype == 'INS':
        # SVim does not output inserted sequence, return a string of SVLEN Ns instead
            endpos=record.stop
            chr_end=''
            altallele='N' * record.info['SVLEN']
        else:
            endpos=record.stop
            chr_end=''
            altallele=record.alts
        return {'chr_start':record.contig,
                'startpos':record.pos,
                'endpos':endpos,
                'alt':altallele,
                'ref':record.ref,
                'svtype':svtype,
                'chr_end':chr_end}
    if caller == 'savana':
        # Infer SV type from the strandedness of the breakends (see https://github.com/cortes-ciriano-lab/savana?tab=readme-ov-file#note-on-sv-types)
        bp_notations={'+-':'DEL', '-+':'DUP', '++':'INV', '--':'INV'}
        svtype=record.info['SVTYPE']
        if svtype == 'INS':
            altallele='N' * record.info['SVLEN']
            chr_end=record.contig
            endpos=record.pos + record.info['SVLEN']
        else:
            svtype=bp_notations[record.info['BP_NOTATION'][0]]
            # Extract end positions from mate
            chr_end,endpos=extract_chr_mate(record)
            altallele=record.alts
        return {'chr_start':record.contig,
                'startpos':record.pos,
                'endpos':endpos,
                'alt':altallele,
                'ref':record.ref,
                'svtype':svtype,
                'chr_end':chr_end}
    else:
        raise Exception('Parameter "caller" needs to be one of "sniffles", "svim" or "savana"')
    
def get_ref_sequence(chr, startpos, final_size, reference):
    """
    Generate a 'reference' sequence for a given SV breakpoint. The aim is to extend the flanks of the reference with the start of the breakpoint in the 'middle' of the sequence. In cases where the SV is too close to the start of the chromosome, the reference will be extended to the left of the breakpoint, if the SV is close to the end of the chromosome, the reference will be extended to the right. In very rare cases where an insertion/deletion occurs next to the end of the chromosome, the reference will also be extended to the right until reaching the desired length

    Parameters:
        - chr: str. Chromosome in starting position of the SV breakpoint
        - startpos: Int. Starting position of the SV breakpoint
        - finale size. Int. TOTAL length of the target sequence 
        - reference. SeqIO dict. Reference genome to use for sequence extraction
        
    Returns:
        - Seq object of the desired size surrounding the target SV
    Raises:
        Exception. If the desired size cannot be reached

    """

    ref_padding = final_size // 2
    # Ensure startpos does not go below 0
    padded_startpos = max((startpos - ref_padding), 0)
    ref_startseq = reference[chr][padded_startpos:startpos].seq
    # If the SV is not near the ends of a chromosome, this will effectively be final_size // 2    
    length_left = final_size - len(ref_startseq)
    # Get reference sequence
    ref_seq = ref_startseq + reference[chr][startpos:(startpos + length_left)].seq
    # The only cases where the reference still won't have the desired length is if the SV happens to be a deletion/insertion near the end of the chromosome (and not including the deletion would then mean the position goes beyond the chromosome length). The if below deals with these very rare cases by extending further to the right of the breakpoint.
    if len(ref_seq) < final_size:
        length_still_left = final_size - len(ref_seq)
        ref_seq = reference[chr][(padded_startpos - length_still_left):padded_startpos].seq + ref_seq
    
    # If length is still shorter than required, check whether the total chromosome is smaller than the padding + SV
    # If this is the case, output a warning
    if len(ref_seq) < final_size and len(reference[chr]) < final_size:
        warn('The SV+padding is larger than the reference chromosome, a shorter sequence has been output')
        return ref_seq

    # If length is still different, something else has gone wrong
    if len(ref_seq) != final_size:
        raise Exception('Desired length for the reference has not been reached, something has gone wrong')

    return ref_seq


def get_sv_sequence(pos_dict, reference, sequence_size):
    """
    Function to extract the sequence surrounding and SV

    Parameters:
        - pos_dict: dict. Dictionary with position information on an SV. {chr_start:chromosome start,startpos:position start,endpos:position end,alt:alt allele,ref:reference allele,svtype:SV type,chr_end:chromosome end}
        - reference. SeqIO dict. Reference genome to use for sequence extraction
        - sequence_size. Int. TOTAL length of the target sequence 
    Returns:
        - Seq object of the desired size (unless the positions are beyonde the chromosome) surrounding the target SV
        - Seq object of the reference sequence around the breakpoint of the SV without the SV. This sequence needs to be of the same size as the input, it will include slightly different flanking regions
    Raises:
        Warning. If the duplicated or inserted sequences are larger than the target length, no sequence will be produced
        Warning. If the svtype field in the input dict is not one of DEL, INS, DUP, INV or BND, no sequence will be produced
    """

    # Get the start and end positions of the variants (needs to be 0 based)
    chr_start = pos_dict['chr_start']
    pos_start = pos_dict['startpos'] - 1
    chr_end = pos_dict['chr_end']
    pos_end = pos_dict['endpos'] - 1

    if pos_dict['svtype'] == 'DEL':
        # For deletions the padding is just half of the final sequence length (no extra sequence is present)
        padding_size = sequence_size // 2
        endseq = reference[chr_start][pos_end:(pos_end + padding_size)].seq
        # Generate N reference sequence (ensure the resulting sequence does not go beyond the expected sequence size, this can happen if the deletion is very large)
        # If the resulting sequence is going to be larger than desired, output only up to desired size limit 
        deletion_size = pos_end - pos_start
        if deletion_size > padding_size:
            endseq_n = Seq('N' * padding_size)
        else:
            endseq_n = Seq('N' * deletion_size) + reference[chr_start][pos_end:(pos_end + (padding_size - deletion_size))].seq

    elif pos_dict['svtype'] == 'DUP':
        # Skip if duplicate is larger than target sequence size
        if (pos_end - pos_start) > sequence_size:
            warn(f'Duplicate {chr_start}:{pos_start} is larger than target size of {sequence_size}, no sequence generated')
            return '','',''
        else:
            # Padding size is half of the total desired size minus double the variant size
            padding_size = (sequence_size - ((pos_end - pos_start) * 2)) // 2
            dup_seq = reference[chr_start][pos_start:pos_end].seq
            endseq = dup_seq + reference[chr_start][pos_start:(pos_end + padding_size)].seq
            endseq_n = Seq('N' * (pos_end - pos_start)) + reference[chr_start][pos_start:(pos_end + padding_size)].seq

    elif pos_dict['svtype'] == 'INV':
        # Skip if inverted sequence is larger than target sequence size
        if (pos_end - pos_start) > sequence_size:
            warn(f'Inversion {chr_start}:{pos_start} is larger than target size of {sequence_size}, no sequence generated')
            return '','',''
        else:
            # Padding size in inversions is simply half of the total size minus size of the variant 
            padding_size = (sequence_size - (pos_end - pos_start)) // 2
            inverted_seq = reference[chr_start][pos_start:pos_end].seq[::-1]
            endseq = inverted_seq + reference[chr_start][pos_end:(pos_end + padding_size)].seq
            endseq_n = Seq('N' * (pos_end - pos_start)) + reference[chr_start][pos_end:(pos_end + padding_size)].seq

    elif pos_dict['svtype'] == 'INS':
        # Skip if inserted sequence is larger than target sequence size
        inserted_seq=Seq(pos_dict['alt'][0])
        if  len(inserted_seq) > sequence_size:
            warn(f'Insertion {chr_start}:{pos_start} is larger than target size of {sequence_size}, no sequence generated')
            return '','',''
        else:
            # Padding size in insertions is half of the total size minus size of the inserted sequence
            padding_size = (sequence_size - len(inserted_seq)) // 2
            endseq = inserted_seq + reference[chr_start][pos_end:(pos_end + padding_size)].seq
            endseq_n = Seq(len(inserted_seq) * 'N') + reference[chr_start][pos_end:(pos_end + padding_size)].seq

    elif pos_dict['svtype'] == 'BND':
        # NOTE, BNDs are treated very simplistically, simply extending either side of the breakpoints
        # This is likely to mis-represent complex SVs, but those should be resolved before arriving to this script
        # Also NOTE, if the VCF contains both breakend mates, they will get a sequence in opposite directions (extended_seq-mate1:mate2-extended_seq and extended_seq-mate2:mate1-extended_seq), to avoid this behaviour, the VCF needs to be filtered beforehand
        padding_size = sequence_size // 2
        endseq = reference[chr_end][pos_end:(pos_end + sequence_size // 2)].seq
        endseq_n = Seq((sequence_size // 2) * 'N')

    else:
        warn('The SVTYPE provided did not match any of DEL, DUP, INS, INV, INS nor BND, skipping')
        return '','',''
    
    # Get SV sequence (ensure start position does not go below 0, by default no positions beyond the max size of the chr will be extracted, so no adjustment needed on the other end)
    startpos_padded = max((pos_start - padding_size), 0)
    sv_seq = reference[chr_start][startpos_padded:pos_start].seq + endseq

    # Ensure reference sequence is as long as the SV sequence
    ref_seq = get_ref_sequence(chr_start, pos_start, len(sv_seq), reference)

    # Get N reference
    n_seq = reference[chr_start][startpos_padded:pos_start].seq + endseq_n

    return sv_seq,ref_seq,n_seq

def main():
    parser = argparse.ArgumentParser(description='Script to call the sequence in the alternative allele of SVs in an input VCF')
    parser.add_argument('-i', '--input_vcf', help = 'Input structural variant VCF', required = True)
    parser.add_argument('-c', '--caller', help = 'Caller used for generating the VCF, one of sniffles, svim or savana', required = True)
    parser.add_argument('-r', '--reference', help = 'Reference genome used in SV calling, fasta format', required = True)
    parser.add_argument('-s', '--sequence_size', help = 'Size of the resulting sequence with the SV', required = True)
    parser.add_argument('-o', '--out_dir', default='.', help = 'Output directory')

    args = parser.parse_args()

    vcf = pysam.VariantFile(args.input_vcf, 'r')
    # Load reference FASTA
    ref_fa = SeqIO.to_dict(SeqIO.parse(args.reference, 'fasta'))
    
    # Output
    out_filename = Path(args.out_dir, 'sv_sequences.csv')

    # Set header writer to True for the first iteration
    write_header = True

    with open(out_filename, mode="w", newline="") as f:
        writer = None
        for record in vcf.fetch():
            positions = extract_sv_startend(record, args.caller)
            sv_sequence,ref_sequence,n_sequence = get_sv_sequence(positions, ref_fa, int(args.sequence_size))
            this_record_dict = {'sv_id':record.id,'chr':record.contig,'pos':record.pos,'caller':args.caller,'sv_seq':str(sv_sequence),'ref_seq':str(ref_sequence),'n_seq':str(n_sequence)}

            # Initialise writer on the first iteration
            if writer is None:  
                writer = csv.DictWriter(f, fieldnames=this_record_dict.keys())
            if write_header:
                writer.writeheader()
                write_header = False

            writer.writerow(this_record_dict)

if __name__ == "__main__":
        main()