import sys
import os
import math
from pathlib import Path
import argparse
import sassy 
from Bio import SeqIO

# Edited from https://github.com/rickbeeloo/barbell-evals/blob/master/assemblies/scripts/contam.py#L536
# But here we do not polish the assembly, and we just focus to find contaminants.

SEARCHER = sassy.Searcher("iupac", alpha = 0.5)

def calculate_mismatch(seq_record):                                                                             
    """
    Calculate allowed mismatches in the contaminant sequence to be considered as contaminant.                                      
    Get the length of the barcode sequence and then use it to plug in a formula in the barbell paper.
    """                                                                                                                                                                                                
    length = len(seq_record)
    mask_length = 24 # native and rapid barcode sequences have 24 bp length
    effective_length = length - mask_length
    theoretical_errors = 0.51 * effective_length
    penalty_constant = 1.7312
    penalty = penalty_constant * math.sqrt(effective_length)
    tetha = math.ceil(theoretical_errors - penalty)
    return max(0, tetha)  

def _collapse_overlap(matches):
    # If two matches overlap, based on text_start, and text_end, choose the one with lowest 'cost'
    matches = sorted(matches, key=lambda x: x[0].text_start)
    i = 0  
    while i < len(matches) - 1:  
        if matches[i][0].text_end > matches[i+1][0].text_start:
            if matches[i][0].cost < matches[i+1][0].cost:
                matches.pop(i+1) 
            else:  
                matches.pop(i) 
        else:  
            i += 1  
    return matches  

def _get_matches(seq, barcode_seqs): 
    seq = seq.encode('utf-8')  
    matches = [] 
    
    # Unpack k directly from the tuple
    for (barcode, barcode_seq, k) in barcode_seqs:
        # Use the pre-calculated k
        search_results = SEARCHER.search(barcode_seq, seq, k=k)
        if len(search_results) > 0:
            for m in search_results: 
                matches.append((m, barcode))  
    return _collapse_overlap(matches) 
 
def search_contam(assembly_fasta, barcode_seqs):
    contam_matches = []  
  
    for record in SeqIO.parse(assembly_fasta, "fasta"):
        matches = _get_matches(str(record.seq), barcode_seqs)
        for m, barcode in matches: 
            start = m.text_start
            end = m.text_end  
            cost = m.cost
            contig = record.id  
            contam_matches.append((barcode, contig, start, end, cost))
    return contam_matches      

def main():
    parser = argparse.ArgumentParser(description="Detect adapter contamination using sassy and auto-calculated mismatches.")
    parser.add_argument("-a", "--assembly", required=True, help="Exact path to the assembly FASTA file")
    parser.add_argument("-b", "--barcodes", required=True, help="Exact path to the barcode FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Exact path for the output TSV file")
    args = parser.parse_args()

    # Load barcode sequences, calculate k, and log it immediately
    barcode_seqs = []
    
    print(f"Searching for {args.barcodes} contamination in: {args.assembly}")
    print(f"Loading barcodes and calculating mismatch thresholds (k):")

    for record in SeqIO.parse(args.barcodes, "fasta"):
        seq_str = str(record.seq)
        
        # Calculate k once per barcode
        k = calculate_mismatch(seq_str.encode('utf-8'))
        
        # Print to your log 
        print(f"  - {record.id}: length={len(seq_str)}, allowed_mismatches={k}")
        
        # Append k to the tuple so we can use it later without global variables
        barcode_seqs.append((record.id, seq_str.encode('utf-8'), k))
 
    
    matches = search_contam(args.assembly, barcode_seqs)
    
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True) 
    assembly_label = Path(args.assembly).name
    
    with open(args.output, 'w') as f:
        for match in matches:
            f.write(f"{assembly_label}\t{match[0]}\t{match[1]}\t{match[2]}\t{match[3]}\t{match[4]}\t{args.barcodes}\n")


if __name__ == "__main__":
    main()
