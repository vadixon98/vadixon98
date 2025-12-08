#!/usr/bin/env python3
"""
MSA Position-Specific Amino Acid Usage Analysis
================================================
This script analyzes amino acid usage at each position in a multiple sequence
alignment (MSA), grouped by a Group variable (e.g., genus or species).

Inputs:
1. Table of selected entries (CSV/TXT) with Group column
2. MSA output file (FASTA format)

Outputs:
- Amino acid usage at each position for each group (counts and percentages)

Author: Victoria Dixon
"""

import sys
import argparse
import csv
from collections import defaultdict, Counter
from Bio import SeqIO


def parse_entry_table(table_file):
    """
    Parse the table of selected entries with Group information.
    
    Expected format: CSV or TXT with at least two columns:
    - One column identifying the sequence (e.g., sequence ID, name)
    - A 'Group' column (or similar) with the group assignment
    
    Returns:
        sequence_to_group: dict mapping sequence identifier -> group name
    """
    sequence_to_group = {}
    
    # Try to detect delimiter
    with open(table_file, 'r') as f:
        sample = f.read(1024)
        f.seek(0)
        sniffer = csv.Sniffer()
        delimiter = sniffer.sniff(sample).delimiter
    
    with open(table_file, 'r') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        
        # Try to find the group column (case-insensitive)
        fieldnames = [f.lower() for f in reader.fieldnames]
        group_col = None
        id_col = None
        
        for i, field in enumerate(reader.fieldnames):
            if 'group' in field.lower():
                group_col = field
            if any(keyword in field.lower() for keyword in ['id', 'name', 'sequence', 'entry', 'accession']):
                id_col = field
        
        if not group_col:
            # Try common variations
            for field in reader.fieldnames:
                if any(keyword in field.lower() for keyword in ['genus', 'species', 'taxon', 'class']):
                    group_col = field
                    break
        
        if not group_col:
            raise ValueError("Could not find 'Group' column in table file. "
                           "Please ensure your table has a column named 'Group' (or similar).")
        
        if not id_col:
            # Use first column as ID
            id_col = reader.fieldnames[0]
        
        for row in reader:
            seq_id = row[id_col].strip()
            group = row[group_col].strip()
            if seq_id and group:
                sequence_to_group[seq_id] = group
    
    return sequence_to_group


def parse_msa(msa_file):
    """
    Parse MSA file (FASTA format).
    
    Returns:
        alignment: list of (sequence_id, sequence) tuples
        max_length: maximum sequence length in alignment
    """
    alignment = []
    max_length = 0
    
    for record in SeqIO.parse(msa_file, 'fasta'):
        seq_id = record.id
        # Also try description in case ID is in description
        seq_str = str(record.seq)
        alignment.append((seq_id, seq_str))
        max_length = max(max_length, len(seq_str))
    
    return alignment, max_length


def normalize_sequence_id(seq_id, sequence_to_group):
    """
    Try to match sequence ID from MSA with entries in the table.
    Handles various ID formats and partial matches.
    """
    # Direct match
    if seq_id in sequence_to_group:
        return seq_id
    
    # Try without version numbers (e.g., "NP_123456.1" -> "NP_123456")
    if '.' in seq_id:
        base_id = seq_id.split('.')[0]
        if base_id in sequence_to_group:
            return base_id
    
    # Try partial matches (check if any table ID is contained in seq_id or vice versa)
    for table_id in sequence_to_group.keys():
        if table_id in seq_id or seq_id in table_id:
            return table_id
    
    # Try case-insensitive match
    seq_id_lower = seq_id.lower()
    for table_id in sequence_to_group.keys():
        if table_id.lower() == seq_id_lower:
            return table_id
    
    return None


def analyze_amino_acid_usage(table_file, msa_file):
    """
    Main function to analyze amino acid usage by position and group.
    """
    print(f"Analyzing amino acid usage from:")
    print(f"  Table file: {table_file}")
    print(f"  MSA file: {msa_file}\n")
    print("=" * 70)
    
    # Parse inputs
    try:
        sequence_to_group = parse_entry_table(table_file)
        print(f"Loaded {len(sequence_to_group)} sequences with group assignments")
    except Exception as e:
        print(f"Error parsing table file: {e}", file=sys.stderr)
        sys.exit(1)
    
    try:
        alignment, max_length = parse_msa(msa_file)
        print(f"Loaded {len(alignment)} sequences from MSA (max length: {max_length})\n")
    except Exception as e:
        print(f"Error parsing MSA file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Group sequences by their group assignment
    group_sequences = defaultdict(list)
    unmatched = []
    
    for seq_id, sequence in alignment:
        matched_id = normalize_sequence_id(seq_id, sequence_to_group)
        if matched_id and matched_id in sequence_to_group:
            group = sequence_to_group[matched_id]
            group_sequences[group].append((seq_id, sequence))
        else:
            unmatched.append(seq_id)
    
    if unmatched:
        print(f"Warning: {len(unmatched)} sequences could not be matched to groups:")
        print(f"  {', '.join(unmatched[:10])}{'...' if len(unmatched) > 10 else ''}\n")
    
    if not group_sequences:
        print("Error: No sequences could be matched to groups. Please check your sequence IDs.")
        sys.exit(1)
    
    print(f"Matched sequences to {len(group_sequences)} group(s):")
    for group, seqs in group_sequences.items():
        print(f"  {group}: {len(seqs)} sequence(s)")
    print()
    
    # Analyze each position
    print("=" * 70)
    print("AMINO ACID USAGE BY POSITION AND GROUP")
    print("=" * 70)
    print()
    
    # Get all groups
    all_groups = sorted(group_sequences.keys())
    
    # Analyze each position in the alignment
    for pos in range(max_length):
        print(f"Position {pos + 1}:")
        print("-" * 70)
        
        # Count amino acids at this position for each group
        for group in all_groups:
            sequences = group_sequences[group]
            aa_at_position = []
            
            for seq_id, sequence in sequences:
                if pos < len(sequence):
                    aa = sequence[pos].upper()
                    # Treat gaps and ambiguous characters
                    if aa == '-' or aa == '.':
                        aa = 'gap'
                    aa_at_position.append(aa)
            
            if not aa_at_position:
                continue
            
            # Count amino acids
            aa_counts = Counter(aa_at_position)
            total = len(aa_at_position)
            
            print(f"  {group} ({total} sequence(s)):")
            
            # Sort by count (descending), then by amino acid name
            sorted_aas = sorted(aa_counts.items(), key=lambda x: (-x[1], x[0]))
            
            for aa, count in sorted_aas:
                percentage = (count / total) * 100
                if aa == 'gap':
                    print(f"    gap: {count}/{total} ({percentage:.1f}%)")
                else:
                    print(f"    {aa}: {count}/{total} ({percentage:.1f}%)")
        
        print()
    
    print("=" * 70)
    print("Analysis complete!")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze position-specific amino acid usage in MSA by group'
    )
    parser.add_argument(
        'table_file',
        type=str,
        help='Path to table file (CSV/TXT) with Group column and sequence identifiers'
    )
    parser.add_argument(
        'msa_file',
        type=str,
        help='Path to MSA file (FASTA format)'
    )
    
    args = parser.parse_args()
    
    try:
        analyze_amino_acid_usage(args.table_file, args.msa_file)
    except FileNotFoundError as e:
        print(f"Error: File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

