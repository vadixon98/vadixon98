#!/usr/bin/env python3
"""
Loop Region Analysis from PDB Files
====================================
This script identifies loop regions in a PDB file and calculates:
1. Average B-factor (temperature factor) for atoms in each loop region
2. Residue names and sequence numbers for each loop region

Loops are identified as regions that are NOT part of helices or sheets
as defined in the HELIX and SHEET records of the PDB file.

Author: Victoria Dixon
"""

import sys
import argparse
from collections import defaultdict


def parse_pdb_secondary_structure(pdb_file):
    """
    Parse HELIX and SHEET records from PDB file to identify structured regions.
    
    Returns:
        structured_regions: set of (chain_id, resnum) tuples that are in helices or sheets
    """
    structured_regions = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            # Parse HELIX records (format: HELIX  serial chainID resSeq1 resSeq2)
            if line.startswith('HELIX'):
                try:
                    chain_id = line[19].strip()
                    resseq1 = int(line[21:25].strip())
                    resseq2 = int(line[33:37].strip())
                    # Add all residues in this helix
                    for resnum in range(resseq1, resseq2 + 1):
                        structured_regions.add((chain_id, resnum))
                except (ValueError, IndexError):
                    continue
            
            # Parse SHEET records (format: SHEET  serial chainID resSeq1 resSeq2)
            elif line.startswith('SHEET'):
                try:
                    chain_id = line[21].strip()
                    resseq1 = int(line[22:26].strip())
                    resseq2 = int(line[33:37].strip())
                    # Add all residues in this sheet
                    for resnum in range(resseq1, resseq2 + 1):
                        structured_regions.add((chain_id, resnum))
                except (ValueError, IndexError):
                    continue
    
    return structured_regions


def parse_pdb_atoms(pdb_file):
    """
    Parse ATOM records from PDB file.
    
    Returns:
        atoms_by_residue: dict mapping (chain_id, resnum) -> list of (resname, b_factor) tuples
    """
    atoms_by_residue = defaultdict(list)
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    chain_id = line[21].strip()
                    resname = line[17:20].strip()
                    resnum = int(line[22:26].strip())
                    b_factor = float(line[60:66].strip())
                    
                    atoms_by_residue[(chain_id, resnum)].append((resname, b_factor))
                except (ValueError, IndexError):
                    continue
    
    return atoms_by_residue


def identify_loops(structured_regions, atoms_by_residue):
    """
    Identify loop regions as residues that are not in helices or sheets.
    
    Returns:
        loops: list of lists, where each inner list contains (chain_id, resnum) tuples
               representing a contiguous loop region
    """
    # Get all residues that have atoms
    all_residues = sorted(set(atoms_by_residue.keys()))
    
    # Identify loop residues (not in structured regions)
    loop_residues = [res for res in all_residues if res not in structured_regions]
    
    if not loop_residues:
        return []
    
    # Group contiguous loop residues
    loops = []
    current_loop = [loop_residues[0]]
    
    for i in range(1, len(loop_residues)):
        prev_chain, prev_resnum = loop_residues[i-1]
        curr_chain, curr_resnum = loop_residues[i]
        
        # Check if residues are contiguous (same chain and consecutive)
        if prev_chain == curr_chain and curr_resnum == prev_resnum + 1:
            current_loop.append(loop_residues[i])
        else:
            # End of current loop, start new one
            if current_loop:
                loops.append(current_loop)
            current_loop = [loop_residues[i]]
    
    # Add the last loop
    if current_loop:
        loops.append(current_loop)
    
    return loops


def analyze_loops(pdb_file):
    """
    Main function to analyze loop regions in a PDB file.
    """
    print(f"Analyzing loops in PDB file: {pdb_file}\n")
    print("=" * 70)
    
    # Parse secondary structure
    structured_regions = parse_pdb_secondary_structure(pdb_file)
    print(f"Found {len(structured_regions)} residues in helices/sheets")
    
    # Parse atoms
    atoms_by_residue = parse_pdb_atoms(pdb_file)
    print(f"Found {len(atoms_by_residue)} residues with atoms\n")
    
    # Identify loops
    loops = identify_loops(structured_regions, atoms_by_residue)
    
    if not loops:
        print("No loop regions identified in this structure.")
        return
    
    print(f"Identified {len(loops)} loop region(s):\n")
    
    # Analyze each loop
    for i, loop in enumerate(loops, 1):
        print(f"Loop {i}:")
        print("-" * 70)
        
        # Collect all B-factors for atoms in this loop
        b_factors = []
        residue_info = []
        
        for chain_id, resnum in loop:
            if (chain_id, resnum) in atoms_by_residue:
                for resname, b_factor in atoms_by_residue[(chain_id, resnum)]:
                    b_factors.append(b_factor)
                    # Store residue info (avoid duplicates)
                    if (resname, chain_id, resnum) not in residue_info:
                        residue_info.append((resname, chain_id, resnum))
        
        # Calculate average B-factor
        if b_factors:
            avg_bfactor = sum(b_factors) / len(b_factors)
            print(f"  Average B-factor: {avg_bfactor:.2f} Å²")
        else:
            print(f"  Average B-factor: N/A (no atoms found)")
        
        # Print residue information
        print(f"  Residues ({len(residue_info)} residues):")
        residue_strings = []
        for resname, chain_id, resnum in sorted(residue_info, key=lambda x: (x[1], x[2])):
            residue_strings.append(f"{resname} {chain_id}:{resnum}")
        print(f"    {', '.join(residue_strings)}")
        print()


def main():
    parser = argparse.ArgumentParser(
        description='Analyze loop regions in a PDB file and calculate average B-factors'
    )
    parser.add_argument(
        'pdb_file',
        type=str,
        help='Path to the PDB file'
    )
    
    args = parser.parse_args()
    
    try:
        analyze_loops(args.pdb_file)
    except FileNotFoundError:
        print(f"Error: File '{args.pdb_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

