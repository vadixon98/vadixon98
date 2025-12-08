#!/usr/bin/env python3
"""
Secondary Structure Analysis from PDB Files
============================================
This script identifies alpha helices and beta sheets in a PDB file and calculates:
1. Distance between alpha carbons of first and last residues in each helix
2. Distance between alpha carbons of first and last residues in each sheet

Distances are calculated using the 3D Euclidean distance formula:
d = sqrt((x2-x1)² + (y2-y1)² + (z2-z1)²)

Author: Victoria Dixon
"""

import sys
import argparse
import math
from collections import defaultdict


def calculate_distance(coord1, coord2):
    """
    Calculate Euclidean distance between two 3D points.
    
    Args:
        coord1: tuple of (x, y, z) coordinates
        coord2: tuple of (x, y, z) coordinates
    
    Returns:
        Distance in Angstroms (Å)
    """
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return distance


def parse_pdb_atoms(pdb_file):
    """
    Parse ATOM records from PDB file, focusing on alpha carbons (CA).
    
    Returns:
        ca_coords: dict mapping (chain_id, resnum) -> (x, y, z) coordinates
    """
    ca_coords = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    atom_name = line[12:16].strip()
                    # Only process alpha carbon atoms
                    if atom_name == 'CA':
                        chain_id = line[21].strip()
                        resnum = int(line[22:26].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        
                        ca_coords[(chain_id, resnum)] = (x, y, z)
                except (ValueError, IndexError):
                    continue
    
    return ca_coords


def parse_helices(pdb_file):
    """
    Parse HELIX records from PDB file.
    
    Returns:
        helices: list of dicts with keys: chain_id, resseq1, resseq2, helix_class
    """
    helices = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HELIX'):
                try:
                    helix_id = line[7:10].strip()
                    chain_id = line[19].strip()
                    resseq1 = int(line[21:25].strip())
                    resseq2 = int(line[33:37].strip())
                    helix_class = line[38:40].strip()  # Class of helix (e.g., 1 = right-handed alpha)
                    
                    helices.append({
                        'helix_id': helix_id,
                        'chain_id': chain_id,
                        'resseq1': resseq1,
                        'resseq2': resseq2,
                        'helix_class': helix_class
                    })
                except (ValueError, IndexError):
                    continue
    
    return helices


def parse_sheets(pdb_file):
    """
    Parse SHEET records from PDB file.
    
    Returns:
        sheets: list of dicts with keys: sheet_id, chain_id, resseq1, resseq2, strand_id
    """
    sheets = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('SHEET'):
                try:
                    sheet_id = line[7:10].strip()
                    strand_id = line[11:14].strip()
                    chain_id = line[21].strip()
                    resseq1 = int(line[22:26].strip())
                    resseq2 = int(line[33:37].strip())
                    
                    sheets.append({
                        'sheet_id': sheet_id,
                        'strand_id': strand_id,
                        'chain_id': chain_id,
                        'resseq1': resseq1,
                        'resseq2': resseq2
                    })
                except (ValueError, IndexError):
                    continue
    
    return sheets


def analyze_secondary_structure(pdb_file):
    """
    Main function to analyze secondary structure in a PDB file.
    """
    print(f"Analyzing secondary structure in PDB file: {pdb_file}\n")
    print("=" * 70)
    
    # Parse alpha carbon coordinates
    ca_coords = parse_pdb_atoms(pdb_file)
    print(f"Found {len(ca_coords)} alpha carbon atoms\n")
    
    # Parse helices
    helices = parse_helices(pdb_file)
    print(f"Found {len(helices)} alpha helix(es)")
    
    # Parse sheets
    sheets = parse_sheets(pdb_file)
    print(f"Found {len(sheets)} beta sheet strand(s)\n")
    
    # Analyze helices
    if helices:
        print("ALPHA HELICES:")
        print("-" * 70)
        for i, helix in enumerate(helices, 1):
            chain_id = helix['chain_id']
            resseq1 = helix['resseq1']
            resseq2 = helix['resseq2']
            
            # Get coordinates for first and last residues
            coord1_key = (chain_id, resseq1)
            coord2_key = (chain_id, resseq2)
            
            if coord1_key in ca_coords and coord2_key in ca_coords:
                coord1 = ca_coords[coord1_key]
                coord2 = ca_coords[coord2_key]
                distance = calculate_distance(coord1, coord2)
                
                print(f"Helix {i} (ID: {helix['helix_id']}, Class: {helix['helix_class']}):")
                print(f"  Chain: {chain_id}")
                print(f"  Residue range: {resseq1} to {resseq2}")
                print(f"  Distance between CA atoms: {distance:.3f} Å")
                print()
            else:
                print(f"Helix {i} (ID: {helix['helix_id']}):")
                print(f"  Warning: Missing CA coordinates for residues {resseq1} or {resseq2}")
                print()
    else:
        print("ALPHA HELICES:")
        print("-" * 70)
        print("No alpha helices found in this structure.\n")
    
    # Analyze sheets
    if sheets:
        print("BETA SHEETS:")
        print("-" * 70)
        for i, sheet in enumerate(sheets, 1):
            chain_id = sheet['chain_id']
            resseq1 = sheet['resseq1']
            resseq2 = sheet['resseq2']
            
            # Get coordinates for first and last residues
            coord1_key = (chain_id, resseq1)
            coord2_key = (chain_id, resseq2)
            
            if coord1_key in ca_coords and coord2_key in ca_coords:
                coord1 = ca_coords[coord1_key]
                coord2 = ca_coords[coord2_key]
                distance = calculate_distance(coord1, coord2)
                
                print(f"Sheet Strand {i} (Sheet ID: {sheet['sheet_id']}, Strand ID: {sheet['strand_id']}):")
                print(f"  Chain: {chain_id}")
                print(f"  Residue range: {resseq1} to {resseq2}")
                print(f"  Distance between CA atoms: {distance:.3f} Å")
                print()
            else:
                print(f"Sheet Strand {i} (Sheet ID: {sheet['sheet_id']}):")
                print(f"  Warning: Missing CA coordinates for residues {resseq1} or {resseq2}")
                print()
    else:
        print("BETA SHEETS:")
        print("-" * 70)
        print("No beta sheets found in this structure.\n")
    
    # Summary message if no secondary structure
    if not helices and not sheets:
        print("=" * 70)
        print("This protein does not contain helices or sheets as defined in the PDB file.")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze secondary structure (helices and sheets) in a PDB file'
    )
    parser.add_argument(
        'pdb_file',
        type=str,
        help='Path to the PDB file'
    )
    
    args = parser.parse_args()
    
    try:
        analyze_secondary_structure(args.pdb_file)
    except FileNotFoundError:
        print(f"Error: File '{args.pdb_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

