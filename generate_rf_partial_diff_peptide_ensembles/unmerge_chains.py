import argparse
import sys
import os
from pathlib import Path
import mdtraj as md
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

BACKBONE_ATOMS = {"N", "CA", "C", "O"}

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
    "ASX": "B",
    "GLX": "Z",
}

def extract_chain_per_residue(original_pdb):
    """
    Extract chain ID and residue index per residue based on CA atoms.
    Returns list of tuples: (chain_id, resseq, icode, resname) in residue order.
    """
    residue_info = []
    seen = set()

    with open(original_pdb, "r") as f:
        
        for line in f:
            if not line.startswith(("ATOM")):
                continue

            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            
            resname = line[17:20].strip()
            if resname not in AA3_TO_1: # skip non-standard residues
                continue

            chain_id = line[21]
            resseq = line[22:26].strip()

            if (chain_id, resseq, resname) in seen:
                print(f"Warning: Duplicate residue found in original PDB: {(chain_id, resseq, resname)}. Skipping duplicate.")
                continue  # skip duplicate residues (e.g. due to alternate locations)
            seen.add((chain_id, resseq, resname))
            residue_info.append((chain_id, resseq, resname))

    return residue_info


def unmerge_backbone(original_pdb, merged_pdb, output_pdb):
    residue_info_list = extract_chain_per_residue(original_pdb)

    if not residue_info_list:
        raise ValueError("No CA atoms found in original PDB. Cannot extract residue information.")

    residue_index = -1
    current_residue_key = None

    with open(merged_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:

            if not line.startswith(("ATOM")):
                outfile.write(line)
                continue
                
            # Identify residue by (resSeq + resName)
            res_key = (line[22:26].strip(), line[17:20].strip())

            if res_key != current_residue_key:
                residue_index += 1
                current_residue_key = res_key
                current_chain, orig_resseq, orig_resname = residue_info_list[residue_index]

                # print(f"merged file res_key: {res_key}, \
                #         original file residue: {(current_chain, orig_resseq, orig_resname)}")

                # chain A is the peptide chain that's partial diffused, 
                # all res in that chain would be GLY in the merged file and won't match the original resname, so we skip the resname check for chain A
                if current_chain != 'A': 
                    assert orig_resname == res_key[1], f"Residue name mismatch: {orig_resname} vs {res_key[1]}"

                if residue_index >= len(residue_info_list):
                    raise ValueError(
                        "More residues in merged file than original file."
                    )

            # Replace chain ID (column 22)
            new_line = line[:21] + current_chain + line[22:]
            outfile.write(new_line)

    if residue_index + 1 != len(residue_info_list):
        print(
            f"residue_index: {residue_index+1}, "
            f"len(residue_info_list): {len(residue_info_list)}, "
            f"last_orig_residue: {residue_info_list[-1]}"
        )
        raise ValueError(
            "Original file has more residues than merged file."
        )

    print(f"Chains successfully restored → {output_pdb}")


def process_pdb_file(pdb_file, subdir, orig_pdb_dir, generated_pdbs_dir):
    """
    Process a single PDB file: unmerge chains and load with mdtraj.
    Returns the trajectory object.
    """
    subdir_path = os.path.join(generated_pdbs_dir, subdir)
    pdb_path = os.path.join(subdir_path, pdb_file)
    
    # Create temporary file for unmerged PDB
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
        temp_pdb_path = tmp.name
    
    try:
        # Call unmerge_backbone
        orig_pdb = os.path.join(orig_pdb_dir, f"{subdir}.pdb")
        unmerge_backbone(orig_pdb, pdb_path, temp_pdb_path)
        
        # Load the unmerged structure with mdtraj
        traj = md.load(temp_pdb_path)
        return traj
        
    finally:
        # Clean up temporary file
        if os.path.exists(temp_pdb_path):
            os.remove(temp_pdb_path)


def process_subdirectory(subdir, orig_pdb_dir, generated_pdbs_dir, output_dir):
    """
    Process all PDB files in a single subdirectory.
    Returns subdir name and success status.
    """
    subdir_path = os.path.join(generated_pdbs_dir, subdir)
    
    # Skip if not a directory
    if not os.path.isdir(subdir_path):
        return subdir, False
    
    print(f"\nProcessing subdirectory: {subdir}")
    
    # Create corresponding subdirectory in output
    output_subdir = os.path.join(output_dir, subdir)
    os.makedirs(output_subdir, exist_ok=True)
    
    # List to store unmerged structures
    unmerged_structures = []
    
    # Loop through all PDB files in the subdirectory
    pdb_files = sorted([f for f in os.listdir(subdir_path) if f.endswith('.pdb')])
    
    for pdb_file in pdb_files:
        try:
            traj = process_pdb_file(pdb_file, subdir, orig_pdb_dir, generated_pdbs_dir)
            unmerged_structures.append(traj)
        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")
    
    # Save combined structures as PDB + XTC
    if unmerged_structures:
        # Concatenate all trajectories
        combined_traj = unmerged_structures[0]
        for traj in unmerged_structures[1:]:
            combined_traj = combined_traj.join(traj)
        
        # Save to output subdirectory
        output_pdb = os.path.join(output_subdir, f"{subdir}.pdb")
        output_xtc = os.path.join(output_subdir, f"{subdir}.xtc")
        
        combined_traj[0].save(output_pdb)  # Save first frame as PDB
        combined_traj.save(output_xtc)     # Save trajectory as XTC
        
        print(f"Saved: {output_pdb}, {output_xtc}")
        return subdir, True
    
    return subdir, False


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--orig_pdb_dir", 
        type=str, 
        help="Directory containing the original PDB file input to RFdiffusion.")
    parser.add_argument(
        "--generated_pdbs_dir", 
        type=str, 
        help="Directory containing subdirectories with the merged PDB files.")
    parser.add_argument(
        "--output_dir", 
        type=str, 
        help="Directory to save the unmerged PDB and XTC files.")
    parser.add_argument(
        "--num_workers",
        type=int,
        default=None,
        help="Number of parallel workers. Default: number of CPU cores.")
    
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Test one
    # subdir = '2xrw'
    # print(f"Testing on subdirectory: {subdir}")
    # process_subdirectory(subdir, args.orig_pdb_dir, args.generated_pdbs_dir, args.output_dir)

    # Get all subdirectories
    subdirs = [
        subdir for subdir in sorted(os.listdir(args.generated_pdbs_dir))
        if os.path.isdir(os.path.join(args.generated_pdbs_dir, subdir))
    ]
    
    # Process each subdirectory in parallel
    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        futures = {
            executor.submit(
                process_subdirectory, subdir, args.orig_pdb_dir, 
                args.generated_pdbs_dir, args.output_dir
            ): subdir for subdir in subdirs
        }
        
        for future in as_completed(futures):
            subdir, success = future.result()
            if not success:
                print(f"Warning: No structures saved for {subdir}")