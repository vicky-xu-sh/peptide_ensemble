import gemmi


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

def build_chain_sequence(chain: gemmi.Chain) -> Tuple[str, float]:
    residues = list(chain.get_polymer())
    if not residues:
        return "", 0.0

    aa_like = 0
    seq_chars: List[str] = []
    for residue in residues:
        resname = residue.name.upper()
        if resname in AA3_TO_1:
            seq_chars.append(AA3_TO_1[resname])
            aa_like += 1
        else:
            seq_chars.append("X")

    aa_fraction = aa_like / max(1, len(residues))
    return "".join(seq_chars), aa_fraction


def get_chain_seq(chain):
    sequence = []
    for res in chain:
        # skip water 
        if res.name in ('HOH', 'WAT'):
            continue

        # look up residue one-letter code, remove spaces
        code = gemmi.find_tabulated_residue(res.name).one_letter_code.strip()
        
        if code and code.isupper():
            sequence.append(code)
        else:
            # print("Unknown residue:", res.name)
            sequence.append('X')
            
    return "".join(sequence)

def read_pdb_ids(csv_path: str) -> List[str]:
    """
    Reads PDB IDs from a CSV file. The CSV must have a header row with a 'pdb_id' column.
    """
    with open(csv_path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {csv_path}")
        field_map = {name.lower(): name for name in reader.fieldnames}
        if "pdb_id" not in field_map:
            raise ValueError(
                f"CSV must contain a 'pdb_id' column. Found: {reader.fieldnames}"
            )
        col = field_map["pdb_id"]
        ids = [row[col].strip() for row in reader if row.get(col)]
    return [pid for pid in ids if pid]