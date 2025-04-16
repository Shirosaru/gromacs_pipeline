import os
from ImmuneBuilder import ABodyBuilder2
from Bio import SeqIO

def extract_hl_sequences(fasta_file):
    sequences = {'H': None, 'L': None}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Try to guess chain type from description or header
        desc = record.description.lower()
        if 'heavy' in desc or record.id.endswith('_H'):
            sequences['H'] = str(record.seq)
        elif 'light' in desc or record.id.endswith('_L'):
            sequences['L'] = str(record.seq)
    if not sequences['H'] or not sequences['L']:
        raise ValueError(f"Could not find both H and L chains in {fasta_file}")
    return sequences

def run_pipeline(input_dir, output_dir):
    predictor = ABodyBuilder2()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith(".fasta") or file_name.endswith(".fa"):
            input_path = os.path.join(input_dir, file_name)
            try:
                sequences = extract_hl_sequences(input_path)
                antibody = predictor.predict(sequences)

                output_file = os.path.join(
                    output_dir,
                    os.path.splitext(file_name)[0] + ".pdb"
                )
                antibody.save(output_file)
                print(f"Saved PDB: {output_file}")
            except Exception as e:
                print(f"Failed to process {file_name}: {e}")

# Example usage:
# run_pipeline("path/to/fasta_files", "path/to/output_pdbs")