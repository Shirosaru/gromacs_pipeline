import os
from ImmuneBuilder import ABodyBuilder2
import os
import pandas as pd
import tempfile
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd

# Extract H and L sequences from a FASTA file
def extract_hl_sequences(fasta_file):
    sequences = {'H': None, 'L': None}
    for record in SeqIO.parse(fasta_file, "fasta"):
        desc = record.description.lower()
        if 'heavy' in desc or record.id.endswith('_H'):
            sequences['H'] = str(record.seq)
        elif 'light' in desc or record.id.endswith('_L'):
            sequences['L'] = str(record.seq)
    if not sequences['H'] or not sequences['L']:
        raise ValueError(f"Could not find both H and L chains in {fasta_file}")
    return sequences

# Load CSV with sequences
def load_sequences(csv_path):
    if csv_path.endswith('.gz'):
        return pd.read_csv(csv_path, compression='gzip')
    else:
        return pd.read_csv(csv_path)

# Run immune model builder from CSV input
def run_pipeline(csv_path, output_dir):
    df = load_sequences(csv_path)
    predictor = ABodyBuilder2()

    os.makedirs(output_dir, exist_ok=True)

    for idx, row in df.iterrows():
        therapeutic = row['ID']
        heavy_seq = row['heavy_sequence']
        light_seq = row['light_sequence']

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as tmp_fasta:
            records = [
                SeqRecord(Seq(heavy_seq), id=f"{therapeutic}_H", description="heavy"),
                SeqRecord(Seq(light_seq), id=f"{therapeutic}_L", description="light")
            ]
            SeqIO.write(records, tmp_fasta, "fasta")
            tmp_fasta_path = tmp_fasta.name

        try:
            sequences = extract_hl_sequences(tmp_fasta_path)
            antibody = predictor.predict(sequences)
            output_file = os.path.join(output_dir, f"{therapeutic}.pdb")
            antibody.save(output_file)
            print(f"[✓] Saved PDB for {therapeutic}: {output_file}")
        except Exception as e:
            print(f"[✗] Failed to process {therapeutic}: {e}")
        finally:
            os.remove(tmp_fasta_path)

# === Example usage ===
csv_input = "/home2/Gromacs/gromacs_pipeline/fasta/paired_sequences.csv"
pdb_output_dir = "/home2/Gromacs/gromacs_pipeline/fasta/output_pdbs"

run_pipeline(csv_input, pdb_output_dir)