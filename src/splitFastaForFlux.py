import os
from Bio import SeqIO

# === USER SETTINGS ===
# TODO Fix this path
input_fasta = "/home/bio/Molecular/flux/cusannot.sorted.fa"  # <-- Change this to your FASTA file
output_dir = "./data/flux/transcripts/" # Folder to hold split files

# === Create output directory ===
os.makedirs(output_dir, exist_ok=True)

# === Split FASTA ===
count = 0
for record in SeqIO.parse(input_fasta, "fasta"):
    count += 1
    out_path = os.path.join(output_dir, f"transcript{count}.fa")
    with open(out_path, "w") as out_handle:
        SeqIO.write(record, out_handle, "fasta")
    if count % 1000 == 0:
        print(f"Written {count} sequences...")

print(f"Done! Total sequences split: {count}")
print(f"Output folder: {output_dir}")