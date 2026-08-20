import pickle
from Bio import SeqIO

# Load PKL
with open('uniprot_arabidopsis.pkl', 'rb') as f:
    pkl_db = pickle.load(f)

# Load original FASTA for comparison
fasta_count = 0
fasta_ids = set()
for record in SeqIO.parse('uniprot_arabidopsis.fasta', 'fasta'):
    fasta_count += 1
    fasta_ids.add(record.id)
    # Also collect simple IDs
    if '|' in record.id:
        parts = record.id.split('|')
        if len(parts) >= 2:
            fasta_ids.add(parts[1])

print("=" * 60)
print("PKL FILE VERIFICATION")
print("=" * 60)
print(f"Sequences in FASTA: {fasta_count}")
print(f"Unique IDs in FASTA: {len(fasta_ids)}")
print(f"Entries in PKL: {len(pkl_db)}")
print()

# Check coverage
pkl_ids = set(pkl_db.keys())
missing_in_pkl = fasta_ids - pkl_ids
extra_in_pkl = pkl_ids - fasta_ids

print(f"IDs in FASTA but missing in PKL: {len(missing_in_pkl)}")
if missing_in_pkl:
    print("  Examples:", list(missing_in_pkl)[:5])

print(f"IDs in PKL but not in FASTA: {len(extra_in_pkl)}")
if extra_in_pkl:
    print("  Examples:", list(extra_in_pkl)[:5])

# Test with actual data
print("\n" + "=" * 60)
print("TESTING WITH YOUR DATASET")
print("=" * 60)

from pathlib import Path
data_file = Path('DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt')
total_pairs = 0
found_pairs = 0
missing_pairs = []

with open(data_file) as f:
    for i, line in enumerate(f):
        if i < 2 or line.startswith('#') or line.startswith('Protein'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        protA, protB, label = parts[0], parts[1], parts[2]
        total_pairs += 1
        
        if protA in pkl_db and protB in pkl_db:
            found_pairs += 1
        else:
            missing_pairs.append((protA, protB, 
                                 protA in pkl_db, protB in pkl_db))

print(f"Total protein pairs: {total_pairs}")
print(f"Pairs with BOTH sequences: {found_pairs} ({found_pairs/total_pairs*100:.1f}%)")
print(f"Pairs with MISSING sequences: {len(missing_pairs)} ({len(missing_pairs)/total_pairs*100:.1f}%)")

if missing_pairs:
    print("\nFirst 10 missing pairs:")
    for protA, protB, found_a, found_b in missing_pairs[:10]:
        print(f"  {protA} (found={found_a}) + {protB} (found={found_b})")

print("\n" + "=" * 60)