# Run this on your HPC to merge UniParC sequences with existing database
# Place this in ~/AraPPLM/ directory

import pickle
import json
import os
from Bio import SeqIO

print("=" * 60)
print("MERGING UniParc SEQUENCES WITH EXISTING DATABASE")
print("=" * 60)

# 1. Load original database
print("\n1. Loading original UniProt database...")
with open('uniprot_arabidopsis.pkl', 'rb') as f:
    seq_db = pickle.load(f)
print(f"   Original entries: {len(seq_db)}")

# 2. Load ID mapping from JSON
print("\n2. Loading ID mapping...")
with open('uniparc_id_mapping.json', 'r') as f:
    id_mapping = json.load(f)
print(f"   Mappings loaded: {len(id_mapping)}")

# 3. Load UniParC FASTA sequences
print("\n3. Loading UniParC sequences...")
uniparc_sequences = {}
for record in SeqIO.parse('uniParC_348.fasta', 'fasta'):
    uniparc_id = record.id
    seq = str(record.seq)
    uniparc_sequences[uniparc_id] = seq
print(f"   UniParC sequences loaded: {len(uniparc_sequences)}")

# 4. Merge: Map original IDs to sequences via UniParC
print("\n4. Merging sequences...")
merged_count = 0
failed_mappings = []

for original_id, uniparc_id in id_mapping.items():
    if uniparc_id in uniparc_sequences:
        seq = uniparc_sequences[uniparc_id]
        seq_db[original_id] = seq  # Map original ID to sequence
        seq_db[uniparc_id] = seq   # Also store by UniParC ID
        merged_count += 1
    else:
        failed_mappings.append((original_id, uniparc_id))

print(f"   Successfully merged: {merged_count}")
if failed_mappings:
    print(f"   Failed mappings: {len(failed_mappings)}")
    print(f"   Examples: {failed_mappings[:5]}")

# 5. Test coverage
print("\n5. Testing coverage with DeepAraPPI dataset...")
with open('DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt') as f:
    total = 0
    found = 0
    missing_pairs = []
    
    for i, line in enumerate(f):
        if i < 2 or line.startswith('#') or line.startswith('Protein'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        total += 1
        protA, protB = parts[0], parts[1]
        
        if protA in seq_db and protB in seq_db:
            found += 1
        else:
            missing_pairs.append((protA, protB))

print(f"   Total pairs: {total}")
print(f"   Pairs with BOTH sequences: {found} ({found/total*100:.1f}%)")
print(f"   Pairs with MISSING sequences: {len(missing_pairs)} ({len(missing_pairs)/total*100:.1f}%)")

if missing_pairs:
    print(f"\n   First 10 still missing:")
    for protA, protB in missing_pairs[:10]:
        found_a = protA in seq_db
        found_b = protB in seq_db
        print(f"     {protA} ({'✓' if found_a else '✗'}) + {protB} ({'✓' if found_b else '✗'})")

# 6. Save final database
print("\n6. Saving final database...")
with open('uniprot_final.pkl', 'wb') as f:
    pickle.dump(seq_db, f)

file_size = os.path.getsize('uniprot_final.pkl') / 1024 / 1024
print(f"   Saved to: uniprot_final.pkl")
print(f"   File size: {file_size:.1f} MB")
print(f"   Total entries: {len(seq_db)}")

# 7. Final verification
print("\n7. Final verification...")
test_ids = ['A0A178VKB9', 'A0A384KCR8', 'Q9LUI9', 'F4IDF1']
print("   Testing sample IDs:")
for tid in test_ids:
    if tid in seq_db:
        print(f"     ✓ {tid}: {len(seq_db[tid])} aa")
    else:
        print(f"     ✗ {tid}: NOT FOUND")

print("\n" + "=" * 60)
print("MERGE COMPLETE!")
print("=" * 60)