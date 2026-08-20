from Bio import SeqIO
import pickle
import re

print("Creating sequence database...")
seq_db = {}

for record in SeqIO.parse('/home/users/nus/e0969134/AraPPLM/data/uniprot_arabidopsis.fasta', 'fasta'):
    seq = str(record.seq)
    full_id = record.id
    
    # Store by full ID
    seq_db[full_id] = seq
    
    # Extract simple ID from format: sp|Q9LUI9|PIAL1_ARATH
    match = re.match(r'(sp|tr)\|([A-Z0-9]+)\|', full_id)
    if match:
        simple_id = match.group(2)
        seq_db[simple_id] = seq  # Now Q9LUI9 maps to the sequence!

print(f"Total mappings: {len(seq_db)}")
print(f"Unique sequences: {len(set(seq_db.values()))}")

# Test with IDs from your data
test_ids = ['Q9LUI9', 'F4IDF1', 'Q9M8S6']
for tid in test_ids:
    if tid in seq_db:
        print(f"✅ {tid}: {len(seq_db[tid])} aa")
    else:
        print(f"❌ {tid}: NOT FOUND - checking FASTA...")
        # Search in FASTA
        for rec in SeqIO.parse('/home/users/nus/e0969134/AraPPLM/data/uniprot_arabidopsis.fasta', 'fasta'):
            if tid in rec.id:
                print(f"   Found as: {rec.id}")
                break

# Save
with open('/home/users/nus/e0969134/AraPPLM/data/uniprot_arabidopsis.pkl', 'wb') as f:
    pickle.dump(seq_db, f)

print(f"\nSaved to uniprot_arabidopsis.pkl")