import pickle
from Bio import SeqIO

# Load original missing IDs
with open('missing_ids.txt') as f:
    missing_ids = [line.strip() for line in f if line.strip()]

# Load UniParc FASTA
uniparc_ids = []
for record in SeqIO.parse('uniparc_348.fasta', 'fasta'):
    uniparc_ids.append(record.id)

print(f"Missing IDs: {len(missing_ids)}")
print(f"UniParc IDs: {len(uniparc_ids)}")

# Create mapping (you need to manually match these!)
# The order should match what UniProt gave you
id_mapping = {}
for i, orig_id in enumerate(missing_ids):
    if i < len(uniparc_ids):
        id_mapping[orig_id] = uniparc_ids[i]

# Save mapping
import json
with open('uniparc_id_mapping.json', 'w') as f:
    json.dump(id_mapping, f, indent=2)

print(f"Created mapping for {len(id_mapping)} IDs")
print("Check data/id_mapping.json to verify the mappings are correct!")