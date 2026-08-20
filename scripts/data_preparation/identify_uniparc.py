import pickle

# Load your PKL
with open('uniprot_arabidopsis.pkl', 'rb') as f:
    pkl_db = pickle.load(f)

# Extract missing IDs from your dataset
missing_ids = set()
with open('DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt') as f:
    for i, line in enumerate(f):
        if i < 2 or line.startswith('#') or line.startswith('Protein'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        protA, protB = parts[0], parts[1]
        if protA not in pkl_db:
            missing_ids.add(protA)
        if protB not in pkl_db:
            missing_ids.add(protB)

print(f"Missing IDs: {len(missing_ids)}")

# Save to file for download
with open('missing_ids.txt', 'w') as f:
    for mid in sorted(missing_ids):
        f.write(mid + '\n')

print(f"Saved to missing_ids.txt")