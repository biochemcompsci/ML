import numpy as np

# lookup table to map 1-hot AA vector encoding to single-letter code
aa_lookup_table = np.array(
    ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'])

# Max distance (A) for residue contact
max_interface_contact_dist = 4.5

# Keys for antigen, antibody chains – based on SabDab naming
pdb_key = 'pdb'
antigen_key = 'antigen_chain'
ab_heavy_chain_key = 'Hchain'
ab_light_chain_key = 'Lchain'