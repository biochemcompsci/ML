from Bio.SeqUtils import IUPACData
import pandas as pd
from Antibody import constants


def remove_non_aa_from_chain(pdb_chain):

    non_aa_res_ids = []

    for curr_residue in pdb_chain:

        if curr_residue.get_resname().capitalize() not in IUPACData.protein_letters_3to1:
            # add non amino acid residues to beginning of list so they can be removed in order
            non_aa_res_ids.insert(0, curr_residue.id)

    for curr_non_aa_id in non_aa_res_ids:
        pdb_chain.detach_child(curr_non_aa_id)


def load_sabdab_metadata(summary_file_path, load_cols):

    # SabDab download summary file
    # Details in "Summary file fields" section of http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/#formats
    pdb_metadata = pd.read_csv(summary_file_path, delimiter='\t', header=0,
                               usecols=load_cols)

    pdb_metadata = pdb_metadata.drop_duplicates(subset=constants.pdb_key, keep='first')

    return pdb_metadata