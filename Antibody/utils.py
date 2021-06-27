from Bio.SeqUtils import IUPACData
from Bio.SeqUtils import seq1
import pandas as pd
import numpy as np
import pickle
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


def update_seq_and_contact_encoding(neighbor_search, residues, contact_query_chains, seq_encoding, contact_encoding):

    seq_encoding_update = seq_encoding
    contact_encoding_update = contact_encoding

    for res_idx, curr_residue in enumerate(residues):

        curr_AA = seq1(curr_residue.get_resname(), undef_code='?')

        if curr_AA != '?':

            seq_encoding_update[res_idx, np.where(constants.aa_lookup_table == curr_AA)] = 1

            # if current residue has an atom in contact with antibody chain(s), update epitope encoding vector
            curr_residue_atoms = curr_residue.get_atoms()
            for curr_atom in curr_residue_atoms:

                curr_atom_contacts = neighbor_search.search(curr_atom.get_coord(),
                                                            constants.max_interface_contact_dist,
                                                            'C')

                if any([contact_chains.get_id() in contact_query_chains
                        for contact_chains in curr_atom_contacts]):
                    contact_encoding_update[res_idx] = 1
                    break

    return seq_encoding_update, contact_encoding_update


def load_and_merge_serialized_dicts(serialized_dict_files):

    merged_dict = {}

    for curr_dict_file in serialized_dict_files:

        input_file = open(curr_dict_file, 'rb')
        merged_dict.update(pickle.load(input_file))
        input_file.close()

    return merged_dict

if __name__ == "__main__":

    data_path = '/home/biochemcompsci/datasets/epitope/sabdab/'

    # ******* sanity check antigen, antibody sequence and epitope, paratope encodings
    serialized_encoding_file = 'light_chain_paratopes_encoding_merged.pickle'
    input_file = open(data_path + serialized_encoding_file, 'rb')

    encoding = pickle.load(input_file)

    print()

    # ******* merge serialized dictionaries ******
    #
    #
    # # antigen_seqs_input_file = 'antigen_seqs_encoding.pickle'
    # # antigen_seqs_input_file_2 = 'antigen_seqs_encoding_2.pickle'
    # # antigen_epitopes_input_file = 'antigen_epitopes_encoding.pickle'
    # # antigen_epitopes_input_file_2 = 'antigen_epitopes_encoding_2.pickle'
    # # heavy_chain_seqs_input_file = 'heavy_chain_seqs_encoding.pickle'
    # # heavy_chain_seqs_input_file_2 = 'heavy_chain_seqs_encoding_2.pickle'
    # # heavy_chain_paratopes_input_file = 'heavy_chain_paratopes_encoding.pickle'
    # # heavy_chain_paratopes_input_file_2 = 'heavy_chain_paratopes_encoding_2.pickle'
    # # light_chain_seqs_input_file = 'light_chain_seqs_encoding.pickle'
    # # light_chain_seqs_input_file_2 = 'light_chain_seqs_encoding_2.pickle'
    # light_chain_paratopes_input_file = 'light_chain_paratopes_encoding.pickle'
    # light_chain_paratopes_input_file_2 = 'light_chain_paratopes_encoding_2.pickle'
    #
    # merged_output_file = 'light_chain_paratopes_encoding_merged.pickle'
    #
    # merged_encoding = load_and_merge_serialized_dicts([data_path + light_chain_paratopes_input_file,
    #                                                    data_path + light_chain_paratopes_input_file_2])
    #
    # # serialize merged encoding dictionaries to a file
    # out_file = open(data_path + merged_output_file, 'wb')
    # pickle.dump(merged_encoding, out_file)
    # out_file.close()
