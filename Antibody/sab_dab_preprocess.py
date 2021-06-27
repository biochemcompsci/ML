import warnings
from Bio.PDB import NeighborSearch, PDBParser
import numpy as np
import pickle
from Antibody import constants
from Antibody import utils

warnings.filterwarnings('ignore')

if __name__ == "__main__":

    data_path = '/home/biochemcompsci/datasets/epitope/sabdab/'
    sabdab_input_file = '20210623_0909829_summary.tsv'
    antigen_seqs_output_file = 'antigen_seqs_encoding_2.pickle'
    antigen_epitopes_output_file = 'antigen_epitopes_encoding_2.pickle'
    heavy_chain_seqs_output_file = 'heavy_chain_seqs_encoding_2.pickle'
    heavy_chain_paratopes_output_file = 'heavy_chain_paratopes_encoding_2.pickle'
    light_chain_seqs_output_file = 'light_chain_seqs_encoding_2.pickle'
    light_chain_paratopes_output_file = 'light_chain_paratopes_encoding_2.pickle'

    pdb_metadata = utils.load_sabdab_metadata(data_path + sabdab_input_file,
                                              [constants.pdb_key,
                                               constants.ab_heavy_chain_key,
                                               constants.ab_light_chain_key,
                                               constants.antigen_key])

    # sequence and contact residue encodings stored in dictionary with PDB ID as key
    antigen_seqs_encoding = {}
    antigen_epitopes_encoding = {}
    heavy_chain_seqs_encoding = {}
    heavy_chain_paratopes_encoding = {}
    light_chain_seqs_encoding = {}
    light_chain_paratopes_encoding = {}

    curr_pdb = PDBParser()

    # load & parse all PDBs in SabDab download summary file
    for pdb_idx, curr_pdb_id in enumerate(pdb_metadata[constants.pdb_key]):

        print(f"Processing PDB: {curr_pdb_id}")

        curr_structure = curr_pdb.get_structure(curr_pdb_id, data_path + curr_pdb_id + '.pdb')

        curr_struct_model = curr_structure[0]

        curr_antigen_chain_id = pdb_metadata[constants.antigen_key].iloc[pdb_idx]
        curr_Hchain_id = pdb_metadata[constants.ab_heavy_chain_key].iloc[pdb_idx]
        curr_Lchain_id = pdb_metadata[constants.ab_light_chain_key].iloc[pdb_idx]

        # ignore PDBs with ambiguous chain IDs
        if len(str(curr_antigen_chain_id)) == 1 and len(str(curr_Hchain_id)) == 1 and len(str(curr_Lchain_id)) == 1:

            # remove non amino acids from antigen, heavy, light chain
            utils.remove_non_aa_from_chain(curr_struct_model[curr_antigen_chain_id])
            utils.remove_non_aa_from_chain(curr_struct_model[curr_Hchain_id])
            utils.remove_non_aa_from_chain(curr_struct_model[curr_Lchain_id])

            # NeighborSearch module instantiated for epitope-paratope contact residue search
            curr_atom_list = [atom for atom in curr_struct_model.get_atoms()]
            curr_epitope_paratope_contacts = NeighborSearch(curr_atom_list)

            # ******** ANTIGEN ENCODING ********************
            # encode antigen sequence in 1-hot n X len(AA codes) matrix, n = seq length
            curr_antigen_seq = np.zeros(
                (len(curr_struct_model[curr_antigen_chain_id].child_list), len(constants.aa_lookup_table)), dtype=int)

            # encode epitope labels in 1-hot n vector, n = seq length
            curr_antigen_epitope = np.zeros(
                len(curr_struct_model[curr_antigen_chain_id].child_list), dtype=int)

            # Get antigen and antibody residues from the structure model
            curr_antigen_residues = curr_struct_model[curr_antigen_chain_id].get_residues()

            # # update antigen sequence 1-hot matrix enconding
            # # update antigen epitope 1-hot vector encoding
            curr_antigen_seq, curr_antigen_epitope = utils.update_seq_and_contact_encoding(curr_epitope_paratope_contacts,
                                                                                           curr_antigen_residues,
                                                                                           [curr_Hchain_id, curr_Lchain_id],
                                                                                           curr_antigen_seq,
                                                                                           curr_antigen_epitope)

            antigen_seqs_encoding[curr_pdb_id] = curr_antigen_seq
            antigen_epitopes_encoding[curr_pdb_id] = curr_antigen_epitope

            # ******** HEAVY CHAIN ENCODING ********************
            # encode heavy chain sequence in 1-hot n X len(AA codes) matrix, n = seq length
            curr_Hchain_seq = np.zeros(
                (len(curr_struct_model[curr_Hchain_id].child_list), len(constants.aa_lookup_table)), dtype=int)

            # encode heavy chain paratope labels in 1-hot n vector, n = seq length
            curr_Hchain_paratope = np.zeros(len(curr_struct_model[curr_Hchain_id].child_list), dtype=int)

            # Get heavy chain residues from the structure model
            curr_Hchain_residues = curr_struct_model[curr_Hchain_id].get_residues()

            # update heavy chain sequence 1-hot matrix enconding
            # update heavy chain paratope 1-hot vector encoding
            curr_Hchain_seq, curr_Hchain_paratope = utils.update_seq_and_contact_encoding(curr_epitope_paratope_contacts,
                                                                                          curr_Hchain_residues,
                                                                                          [curr_antigen_chain_id],
                                                                                          curr_Hchain_seq,
                                                                                          curr_Hchain_paratope)

            heavy_chain_seqs_encoding[curr_pdb_id] = curr_Hchain_seq
            heavy_chain_paratopes_encoding[curr_pdb_id] = curr_Hchain_paratope

            # ******** LIGHT CHAIN ENCODING ********************
            # encode light chain sequence in 1-hot n X len(AA codes) matrix, n = seq length
            curr_Lchain_seq = np.zeros(
                (len(curr_struct_model[curr_Lchain_id].child_list), len(constants.aa_lookup_table)), dtype=int)

            # encode light chain paratope labels in 1-hot n vector, n = seq length
            curr_Lchain_paratope = np.zeros(len(curr_struct_model[curr_Lchain_id].child_list), dtype=int)

            # Get light chain residues from the structure model
            curr_Lchain_residues = curr_struct_model[curr_Lchain_id].get_residues()

            # update light chain sequence 1-hot matrix enconding
            # update light chain paratope 1-hot vector encoding
            curr_Lchain_seq, curr_Lchain_paratope = utils.update_seq_and_contact_encoding(curr_epitope_paratope_contacts,
                                                                                          curr_Lchain_residues,
                                                                                          [curr_antigen_chain_id],
                                                                                          curr_Lchain_seq,
                                                                                          curr_Lchain_paratope)

            light_chain_seqs_encoding[curr_pdb_id] = curr_Lchain_seq
            light_chain_paratopes_encoding[curr_pdb_id] = curr_Lchain_paratope

            print(f"Finished processing PDB: {curr_pdb_id}")

    # serialize encoding dictionaries to a file
    out_file = open(data_path + antigen_seqs_output_file, 'wb')
    pickle.dump(antigen_seqs_encoding, out_file)
    out_file.close()

    out_file = open(data_path + antigen_epitopes_output_file, 'wb')
    pickle.dump(antigen_epitopes_encoding, out_file)
    out_file.close()

    out_file = open(data_path + heavy_chain_seqs_output_file, 'wb')
    pickle.dump(heavy_chain_seqs_encoding, out_file)
    out_file.close()

    out_file = open(data_path + heavy_chain_paratopes_output_file, 'wb')
    pickle.dump(heavy_chain_paratopes_encoding, out_file)
    out_file.close()

    out_file = open(data_path + light_chain_seqs_output_file, 'wb')
    pickle.dump(light_chain_seqs_encoding, out_file)
    out_file.close()

    out_file = open(data_path + light_chain_paratopes_output_file, 'wb')
    pickle.dump(light_chain_paratopes_encoding, out_file)
    out_file.close()
