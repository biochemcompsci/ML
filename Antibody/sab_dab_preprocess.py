from Bio.PDB import NeighborSearch, PDBParser
from Bio.SeqUtils import seq1
import numpy as np
from Antibody import constants
from Antibody import utils

if __name__ == "__main__":

    data_path = '/home/biochemcompsci/datasets/epitope/sabdab/'

    pdb_metadata = utils.load_sabdab_metadata(data_path + '20210626_0570866_summary.tsv',
                                              [constants.pdb_key,
                                               constants.ab_heavy_chain_key,
                                               constants.ab_light_chain_key,
                                               constants.antigen_key])

    curr_pdb = PDBParser()

    # load & parse all PDBs in SabDab download summary file
    for pdb_idx, curr_pdb_id in enumerate(pdb_metadata[constants.pdb_key]):

        curr_structure = curr_pdb.get_structure(curr_pdb_id, data_path + curr_pdb_id + '.pdb')
        curr_struct_model = curr_structure[0]

        curr_antigen_chain_id = pdb_metadata[constants.antigen_key].iloc[pdb_idx]
        curr_Hchain_id = pdb_metadata[constants.ab_heavy_chain_key].iloc[pdb_idx]
        curr_Lchain_id = pdb_metadata[constants.ab_light_chain_key].iloc[pdb_idx]

        # remove non amino acids from antigen, heavy, light chain
        utils.remove_non_aa_from_chain(curr_struct_model[curr_antigen_chain_id])
        utils.remove_non_aa_from_chain(curr_struct_model[curr_Hchain_id])
        utils.remove_non_aa_from_chain(curr_struct_model[curr_Lchain_id])

        # Get antigen and antibody residues from the structure model
        curr_antigen_residues = curr_struct_model[curr_antigen_chain_id].get_residues()
        curr_Hchain_residues = curr_struct_model[curr_Hchain_id].get_residues()
        curr_Lchain_residues = curr_struct_model[curr_Lchain_id].get_residues()

        # NeighborSearch module instantiated for epitope-paratope contact residue search
        curr_atom_list = [atom for atom in curr_struct_model.get_atoms()]
        curr_epitope_paratope_contacts = NeighborSearch(curr_atom_list)

        # encode sequence in 1-hot n X len(AA codes) matrix, n = seq length
        curr_antigen_seq = np.zeros((len(curr_struct_model[curr_antigen_chain_id].child_list), len(constants.aa_lookup_table)), dtype=int)

        # encode epitope labels in 1-hot n vector, n = seq length
        curr_antigen_epitope = np.zeros(len(curr_struct_model[curr_antigen_chain_id].child_list), dtype=int)

        # update antigen sequence 1-hot matrix enconding
        # update antigen epitope 1-hot vector encoding
        for res_idx, curr_residue in enumerate(curr_antigen_residues):

            curr_AA = seq1(curr_residue.get_resname(), undef_code='?')

            if curr_AA != '?':

                curr_antigen_seq[res_idx, np.where(constants.aa_lookup_table == curr_AA)] = 1

                # if current residue has an atom in contact with antibody chain(s), update epitope encoding vector
                curr_residue_atoms = curr_residue.get_atoms()
                for curr_atom in curr_residue_atoms:

                    curr_atom_contacts = curr_epitope_paratope_contacts.search(curr_atom.get_coord(), constants.max_interface_contact_dist, 'C')
                    print(curr_residue)
                    if any([contact_chains.get_id() in [curr_Hchain_id, curr_Lchain_id]
                            for contact_chains in curr_atom_contacts]):

                        curr_antigen_epitope[res_idx] = 1
                        print(curr_antigen_epitope)
                        break

