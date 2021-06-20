import numpy as np
import _pickle as cPickle
import dill

if __name__ == "__main__":

    # required to read Python 2 pickled data into Python 3
    dill._dill._reverse_typemap['ObjectType'] = object

    epitope_data = ["/home/biochemcompsci/datasets/epitope/pittala_context_aware_structure_ag_ab/data_epipred/train.cpkl",
                    "/home/biochemcompsci/datasets/epitope/pittala_context_aware_structure_ag_ab/data_epipred/test.cpkl",
                    "/home/biochemcompsci/datasets/epitope/pittala_context_aware_structure_ag_ab/data_epipred/val.cpkl",
                    "/home/biochemcompsci/datasets/epitope/pittala_context_aware_structure_ag_ab/data_paratope/train.cpkl",
                    "/home/biochemcompsci/datasets/epitope/pittala_context_aware_structure_ag_ab/data_paratope/test.cpkl",
                    "/home/biochemcompsci/datasets/epitope/pittala_context_aware_structure_ag_ab/data_paratope/val.cpkl"]

    # lookup table to map 1-hot AA vector encoding to single-letter code
    aa_lookup_table = np.array(
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'])

    # load labeled antigens into array
    antigen_sequences = []
    labeled_epitopes = []

    for currFile in epitope_data:

        print("Loading data from file: " + str(currFile))

        # encoding=bytes required to read Python 2 pickled data into Python 3
        # data is a list of dictionaries describing Ab-Ag interactions with format described at https://github.com/vamships/PECAN
        epitopes_processed = {}
        with open(currFile, "rb") as epitope_file:
            epitopes_processed = cPickle.load(epitope_file, encoding="bytes")

        for idx in range(0, len(epitopes_processed) - 1):

            if idx % 2 == 0:

                # need str.encode for Python 2 to Python 3 interop
                curr_antigen = epitopes_processed[idx].get(str.encode('l_vertex'))

                # loop through all sequence positions
                curr_sequence = ""
                for idx_b in range(0, len(curr_antigen) - 1):

                    one_hot_aa_vec = np.array(curr_antigen[idx_b][0:21])
                    curr_sequence = curr_sequence + aa_lookup_table[one_hot_aa_vec.nonzero()][0]

                antigen_sequences.append(curr_sequence)
                labeled_epitopes.append(epitopes_processed[idx].get(str.encode('label'))[:, 1])
                # print(curr_sequence)
                # print(labeled_epitopes[len(labeled_epitopes) - 1])

    print(len(antigen_sequences))
    print(len(labeled_epitopes))



