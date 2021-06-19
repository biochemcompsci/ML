import numpy as np
import _pickle as cPickle
import dill

if __name__ == "__main__":

    # required to read Python 2 pickled data into Python 3
    dill._dill._reverse_typemap['ObjectType'] = object

    epitope_data = "/home/biochemcompsci/datasets/epitope/pittala_context_aware_structure_ag_ab/data_epipred/train.cpkl"

    # encoding=bytes required to read Python 2 pickled data into Python 3
    # data is a list of dictionaries describing Ab-Ag interactions with format described at https://github.com/vamships/PECAN
    with open(epitope_data, "rb") as epitope_file:
        epitopes_processed = cPickle.load(epitope_file, encoding="bytes")

    # lookup table to map 1-hot AA vector encoding to single-letter code
    aa_lookup_table = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'])

    print(epitopes_processed[0])

