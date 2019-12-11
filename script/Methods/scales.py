import pickle
import numpy as np
def stScales(seq, stscale_data_path = "data/stScales.pkl"):
    """
        Args: amino acid sequence, file path to stscale data
        Returns: Average of the ST-Scales of each amino acid in seq
    """

    # Open stscale data file (contains stscale for each amino acid)
    with open(stscale_data_path, 'rb') as f:
        data = pickle.load(f)

    # Calculate average
    total_stscale = np.zeros(8)
    for AA in list(seq):
        total_stscale = np.add(total_stscale, np.array(data[AA]))
    return total_stscale/len(seq)