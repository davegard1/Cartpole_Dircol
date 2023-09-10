import numpy as np

# x = np.arange(0.0, 2.0, 0.01)
# y = 1 + np.sin(2*np.pi*x)
# z = 1 + np.cos(2*np.pi*x)

def load_data(data_file,delim=","):

    data = np.loadtxt(data_file, comments="#", delimiter=delim, unpack=False)

    n = len(data) # rows of data
    m = len(data[0]) # cols of data

    return data


