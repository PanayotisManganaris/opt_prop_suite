# Import libraries
import numpy as np
import logging
import math

logfmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(filename='run.log', level=logging.DEBUG, datefmt="%Y-%m-%d %H:%M:%S", format=logfmt)

# Extract the pressures second derivates of energy
def check_vcr_convergence(vcr_output_name):

    # Read energy values from vc-relax output file
    vcr_output_file = open(vcr_output_name,"r")
    lines = vcr_output_file.readlines()
    E = []
    for line in lines:
        if '!' in line:
            data = line.split()
            E.append(data[4])
    vcr_output_file.close()
    E = [float(e) for e in E]

    # Get second derivative of energy values
    dE = np.diff(E, 2)

    # Read pressure values
    stresses = []
    for i, line in enumerate(lines):
        if '(kbar)' in line:
            stress = float(line.split()[-1])
            stresses.append(stress)
    return dE, stresses

# Check for negative frequencies in the spectra df
def check_spectra_convergence(spectradf, threshold = 0):
    negative_frequencies = spectradf[spectradf['[cm-1]'] < threshold]
    return negative_frequencies

# Run all convergence checks simultaneously, return convergence boolean
def check_convergence(vcr_output_name, spectradf):
    dE, stresses = check_vcr_convergence(vcr_output_name)
    negative_frequencies = check_spectra_convergence(spectradf)

    if (dE[-1] > 1e-5):
        logging.info('The energy does not seem to be well converged')
        print('a')
        return False
    if (abs(stresses[-1]) > 10):
        logging.info('The stresses are larger than 10 kbar')
        print('b')
        return False
    if (len(negative_frequencies) > 3):
        logging.info('There are more than 3 negative frequencies')
        print('c')
        return False

    return True

# Estimate the number of required num_CPUs and walltime for a given structure
def get_queue_settings(struct_dict):
    # Check number of atoms
    n_atoms = len(struct_dict['sites'])

    # Try to scale up to 4 num_CPUs in 1:00:00 walltime, then increase walltime if more than 4 num_CPUs are needed
    f = 4
    num_CPUs = math.ceil(n_atoms/f)
    walltime = "01:00:00"
    if num_CPUs > 4:
        n = math.ceil(num_CPUs/4)
        num_CPUs = 4
        if n < 10:
            walltime = f"0{n}:00:00"
        else:
            walltime = f"{n}:00:00"

    return num_CPUs, walltime
