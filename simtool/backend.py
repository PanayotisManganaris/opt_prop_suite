"""
Library of QuantumEspresso simulator IO tools designed for nanohub
simtools to be published using the Anaconda2020.11 jupyter tool
default kernel.
"""
#logging edict
import logging
logfmt = '[%(levelname)s] %(asctime)s - %(message)s'
logging.basicConfig(filename='run.log', level=logging.DEBUG, datefmt="%Y-%m-%d %H:%M:%S", format=logfmt)

#data tools
import pandas as pd
import numpy as np

#cli utilities
import io
import shutil
import subprocess

#nanohub utilities
import hublib.use
import fileinput
from simtool import DB, parse

#data handling
from pymatgen.core import Structure, Lattice, Element, Composition

#misc tools
from typing import Any, NoReturn, Tuple, Union
import os
import re
import math
import itertools
from functools import partial
np.set_printoptions(formatter={'float': '{: 0.16f}'.format}, suppress=False)

#summary plots
import plotly.express as px

def copyAndSaveFileAsOutput(outputVariableName,inputPath):
    """ saves an ouput variable as a file at inputpath """
    if inputPath:
        if inputPath.startswith('file://'):
            resultFile = os.path.basename(inputPath[7:])
            if resultFile != inputPath[7:]:
                shutil.copy2(inputPath[7:],resultFile)
        else:
            resultFile = os.path.basename(inputPath)
            if resultFile != inputPath:
                shutil.copy2(inputPath,resultFile)
        DB.save(outputVariableName,file=resultFile)

def copyAndSaveFileAsBogusOutput(outputVariableName,inputPath):
    """
    saves an output variable as a placeholder file if file does not
    exist
    """
    if inputPath:
        if inputPath.startswith('file://'):
            resultFile = os.path.basename(inputPath[7:])
            if resultFile != inputPath[7:]:
                shutil.copy2(inputPath[7:],resultFile)
        else:
            resultFile = os.path.basename(inputPath)
            if resultFile != inputPath:
                shutil.copy2(inputPath,resultFile)
        try:
            DB.save(outputVariableName,file=resultFile + "_bogus")
        except FileNotFoundError as e:
            print("%s" % (e.args[0]))

## input verification
def verify_composition(struct:Structure)->Union[Composition,None]:
    """
    Validates a formula string.
    Produces standard composition object.
    Logs particularities of the transcription
    """
    compound = struct.formula
    logging.info(f"User Entered: {compound}")
    vcomp = Composition.ranked_compositions_from_indeterminate_formula(compound, lock_if_strict=True)
    logging.info(f"User likely refering to one of these compositions: {vcomp}")
    if isinstance(vcomp, list) and vcomp:
        vcomp = vcomp[0] #ordered by likelihood so pick best candidate
        logging.debug(f"vcomp is type {type(vcomp)}")
        if not vcomp.valid:
            logging.info(f"the input composition contains dummy species")
        return vcomp
    else:
        raise ValueError(f"The compound: {compound} is not recognizable as a valid chemical formula")

def verify_site_fidelity(struct:Structure)->NoReturn:
    """
    test resulting periodicsites objects for significant digits.
    Warn user of insufficient precision and advise set nosym TRUE.
    """
    pass

# Attempt to generate pseudopotentials on the spot?
# something like: pmg potcar --symbols Li_sv O --functional PBE
# Can potcar be converted to UPF format?

## simulation transcription
def get_constituents(struct:Structure)->list:
    constituent_elements = [element for element in itertools.chain(
        *[site.species.elements for site in struct.sites]
    )
                            ]
    return constituent_elements

def order_constituents(struct:Structure)->np.ndarray:
    """ sort formula elements and atomic mass by descending atomic mass """
    uniquel = np.unique(get_constituents(struct))
    amass_ar = np.array([el.atomic_mass for el in uniquel])
    massive_first = np.argsort(amass_ar).tolist()[::-1]
    amass_falling = amass_ar[massive_first]
    el_falling = np.array([el.symbol for el in uniquel])[massive_first]
    logging.debug(f"Elements ordered by descending mass: {el_falling}. Descending mass: {amass_falling}")
    return el_falling

def get_pps(struct:Structure, pp_class:str)->Tuple[list,list]:
    """ identifies pseudopotential files to use in simulation """
    symbols = order_constituents(struct)
    pp_files = os.listdir(f"./pseudo/pseudo_{pp_class}/")

    logging.info(f'retrieving pseudo-potentials of type {pp_class} for elements: {symbols}')
    # find corresponding .upf files from pseudo list
    pps = []
    for symbol in symbols:
        for pp in pp_files:
            identifier = re.split('\.|\_',pp)[0]
            if symbol == identifier:
                pps.append(pp)
    logging.debug(f"Pseudopotentials in order of descending mass: {pps}")
    if len(pps) != len(symbols):
        logging.warning(f"There are not as many Pseudopotentials as Species. This may not work")
    return pps #should share order of symbols

def get_functional_from_pp(pps:list)->str:
    """
    Frees the user from the need to consider valid combinations of
    pseudo-potentials and xc functionals.
    Reads the pp contents for the functional
    """
    return 'pbe'

## extract input data


## 
