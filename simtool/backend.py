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
from pymatgen.core import Structure, Lattice, Element, Composition, sites
from pymatgen.io.pwscf import PWInput, PWOutput

#misc tools
from typing import Any, NoReturn, Tuple, Union, Callable
import os
import re
import math
import itertools
from functools import partial
np.set_printoptions(formatter={'float': '{: 0.16f}'.format}, suppress=False)

#summary plots
import plotly.express as px

#db and ppdb instantiated in the simtool global scope.

def copyAndSaveFileAsOutput(outputVariableName,inputPath):
    """ saves an output variable as a file at inputpath """
    if inputPath:
        if inputPath.startswith('file://'):
            resultFile = os.path.basename(inputPath[7:])
            if resultFile != inputPath[7:]:
                shutil.copy2(inputPath[7:],resultFile)
        else:
            resultFile = os.path.basename(inputPath)
            if resultFile != inputPath:
                shutil.copy2(inputPath,resultFile)
        db.save(outputVariableName,file=resultFile)

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
            db.save(outputVariableName,file=resultFile + "_bogus")
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

## Parse and Prepare Data for Transcription
#outdated utility?
def get_constituents(struct:Structure)->list:
    constituent_elements = [element for element in itertools.chain(
        *[site.species.elements for site in struct.sites]
    )
                            ]
    return constituent_elements

def order_constituents(struct:Structure)->Tuple[np.ndarray, np.ndarray]:
    """ sort constituent elements and atomic mass by descending atomic mass """
    uniquel = np.unique(get_constituents(struct))
    amass_ar = np.array([el.atomic_mass for el in uniquel])
    massive_first = np.argsort(amass_ar).tolist()[::-1]
    amass_falling = amass_ar[massive_first]
    el_falling = np.array([el.symbol for el in uniquel])[massive_first]
    logging.debug(f"Elements ordered by descending mass: {el_falling}. Descending mass: {amass_falling}")
    return el_falling, amass_falling

def get_pps(struct:Structure, ppdb:dict)->list:
    """ identifies pseudopotential files to use in simulation """
    symbols, _ = order_constituents(struct)
    pps = _lookup_pp(symbols, ppdb)
    logging.debug(f"Pseudopotentials in order of descending mass: {pps}")
    if len(pps) != len(symbols):
        logging.warning(f"There are not as many Pseudopotentials as Species. This may not work")
    return pps #should share order of symbols

def _lookup_pp(symbols:list, ppdb:dict)->list:
    pps = []
    for symbol in symbols:
        pps.append(ppdb[symbol])
    return pps


## assign information to sites using setters 
SiteObj = Union[sites.Site, sites.PeriodicSite]

def set_site_properties(struct:Structure,
                        key:str,
                        setter:Callable[[SiteObj], Any])->Structure:
    """
    Sets a structure's site objects property key to a value obtained by a setter function.
    
    """
    #notice, for multiple elements in a single site (disordered structures)
    #the appropriate setter will have to make sense of the listed members
    for site in struct.sites:
        try:
            site.properties[key] = setter(site)
        except AttributeError as e:
            logging.error(e.args)
    return struct

## setter functions obtain an external data value for any possible value of a given site attribute
def pseudo_setter(site:SiteObj)->str:
    symbols = [element.symbol for element in site.species.elements]
    
    return symbols

## input printers
class BlockPrinter():
    def __init__(self, struct:Structure, pp_class:str):
        self.struct = struct
        self.ppc = pp_class

    def print_functional(self)->str:
        """
        Frees the user from the need to consider valid combinations of
        pseudo-potentials and xc functionals.
        Reads the pp contents for the functional
        """
        pps = get_pps(self.struct, self.ppc)
        return 'pbe'
    
    def print_matrix(self)->str:
        """
        Create cell parameters block.
        Control structures' lattice matrix units.
        Defaults to Å
        maybe Bohr radii
        to much precision in the input files causes convergence issues. to little causes symmetry issues.
        """
        bio = io.BytesIO()
        np.savetxt(bio, self.struct.lattice.matrix*1.88973, fmt="%0.8f")
        cell_parameters_block = bio.getvalue().decode('latin1')
        logging.info(
            f"""The standardized conventional lattice parameters are obtained:
            {cell_parameters_block}"""
        )
        return cell_parameters_block
    
    def print_species(self)->str:
        """ pw's atomic Species block and dm's amass_block """
        el_falling, amass_falling = order_constituents(self.struct)
        pps = get_pps(self.struct, self.ppc)
        atomic_species_block = """"""
        for pp, amass, species in zip(pps, amass_falling, el_falling):
            block_line = f"{species} {amass} {pp}\n"
            atomic_species_block += block_line
        logging.info(
            f"""The atomic_species_block:
            {atomic_species_block}"""
        )
        return atomic_species_block
    
    def print_amass(self)->str:
        _, amass_falling = order_constituents(self.struct)
        amass_block = """"""
        for ind, amass in enumerate(amass_falling):
            block_line = f"amass({ind+1})={amass},\n  "
            amass_block += block_line
        logging.info(
            f"""The amass_block:
            {amass_block}"""
        )
        return amass_block
    
    @staticmethod
    def get_fractional_coords(periodicsitesobj)->np.ndarray:
        at_site = np.array(
            [str(*periodicsitesobj.species.to_reduced_dict.keys()),
             periodicsitesobj.a, periodicsitesobj.b, periodicsitesobj.c]
        )
        logging.debug(f"atomic sites array: {at_site}")
        # str might break if the species has whitespaces in it... I think. Or maybe that's just for unpacked lists
        return at_site
    
    def print_frac_coords(self)->str:
        atomic_positions_block = np.array(list(map(lambda x: self.get_fractional_coords(x), self.struct.sites)))
        bio = io.BytesIO()
        np.savetxt(bio, atomic_positions_block, fmt="%s", encoding="latin1")
        atomic_positions_block = bio.getvalue().decode('latin1')
        logging.info(
            f"""The atomic sites are obtained:
            {atomic_positions_block}"""
        )
        return atomic_positions_block
    
    def print_kpoints(self)->str:
        """ generate a k-point grid from the provided kpoint magnitude and lattice vectors """
        b = np.array([1/a for a in self.struct.lattice.abc])
        b = b*(1/max(b))
        k = np.ceil(kpoints*b)
        kpoints = f"""{k[0]} {k[1]} {k[2]} 0 0 0"""
        logging.info(
            f"""The kpoints are obtained:
            {kpoints}"""
        )
        return kpoints

## Simulation Pipeline Components
class NanoHUB_PWInput(PWInput):
    """
    https://pymatgen.org/pymatgen.io.pwscf.html

    Handles setting up pw.x input files. It is equipped with a submit
    method for running the resulting job on nanoHUB.
    """
    def __init__(self,
                 struct:Structure,
                 pp_class:str,                 
                 control:dict=None,
                 system:dict=None,
                 electrons:dict=None,
                 ions:dict=None,
                 cell:dict=None,
                 kpoints_mode:str="automatic",
                 kpoints_grid:dict=None,
                 kpoints_shift:dict=None,
                 ):
        self.vcomp = verify_composition(struct)
        self.nat_num = len(struct.sites)
        self.ntyp_num = len(np.unique(get_constituents(struct)))
        logging.debug(f"ntyp = {self.ntyp_num}")

    def _input(self)->str:
        return f"""
         &control
            calculation='vc-relax',
            restart_mode='from_scratch',
            prefix='{self.vcomp.reduced_formula}',
            outdir='./',
            pseudo_dir = './',
            etot_conv_thr=1.0d-6,
            forc_conv_thr=1.0d-6,
        /
         &system    
            ibrav= 0, celldm(1)=1 ,nat= {self.nat_num}, ntyp= {self.ntyp_num},
            ecutwfc = {self.ecutwfc}, ecutrho = {self.ecutrho},
            occupations={self.smearing}, {"smearing='mp', degauss=0.06," if self.smearing == 'smearing' else ""}
            input_dft = {self.bp.print_functional()},
        /
         &electrons
            mixing_beta =0.7,
            conv_thr =1.0d-8,
        /
         &ions
            ion_dynamics='bfgs'
        /
         &cell
            cell_dynamics='bfgs',
            press=0.0,
            cell_factor=2.0,
            press_conv_thr=0.5,
        /
         CELL_PARAMETERS (alat= 1.00000000)
         {self.bp.print_matrix()}
         ATOMIC_SPECIES
         {self.bp.print_species()}
         ATOMIC_POSITIONS (crystal)
         {self.bp.print_frac_coords()}
         K_POINTS (automatic)
         {self.bp.print_kpoints()}
         """

    def write_input_file(self)->NoReturn:
        with open(f"{self.vcomp.reduced_formula}.vc-relax.in", "w") as input_file:
            input_file.write(self._input())
        
    def submit_job(self)->NoReturn:
        pass
    
    def write_output_file(self)->NoReturn:
        pass
