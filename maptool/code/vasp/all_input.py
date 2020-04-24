#!/usr/bin/env python
import os
import re
import itertools
import six

from math import ceil
from pymatgen.io.vasp.sets import MITRelaxSet,MITNEBSet,MITMDSet,MPRelaxSet,MPHSERelaxSet,MPStaticSet,\
              MPHSEBSSet,MPNonSCFSet,MPSOCSet,MVLElasticSet,MVLGWSet,MVLSlabSet,MVLGBSet, MVLNPTMDSet
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar,Kpoints,Kpoints,Potcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.io import zopen
from tabulate import tabulate
from monty.json import MSONable
from pymatgen.util.string import str_delimited
from pymatgen.util.io_utils import clean_lines
from pymatgen.electronic_structure.core import Magmom
from maptool.util.utils import *
from maptool import NAME,mlog
from maptool.code.vasp.templates import *


def generate_all_input(struct):
    # select one below

    sepline(ch=' generate input files ', sp='-')
    print("your choice?")
    inputsets = ['MITRelaxSet', 'MITNEBSet', 'MITMDSet', 'MPRelaxSet', 'MPHSERelaxSet', 'MPStaticSet',
                 'MPHSEBSSet', 'MPNonSCFSet', 'MPSOCSet', 'MVLElasticSet',
                 'MVLGWSet', 'MVLSlabSet', 'MVLGBSet',
                 'MVLNPTMDSet']
    print('{} >>> {}'.format('1 ', 'MIT Relax Set'))
    print('{} >>> {}'.format('2 ', 'MIT NEB Set'))
    print('{} >>> {}'.format('3 ', 'MIT MD Set'))
    print('{} >>> {}'.format('4 ', 'MP Relax Set'))
    print('{} >>> {}'.format('5 ', 'MP HSE Relax Set'))
    print('{} >>> {}'.format('6 ', 'MP None SCF Set'))
    print('{} >>> {}'.format('7 ', 'MP SOC Set'))
    print('{} >>> {}'.format('8 ', 'MVL Elastic Set'))
    print('{} >>> {}'.format('9 ', 'MVL GW Set'))
    print('{} >>> {}'.format('10', 'MVL Slab Set'))
    print('{} >>> {}'.format('11', 'MVL GB Set'))
    print('{} >>> {}'.format('12', 'MVL NPT Set'))
    wait_sep()
    in_str = wait()
    choice = int(in_str)
    selected_set = inputsets[choice - 1]
    cmd = selected_set + '(struct)'
    print(cmd)
    outset = eval(cmd)
    outset.incar['NSW'] = 100
    print('writting the input files !')
    outset.write_input('./input', include_cif=True)
    