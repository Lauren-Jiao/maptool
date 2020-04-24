#!/usr/bin/env python
import os
import re
import itertools
import six

from pymatgen.io.vasp.inputs import Poscar,Kpoints,Kpoints,Potcar
from maptool.util.utils import *
from maptool import NAME,mlog
from maptool.code.vasp.templates import *

def generate_potcar(struct,dirname='.'):
    """
    Generate POTCAR according to input structure via pymatgen's POTCAR module
    """
    avail_pot =" ".join(Potcar.FUNCTIONAL_CHOICES)
    tip="""
    Available Pseudo-potentials are:
    """
    tip+=avail_pot+'\n'
    warn_tip(1,tip)
    print('your choice ?')
    wait_sep()
    in_str=wait()
    assert in_str.upper() in avail_pot

    potcar=Potcar([el.value for el in struct.types_of_specie],functional=in_str)
    potcar.write_file(os.path.join(dirname, "POTCAR"))