#!/usr/bin/env python
import os
import re
import itertools
import six

from math import ceil
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.sets import MITRelaxSet,MITNEBSet,MITMDSet,MPRelaxSet,MPHSERelaxSet,MPStaticSet,\
              MPHSEBSSet,MPNonSCFSet,MPSOCSet,MVLElasticSet,MVLGWSet,MVLSlabSet,MVLGBSet, MVLNPTMDSet
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Poscar,Kpoints,Potcar
from maptool.util.utils import *
from maptool import mlog
from maptool.code.vasp.templates import *

def generate_kpoint(struct,dirname='.'):
    '''
    Generate KPOINTS file according to user's choice. Currently supports:
      1 : automatic k-grid
      2 : Band structure k-path
      3 : HSE06 band structure k-path
      4 : 3D plot k-grid
    '''
    sepline(ch=' generate KPOINTS file ',sp='-')
    print("your choice?")
    print('''\
    1 >>> automatic k-grid 
    2 >>> Band structure k-path
    3 >>> HSE06 band structure
    4 >>> 3D plot k-grid ''')

    wait_sep()
    in_str = wait()
    choice=int(in_str)
    assert choice in [1,2,3,4]

    if choice==1:
       auto_kgrid(struct,dirname)
    elif choice==2:
       band_structure_kpath(struct,dirname)
    elif choice==3:
       pass
    else:
      generate_3Dkpoints(struct,dirname)

def auto_kgrid(struct,dirname):
    '''
    Generate KPOINTS according given mesh density.
      input the dimensionality and mesh grid density
      dimensionality can be 0D 1D 2D 3D
      500 for low grid density
      1000 for medium grid density
      2000 for high grid density
      3000 for accurate density
      input format: 1 1000

      note: in 1D system, mesh grids for x & y direction = 1
            in 2D system, mesh grids for z direction = 1
    '''
    print(" input the dimensionality and mesh grid density ")
    print(" dimensionality can be 0D 1D 2D 3D")
    print(" 500 for low grid density")
    print(" 1000 for medium grid density")
    print(" 2000 for high grid density")
    print(" 3000 for accurate density")
    print(" input format: 1 1000")
    wait_sep()
    in_str = wait()
    data = [int(x) for x in in_str]
    dim = data[0]
    grid_density = data[1]
    kps = Kpoints.automatic_density(struct, grid_density)
    if dim==0:
       kps.kpts=[[1,1,1]]
    if dim==1:
      # 1 for x and y direction
      # e.g. :
      # ""
      # comment
      # M or G
      # 1 1 <density>
      # ""
       kps.kpts[0][0]=1
       kps.kpts[0][1]=1
    if dim == 2:
        # 1 for z direction
        # e.g. :
        # ""
        # comment
        # M or G
        # <dens1> <dens2> 1
        # ""
        kps.kpts[0][2] = 1
    kps.write_file(os.path.join(dirname, "KPOINTS"))

def generate_3Dkpoints(struct,dirname):
    """
    Generate 3D system KPOINTS file according to given accuracy
    User should input a float number to set the accuracy.
    """
    #print(" input the dimensionality and mesh grid density ")
    tip='''
    Accuracy Levels: (1) Low:    0.04~0.03;
                     (2) Medium: 0.03~0.02; 
                     (2) Fine:   0.02~0.01; 
    '''
    warn_tip(1,tip)
    print("Input KP-Resolved Value (unit: 2*PI/Ang):")
    wait_sep()
    in_str = wait()
    delta_k = float(in_str)
    (lka, lkb, lkc) = struct.lattice.reciprocal_lattice.abc
    ka = ceil(lka / (2 * np.pi) / delta_k)
    kb = ceil(lkb / (2 * np.pi) / delta_k)
    ka_dis = np.linspace(-0.5, 0.5, ka)
    kb_dis = np.linspace(-0.5, 0.5, kb)
    kxx, kyy = np.meshgrid(ka_dis, kb_dis)
    kpts = [[i[0], i[1], 0.0] for i in zip(kxx.flat, kyy.flat)]

    tmp_K="3D K-meshs by maptool with : "+ str(ka)+"x"+ str(kb)+ "\n 1\nReciprocal\n 0 0 0 1"
    # initialize a Kpionts instance from template string
    reciprocal_kpoints=Kpoints.from_string(tmp_K)
    reciprocal_kpoints.labels=None
    reciprocal_kpoints.kpts=kpts
    reciprocal_kpoints.num_kpts=ka*kb
    reciprocal_kpoints.kpts_weights=[1.0]*(ka*kb)
    reciprocal_kpoints.write_file(os.path.join(dirname, "KPOINTS"))
#    print(spk)

def band_structure_kpath(struct,dirname,nkpts=30):
    """
    Generate KPOINTS file for band structure calculation via pymatgen's symmetry
    analyzing system.
    """
    #struct=Structure.from_file('POSCAR')
    #ana_struct=SpacegroupAnalyzer(struct)
    #pst=ana_struct.find_primitive()
    # First brillouin zone
    ibz = HighSymmKpath(struct)
    linemode_kpoints = Kpoints.automatic_linemode(nkpts,ibz)
    linemode_kpoints.write_file(os.path.join(dirname, "KPOINTS"))

def hse06_bandstructure_kpoints(struct, nkpts=20):
    '''
    Generate HSE06 bandstructure KPOINTS
    Append high-symmetry path points to the IBZKPT file and set weight of
    all the high-symmetry path points to zero and then write to "KPOINTS"

    High-symmetry path kpoints is saved as a backup file named 'KPOINTS_bak'

    Note: We asssert the IBZKPT file is valid
    '''
    def chunks(lst, n):
      for i in range(0, len(lst), n):
        yield lst[i: i+n]

    hsk = HighSymmKpath(struct)
    sym_kpts = Kpoints.automatic_linemode(nkpts, hsk)
    sym_kpts.write_file("KPOINTS_bak")

    kpts = sym_kpts.kpts
    nsegs = sym_kpts.num_kpts

    kpoints_result = []
    for rng in chunks(kpts, 2):
      start, end = rng
      kpoints_result.append(np.linspace(start, end, nsegs))
    kpoints_result = np.array(kpoints_result).reshape((-1, 3))

    KPOINTS = open('IBZKPT').readlines()
    for i in range(kpoints_result.shape[0]):
      x, y, z = kpoints_result[i, :]
      KPOINTS.append("{:20.14f}{:20.14f}{:20.14f}{:14}\n".format(x, y, z, 0))
    KPOINTS[1] = "{:8}\n".format(len(KPOINTS) - 3)
    with open("KPOINTS", 'w') as f:
      print("".join(KPOINTS), file=f)
    pass
