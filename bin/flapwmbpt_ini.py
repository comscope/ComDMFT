#!/usr/bin/env python
# A Python script to take .CIF files and generate inputs for comsuite.
#
# - We use the CifFile library
#
import re
import os,sys
import itertools
import numpy as np
from math import ceil, floor, acos, sqrt
import pymatgen as mg
import pymatgen.symmetry.analyzer as pgsa

# sys.path.insert(0, '/global/homes/s/sang/usr/dmft_matdelab_sang_dev/bin/flapwmbpt_input')
import flapwmbpt_input.crystalstructure as fm
from flapwmbpt_input import atom_info
from flapwmbpt_input import flapwmbpt_input_sub
from random import random
from random import seed
import flapwmbpt_input.common as fmc
import json

bohr=0.529177210903

#versionstr = '%(prog)s version 0.0' #\npymatgen version '+mg.__version__
# In this script we convert all quantities to Rydberg atomic units
# (rather than Hartree atomic units).
#
rmt_reduce    = 0.97 # scale the muffin tin radius down to avoid overlaps
regname = re.compile("[a-z,A-Z]+")
error = None

chemical_name =atom_info.chemical_name
chemical_symb =atom_info.chemical_symb
chemical_chrg =atom_info.chemical_chrg
element_rmt   =atom_info.element_rmt
element_lmb   =atom_info.element_lmb
element_ntle  =atom_info.element_ntle
element_augm  =atom_info.element_augm
element_atocc =atom_info.element_atocc
element_ptnl  =atom_info.element_ptnl
element_idmd  =atom_info.element_idmd


def supercell(xtalstr_in, ini0):

    species_primitive=xtalstr_in.species
    islist_primitive=xtalstr_in.equivalent_atoms
    real_lattice_primitive=xtalstr_in.real_lattice
    # print(np.transpose(real_lattice_primitive))
    
    str_supercell=xtalstr_in.str*ini0["supercell"]

    seed(1)

    nspecies=len(str_supercell.species)
    for ispecies in range(nspecies):
        orig_species=str_supercell.species[ispecies]        
        str_supercell.replace(ispecies,mg.Species(orig_species,oxidation_state=random()-0.5))
    str_supercell=pgsa.SpacegroupAnalyzer(str_supercell).get_primitive_standard_structure()
    sga_supercell=pgsa.SpacegroupAnalyzer(str_supercell)

    xtalstr_supercell=fm.crystalstructure(str_supercell)
    # print(xtalstr_supercell.frac_coords)
    # print(np.dot(xtalstr_supercell.real_lattice,xtalstr_supercell.frac_coords))        
    
    if not ini0['magn']:
        if (xtalstr_supercell.nsymop*ini0["supercell"][0]*ini0["supercell"][1]*ini0["supercell"][2]>48):
            print('supercell construction failed')
            sys.exit()
        else:
            islist=[]
            specieslist=[]
            for ii in range(len(species_primitive)):
                for jj in range(ini0["supercell"][0]*ini0["supercell"][1]*ini0["supercell"][2]):
                    islist.append(islist_primitive[ii])
                    specieslist.append(species_primitive[ii])
                
            xtalstr_supercell.equivalent_atoms=np.array(islist)
            xtalstr_supercell.species=specieslist
            nsymop_old=xtalstr_supercell.nsymop
            xtalstr_supercell.nsymop=nsymop_old*ini0["supercell"][0]*ini0["supercell"][1]*ini0["supercell"][2]            
            new_rotation=np.zeros((xtalstr_supercell.nsymop, 3, 3))
            new_translation=np.zeros((xtalstr_supercell.nsymop, 3))
            cnt=0
            for ii in range(ini0["supercell"][0]*ini0["supercell"][1]*ini0["supercell"][2]):
                for isym in range(nsymop_old):                            
                    new_rotation[cnt,:,:]=xtalstr_supercell.rotation[isym,:,:]
                    shiftvec=xtalstr_supercell.frac_coords[:,ii]-xtalstr_supercell.frac_coords[:,0]
                    new_translation[cnt,:]=xtalstr_supercell.translation[isym,:]+shiftvec
                    cnt=cnt+1
            xtalstr_supercell.rotation=new_rotation
            xtalstr_supercell.translation=new_translation
            # print(cnt)

            xtalstr_supercell.write_symmetry_operation_input()
            
        
    return xtalstr_supercell


    
    # self.equivalent_atoms
    # self.species
    # self.write_symmetry_operation_input()

    # self.rotation=np.array(self.sga.get_symmetry_dataset()['rotations'])
    # self.translation=np.array(self.sga.get_symmetry_dataset()['translations'])

    
def rad_inscribed_sphere(vec_mat):
    v_temp=np.cross(vec_mat[:,0], vec_mat[:,1])
    r=abs(np.dot(v_temp,vec_mat[:,2]))/np.sqrt(sum(v_temp**2))
    v_temp=np.cross(vec_mat[:,1], vec_mat[:,2])
    r=min(r,abs(np.dot(v_temp,vec_mat[:,0]))/np.sqrt(sum(v_temp**2)))    
    v_temp=np.cross(vec_mat[:,2], vec_mat[:,0])
    r=min(r,abs(np.dot(v_temp,vec_mat[:,1]))/np.sqrt(sum(v_temp**2)))
    r=r/2.0
    return r

def cif2float(cifnum):
    """
    Convert a cif-floating point number that may include an uncertainty
    indication to a proper floating point number. 
    In a .cif file the value "0.4254(4)" is a floating point number where
    the digit in brackets gives the uncertainty. To convert this number to
    a regular Python floating point number the uncertainty needs to be 
    eliminated and the resulting string converted.
    """
    ii = cifnum.find("(")
    if ii >= 0:
      pnum = float(cifnum[:ii])
    else:
      pnum = float(cifnum)
    return pnum

def cif2element(label):
    """
    Convert a label for an atom to the corresponding chemical symbol
    of the element. Examples of labels are "Al1", "Al3+", "boron2a",
    "H*251", etc.
    This algorithm could certainly be implemented more efficiently.
    """
    tag = label.lower()
    tag = re.match(regname,tag)
    symbol = None
    if tag:
        tag = tag.string[:tag.end()]
        for ii in range(1,119):
            if tag == chemical_name[ii] or tag == chemical_symb[ii]:
                symbol = chemical_symb[ii]
                break
    if symbol == None:
        error = "Unknown element: "+tag
        print(error)
    return symbol

def translate_elements(ini0):
    """
    CIF files may specify elements in funny ways. The structure 
    component just extracts whatever the CIF file contains but this
    may not be suitable for any program. This routine translates the
    elements from the way the CIF file specifies them into chemical
    symbols for the elements.
    The updated dictionary is returned.
    """
    elem_in = ini0["symbol"]
    elem_out = []
    for elem in elem_in:
        elem_out.append(cif2element(elem.value))
    ini0["symbol"] = elem_out
    return ini0

def rkm_fact(ielm):
    """
    Following the Wien2K scheme for the muffin-tin radius adjustment.
    These factors are used to compute the relative size of two atoms.
    I.e. for a pair of atoms of ielm and jelm separated by a distance
    D (and Radius(ielm)+Radius(jelm) > D) then the new radii are 
    computed as

      Radius(ielm) = 0.5*(1+(rkm_fact(ielm)-rkm_fact(jelm))*0.005)*D
      Radius(jelm) = 0.5*(1-(rkm_fact(ielm)-rkm_fact(jelm))*0.005)*D

    This function returns the rkm_fact factors that Wien2K uses.
    See Wien2K setrmt_lapw for details.
    """
    if ielm == 3 or ielm == 13 or ielm == 14:
      # Li, Al, Si
      return 45.0
    elif ielm == 4 or ielm == 5:
      # Be, B
      return 50.0
    elif ielm == 6 or ielm == 15:
      # C, P
      return 55.0
    elif ielm == 7 or ielm == 16:
      # N, S
      return 60.0
    elif (ielm == 8 or (ielm >= 11 and ielm <= 13) or ielm == 17 or 
          ielm == 19 or ielm == 20 or ielm == 37 or ielm == 38 or
          ielm == 55 or ielm == 56):
      # O, Na, Mg, Cl, K, Ca, Rb, Sr, Cs, Ba
      return 65.0
    elif ielm == 9:
      # F
      return 70.0
    elif ((ielm >= 21 and ielm <= 24) or (ielm >= 39 and ielm <= 42) or
          (ielm >= 31 and ielm <= 35)):
      # Sc-Cr, Ga-Br, Y-Mo
      return 75.0
    elif ((ielm >= 25 and ielm <= 30) or (ielm >= 44 and ielm <= 53) or 
          ielm == 57 or ielm == 58 or (ielm >= 72 and ielm <= 75)):
      # Mn-Zn, Ru-I, La, Ce, Hf-Re
      return 80.0
    elif ((ielm >= 59 and ielm <= 71) or (ielm >= 76 and ielm <= 85) or
          (ielm >= 87 and ielm <= 118)):
      # Pr-Lu, Os-At, Fr-Og
      return 85.0
    else:
      return 60.0

# def establish_mt_radii(ini_struct):
#     """
#     Takes the elements of the sites, the radii stored in element_rmt,
#     and the distance matrix and generates muffin tin radii that are
#     compatible with the current material structure.
#     """
#     element_rad = {}
#     elements = ini_struct["symbol"]
#     distmat  = ini_struct["dist_matrix"]
#     a        = ini_struct["a"]
#     b        = ini_struct["b"]
#     c        = ini_struct["c"]
#     lena     = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#     lenb     = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
#     lenc     = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])
#     minlen   = min(lena,lenb,lenc)
#     #
#     # Initialize the muffin tin radii taking the lattice vectors into account
#     #
#     for ii in range(0,len(elements)):
#         eli = elements[ii]
#         chg = chemical_chrg[eli]
#         rmi = element_rmt[chg]
#         # the muffin tin radius must be smaller than half the shortest lattice vector
#         element_rad[chg] = min(rmi,0.5*minlen)
#     #
#     # Adjust the muffin tin radii based on the interatomic distances
#     #
#     for ii in range(0,len(elements)):
#         row = distmat[ii,:]
#         eli = elements[ii]
#         chg = chemical_chrg[eli]
#         for jj in range(0,len(elements)):
#             if ii < jj and row[jj] < 1.0e-6:
#                 print("ERROR: atoms ",ii+1, " and ",jj+1, " are in the same position!")
#                 print("Atom ",ii+1,":")
#                 print("  Element    : ",eli)
#                 print("  Coordinates: ",ini_struct["xcoord"][ii],
#                                         ini_struct["ycoord"][ii],
#                                         ini_struct["zcoord"][ii])
#                 print("Atom ",jj+1,":")
#                 print("  Element    : ",elj)
#                 print("  Coordinates: ",ini_struct["xcoord"][jj],
#                                         ini_struct["ycoord"][jj],
#                                         ini_struct["zcoord"][jj])
#             if ii != jj and row[jj] > 0.0:
#                 rmi = element_rmt[chg]
#                 rr  = row[jj]
#                 elj = elements[jj]
#                 nmj = chemical_chrg[elj]
#                 rmj = element_rmt[nmj]
#                 if rmi+rmj > rr:
#                     fi = rkm_fact(chg)
#                     fj = rkm_fact(nmj)
#                     rmi = 0.5*(1.0+(fi-fj)*0.005)*rr
#                     rmj = 0.5*(1.0-(fi-fj)*0.005)*rr
#                     scale = rr/(rmi+rmj)
#                 if chg in element_rad:
#                     element_rad[chg] = min(element_rad[chg],rmi)
#                 else:
#                     element_rad[chg] = rmi
#                 if nmj in element_rad:
#                     element_rad[nmj] = min(element_rad[nmj],rmj)
#                 else:
#                     element_rad[nmj] = rmj


#     # # scale up

#     # scalefac=0
#     # for ii in range(len(elements)):
#     #     eli = elements[ii]
#     #     chg = chemical_chrg[eli]        
#     #     for jj in range(len(elements)):
#     #         elj = elements[jj]            
#     #         nmj = chemical_chrg[elj]            
#     #         scaletemp=distmat[ii,jj]/(element_rad[chg]+element_rad[nmj])
#     #         scalefac=max(scaletemp
            
#     ini_struct["element_rad"] = element_rad
#     # print('rad', eli,element_rad) 
#     return ini_struct

def establish_mt_radii(ini0):
    """
    Takes the elements of the sites, the radii stored in element_rmt,
    and the distance matrix and generates muffin tin radii that are
    compatible with the current material structure.
    """
    element_rad = {}
    elements = ini0["symbol"]
    distmat  = ini0["dist_matrix"]
    a        = ini0["a"]
    b        = ini0["b"]
    c        = ini0["c"]
    lena     = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
    lenb     = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])
    lenc     = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])
    minlen   = min(lena,lenb,lenc)
    #
    # Initialize the muffin tin radii taking the lattice vectors into account
    #
    # print(element_rmt)
    # print(distmat)
    for ii in range(0,len(elements)):

        eli = elements[ii]
        chg = chemical_chrg[eli]
        # print(element_rmt[chg], distmat[ii,ii]/2)        
        element_rad[chg] = max(element_rmt[chg], distmat[ii,ii]/2)        
    #
    # Adjust the muffin tin radii based on the interatomic distances
    #
    # print(element_rad)
    for ii in range(len(elements)):
        eli = elements[ii]
        chg = chemical_chrg[eli]
        rmi = element_rad[chg]
        for jj in range(len(elements)):
            elj = elements[jj]
            nmj = chemical_chrg[elj]
            rmj = element_rad[nmj]
            rr  = distmat[ii,jj]                      
            if ((rmi+rmj > rr) & (ii!=jj)):
                fi = rkm_fact(chg)
                fj = rkm_fact(nmj)
                rmi = 0.5*(1.0+(fi-fj)*0.005)*rr
                rmj = 0.5*(1.0-(fi-fj)*0.005)*rr
                scale = rr/(rmi+rmj)
                element_rad[chg] = min(element_rad[chg],rmi)
                element_rad[nmj] = min(element_rad[nmj],rmj)
                # print(ii,jj,element_rad)
    ini0["element_rad"] = element_rad

    return ini0

def volume_sphere(radius):
    """
    Return the volume of a sphere given its radius

    The volume of a sphere is 4/3*pi*(r^3)
    """
    pi = acos(-1.0)
    vol = 4.0/3.0*pi*(radius**3)
    return vol

def establish_atoms_volume(ini0):
    """
    Calculate the total volume of the atoms within the muffin tin radii

    For all atoms in the unit cell calculate the volume of each atom
    and add all volumes up. Add the total volume to the ini0 
    dictionary.
    """
    elements = ini0["symbol"]
    radii    = ini0["element_rad"]
    atoms_volume = 0.0
    for element in elements:
        number = chemical_chrg[element]
        radius = radii[number]
        atoms_volume += volume_sphere(radius)
    ini0["atoms_volume"] = atoms_volume
    return ini0

def establish_Kmin(ini0):
    """
    Establish the minimum plane wave cut-off

    From the relation RK=2*l where R is the muffin tin radius, K
    is the plane wave cut-off and l is the highest angular momentum
    of the occupied orbitals of an atom we can find K=2*l/R and use
    that to establish the smallest plane wave cut-off that still
    provides a qualitatively correct description of the system.

    # In addition Klim is established which derives from the fact that 
    # the GW code currently cannot handle l quantum numbers larger than
    # 10. Hence Klim is the minimum value of 10/R across all elements.

    # The values of Kmin and Klim are returned in the ini0
    # dictionary.

    #References#

    [1] D.J. Singh, "Planewaves, pseudo-potentials and the LAPW method",
        Springer (1994), ISBN: 978-1-4757-2314-4, 
        DOI: 10.1007/978-1-4757-2312-0, pp. 62-63.
    """
    element_rad = ini0["element_rad"]
    elements    = ini0["symbol"]
    Kmin        = 0.00
    for ii in range(0,len(elements)):
        eli  = elements[ii]
        chg  = chemical_chrg[eli]
        # print(lmax_ef(chg))
        # ll   = (lmax_ef(chg))*2+1
        # ll   = max(lmax_ef(chg)*2, 4)
        ll   = lmax_ef(chg)*2+ini0['mt_orb_l']
        rr   = ini0["element_rad"][chg]
        Kmin = max(Kmin,ll/rr)
    ini0["Kmin"] = Kmin
    return ini0


def establish_Kmax(ini0,Kmax=0):

    element_rad = ini0["element_rad"]
    elements    = ini0["symbol"]
    # print(ini0["l_max"][83])
    for ii in range(0,len(elements)):
        eli  = elements[ii]
        chg  = chemical_chrg[eli]
        ll   = ini0["l_max"][chg]                
        rr   = ini0["element_rad"][chg]
        Kmax = max(Kmax,ll/rr)
    ini0["Kmax"] = Kmax
    return ini0


def establish_r_grid(ini0):
    """
    Establish the real space grid

    Given the K_max and the lattice vectors work out how many grid
    points are needed in each dimension.

    The number of grid points in all dimensions is returned in
    ini0 under key "nrdiv". Also returned is "mdiv" which
    is about 4/3*nrdiv.

    In addition this function takes the packing factor of the
    material into account. For packing factors larger than 0.5
    the regular approach to calculating the grid sizes is sufficient.
    For packing factors smaller than 0.3 about twice as many plane
    waves are needed to represent the interstitial region well. 
    Between 0.3 and 0.5 the grid sizes are scaled linearly with
    a factor 2 to 1.
    """
    Vcel  = ini0["cell_volume"]
    Vatm  = ini0["atoms_volume"]
    SpGr  = ini0["spacegroup"]
    pfac  = Vatm/Vcel  # packing factor
    # print(Vcel, Vatm,SpGr, pfac)

    scale_fcc = 1.9
    scale_03 = 2.2

    # scale_fcc = 1.8
    # scale_03 = 2.0
    
    scale = (scale_fcc+(0.74-pfac)/0.44*(scale_03-scale_fcc))

    print('pfac scale', pfac, scale)        
    pi    = acos(-1.0)
    r43   = 4.0/3.0
    Kmax  = ini0["Kmax"]

    a     = np.array(ini0["a*"])
    b     = np.array(ini0["b*"])
    c     = np.array(ini0["c*"])


    ii_max, jj_max, kk_max=flapwmbpt_input_sub.real_grid(a,b,c,Kmax*scale,ini0["cut_lapw_ratio"])

    maxvalue=max(ii_max,jj_max,kk_max)
    nra   = maxvalue
    nrb   = maxvalue
    nrc   = maxvalue

    tempmat=np.zeros((3,3))
    tempmat[:,0]=a*nra
    tempmat[:,1]=b*nrb
    tempmat[:,2]=c*nrc    
    cutoff=rad_inscribed_sphere(tempmat)*ini0["cut_lapw_ratio"]



    element_rad = ini0["element_rad"]
    elements    = ini0["symbol"]

    rkmax_ratio=0.0    
    for ii in range(0,len(elements)):
        eli  = elements[ii]
        chg  = chemical_chrg[eli]
        ll   = ini0["l_max"][chg]                
        rr   = ini0["element_rad"][chg]
        rkmax_ratio=max(rkmax_ratio, rr*cutoff/ll)

    if (rkmax_ratio > ini0['rklmax']):
        reducefac=ini0['rklmax']/rkmax_ratio
    else:
        reducefac=1.0

        
    nra   = int(round(maxvalue*reducefac))
    nrb   = int(round(maxvalue*reducefac))
    nrc   = int(round(maxvalue*reducefac))

    mra   = int(round(maxvalue*reducefac*2.0))
    mrb   = int(round(maxvalue*reducefac*2.0))
    mrc   = int(round(maxvalue*reducefac*2.0))
    

    tempmat=np.zeros((3,3))
    tempmat[:,0]=a*nra
    tempmat[:,1]=b*nrb
    tempmat[:,2]=c*nrc    
    cutoff=rad_inscribed_sphere(tempmat)*ini0["cut_lapw_ratio"]
    tempmat=np.zeros((3,3))
    tempmat[:,0]=a*mra
    tempmat[:,1]=b*mrb
    tempmat[:,2]=c*mrc        
    cutoffro=rad_inscribed_sphere(tempmat)*ini0["cut_pb_ratio"]*2.0        
        
    mdiv  = []
    nrdiv = []

    numplw=4.0*pi/3.0*(cutoffro+2*cutoff)**3/abs(np.dot(np.cross(a,b), c))

    mdiv.append(mra)
    mdiv.append(mrb)
    mdiv.append(mrc)
    nrdiv.append(nra)
    nrdiv.append(nrb)
    nrdiv.append(nrc)
    ini0["mdiv"]  = mdiv
    ini0["nrdiv"] = nrdiv
    # rkmax==l test



    for ii in range(0,len(elements)):
        eli  = elements[ii]
        chg  = chemical_chrg[eli]
        ll   = ini0["l_max"][chg]                
        rr   = ini0["element_rad"][chg]
        print(eli, 'rkmax/l=', rr*cutoff/ll)
    
    return ini0

def establish_lmax(ini0):
    """
    For every element establish the maximum angular momentum to be used

    From the relation RK=l_max [1] where R is the muffin tin radius, K
    is the plane wave cut-off and l_max is the highest angular momentum
    of the orbitals of an atom we can find l_max for each element.
    The l_max values are stored in a dictionary.

    The values of l_max are returned in the ini0 dictionary.

    #References#

    [1] D.J. Singh, L. Nordstrom, "Planewaves, pseudo-potentials and
        the LAPW method", Springer (1994), ISBN: 978-0-387-29684-5,
        pp. 62-63.
    """
    element_rad = ini0["element_rad"]
    elements    = ini0["symbol"]
    Kmin        = ini0["Kmin"]
    l_max       = {}

    for ii in range(0,len(elements)):
        eli        = elements[ii]
        chg        = chemical_chrg[eli]
        rr         = element_rad[chg]
        l_max[chg] = int(round(rr*Kmin))
    ini0["l_max"] = l_max
    # print('l_max', l_max)
    return ini0


def mt_expanding_basis_n(chg,l):
    if ((chg>=81) & (chg<=88)):
        if (l==3):
            return l+2 # due to too deep 4f states
        else:
            return l+1
    else:
        return l+1

    

def expand_atomic_basis(ini0):
    """
    For every element (atomic type) expand the basis up to the given l_max

    For every atomic type, using the valence basis set information and the
    maximum angular momentum, expand the stored basis set information to
    the full basis set. This includes adding missing angular momenta but
    also for the lower angular momenta add more radial functions up to 
    the appropriate N-quantum number.

    The corresponding dictionaries are returned in ini0. The names
    of the dictionaries are elmnt_exp_* corresponding to the global names
    element_* from which they were derived.
    """

    elements        = ini0["symbol"]
    l_max           = ini0["l_max"]

    elmnt_exp_lmb={}
    elmnt_exp_ntle={}
    elmnt_exp_augm={}
    elmnt_exp_atocc={}
    elmnt_exp_ptnl={}
    elmnt_exp_idmd={}


    
    for ii in range(0,len(elements)):

        augm       = []
        atocc      = []
        ptnl       = []
        idmd       = []
        ntle       = []

        eli        = elements[ii]
        chg        = chemical_chrg[eli]
        lmax       = l_max[chg]
        # print('lmax_info', lmax, element_lmb[chg])

        if (lmax > element_lmb[chg]):
            nmax=np.max([np.max([int(i) for i in element_ptnl[chg]]), mt_expanding_basis_n(chg,lmax)+1])
            if (lmax <= 2):
                nmax=nmax+1
            nmax=nmax+ini0['mt_orb_n']

            # adding more radial function
            for ll in range(element_lmb[chg]+1):
                kstart=int(np.sum(element_ntle[chg][:ll]))
                kend=kstart+element_ntle[chg][ll]-1
                nmax_l=int(element_ptnl[chg][kend])
                ntle.append(element_ntle[chg][ll]+nmax-nmax_l)            
                
                for kk in range(kstart, kend+1):
                    augm.append(element_augm[chg][kk])
                    atocc.append(element_atocc[chg][kk])
                    ptnl.append(element_ptnl[chg][kk])
                    idmd.append(element_idmd[chg][kk])
                        
                for nind in range(nmax-nmax_l):
                    augm.append("LOC")
                    atocc.append(0.0)
                    ptnl.append(nmax_l+nind+1+0.8)
                    idmd.append(1)

            # adding higher angular momentum function
            for ll in range(element_lmb[chg]+1, lmax+1):
                ntle.append(nmax-mt_expanding_basis_n(chg,ll)+1)
                augm.append("APW")
                atocc.append(0.0)
                ptnl.append(mt_expanding_basis_n(chg,ll)+0.8)
                idmd.append(0)
                
                for nind in range(mt_expanding_basis_n(chg,ll)+1, nmax+1):
                    augm.append("LOC")
                    atocc.append(0.0)
                    ptnl.append(nind+0.8)
                    idmd.append(1)

        else:
            nmax=0
            for ll in range(lmax+1):
                kstart=int(np.sum(element_ntle[chg][:ll]))
                kend=kstart+element_ntle[chg][ll]-1
                if (element_augm[chg][kend]=='APW'):
                    nmax_l=int(element_ptnl[chg][kend])+1
                else:
                    nmax_l=int(element_ptnl[chg][kend])                    
                nmax=np.max([nmax, nmax_l])

            if (lmax <= 2):
                nmax=nmax+1

            nmax=nmax+ini0['mt_orb_n']
            
            # print('nmax',nmax)
            # adding more radial function
            for ll in range(lmax+1):
                kstart=int(np.sum(element_ntle[chg][:ll]))
                kend=kstart+element_ntle[chg][ll]-1
                nmax_l=int(element_ptnl[chg][kend])
                nradf=element_ntle[chg][ll]+nmax-nmax_l                
                ntle.append(nradf)
                # print('nradf',ll, nradf)
                # print('nmax_l',ll, nmax,nmax_l)

                if (nmax<nmax_l):
                    for kk in range(kstart, kend+1+nmax-nmax_l):
                        augm.append(element_augm[chg][kk])
                        atocc.append(element_atocc[chg][kk])
                        ptnl.append(element_ptnl[chg][kk])
                        idmd.append(element_idmd[chg][kk])
                else:
                    for kk in range(kstart, kend+1):
                        # print(int(element_ptnl[chg][kk]))
                        augm.append(element_augm[chg][kk])
                        atocc.append(element_atocc[chg][kk])
                        ptnl.append(element_ptnl[chg][kk])
                        idmd.append(element_idmd[chg][kk])
                    for nind in range(nmax-nmax_l):
                        # print(nmax_l+nind+1)
                        augm.append("LOC")
                        atocc.append(0.0)
                        ptnl.append(nmax_l+nind+1+0.8)
                        idmd.append(1)            
                    

        # print(lmax)
        # print(ntle)
        # print(augm)
        # print(atocc)
        # print(ptnl)
        # print(idmd)
        elmnt_exp_lmb[chg]   = lmax
        elmnt_exp_ntle[chg]  = ntle
        elmnt_exp_augm[chg]  = augm
        elmnt_exp_atocc[chg] = atocc
        elmnt_exp_ptnl[chg]  = ptnl
        elmnt_exp_idmd[chg]  = idmd
    ini0["element_lmb"]   = elmnt_exp_lmb
    ini0["element_ntle"]  = elmnt_exp_ntle
    ini0["element_augm"]  = elmnt_exp_augm
    ini0["element_atocc"] = elmnt_exp_atocc
    ini0["element_ptnl"]  = elmnt_exp_ptnl
    ini0["element_idmd"]  = elmnt_exp_idmd
    return ini0

def write_inifile(ini0,inifile):
    """
    Take a dictionary with all the relevant ini settings and write an
    input file for the GW+DMFT code.
    """
    inifile.write("TEXT band structure calculation\n")
    inifile.write("CONTROL   iter_dft=%3i  iter_hf=%3i  iter_gw=%3i iter_qp=%3i\n" % 
                  (ini0["iter_dft"],ini0["iter_hf"],
                   ini0["iter_gw"],ini0["iter_qp"]))
    inifile.write("          admix=%5.3f  adspin=%5.3f adm_gw=%5.3f acc_it_gw=%5.3f\n" %(ini0['admix'],ini0['adspin'], ini0['adm_gw'], ini0['acc_it_gw']))
    inifile.write("          iexch=005 scal_spin= 1.0000\n")
    inifile.write("          nproc_tau= %3i nproc_k= %3i\n" %(ini0['nproc_tau_flapwmbpt'], ini0['nproc_k_flapwmbpt']))
    inifile.write("          irel=%1i clight=274.074e+00 rel_interst=F irel_core=%1i\n" %(ini0['irel'], ini0['irel_core']))
    if ini0["restart"]:
        restart="T"
    else:
        restart="F"
    inifile.write("          temperature=%10.2f  restart=%s\n" % (ini0["temperature"],restart))
    inifile.write("FILES\n")
    inifile.write("  allfile=com\n")
    inifile.write("SYM symgen=input\n")
    nsort = np.max(ini0["islist"])
    inifile.write("STRUCTURE  par=%11.7f  natom=%3d nsort=%3d istruct=%3d\n" %
                  (ini0["par"],ini0["natom"],nsort,ini0["istruc"]))
    inifile.write("      is=")
    islist = ini0["islist"]
    for isnum in islist:
        inifile.write("%3d" % isnum)
    inifile.write("\n")
    inifile.write("      b/a=%9.6f  c/a=%9.6f\n" % (ini0["b_a"],ini0["c_a"]))
    inifile.write("      a=%21.16f%21.16f%21.16f\n" %(ini0["a"][0],ini0["a"][1],ini0["a"][2]))
    inifile.write("      b=%21.16f%21.16f%21.16f\n" %(ini0["b"][0],ini0["b"][1],ini0["b"][2]))
    inifile.write("      c=%21.16f%21.16f%21.16f\n" %(ini0["c"][0],ini0["c"][1],ini0["c"][2]))
    natom = ini0["natom"]
    for ii in range(0,natom):
        inifile.write("    tau=%21.16f%21.16f%21.16f\n" %
                      (ini0["a_coeff"][ii],ini0["b_coeff"][ii],ini0["c_coeff"][ii]))
    mdiv  = ini0["mdiv"]
    nrdiv = ini0["nrdiv"]
    inifile.write("REAL SPACE MESHES mdiv=%4d %4d %4d\n" % (mdiv[0], mdiv[1], mdiv[2]))
    inifile.write("                 nrdiv=%4d %4d %4d\n" % (nrdiv[0],nrdiv[1],nrdiv[2]))
    inifile.write("BASIS  cut_lapw_ratio=%4.3f cut_pb_ratio=%4.3f\n" %(ini0["cut_lapw_ratio"], ini0["cut_pb_ratio"]))
    inifile.write("       eps_pb=1.e-03\n")
    inifile.write("ZONES nbndf=   0\n")
    if ini0["band"]:
        band = ' T'
    else:
        band = ' F'
    if ini0["dos"]:
        dos = ' T'
    else:
        dos = ' F'
    inifile.write("DOS   emindos=-15.000  emaxdos= 15.000   ndos= 3000\n")
    inifile.write("      n_cont_frac=  30 e_small=2.e-02\n")
    inifile.write("      dos=%s           bandstructure=%s\n" % (dos,band))
    inifile.write("K_POINT  ndiv=%3d %3d %3d  metal=T n_k_div= 27 k_line=010\n" %(ini0['kmesh'][0], ini0['kmesh'][1], ini0['kmesh'][2]))
    inifile.write("MULTI_SCF vv0=  1.00\n")
    inifile.write("MAGNET  b_extval=   0.000000 iter_h_ext=%s\n" %(ini0['iter_h_ext']))
    inifile.write("        b_ext=  0.000  0.000  1.000\n")
    inifile.write("TAU MESH n_tau=   46 n_tau_int= 1200\n")
    inifile.write("OMEGA MESH n_omega_exa=   29 n_omega_asy=   18 omega_max=  200.00 \n")
    inifile.write("           interp_omega_d= 2\n")
    inifile.write("NU MESH n_nu_exa=   29 n_nu_asy=   18 nu_max=  200.00\n")
    inifile.write("        interp_nu_d= 2\n")
    inifile.write("ATOMIC DATA --------------------------------------------------------\n")
    element_rad = ini0["element_rad"]
    isold = 0
    for ii in range(0,natom):
        isnew = ini0["islist"][ii]
        if isold != isnew:
            symbol = ini0["symbol"][ii]
            if len(symbol) == 1:
                symb2 = symbol+"  "
            elif len(symbol) == 2:
                symb2 = symbol+" "
            elif len(symbol) == 3:
                symb2 = symbol
            else:
                error = "Strange chemical symbol:"+symbol
                print(error)
            number = chemical_chrg[symbol]
            smt    = element_rad[number]*rmt_reduce
            inifile.write("  txtel=%s   z=%5.1f magn_shift= %5.3f\n" %
                          (symb2,float(number), ini0['magn_shift']))
            inifile.write("  smt=%8.5f h= 0.0120 nrad= 1216 z_dop=0.000\n" % 
                          smt)
            lmb = ini0["element_lmb"][number]
            # lmpb is the maximum l-quantum number for the product basis
            # we set this to the same value as lmb for now...
            lmpb = min(ini0["element_lmb"][number]+ini0['mt_pb_l'], ini0["mt_pb_l_max"])
            inifile.write("  lmb=%2d  lmpb=%2d\n" % (lmb,lmpb))
            ntle = ini0["element_ntle"][number]
            inifile.write("  lim_pb_mt=")
            for ii in range(0,lmpb+1):
                inifile.write("%3d" % 30)
            inifile.write("\n")
            inifile.write("  ntle=")
            for ii in range(0,len(ntle)):
                inifile.write("%3d" % ntle[ii])
            inifile.write("\n")
            inifile.write("  l  augm  atocc   ptnl  corr idmd\n")
            kk = 0
            for ii in range(0,len(ntle)):
                ntlen = ntle[ii]
                for jj in range(0,ntlen):
                    inifile.write("%3d   %s%7.3f%7.3f     %s    %1d\n" %
                                  (ii,ini0["element_augm"][number][kk],
                                   ini0["element_atocc"][number][kk],
                                   ini0["element_ptnl"][number][kk],"N",
                                   ini0["element_idmd"][number][kk]))
                    kk+=1
        isold = isnew

def write_kpathfile(ini0):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpath, and write the data to the
    kpathfile.
    """

    kpathfile = open("kpath",'w')

    length = len(ini0["kpath_label"])
    print('kpath info')
    for ii in range(length-1):
        kpathfile.write("%5s %12.8f %12.8f %12.8f       " %(ini0["kpath_label"][ii],ini0["kpath"][ii][0],ini0["kpath"][ii][1],ini0["kpath"][ii][2]))
        kpathfile.write("%5s %12.8f %12.8f %12.8f\n" %(ini0["kpath_label"][ii+1],ini0["kpath"][ii+1][0],ini0["kpath"][ii+1][1],ini0["kpath"][ii+1][2]))
        print(ini0["kpath_label"][ii], ini0["kpath_label"][ii+1], np.sqrt(np.sum((ini0["a*"]*(ini0["kpath"][ii+1][0]-ini0["kpath"][ii][0])+ini0["b*"]*(ini0["kpath"][ii+1][1]-ini0["kpath"][ii][1])+ini0["c*"]*(ini0["kpath"][ii+1][2]-ini0["kpath"][ii][2]))**2)))
    kpathfile.close()
    return None

def write_plot1dfile(ini0,plot1dfile):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpath, and write the data to the
    plot1dfile for the Elk code.
    """
    length = len(ini0["kpath_label"])
    plot1dfile.write("plot1d\n")
    plot1dfile.write("  %d\n" % length)
    plot1dfile.write("  200\n")
    for ii in range(0,length):
        point = ini0['kpath_label'][ii]
        kpoint = ini0['kpath'][ii]
        plot1dfile.write("%12.8f %12.8f %12.8f\n" %
                        (kpoint[0],kpoint[1],kpoint[2]))

def write_kpointfile(ini0):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpoints, and write the data to the
    kpointfile.
    """

    kpointfile=open("kpoints",'w')


    length = len(ini0["kpoint_label"])

    kpointfile.write("\n")
    kpointfile.write("frac\n")
    kpointfile.write("%5i\n" %(length))    
    for ii in range(length):
        kpointfile.write("%5i  %12.8f %12.8f %12.8f   %s\n" %
                         (ii+1, ini0['kpoint'][ii][0], ini0['kpoint'][ii][1], ini0['kpoint'][ii][2], ini0['kpoint_label'][ii]))
    kpointfile.close()
    return None


def write_klistfile(ini0,kpointfile):
    """
    Take a dictionary with all the relevant information for the
    structure, extract the Kpoints, and write the data to the
    Wien2k klistfile.

    The Wien2K klist format is special in that the k-points are
    given in terms of integers instead of fractional coordinates.
    So we need to convert all fractional coordinates of each
    k-point to integers first.
    """
    import fractions

    case_in1 = case+".klist_band"
    klist_file = open(case_in1,"w")
    
    if not ini0["kpoint"]:
        raise ValueError("None is an invalid value")
    kpoint=ini0["kpoint"]
    kpath=ini0["kpoint_label"]
    length = len(kpath)
    for ii in range(0,length):
        point  = str(kpath[ii])
        kpoint = kpoint[ii]
        fracx  = fractions.Fraction.from_float(kpoint[0]).limit_denominator(10000)
        fracy  = fractions.Fraction.from_float(kpoint[1]).limit_denominator(10000)
        fracz  = fractions.Fraction.from_float(kpoint[2]).limit_denominator(10000)
        idv    = fracx.denominator
        if fracy.denominator != fracx.denominator:
            idv *= fracy.denominator
        if (fracz.denominator != fracy.denominator and 
            fracz.denominator != fracx.denominator):
            idv *= fracz.denominator
        isx   = int(kpoint[0]*idv+0.5)
        isy   = int(kpoint[1]*idv+0.5)
        isz   = int(kpoint[2]*idv+0.5)
        kpointfile.write("%-10s%5d%5d%5d%5d%5.2f%5.2f%5.2f%3s\n" %
                         (point,isx,isy,isz,idv,2.0,0.0,0.0,"   "))
    kpointfile.write("END\n")
    klist_file.close()

def any2utf8(infile):
    """
    Convert the encoding of a text file to UTF8 encoding

    CIF files are supposed to contain only ASCII characters but some sources
    generate CIF files with other characters. This may cause problems if 
    these files are encoded in anything else but UTF8. Python and therefore
    Pymatgen expect CIF files to be UTF8 encoded. If they are not the code
    will throw an exception. 

    To avoid that problem this function takes the name of a CIF file. It
    detects the encoding of that file, reads it in using that encoding,
    writes out the contents to a temporary file using UTF8 encoding,
    and returns the name of the temporary file.
    """
    import codecs
    import os
    from chardet.universaldetector import UniversalDetector
    #
    # Generate a unique filename
    #
    pid = os.getpid()
    outfile = "/tmp/tmp"+str(pid)+".cif"
    #
    # Detect the encoding of the input file
    #
    detector = UniversalDetector()
    detector.reset()
    file = open(infile,'rb')
    for line in file:
        detector.feed(line)
        if detector.done:
            break
    detector.close()
    file.close()
    #
    # Read the input file and decode the contents
    #
    codec = detector.result['encoding']
    f = codecs.open(infile,mode='r',encoding=codec)
    contents_in = f.readlines()
    f.close()
    #
    # Map all characters to ASCII characters (non ASCII characters are skipped)
    #
    contents_out = []
    for line in contents_in:
        line_out = ""
        for char in line:
            if ord(char) <= 127:
                line_out = line_out + char
        contents_out.append(line_out)
    #
    # Write the contents to the temporary file UTF8 encoded
    #
    f = codecs.open(outfile,mode='w',encoding='utf8')
    f.writelines(contents_out)
    f.close()
    #
    # Return the temporary filename
    #
    return outfile

def retr_lattice_vecs(xtalstr,ini0):
    """
    Retrieve the lattice vectors.

    The lattice vectors are extracted from the structure instance.
    They are added to the dictionary passed in. The updated dictionary
    is returned. The lattice vectors are given in Bohr.
    """
    # print(xtalstr.real_lattice)
    ini0["a"] = xtalstr.real_lattice[:,0]
    ini0["b"] = xtalstr.real_lattice[:,1]
    ini0["c"] = xtalstr.real_lattice[:,2]
    return ini0

def retr_recip_vecs(xtalstr,ini0):
    """
    Retrieve the lattice vectors.

    The lattice vectors are extracted from the structure instance.
    They are added to the dictionary passed in. The updated dictionary
    is returned. The lattice vectors are given in Bohr.
    """
    ini0["a*"] = xtalstr.recip_lattice[:,0]
    ini0["b*"] = xtalstr.recip_lattice[:,1]
    ini0["c*"] = xtalstr.recip_lattice[:,2]
    return ini0

def retr_sites(xtalstr,ini0):
    """
    Retrieve the atomic positions and elements in fractional coordinates.

    The atomic positions and elements are extracted from the structures
    instance. They are added to the dictionary passed in. The updated dictionary
    is returned.
    """
    ini0["natom"]  = xtalstr.natom
    # print(xtalstr.frac_coords)    
    ini0["a_coeff"] = xtalstr.frac_coords[0,:]
    ini0["b_coeff"] = xtalstr.frac_coords[1,:]
    ini0["c_coeff"] = xtalstr.frac_coords[2,:]
    # print(ini0["a_coeff"])
    # print(ini0["b_coeff"])
    # print(ini0["c_coeff"])
    ini0["islist"] = xtalstr.equivalent_atoms
    ini0["symbol"] = xtalstr.species
    return ini0

def retr_distance_matrix(xtalstr,ini0):
    """
    Retrieve the distances between sites.
    """

    # raw_mat = xtalstr.dist_mat
    # dist_mat = []
    # for ii in range(0,len(raw_mat)):
    #     raw_row = raw_mat[ii]
    #     dist_row = []
    #     for jj in range(0,len(raw_row)):
    #         dist_row.append(raw_row[jj])
    #     dist_mat.append(dist_row)    
    
    ini0["dist_matrix"] = xtalstr.dist_mat
    # print(type(raw_mat))
    # print(type(dist_mat))
    return ini0

def retr_lattice_type(xtalstr,ini0):
    """
    Retrieve lattice type.

    Lookup the crystal_system, find out whether the structure is
    Bcc, Fcc or something else and lookup the corresponding lattice
    type. Add the lattice type to the dictionary and return the
    updated dictionary.

    The logic is copied from pymatgen.symmetry.bandstructure
    HighSymmKpath. That routine acknowledges:

    W. Setyawan, S. Curtarolo (2010), "High-throughput electronic
    band structure calculations: Challenges and tools".
    Computational Materials Science, 49(2), 299-312. 
    DOI: 10.1016/j.commatsci.2010.05.010
    """
    lattice    = xtalstr.crystal_system
    hall       = xtalstr.spacegroup_hall
    lattice_tp = None
    if   lattice == "cubic":
        if   "P" in hall:
            lattice_tp = "cubic"
        elif "F" in hall:
            lattice_tp = "cubic_fcc"
        elif "I" in hall:
            lattice_tp = "cubic_bcc"
        else:
            print("Unexpected lattice and Hall symbol combination: %s, %s" % (lattice,hall))
            exit(1)
    elif lattice == "tetragonal":
        if   "P" in hall:
            lattice_tp = "tetragonal"
        elif "I" in hall:
            lattice_tp = "tetragonal_bcc"
        else:
            print("Unexpected lattice and Hall symbol combination: %s, %s" % (lattice,hall))
            exit(1)
    elif lattice == "orthorhombic":
        if   "P" in hall:
            lattice_tp = "orthorhombic"
        elif "A" in hall:  # Base centered, no separately recognized
            lattice_tp = "orthorhombic"
        elif "B" in hall:  # Base centered, no separately recognized
            lattice_tp = "orthorhombic"
        elif "C" in hall:  # Base centered, no separately recognized
            lattice_tp = "orthorhombic"
        elif "F" in hall:
            lattice_tp = "orthorhombic_fcc"
        elif "I" in hall:
            lattice_tp = "orthorhombic_bcc"
        else:
            print("Unexpected lattice and Hall symbol combination: %s, %s" % (lattice,hall))
            exit(1)
    elif lattice == "hexagonal":
        lattice_tp = "hexagonal"
    elif lattice == "monoclinic":
        lattice_tp = "monoclinic"
    elif lattice == "trigonal" or lattice == "rhombohedral":
        lattice_tp = "rhombohedral"
    elif lattice == "triclinic":
        lattice_tp = "triclinic"
    else:
        print("Unexpected lattice: %s" % lattice)
        exit(1)
    ini0["hall"]         = hall
    if lattice_tp:
        if (lattice_tp == 'cubic'):
            ini0["istruc"]     = -1
        elif (lattice_tp == 'cubic_bcc'):
            ini0["istruc"]     = -2
        elif (lattice_tp == 'cubic_fcc'):
            ini0["istruc"]     = -3
        elif (lattice_tp == 'hexagonal'):
            ini0["istruc"]     = -4
        elif (lattice_tp == 'tetragonal'):
            ini0["istruc"]     = -5
        elif (lattice_tp == 'tetragonal_bcc'):
            ini0["istruc"]     = -6
        elif (lattice_tp == 'orthorhombic'):
            ini0["istruc"]     = -7
        elif (lattice_tp == 'orthorhombic_bcc'):
            ini0["istruc"]     = -8
        elif (lattice_tp == 'orthorhombic_fcc'):
            ini0["istruc"]     = -9
        elif (lattice_tp == 'monoclinic'):
            ini0["istruc"]     = -10
        elif (lattice_tp == 'rhombohedral'):
            ini0["istruc"]     = -11
        elif (lattice_tp == 'triclinic'):
            ini0["istruc"]     = -13
    else:
      ini0["istruc"]     = 0
    ini0["lattice_type"] = lattice_tp
    return ini0

def retr_cell_volume(xtalstr,ini0):
    """
    Retrieve the cell volume
    """
    ini0["cell_volume"] = xtalstr.volume
    return ini0

def retr_kpath(xtalstr,ini0):
    """
    Retrieve the kpath
    """
    ini0["kpath"] = xtalstr.kpath
    ini0["kpath_label"] = xtalstr.kpath_label    
    return ini0

def retr_kpoint(xtalstr,ini0):
    """
    Retrieve the kpoints
    """
    ini0["kpoint"] = xtalstr.kpoint
    ini0["kpoint_label"] = xtalstr.kpoint_label    
                
    return ini0

def retr_spacegroup_number(xtalstr,ini0):
    """
    Retrieve the space group number.

    When choosing real space grids some space groups imply restrictions
    on the number of points along the different lattice vectors.
    To impose the correct restrictions we need to know the space group
    number of the current system. This function add this number to
    the dictionary.
    """
    ini0["spacegroup"] = xtalstr.spacegroup_number
    return ini0


def write_symmetry_operation_input(nsym, rotation, translation):
    f=open('symmetry_operations.dat', 'w')
    f.write("   number of symmetry operations= %3d\n" % (nsym))
    for ii in range(nsym):    
        f.write("    symmetry operation %3d\n" % (ii+1))
        for jj in range(0,3):
            line = (rotation[ii][jj][0],rotation[ii][jj][1],rotation[ii][jj][2],translation[ii][jj])
            f.write("    (%14.10f,%14.10f,%14.10f)    (%14.10f)\n" % line)    
    f.close()


def write_comsuite(ini):
    """
    Take the input data and generate the input files for the Comsuite of
    programs.
    """
    filename = "ini"
    inifile = open(filename,'w')
    write_inifile(ini,inifile)
    inifile.close()    
    write_kpathfile(ini)
    write_kpointfile(ini)

def unique_species(ini0):
    """
    Return the list of different chemical elements there are in the
    current structure.
    """
    natom = ini0["natom"]
    elmlist = []
    for ii in range(0,natom):
        symbol = ini0["symbol"][ii]
        if not symbol in elmlist:
            elmlist.append(symbol)
    return elmlist

def lmax_ef(chg):
    if (chg <= 20):
        lmax=1
    elif (chg <= 30):
        lmax=2
    elif (chg <= 38):
        lmax=1
    elif (chg <= 48):
        lmax=2
    elif (chg <= 56):
        lmax=1
    elif (chg <= 80):
        lmax=3
    elif (chg <= 88):
        lmax=1
    elif (chg <= 103):
        lmax=3
    return lmax



def check_key_in_string(key,dictionary):
    if (key not in dictionary):
        print('missing \''+key+'\' in '+dictionary['name'], flush=True)
        sys.exit()
    return None


def read_comdmft_ini_control():
    vglobl={}
    vlocal={}
    with open('comdmft.ini') as f_ini:
        code = compile(f_ini.read(), "comdmft.ini", 'exec')
        exec(code, vglobl, vlocal)
        f_ini.close()
    control=vlocal['control']
    control['name']='control'

    if (('mpi_prefix' in control) | ('mpi_prefix_flapwmbpt' in control)):
        control['mpi_prefix_flapwmbpt']=control.get('mpi_prefix_flapwmbpt', control['mpi_prefix'])
    else:
        print('no mpi_prefix for flapwmbpt')
        sys.exit()
    
    return control


def read_comdmft_ini_wan():
    vglobl={}
    vlocal={}
    with open('comdmft.ini') as f_ini:
        code = compile(f_ini.read(), "comdmft.ini", 'exec')
        exec(code, vglobl, vlocal)
        f_ini.close()
        
    if ('wan_hmat' in vlocal):
        wan_hmat=vlocal['wan_hmat']
        wan_hmat['froz_win_min']=wan_hmat.get('froz_win_min', -10.0)
        wan_hmat['froz_win_max']=wan_hmat.get('froz_win_max', 10.0)

        wan_hmat['dis_win_min']=wan_hmat.get('dis_win_min', wan_hmat['froz_win_min']-40.0)
        wan_hmat['dis_win_max']=wan_hmat.get('dis_win_max', wan_hmat['froz_win_max']+40.0)
    
        wan_hmat['num_iter']=wan_hmat.get('num_iter', 0)
        wan_hmat['dis_num_iter']=wan_hmat.get('dis_num_iter', 100)

        wan_hmat['cut_low']=wan_hmat.get('cut_low', 0.4)
        wan_hmat['cut_froz']=wan_hmat.get('cut_froz', 0.10)
        wan_hmat['cut_total']=wan_hmat.get('cut_total', 0.0)
        wan_hmat['write_wan']=True

        if (vlocal['control']['method']=='lqsgw'):
            wan_hmat['rmode']=wan_hmat.get('rmode', 0)
            wan_hmat['radfac']=wan_hmat.get('radfac', 1.0)        
        if (vlocal['control']['method']=='dft'):
            wan_hmat['rmode']=wan_hmat.get('rmode', 0)                
            wan_hmat['radfac']=wan_hmat.get('radfac', 1.0)        
        wan_hmat['name']='wan_hmat'
    
        return wan_hmat
    else:
        return None


def read_comdmft_ini_fm():
    vglobl={}
    vlocal={}
    with open('comdmft.ini') as f_ini:
        code = compile(f_ini.read(), "comdmft.ini", 'exec')
        exec(code, vglobl, vlocal)
        f_ini.close()
    # print vglobl
    # print 'here'
    fm_dict=vlocal['flapwmbpt']
    fm_dict['name']='flapwmbpt'
    control=read_comdmft_ini_control()



    if (control['method']=='dft'):
        check_key_in_string('iter_dft', fm_dict) #
        fm_dict['iter_hf']=0
        fm_dict['iter_qp']=0
        fm_dict['iter_gw']=0
    elif (control['method']=='hf'):
        check_key_in_string('iter_hf', fm_dict) #
        fm_dict['iter_dft']=fm_dict.get('iter_dft', 0)
        fm_dict['iter_qp']=0
        fm_dict['iter_gw']=0        
    if (control['method']=='lqsgw'):
        check_key_in_string('iter_dft', fm_dict) #
        check_key_in_string('iter_lqsgw', fm_dict) #
        fm_dict['iter_qp']=fm_dict['iter_lqsgw']
        fm_dict['iter_hf']=fm_dict.get('iter_hf', 0)
        fm_dict['iter_gw']=0                
    if (control['method']=='gw'):
        check_key_in_string('iter_dft', fm_dict) #
        check_key_in_string('iter_gw', fm_dict) #
        fm_dict['iter_hf']=fm_dict.get('iter_hf', 0)
        fm_dict['iter_qp']=fm_dict.get('iter_lqsgw', 0)

                
    control['restart']=control.get('restart', False)
    check_key_in_string('nproc_k_flapwmbpt', control) #
    
    if ((control['method']=='lqsgw') | (control['method']=='gw')):    
        check_key_in_string('nproc_tau_flapwmbpt', control) #
    else:
        control['nproc_tau_flapwmbpt']=1
    fm_dict['restart']=control['restart']
    fm_dict['nproc_k_flapwmbpt']=control['nproc_k_flapwmbpt']
    fm_dict['nproc_tau_flapwmbpt']=control['nproc_tau_flapwmbpt']    
    
    check_key_in_string('rel', fm_dict) #
    check_key_in_string('magn', fm_dict) #
    check_key_in_string('cif', fm_dict) #
    check_key_in_string('kmesh', fm_dict) #    


    fm_dict['rklmax']=fm_dict.get('rklmax', 2.3) #    
    fm_dict['mt_orb_l']=fm_dict.get('mt_orb_l', 0) #
    fm_dict['mt_orb_n']=fm_dict.get('mt_orb_n', 0) #
    fm_dict['mt_pb_l']=fm_dict.get('mt_pb_l', 2) #
    fm_dict['mt_pb_l_max']=fm_dict.get('mt_pb_l_max', 10) #        
    fm_dict['dft_mix']=fm_dict.get('dft_mix', 0.1) #
    fm_dict['supercell']=fm_dict.get('supercell', [1,1,1]) #    
    fm_dict['admix']=fm_dict.get('admix', fm_dict['dft_mix']) #
    fm_dict['adspin']=fm_dict.get('admix', 0.6) #

    fm_dict['gw_mix']=fm_dict.get('gw_mix', 0.1) #
    fm_dict['adm_gw']=fm_dict.get('adm_gw', fm_dict['gw_mix']) #
    fm_dict['acc_it_gw']=fm_dict.get('acc_it_gw', fm_dict['gw_mix']) #

    fm_dict['irel']=fm_dict.get('irel', fm_dict['rel']) #
    fm_dict['irel_core']=fm_dict.get('irel_core', fm_dict['rel']) #

    if (fm_dict['magn']):
        fm_dict['iter_h_ext']='1000000' #
        fm_dict['magn_shift']=0.1 #
    else:
        fm_dict['iter_h_ext']='0000000' #
        fm_dict['magn_shift']=0.0 #
#    fm_dict['version']=fm_dict.get('version', versionstr)
    fm_dict['code']='comsuite'
    fm_dict['cell']=fm_dict.get('cell', 'primitive')
    fm_dict['band']=fm_dict.get('band', True)
    fm_dict['dos']=fm_dict.get('dos', True)        
    fm_dict['temperature']=fm_dict.get('temperature', 1000.0)

    return fm_dict

def write_elkfile(ini0,elkfile):
    """
    Take a dictionary with all the relevant input settings and write an
    input file for the Elk program.
    """
    elkfile.write("tasks\n")
    elkfile.write("  0\n")
    if ini0["dos"]:
        elkfile.write("  10\n")
    if ini0["band"]:
        elkfile.write("  20\n")
        elkfile.write("  21\n")
    elkfile.write("\n")
    elkfile.write("isgkmax\n")
    elkfile.write("  -3\n")
    elkfile.write("\n")
    elkfile.write("spinpol\n")
    elkfile.write("  .true.\n")
    elkfile.write("\n")
    if ini0["dos"] or ini0["band"]:
        # vhighq seems rather expensive to run, maybe highq is good enough
        elkfile.write("highq\n")
        elkfile.write("  .true.\n")
        elkfile.write("\n")
    elkfile.write("tempk\n")
    elkfile.write("  %s\n" % ini0["temperature"])
    elkfile.write("\n")
    elkfile.write("scale\n")
    elkfile.write("  %11.7f\n" % ini0["par"])
    elkfile.write("\n")
    elkfile.write("avec\n")
    elkfile.write("  %21.16f %21.16f %21.16f\n" % ini0["a"])
    elkfile.write("  %21.16f %21.16f %21.16f\n" % ini0["b"])
    elkfile.write("  %21.16f %21.16f %21.16f\n" % ini0["c"])
    elk_species_path = os.environ.get('ELK_SPECIES_PATH')
    if not elk_species_path:
        error = "Environment variable ELK_SPECIES_PATH not set"
        elk_species_path = "."
        print(error)
    elkfile.write("\n")
    elkfile.write("atoms\n")
    natom = ini0["natom"]
    elmlist = unique_species(ini0)
    nelm = len(elmlist)
    elkfile.write("  %d\n" % nelm)
    for element in elmlist:
        elmname = element.strip()
        elmname = elmname.capitalize()
        elkfile.write("  '%s.in'\n" % (elk_species_path + "/" + elmname) )
        elkfile.write("  %d\n" % ini0["symbol"].count(element) )
        for ii in range(0,natom):
            symbol = ini0["symbol"][ii]
            if element == symbol:
                elkfile.write("  %21.16f %21.16f %21.16f 0.0 0.0 0.0\n" % (ini0["a_coeff"][ii],ini0["b_coeff"][ii],ini0["c_coeff"][ii]))
    elkfile.write("\n")
    elkfile.write("nempty\n")
    if ini0["dos"] or ini0["band"]:
        elkfile.write("  30\n")
    else:
        elkfile.write("  5\n")
    elkfile.write("\n")
    elkfile.write("ngridk\n")
    elkfile.write("  %d %d %d\n" % (ini0["mdiv"][0],ini0["mdiv"][1],ini0["mdiv"][2]))
    elkfile.write("\n")
    try:
        write_plot1dfile(ini0,elkfile)
    except ValueError:
        pass

def write_elk(ini):
    """
    Take the input data and generate the input files for the Elk
    program.
    """
    #
    # Create an .ini file for data of <key> structure
    #
    filename = "elk.in"
    elkfile = open(filename,'w')
    #
    # Write elk.in
    #
    write_elkfile(ini,elkfile)
    #
    # Close elk.in file
    #
    elkfile.close()    
    #
    # Create an .kpoint file for data of <key> structure
    #
    filename = "kpoints"
    kpointfile = open(filename,'w')
    #
    # Write MatDeLab.kpath
    #
    try:
        write_kpointfile(ini,kpointfile)
        kpointfile.close()
    except ValueError:
        kpointfile.close()
        os.remove(filename)

def write_wien2k(ini):
    """
    Take the input data and generate the input files for Wien2K 
    package.
    """
    import sys
    import subprocess
    case = ini["cif"].split(".")[0]
    case_st2 = case+".struct"
    case_st1 = case_st2+"_st"
    case_st3 = "/tmp/"+case_st2
    make_struct = []
    make_struct.append("cif2struct")
    make_struct.append("/tmp/"+os.path.basename(ini["cif"]))
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_st3)
    make_struct.append(case_st2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("x")
    make_struct.append("symmetry")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_st1)
    make_struct.append(case_st2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    # Run setrmt_lapw to choose the muffin tin radii.
    # The -r flag specifies the percentage reduction of the radii from
    # the just-touching radii. This is a REQUIRED flag because 
    # setrmt_lapw rounds the radii after calculating them to 2 decimal
    # places. The test whether the spheres are overlapping uses at
    # least 5 decimal places. Hence, if the spheres are not reduced
    # the rounding of the radii may cause the non-overlapping requirement
    # to be violated, and the calculation will abort!
    #
    # - setrmt_lapw case -r 3
    #
    case_st2 = case+".struct"
    case_st1 = case_st2+"_setrmt"
    make_struct = []
    make_struct.append("setrmt_lapw")
    make_struct.append(case)
    make_struct.append("-r")
    make_struct.append("3")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_st1)
    make_struct.append(case_st2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()

    case_ins = case+".inst"
    make_struct = []
    make_struct.append("rm")
    make_struct.append("-f")
    make_struct.append(case_ins)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("instgen_lapw")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    make_struct = []
    make_struct.append("x")
    make_struct.append("lstart")
    process = subprocess.Popen(make_struct,stdin=subprocess.PIPE)
    outs,errs = process.communicate(bytes("5\n-6.0\n","utf-8"))
    case_in2 = case+".in0"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in2 = case+".in1"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in2 = case+".vsp"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in3 = case+".in2"
    case_in1 = case_in3+"_ls"
    case_in2 = case_in3+"_sy"
    #make_struct = []
    #make_struct.append("cat")
    #make_struct.append(case_in1)
    #make_struct.append(case_in2)
    #make_struct.append(">")
    #make_struct.append(case_in3)
    line = "cat "+str(case_in1)+" "+str(case_in2)+" > "+str(case_in3)
    result = subprocess.run(line,shell=True)
    #if sys.version_info.major==3:
    #    result = result.decode()
    #
    # If Wien2K thinks there is no inversion symmetry the code needs 
    # .in1c and .in2c files instead of .in1 and .in2 files. 
    # To generate the former files just copy the latter.
    #
    # - cp case.in1 case.in1c
    #
    case_in1 = case+".in1"
    case_in2 = case_in1+"c"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    # - cp case.in2 case.in2c
    #
    case_in1 = case+".in2"
    case_in2 = case_in1+"c"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    make_struct = []
    make_struct.append("x")
    make_struct.append("kgen")
    process = subprocess.Popen(make_struct,stdin=subprocess.PIPE)
    line = "0\n"+str(ini["mdiv"][0])+" "+str(ini["mdiv"][1])+" "+str(ini["mdiv"][2])+"\n0"
    print("line = %s\n" % line)
    outs,errs = process.communicate(bytes(line,"utf-8"))
    make_struct = []
    make_struct.append("x")
    make_struct.append("dstart")
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    case_in2 = case+".inm"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    case_in2 = case+".inc"
    case_in1 = case_in2+"_st"
    make_struct = []
    make_struct.append("cp")
    make_struct.append(case_in1)
    make_struct.append(case_in2)
    result = subprocess.check_output(make_struct)
    if sys.version_info.major==3:
        result = result.decode()
    #
    # Create the .klist_band file
    #
    write_klistfile(ini)


    
def main():
    """
    Run the input generator for the method and ciffile as specified by the 
    command line arguments.
    """

    ini=read_comdmft_ini_fm()
    tmpfile = any2utf8(ini['cif'])
    
    str_in=pgsa.SpacegroupAnalyzer(pgsa.Structure.from_file(any2utf8(ini['cif']))).get_primitive_standard_structure()
    # print(str_in.lattice.matrix/bohr)

    
    xtalstr  = fm.crystalstructure(str_in)
    if (ini['supercell'] != [1,1,1]):
        xtalstr=supercell(xtalstr,ini)

    xtalstr.str.to(fmt='json', filename='crystal_structure.json')
    xtalstr.str.to(fmt='xsf', filename='crystal_structure.xsf')
    print(pgsa.SpacegroupAnalyzer(xtalstr.str).get_symmetry_dataset())
    # json.dump(mg.symmetry.analyzer.SpacegroupAnalyzer(mg.symmetry.analyzer.SpacegroupAnalyzer(xtalstr.str).get_conventional_standard_structure()).get_symmetry_dataset(),'xtal_str_conv.json')    
    
    # xtalstr.str.to(fmt="cif",filename="/tmp/"+os.path.basename(ini['cif']))
    # xtalstr.prints()
    Kmax   = 0.0
    code   = ini['code']
    ini["par"] = 1.0
    ini["b_a"] = 1.0
    ini["c_a"] = 1.0
    ini["cut_lapw_ratio"]=0.61
    ini["cut_pb_ratio"]=0.98
    ini = retr_cell_volume(xtalstr,ini)    
    ini = retr_lattice_vecs(xtalstr,ini)
    ini = retr_recip_vecs(xtalstr,ini)
    ini = retr_sites(xtalstr,ini)
    ini = retr_distance_matrix(xtalstr,ini)
    ini = retr_lattice_type(xtalstr,ini)
    ini = retr_kpath(xtalstr,ini)
    ini = retr_kpoint(xtalstr,ini)
    
    ini = retr_spacegroup_number(xtalstr,ini)

    ini = translate_elements(ini)
    ini = establish_mt_radii(ini)
    ini = establish_atoms_volume(ini)
    ini = establish_Kmin(ini)
    ini = establish_lmax(ini)    
    ini = establish_Kmax(ini,Kmax=Kmax)
    ini = establish_r_grid(ini)

    ini = expand_atomic_basis(ini)

    if code == "comsuite":
        write_comsuite(ini)
    elif code == "elk":
        write_elk(ini)
    elif code == "wien2k":
        write_wien2k(ini)
    else:
        error = "Unknown code suite: "+code+" No input files generated"
        print(error)
        
if __name__ == "__main__":
    main()
