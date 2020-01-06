"""
A class to query and manipulate crystal structure from a CIF file.

This class uses pymatgen underneath to implement most of the
functionality. The idea is provide all the usual data about
a material as well as some of more idiosyncratic data specific to
the programs in the MatDeLab package.

Dependencies:
- pymatgen
- spglib
"""
import array
import numpy
import pymatgen as mg
import math

#print mg.__version__

#
# Angstrom to Bohr conversion factor
#
#    "CODATA recommended values of the fundamental physical constants: 2014", 
#    P.J. Mohr, D. B. Newell, B. N. Taylor,
#    the Committee on Data for Science and Technology (CODATA),
#    August, 2015, <http://dx.doi.org/10.5281/zenodo.22826>
#
# See also:
#
#    "Fundamental Physical Constants", NIST
#    (<http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0>)
#
angstrom2bohr = 1.8897261254535

#
# Dictionary of lattice types
#
# The GW code needs to know the lattice type to select the k-point
# path for bandstructures. For this purpose a set of constants are
# defined and we need to select the right one.
#
lattice_type                     = {}
lattice_type["cubic"]            = -1
lattice_type["cubic_bcc"]        = -2
lattice_type["cubic_fcc"]        = -3
lattice_type["hexagonal"]        = -4
lattice_type["tetragonal"]       = -5
lattice_type["tetragonal_bcc"]   = -6
lattice_type["orthorhombic"]     = -7
lattice_type["orthorhombic_bcc"] = -8
lattice_type["orthorhombic_fcc"] = -9
lattice_type["monoclinic"]       = -10
lattice_type["rhombohedral"]     = -11
lattice_type["triclinic"]        = -13

def map_to_conventional_cell(InKPath):
    """
    Map the fractional coordinates of the high symmetry path to the range 0 to 1

    For some bizarre reason Pymatgen will generate negative fractional coordinates
    of the high symmetry points for some cells. This situation causes problems
    when trying to compare the band structures from different codes. Hence this
    function will generate a sane path will all the high symmetry points mapped
    into the conventional unit cell.
    """
    InKPoints = InKPath['kpoints']
    OutKPoints = {}
    for key,inval in InKPoints.items():
        (kx,ky,kz) = inval
        if kx < 0.0:
            kx += 1.0
        if ky < 0.0:
            ky += 1.0
        if kz < 0.0:
            kz += 1.0
        outval = numpy.array([kx,ky,kz])
        OutKPoints[key] = outval
    OutKPath = {}
    OutKPath['kpoints'] = OutKPoints
    OutKPath['path'] = InKPath['path']
    return OutKPath

class structure: 

    def __init__(self,ciffilename,cellkind):
        """
        Create a structure instance based on the contents of a CIF file.
        Due to the implementation of pymatgen we first have to setup a dummy
        structure to create an instance whose from_file method we can invoke.
        """
        if   cellkind == "primitive":
            self.struct = mg.Structure.from_file(ciffilename,primitive=True)
        elif cellkind == "conventional":
            self.struct = mg.Structure.from_file(ciffilename,primitive=False)
        else:
            print("Unknown cellkind: %s" % cellkind)
            print("Valid options are \"primitive\" or \"conventional\"")
        #if self.struct.num_sites > 1:
        #  self.struct.merge_sites(mode="delete") # remove any duplicate atoms
        self.sga    = mg.symmetry.analyzer.SpacegroupAnalyzer(self.struct)
        self.struct = self.sga.get_refined_structure()
        if cellkind == "primitive":
            self.struct = mg.symmetry.analyzer.SpacegroupAnalyzer(self.struct).find_primitive()
        self.sga    = mg.symmetry.analyzer.SpacegroupAnalyzer(self.struct)
        self.kpath  = mg.symmetry.bandstructure.HighSymmKpath(self.struct)
        self.kpath._kpath = map_to_conventional_cell(self.kpath._kpath)

    def lattice(self):
        """
        Return the lattice of structure.
        """
        return self.struct.lattice

    def frac_coords(self):
        """
        Return the factional coordinates of the atoms.
        """
        return self.struct.frac_coords

    def sites(self):
        """
        Return an iterator over sites.
        """
        return self.struct.sites

    def lattice_type(self):
        """
        Return the lattice type as a string.
        """
        return self.sga.get_lattice_type()

    def spacegroup_hall(self):
        """
        Return the Hall spacegroup symbol.
        """
        hall = self.sga.get_hall()
        return hall

    def norm2(self,vec):
        """
        Calculate the 2-norm of a vector.
        """
        d = 0.0
        for xx in vec:
            d += xx*xx
        d = math.sqrt(d)
        return d

    def get_symmetry_operations(self):
        """
        Return the symmetry operations.
        """
        switch_cartesian=False
        lattice = self.struct.lattice.matrix
        symmops = self.sga.get_symmetry_operations(cartesian=switch_cartesian)
        # If the structure is given in the primitive cell it looks like
        # some of the symmetry operations returned are actually dodgy.
        # So we have to filter the dodgy operators out and return only the
        # sane ones.
        symmops_out = []
        for symmop in symmops:
           rotation = symmop.rotation_matrix
           ok = True
           #for ii in range(0,3):
           #   if ((rotation[0][ii]==rotation[1][ii]) and
           #       (rotation[0][ii]==rotation[2][ii])):
           #     ok = False
           if ok and switch_cartesian:
               vec     = [0.0,0.0,0.0]
               vec1    = [0.0,0.0,0.0]
               vect    = [0.0,0.0,0.0]
               vec[0]  = symmop.translation_vector[0]
               vec[1]  = symmop.translation_vector[1]
               vec[2]  = symmop.translation_vector[2]
               dd      = self.norm2(vec)
               vec1[0] = vec[0]
               vec1[1] = vec[1]
               vec1[2] = vec[2]
               for ii in range(3,-4,-1):
                   for jj in range(3,-4,-1):
                       for kk in range(3,-4,-1):
                           for ix in range(0,3):
                               vect[ix] = vec1[ix] + ii*lattice[0][ix] + jj*lattice[1][ix] + kk*lattice[2][ix]
                           d0 = self.norm2(vect)
                           if d0+1.0e-6 < dd:
                               dd  = d0
                               vec[0] = vect[0]
                               vec[1] = vect[1]
                               vec[2] = vect[2]
               symmop.translation_vector[0] = vec[0]
               symmop.translation_vector[1] = vec[1]
               symmop.translation_vector[2] = vec[2]
               symmop.translation_vector[0] = symmop.translation_vector[0]*angstrom2bohr
               symmop.translation_vector[1] = symmop.translation_vector[1]*angstrom2bohr
               symmop.translation_vector[2] = symmop.translation_vector[2]*angstrom2bohr
           elif ok:
               for ii in range(0,3):
                   xx = symmop.translation_vector[ii]
                   xx = abs(xx-int(xx+0.5))
                   if xx < 1.0e-6:
                       symmop.translation_vector[ii] = 0.0
           if ok:
             symmops_out.append(symmop)
        return symmops_out

    def get_symmetry_generators(self):
        """
        Return the symmetry generators.

        Pymatgen produces the symmetry operations but some programs rather
        work with the group generators. This routine establishes the list
        of group generators from the list of symmetry operations. The
        generators of interest are:
        - E  : The identity
        - RnV: n-Fold rotations axes V
        - SnV: n-Fold improper rotations axes V
        - MV : Mirror planes perpendicular to V
        - I  : Inversion
        The algorithm makes use of the fact that a symmetry operation transforms
        a vector into another vector. Depending on the particulars of the 
        symmetry operation this transformation has the following properties:
        - E  : leaves all vector unchanged
        - RnV: leaves only the vector V unchanged, all other vectors W are
               transformed into vectors W' so that 
               n=360/(angle between W' and W)
        - SnV: Same as RnV*I
        - MV : flips the sign on V, and leaves at least two vectors spanning
               the mirror plane unchanged
        - I  : flips the sign on all vectors.
        The list of symmetry operations also contains products of the generators
        which cannot be resolved this way. Those operations will be discarded
        as only the generators themselves are needed. 

        In this routine we need to be working in an orthonormal basis. Hence
        the symmetry operations will be obtained in Cartesian coordinates.
        """
        switch_cartesian=True
        lattice = self.struct.lattice.matrix
        symmops = self.sga.get_symmetry_operations(cartesian=switch_cartesian)
        # If the structure is given in the primitive cell it looks like
        # some of the symmetry operations returned are actually dodgy.
        # So we have to filter the dodgy operators out and return only the
        # sane ones.
        symmops_out = []
        symmops_tmp = []
        axes=[[ 1.0, 0.0, 0.0], # x-axis
              [ 0.0, 1.0, 0.0], # y-axis
              [ 0.0, 0.0, 1.0], # z-axis
              [ 1.0, 1.0, 0.0], # x,y-mid edge
              [ 1.0, 0.0, 1.0], # x,z-mid edge
              [ 0.0, 1.0, 1.0], # y,z-mid edge
              [-1.0, 1.0, 0.0], # -x,y-mid edge
              [-1.0, 0.0, 1.0], # -x,z-mid edge
              [ 0.0,-1.0, 1.0], # -y,z-mid edge
              [ 1.0, 1.0, 1.0], # x,y,z-diagonal
              [-1.0, 1.0, 1.0], # -x,y,z-diagonal
              [ 1.0,-1.0, 1.0], # x,-y,z-diagonal
              [ 1.0, 1.0,-1.0], # x,-y,z-diagonal
              [ 2.0,-1.0,-1.0], # x,-y,z-diagonal
              [-1.0, 2.0,-1.0], # x,-y,z-diagonal
              [-1.0,-1.0, 2.0], # x,-y,z-diagonal
              [ 2.0, 1.0, 1.0], # x,-y,z-diagonal
              [ 1.0, 2.0, 1.0], # x,-y,z-diagonal
              [ 1.0, 1.0, 2.0]] # x,-y,z-diagonal
        inversion=[[-1.0, 0.0, 0.0],
                   [ 0.0,-1.0, 0.0],
                   [ 0.0, 0.0,-1.0]]
        for symmop in symmops:
            rotation = symmop.rotation_matrix
            tran=numpy.matmul(axes,rotation)
            ixtran=numpy.matmul(tran,inversion)
            (nvec,vlen)=tran.shape
            nsame   = 0 # number of tran vectors that equal axes vectors
            nsign   = 0 # number of tran vectors with opposite sign as axes
                        # vectors
            nixsame = 0 # number of ixtran vectors that equal axes vectors
            nixsign = 0 # number of ixtran vectors with opposite sign as axes
                        # vectors
            isame   =-1 # number of the vector that is the same
            isign   =-1 # number of the vector that flipped sign
            ixsame  =-1 # number of the vector that is the same in ixtran
            ixsign  =-1 # number of the vector that flipped sign in ixtran
            iorth   =-1 # number of the vector orthogonal to isame
            for ii in range(0,nvec):
                if abs((numpy.dot(tran[ii],axes[ii])/
                        numpy.dot(axes[ii],axes[ii]))-1.0)<1.e-6:
                    nsame += 1
                    if isame == -1:
                        isame  = ii
                if abs((numpy.dot(tran[ii],axes[ii])/
                        numpy.dot(axes[ii],axes[ii]))+1.0)<1.e-6:
                    nsign += 1
                    if isign == -1:
                        isign  = ii
                if abs((numpy.dot(ixtran[ii],axes[ii])/
                        numpy.dot(axes[ii],axes[ii]))-1.0)<1.e-6:
                    nixsame += 1
                    if ixsame == -1:
                        ixsame   = ii
                if abs((numpy.dot(ixtran[ii],axes[ii])/
                        numpy.dot(axes[ii],axes[ii]))+1.0)<1.e-6:
                    nixsign += 1
                    if ixsign == -1:
                        ixsign   = ii
            # Workout transformation kind E,I,R,M,S, or U (unidentified)
            if nvec == nsame:
                symmop.kind="E"
                symmop.axis=None
                symmop.n   =None
            elif nvec == nsign:
                symmop.kind="I"
                symmop.axis=None
                symmop.n   =None
            elif 1 == nsame:
                symmop.kind="R"
                symmop.axis=axes[isame]
                iorth=-1
                for ii in range(0,nvec):
                    if abs(numpy.dot(axes[ii],axes[isame]))<1.e-6 and iorth==-1:
                        iorth = ii
                veca = numpy.array(axes[iorth])
                vect = numpy.array(tran[iorth])
                symmop.n=(2*math.acos(-1.0)/
                          math.acos(numpy.dot(vect,veca)/
                               (math.sqrt(numpy.dot(veca,veca))*
                                math.sqrt(numpy.dot(vect,vect))
                               )))
                symmop.n=int(symmop.n+0.1)
            elif 1 == nsign and 2 <= nsame:
                symmop.kind="M"
                symmop.axis=axes[isign]
                symmop.n   =None
            elif 1 == nixsame:
                symmop.kind="S"
                symmop.axis=axes[ixsame]
                iorth=-1
                for ii in range(0,nvec):
                    if abs(numpy.dot(axes[ii],axes[ixsame]))<1.e-6 and iorth==-1:
                        iorth = ii
                veca = numpy.array(axes[iorth])
                vect = numpy.array(tran[iorth])
                symmop.n=(2*math.acos(-1.0)/
                          math.acos(numpy.dot(vect,veca)/
                               (math.sqrt(numpy.dot(veca,veca))*
                                math.sqrt(numpy.dot(vect,vect))
                               )))
                symmop.n=int(symmop.n+0.1)
            else:
                symmop.kind="U"
                symmop.axis=None
                symmop.n   =None
            ok = True
            #for ii in range(0,3):
            #   if ((rotation[0][ii]==rotation[1][ii]) and
            #       (rotation[0][ii]==rotation[2][ii])):
            #     ok = False
            if ok and switch_cartesian:
                #DEBUG
                #print(symmop.translation_vector)
                #DEBUG
                vec     = [0.0,0.0,0.0]
                vec1    = [0.0,0.0,0.0]
                vect    = [0.0,0.0,0.0]
                vec[0]  = symmop.translation_vector[0]
                vec[1]  = symmop.translation_vector[1]
                vec[2]  = symmop.translation_vector[2]
                dd      = self.norm2(vec)
                vec1[0] = vec[0]
                vec1[1] = vec[1]
                vec1[2] = vec[2]
                for ii in range(3,-4,-1):
                    for jj in range(3,-4,-1):
                        for kk in range(3,-4,-1):
                            for ix in range(0,3):
                                vect[ix] = vec1[ix] + ii*lattice[0][ix] + jj*lattice[1][ix] + kk*lattice[2][ix]
                            d0 = self.norm2(vect)
                            if d0+1.0e-6 < dd:
                                dd  = d0
                                vec[0] = vect[0]
                                vec[1] = vect[1]
                                vec[2] = vect[2]
                if abs(vec[0]) < 1.0e-12:
                    vec[0] = 0.0e0
                if abs(vec[1]) < 1.0e-12:
                    vec[1] = 0.0e0
                if abs(vec[2]) < 1.0e-12:
                    vec[2] = 0.0e0
                symmop.translation_vector[0] = vec[0]
                symmop.translation_vector[1] = vec[1]
                symmop.translation_vector[2] = vec[2]
                symmop.translation_vector[0] = symmop.translation_vector[0]*angstrom2bohr
                symmop.translation_vector[1] = symmop.translation_vector[1]*angstrom2bohr
                symmop.translation_vector[2] = symmop.translation_vector[2]*angstrom2bohr
            elif ok:
                for ii in range(0,3):
                    xx = symmop.translation_vector[ii]
                    xx = abs(xx-int(xx+0.4))
                    if xx < 1.0e-6:
                        symmop.translation_vector[ii] = 0.0
            
            symmops_tmp.append(symmop)
        #DEBUG
        #for symmop in symmops_tmp:
        #    print(symmop.kind, symmop.n, symmop.axis,symmop.translation_vector)
        mult=1
        #DEBUG
        #
        # Filter the generators out from the other operations
        #
        icur = 0
        while icur < len(symmops_tmp):
            #DEBUG
            #print("entry: ",icur,mult)
            #DEBUG
            symmop = symmops_tmp[icur]
            mult = 1
            if symmop.kind == "R" or symmop.kind == "S":
                mult = symmop.n-1
            if symmop.kind != "U":
                symmops_out.append(symmop)
            icur += max(icur,1)*mult
        return symmops_out

    def get_equivalent_sites(self):
        """
        Return a list of lists of symmetry equivalent sites in the structure
        """
        symm = self.sga.get_symmetrized_structure()
        sites = symm.equivalent_sites
        return sites

    def get_distance_matrix(self):
        """
        Return the distance matrix in units of Bohr.
        """
        raw_mat = self.struct.distance_matrix
        dist_mat = []
        for ii in range(0,len(raw_mat)):
            raw_row = raw_mat[ii]
            dist_row = []
            for jj in range(0,len(raw_row)):
                dist_row.append(raw_row[jj]*angstrom2bohr)
            dist_mat.append(dist_row)
        return dist_mat

    def get_kpath(self):
        """
        Return the symmetry line path in reciprocal space
        """
        return self.kpath.kpath

    def get_kpoints(self):
        """
        Return the kpoints along the symmetry line path in reciprocal space
        """
        if self.kpath.kpath:
            return self.kpath.get_kpoints(coords_are_cartesian=False)
        else:
            return None

    def get_space_group_number(self):
        """
        Return the space group number
        """
        spgr = int(self.sga.get_space_group_number())
        return spgr

    def crystal_system(self):
        """
        Return the crystal system.
        """
        sys = self.sga.get_crystal_system()
        return sys

    def cell_volume(self):
        """
        Return the volume of the unit cell
        """
        volume = self.struct.volume*(angstrom2bohr**3)
        return volume

    def prints(self):
        """
        Print the contents of the current instance."
        """
        print(self.struct)

def retr_lattice_vecs(struct,ini_struct):
    """
    Retrieve the lattice vectors.

    The lattice vectors are extracted from the structure instance.
    They are added to the dictionary passed in. The updated dictionary
    is returned. The lattice vectors are given in Bohr.
    """
    lattice = struct.lattice()
    vec_a = lattice.matrix[0]
    vec_b = lattice.matrix[1]
    vec_c = lattice.matrix[2]
    vec_a = (vec_a[0]*angstrom2bohr,
             vec_a[1]*angstrom2bohr,
             vec_a[2]*angstrom2bohr)
    vec_b = (vec_b[0]*angstrom2bohr,
             vec_b[1]*angstrom2bohr,
             vec_b[2]*angstrom2bohr)
    vec_c = (vec_c[0]*angstrom2bohr,
             vec_c[1]*angstrom2bohr,
             vec_c[2]*angstrom2bohr)
    ini_struct["a"] = vec_a
    ini_struct["b"] = vec_b
    ini_struct["c"] = vec_c
    return ini_struct

def retr_sites(struct,ini_struct):
    """
    Retrieve the atomic positions and elements in fractional coordinates.

    The atomic positions and elements are extracted from the structures
    instance. They are added to the dictionary passed in. The updated dictionary
    is returned.
    """
    xcoords = []
    ycoords = []
    zcoords = []
    element = []
    islist  = []
    symeqv = struct.get_equivalent_sites()
    ii = 0
    for sites in symeqv:
        ii += 1
        for site in sites:
            xcoords.append(site.a)
            ycoords.append(site.b)
            zcoords.append(site.c)
            islist.append(ii)
            str = site.species_string
            str = str.lower()
            element.append(str)
    natom = len(xcoords)
    ini_struct["natom"]  = natom
    ini_struct["xcoord"] = xcoords
    ini_struct["ycoord"] = ycoords
    ini_struct["zcoord"] = zcoords
    ini_struct["islist"] = islist
    ini_struct["symbol"] = element
    ini_struct["symeqv"] = symeqv
    return ini_struct

def retr_distance_matrix(struct,ini_struct):
    """
    Retrieve the distances between sites.
    """
    dist_mat = struct.get_distance_matrix()
    ini_struct["dist_matrix"] = dist_mat
    return ini_struct

def retr_lattice_type(struct,ini):
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
    lattice    = struct.crystal_system()
    hall       = struct.spacegroup_hall()
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
    ini["hall"]         = hall
    if lattice_tp:
      ini["istruc"]     = lattice_type[lattice_tp]
    else:
      ini["istruc"]     = 0
    ini["lattice_type"] = lattice_tp
    return ini

def retr_cell_volume(struct,ini):
    """
    Retrieve the cell volume
    """
    ini["cell_volume"] = struct.cell_volume()
    return ini

def retr_kpath(struct,ini):
    """
    Retrieve the kpath
    """
    kpath = struct.get_kpath()
    ini["kpath"] = kpath
    return ini

def retr_kpoints(struct,ini):
    """
    Retrieve the kpoints
    """
    kpoints = struct.get_kpoints()
    ini["kpoints"] = kpoints
    return ini
