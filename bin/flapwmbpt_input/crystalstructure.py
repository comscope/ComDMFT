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
import numpy as np
import pymatgen as mg
import math
import pprint
# from common import vec_norm

bohr=0.529177210903

def vec_norm(basis, vector):
    tempmat=np.zeros((3,3))
    for ii in range(3):
        for jj in range(3):
            tempmat[ii,jj]=np.dot(basis[:,ii], basis[:,jj])
    return np.dot(np.dot(vector, tempmat), vector)


class crystalstructure: 

    def __init__(self,str_in):
        """
        Create a structure instance based on the contents of a CIF file.
        Due to the implementation of pymatgen we first have to setup a dummy
        structure to create an instance whose from_file method we can invoke.
        """
        # if   cellkind == "primitive":
        #     self.str=mg.symmetry.analyzer.SpacegroupAnalyzer(mg.Structure.from_file(ciffile)).get_primitive_standard_structure()
        # elif cellkind == "conventional":
        #     self.str=mg.symmetry.analyzer.SpacegroupAnalyzer(mg.Structure.from_file(ciffile)).get_conventional_standard_structure()            
        # else:
        #     print("Unknown cellkind: %s" % cellkind)
        #     print("Valid options are \"primitive\" or \"conventional\"")
        #if self.str.num_sites > 1:
        #  self.str.merge_sites(mode="delete") # remove any duplicate atoms


        self.str=str_in
        self.sga    = mg.symmetry.analyzer.SpacegroupAnalyzer(self.str)
        # if cellkind == "primitive":
        #     self.str = mg.symmetry.analyzer.SpacegroupAnalyzer(self.str).find_primitive()
        # self.sga    = mg.symmetry.analyzer.SpacegroupAnalyzer(self.str)
        # self.kpoints  = mg.symmetry.kpath.KPathSetyawanCurtarolo(self.str).get_kpoints(line_density=150, coords_are_cartesian=False)

        self.nsymop=np.shape(np.array(self.sga.get_symmetry_dataset()['rotations']))[0]
        self.rotation=np.array(self.sga.get_symmetry_dataset()['rotations'])
        self.translation=np.array(self.sga.get_symmetry_dataset()['translations'])
        self.real_lattice=np.transpose(self.str.lattice.matrix)/bohr
        self.recip_lattice=np.transpose(self.str.lattice.reciprocal_lattice.matrix)*bohr
        self.frac_coords=np.transpose(self.str.frac_coords)
        self.atoms=self.str.sites
        self.lattice_type=self.sga.get_lattice_type()
        self.spacegroup_hall=self.sga.get_hall()
        self.spacegroup_number=self.sga.get_space_group_number()

        islist0=list(self.sga.get_symmetry_dataset()['equivalent_atoms'])
        islist=[]
        for ii in range(len(islist0)):
            if (ii ==0):
                islist.append(1)
            else:
                if (islist0[ii] in islist0[:ii]):
                    islist.append(islist[islist0[:ii].index(islist0[ii])])
                else:
                    islist.append(islist[ii-1]+1)
        self.equivalent_atoms=np.array(islist)
        
        self.crystal_system=self.sga.get_crystal_system()
        self.volume=self.str.volume/(bohr**3)
        self.natom=len(self.atoms)
        self.species=self.str.species
        self.dist_mat=self.str.distance_matrix/bohr

        dist_diag=100000000000000
        for ii in range(-3,3):
            for jj in range(-3,3):
                for kk in range(-3,3):
                    if (not((ii==0) & (jj==0) & (kk==0))):
                        if (dist_diag>np.sqrt(vec_norm(self.real_lattice, np.array([ii,jj,kk])))):
                            dist_diag=np.sqrt(vec_norm(self.real_lattice, np.array([ii,jj,kk])))
        for ii in range(self.natom):
            self.dist_mat[ii,ii]=dist_diag
        # print(type(self.dist_mat))                                    
        # print(self.dist_mat)                        
            
        
        # print(np.dot(np.transpose(self.recip_lattice), self.real_lattice))

        # self.test_symm_op()
        self.write_symmetry_operation_input()

        sp = mg.symmetry.bandstructure.HighSymmKpath(self.str, path_type='sc')
        kpath_label = []
        for labellist in sp.kpath["path"]:
            kpath_label.extend(labellist)

            
        kpath_bound = []
        for label in kpath_label:
            red_coord = sp.kpath["kpoints"][label]
            kpath_bound.append(red_coord)

        kpath_label_old=kpath_label
        kpath_label=[]
        for string in  kpath_label_old:
            new_string = string.replace("\\Gamma", "G")
            new_string = new_string.replace("\\Sigma", "s")
            new_string = new_string.replace("_", "")                        
            kpath_label.append(new_string)
            
        self.kpath_label  = kpath_label
        self.kpath  = kpath_bound


        den=150

        kpoint=[]
        kpoint_label=[]        
        for ii in range(len(kpath_bound)-1):
            vec1=kpath_bound[ii]
            vec2=kpath_bound[ii+1]            
            vec3=np.dot(self.recip_lattice, vec1-vec2)
            length=np.sqrt(np.sum(vec3**2))
            num=np.int(np.round(length*den))
            for jj in range(num):
                if (jj == 0):
                    kpoint_label.append(kpath_label[ii])
                else:
                    kpoint_label.append(" ")                    
                vec4=vec1*(num-jj)/num+vec2*(jj)/num
                kpoint.append(vec4)
        kpoint.append(kpath_bound[-1])
        kpoint_label.append(kpath_label[-1])        
        self.kpoint  = kpoint
        self.kpoint_label  = kpoint_label
        

    # def test_symm_op(self):

    #     print(self.atoms)
    #     print(self.equivalent_atoms)
        
    #     for ioper in range(self.nsymop):
    #         rotmat=np.dot(self.real_lattice, np.dot(self.rotation[ioper, :,:], np.transpose(dual_vector(self.real_lattice))))
    #         # print(ioper, rotmat)
    #         # print(self.translation[ioper, :])
    #         # print('\n')
    #         # print(ioper, rotmat)
    #         if (abs(abs(np.linalg.det(rotmat))-1.0)>1e-6):
    #             print(ioper, ' rotation_matrix wrong\n')
    #             sys.exit()

    #     for ioper in range(self.nsymop):
    #         for atom1 in (range(self.natom)):
                
    #             coord1=self.frac_coords[:,atom1]
    #             new_coord1=np.dot(self.rotation[ioper, :,:], coord1)+self.translation[ioper,:]
    #             # rotmat=self.rotation[ioper, :,:]
    #             # print(rotmat)
    #             # new_coord1=np.dot(np.ones((3,3)), coord1)
    #             # print(ioper, np.shape(coord1))
    #             # print(ioper, np.shape(self.rotation[ioper, :,:]))
    #             # print(ioper, key1["species"][0]['element'],coord1,new_coord1)
    #             cnt=0
    #             for atom2 in (range(self.natom)):
    #                 coord2=self.frac_coords[:,atom1]                
    #                 for ii in range(-3,3):
    #                     for jj in range(-3,3):
    #                         for kk in range(-3,3):
    #                             if (self.equivalent_atoms[
    #                                 cnt=cnt+1
    #                                 # print(ioper, key1["species"][0]['element'],key2["species"][0]['element'],ii,jj,kk,coord1,new_coord1,vec_norm(self.real_lattice, coord2+np.array([ii,jj,kk])-new_coord1))
    #             if (cnt!=1):
    #                 print(ioper, key1, cnt, 'sym_operation wrong')
        

# def sym_generator(int_symbol):
#     gen_matrices = symmetry.groups._get_symm_data("generator_matrices")
#     sgencoding = symmetry.groups._get_symm_data("space_group_encoding")
#     abbrev_sg_mapping = symmetry.groups._get_symm_data("abbreviated_spacegroup_symbols")
#     translations = {k: Fraction(v) for k, v in symmetry.groups._get_symm_data("translations").items()}
#     full_sg_mapping = {v["full_symbol"]: k for k, v in  symmetry.groups._get_symm_data("space_group_encoding").items()}

#     # int_symbol=structure_supercell_primitive.get_space_group_info()[0]
#     ii=sgencoding[int_symbol]['int_number']
#     data = sgencoding[int_symbol]
#     # TODO: Support different origin choices.
#     enc = list(data["enc"])
#     inversion = int(enc.pop(0))
#     ngen = int(enc.pop(0))
#     symm_ops = [np.eye(4)]
#     if inversion:
#         symm_ops.append(np.array(
#             [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0],
#              [0, 0, 0, 1]]))
#     for i in range(ngen):
#         m = np.eye(4)
#         m[:3, :3] = gen_matrices[enc.pop(0)]
#         m[0, 3] = translations[enc.pop(0)]
#         m[1, 3] = translations[enc.pop(0)]
#         m[2, 3] = translations[enc.pop(0)]
#         symm_ops.append(m)
#     return sym_ops
    def write_symmetry_operation_input(self):
        f=open('symmetry_operations', 'w')
        f.write("   number of symmetry operations= %3d\n" % (self.nsymop))
        for ii in range(self.nsymop):    
            f.write("    symmetry operation %3d\n" % (ii+1))
            for jj in range(0,3):
                line = (self.rotation[ii][jj][0],self.rotation[ii][jj][1],self.rotation[ii][jj][2],self.translation[ii][jj])
                f.write("    (%14.10f,%14.10f,%14.10f)    (%14.10f)\n" % line)    
        f.close()


    # def prints(self):
    #     """
    #     Print the contents of the current instance."
    #     """
    #     pprint.pprint(self)
