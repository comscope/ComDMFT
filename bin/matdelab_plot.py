#!/usr/bin/env python
#
import sys
import glob
import matplotlib.pyplot as plt
import matdelab_plot.figure_bd as fbd

versionstr = '%(prog)s version 0.0'

#
# Note on units:
#
# The energy scales are in eV
#
def parse_arguments():
    """
    Parse command line arguments for the matdelab_plot.py script when run from
    the command line.
    """
    from argparse import ArgumentParser
    prs = ArgumentParser()
    prs.add_argument("--version",action='version',
                     version=versionstr)
    opt = prs.add_mutually_exclusive_group()
    opt.add_argument('--band',action='store_const',dest='plotType',
                     const='band',default='bd',
                     help='Plot the band structure only')
    opt.add_argument('--dos',action='store_const',dest='plotType',
                     const='dos',default='bd',
                     help='Plot the density of states only')
    opt.add_argument('--bd',action='store_const',dest='plotType',
                     const='bd',default='bd',
                     help='Plot both the band structure and the density of states')
    opt = prs.add_mutually_exclusive_group()
    opt.add_argument('--dft',action='store_const',dest='statesType',
                     const='dft',default='dft',
                     help='Use DFT states')
    opt.add_argument('--qp',action='store_const',dest='statesType',
                     const='qp',default='dft',
                     help='Use Quasi-Particle states')
    prs.add_argument('--e-scale',dest='energyScale',
                     help='Set the energy scale (--e-scale="Elow:Ehigh") where Elow and Ehigh are in eV')
    prs.add_argument('--nstates-scale',dest='nstatesScale',
                     help='Set the number of states scale (--nstates-scale=Nhigh)')
    prs.add_argument('--path',dest='bandPath',
                     help='Select the band structure path (--path="G-H-N")')
    prs.add_argument('--width-ratios',dest='widthRatios',
                     help='Set the relative sizes of the band structure and \
                     density of states plots (--width-ratios=Nband-Ndos)')
    args = prs.parse_args()
    return args

def find_file(file_query):
    """
    Use the glob module to find files with names that match the query.
    """
    return glob.glob(file_query)

def find_sympts(lines):
    """
    Given the contents of a Gnuplot script find the symmetry points

    To draw the band structure the symmetry points along the 
    path of k-points is needed. The GW code stores this information 
    in the Gnuplot script that draws the band structure. Hence 
    this routine parses the Gnuplot script and returns the symmetry
    points as a list of tuples, where every tuple contains 
    ( symmetry_point_label, x-coord) where symmetry_point_label
    is a string representing the name of the symmetry point, and
    x-coord is a floating point value representing a distance along
    the k-point path in arbitrary units.

    Note that in the Gnuplot script the symmetry points are stored as
    '" G " 0.0000', which means that the split command on a single point
    returns ['"','G','"','0.00000'] so we want to extract elements 1 and 3.
    """
    path = []
    for line in lines:
        if 'xtics' in line:
            begin = line.index('(')
            end   = line.index(')')
            points = line[begin+1:end]
            points = points.split(',')
            for point in points:
                point = point.split()
                path.append((point[1],float(point[3])))
    return path

def find_subpath(path,subpath):
    """
    Find a subpath in the band structure path

    Given the overall band structure as a list of tuples of labels and x-values
    and given a list of symmetry point labels, return the section of the path
    that matches the subpath. For example, if the band structure path is

        [("G",0.0),("H",0.18),("P",0.34),("G",0.50),("N",0.63)]

    and subpath is 

        ["H","P","G"]

    then return

        [("H",0.18),("P",0.34),("G",0.50)]

    If subpath does not correspond to a consecutive subsequence of path return
    None. For example is subpath is

        ["H","G"]
    """
    lpath = len(path)
    lsubpath = len(subpath)
    ldiff = lpath - lsubpath
    for ii in range(0,ldiff+1):
        jj = 0
        retpath = []
        match = True
        while jj < lsubpath and match:
            tuple = path[ii+jj]
            (pathlabel,pathval) = tuple
            subpathlabel = subpath[jj]
            match = match and (pathlabel == subpathlabel)
            retpath.append(tuple)
            jj = jj+1
        if match:
            return retpath
    return None

def find_energy_range(lines):
    """
    Given the contents of a Gnuplot script find the energy range

    To draw the band structure an initial energy range is handy. 
    The energy range is given in a Gnuplot script so here the
    script is parsed and the energy range returned as a tuple
    (emin,emax) where emin and emax are given in eV.
    """
    for line in lines:
        if 'yrange' in line:
            begin = line.index('[')
            end   = line.index(']')
            energies = line[begin+1:end]
            energies = energies.split(':')
    return (float(energies[0]),float(energies[1]))

def execute_with_arguments(args):
    print(args)
    properties = {}
    if args.widthRatios:
        data = args.widthRatios.split("-")
        properties["width_ratios"] = [int(data[0]),int(data[1])]
    fig = fbd(args.plotType,properties)
    if args.plotType == "band" or args.plotType == "bd":
        files = find_file("*.gnu")
        if len(files) == 0:
            print("Band structure plot requested but no data was found.")
            sys.exit(10)
        f = open(files[0],'r')
        lines = f.readlines()
        f.close()
        (emin,emax) = find_energy_range(lines)
        if args.energyScale:
            elist = args.energyScale.split(":")
            emin = float(elist[0])
            emax = float(elist[1])
        fig.set_energy_range(emin,emax)
        path = find_sympts(lines)
        if args.bandPath:
            subpath = args.bandPath.split("-")
            path = find_subpath(path,subpath)
            if not path:
                print("Invalid subpath specified")
                sys.exit(20)
        fig.set_path(path)
        if args.statesType == "dft":
            files = find_file("*dft_band*.dat")
        elif args.statesType == "qp":
            files = find_file("*qp_band*.dat")
        if len(files) == 0:
            print("Band structure plot requested but no data was found.")
            sys.exit(10)
        f = open(files[0],'r')
        xlist = []
        ylist = []
        for line in f:
            xy = line.split()
            if len(xy) != 0:
                xlist.append(float(xy[0]))
                ylist.append(float(xy[1]))
            else:
                fig.add_band(xlist,ylist)
                xlist = []
                ylist = []
        f.close()
    if args.plotType == "dos" or args.plotType == "bd":
        if args.statesType == "dft":
            files = find_file("*dft.dos")
        elif args.statesType == "qp":
            files = find_file("*qp.dos")
        if len(files) == 0:
            print("Density of states plot requested but no data was found.")
            sys.exit(20)
        f = open(files[0],'r')
        elist = [] # energies
        alist = [] # alpha-density of states
        blist = [] # beta-density of states
        tlist = [] # total-density of states
        for line in f:
            dat = line.split()
            elist.append(float(dat[0]))
            alist.append(float(dat[1]))
            blist.append(float(dat[2]))
            tlist.append(float(dat[3]))
        f.close()
        emin = elist[0]
        emax = elist[-1]
        if args.energyScale:
            list = args.energyScale.split(":")
            emin = float(list[0])
            emax = float(list[1])
        if args.plotType == "dos":
            fig.set_energy_range(emin,emax)
        if args.nstatesScale:
            nmax = float(args.nstatesScale)
            nmin = 0.0
            fig.set_number_range(nmin,nmax)
        else:
            ilist = [i for i,e in enumerate(elist) if e > emin and e < emax]
            nmin = 0.0
            nmax = 0.0
            for i in ilist:
                nmax = max(nmax,tlist[i])
            fig.set_number_range(nmin,nmax)
        fig.add_dos(elist,alist)
        fig.add_dos(elist,blist)
        fig.add_dos(elist,tlist)
    plt.show()

def main():
    """
    Launch the application
    """
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
