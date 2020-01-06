#!/usr/bin/env python
# A Python script that transforms a Wien2K generated band structure .energy file to 
# a format that matches the Comsuite band structure format.
#
versionstr = '%(prog)s version 0.0'

import sys

def parse_arguments():
    """
    Parse command line arguments for the wien2k2comsuite_band.py script when run
    from the command line.
    """
    import os
    from argparse import ArgumentParser
    cwd      = os.getcwd()
    folder   = os.path.basename(cwd)
    # In Wien2K the name of the CASE has to equal the folder name
    case     = folder 
    case_in  = case+".energy"
    case_scf = case+".scf"
    case_lat = case+".outputkgen"
    prs = ArgumentParser()
    prs.add_argument("--version",action='version',version=versionstr)
    prs.add_argument("-w","--wien2k",dest="input",default=case_in,
                     help="The name of the file containing the Wien2K band structure data")
    prs.add_argument("-f","--fermi",dest="fermi",default=case_scf,
                     help="The name of the file containing the Wien2K Fermi energy")
    prs.add_argument("-l","--latvec",dest="latvec",default=case_lat,
                     help="The name of the file containing the Wien2K lattice vectors")
    prs.add_argument("-c","--comsuite",dest="output",default="band_wien2k.dat",
                     help="The file name for the Comsuite format band structure output data")
    prs.add_argument("--xmax",dest="xmax_cms",
                     help="The maximum X-coord in the Comsuite band structure")
    prs.add_argument("-m","--match",dest="match",default=None,
                     help="The file name for the Comsuite band structure file to match the Wien2K band structure to")
    args = prs.parse_args()
    return args

def read_wien2k_lattice_vectors(filename):
    """
    Read the Wien2K lattice vectors and return
    the vectors
    - a
    - b
    - c
    In Wien2K the lattice vectors are stored in a file with the extension
    .outputkgen, and are listed as
    R1 =
    R2 = 
    R3 =
    so we need to be looking for lines that start with that content.
    """
    a = []
    b = []
    c = []
    line = ""
    with open(filename,"r") as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            if line.find("R1 =") > -1:
                point = line.split()
                a = [ float(point[2]), float(point[3]), float(point[4]) ]
            if line.find("R2 =") > -1:
                point = line.split()
                b = [ float(point[2]), float(point[3]), float(point[4]) ]
            if line.find("R3 =") > -1:
                point = line.split()
                c = [ float(point[2]), float(point[3]), float(point[4]) ]
    return (a,b,c)

def read_wien2k_G_vectors(filename):
    """
    Read the Wien2K G-vectors and return
    the vectors
    - G1
    - G2
    - G3
    In Wien2K the G-vectors are stored in a file with the extension
    .outputkgen, and are listed under
    G1 G2 G3
    so we need to be looking for lines that start with that content.
    """
    a = []
    b = []
    c = []
    line = ""
    with open(filename,"r") as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            if line.find("G1        G2        G3") > -1:
                for ii in range(0,3):
                    line = fp.readline()
                    point = line.split()
                    a.append(float(point[0]))
                    b.append(float(point[1]))
                    c.append(float(point[2]))
                break
    return (a,b,c)

def read_wien2k_bands(filename,a,b,c):
    """
    Read the Wien2K band structure data and return
    - List of 1D positions along the path
    - A dictionary with the data. The keys consist of tuples of
      the point number and the band number
    - The minimum number of bands at any point
    - The number of k-points
    """
    import math
    numpoints = 0
    minbands  = 1000000
    points    = []
    data      = {}
    line      = ""
    with open(filename,"r") as fp:
        #
        # Skip over initial stuff, looking for the Gamma point
        #
        while line.find("G") == -1:
            line = fp.readline()
            if not line:
                raise Exception("ERROR: Not a valid band structure file: {}".format(filename))
                
        #
        # Read bands
        #
        while True:
            #
            # Read k-point line
            #
            numpoints += 1
            nbands = int(line[74:80])
            minbands = min(minbands,nbands)
            sx = float(line[ 0:19])
            sy = float(line[19:38])
            sz = float(line[38:57])
            points.append((sx,sy,sz))
            #
            # Read band energies
            #
            for ii in range(0,nbands):
                line = fp.readline()
                items = line.split()
                iband = int(items[0])
                energy = float(items[1])
                data[(numpoints,iband)] = energy
            line = fp.readline()
            if not line:
                break
    prev = (0.0,0.0,0.0)
    xcoords = []
    r = 0.0
    for current in points:
        # Current is a point in fractional coordinates (xf,yf,zf)
        # this needs to be converted to Cartesian coordinates (xc,yc,zc)
        (xf,yf,zf) = current
        xc         = xf*a[0]+yf*b[0]+zf*c[0]
        yc         = xf*a[1]+yf*b[1]+zf*c[1]
        zc         = xf*a[2]+yf*b[2]+zf*c[2]
        (xp,yp,zp) = prev
        current_c  = (xc,yc,zc)
        d=(xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp)
        d=math.sqrt(d)
        r+=d
        xcoords.append(r)
        prev = current_c
    return (xcoords,data,minbands,numpoints)

def read_wien2k_fermi(filename):
    """
    Read the Wien2K .scf file and return the Fermi energy in units of Rydberg
    The Fermi energy is given on lines starting with ":FER  :". The .scf file
    may contain multiple of such lines. We want to keep the last one and 
    extract the Fermi energy from that.
    """
    line       = ""
    line_fermi = ""
    with open(filename,"r") as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            if line.find(":FER  :") > -1:
                line_fermi = line
    data = line_fermi.split("=")
    fermi_energy = float(data[1])
    return fermi_energy

def extract_xmax(mfile):
    """
    Extract the maximum X-coordinates of the band structure and return
    that value.
    """
    xmax = 0.0
    for line in mfile:
        words = line.split()
        if len(words) == 0:
            return xmax
        else:
            xmax = float(words[0])
    return xmax

def write_comsuite_bands(filename,xmax_cms,xcoords,data,enfermi,minbands,numpoints):
    """
    Write the data to the comsuite file
    """
    fp = open(filename,"w")
    xmax = xcoords[numpoints-1]
    for iband in range(1,minbands+1):
        for ipoint in range(1,numpoints+1):
            energy = (data[(ipoint,iband)]-enfermi)*13.6056980659 
            xcoord = xcoords[ipoint-1]*xmax_cms/xmax
            fp.write("%14.8e %14.8e\n" % (xcoord,energy))
        if iband != minbands:
            fp.write("\n")
    fp.close()

def execute_with_arguments(args):
    """
    Convert the Wien2K .energy file to the Comsuite band structure format.

    The conversion involves a number of steps:
    1. The energies in the .energy file are given in Rydberg and
       in the comsuite band structure we need them in eV.
    2. The .energy file must be read and it lists
       a. the k-point, the k-point name, an integer, the number of bands, the weight
       b. for each band: the band number, the band energy
    3. Find the minimum number of bands
    4. For each band: 
          For each point write: the x-coordinate, the band energy
    """
    infilename    = str(args.input)
    outfilename   = str(args.output)
    fermifilename = str(args.fermi)
    latfilename   = str(args.latvec)
    if args.xmax_cms:
        xmax_cms = float(args.xmax_cms)
    elif args.match:
        matchfilename = str(args.match)
        mfile = open(matchfilename,'rb')
        xmax_cms = extract_xmax(mfile)
    else:
        error = "No information about the new Xmax.\nTry --xmax or --match"
        print(error)
        return
    (a,b,c) = read_wien2k_G_vectors(latfilename)
    (xcoords,data,minbands,numpoints) = read_wien2k_bands(infilename,a,b,c)
    enfermi = read_wien2k_fermi(fermifilename)
    write_comsuite_bands(outfilename,xmax_cms,xcoords,data,enfermi,minbands,numpoints)

def main():
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
