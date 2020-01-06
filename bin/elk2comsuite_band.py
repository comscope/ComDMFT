#!/usr/bin/env python
# A Python script that transforms an Elk generated band structure file to 
# a format that matches the Comsuite band structure format.
#
versionstr = '%(prog)s version 0.0'

import sys

def parse_arguments():
    """
    Parse command line arguments for the elk2comsuite_band.py script when run
    from the command line.
    """
    from argparse import ArgumentParser
    prs = ArgumentParser()
    prs.add_argument("--version",action='version',version=versionstr)
    prs.add_argument("-e","--elk",dest="input",default="BAND.OUT",
                     help="The name of the file containing the Elk band structure data")
    prs.add_argument("-c","--comsuite",dest="output",default="band_elk.dat",
                     help="The file name for the Comsuite format band structure output data")
    prs.add_argument("--xmax",dest="xmax_cms",
                     help="The maximum X-coord in the Comsuite band structure")
    prs.add_argument("-m","--match",dest="match",default=None,
                     help="The file name for the Comsuite band structure file to match the Elk band structure to")
    args = prs.parse_args()
    return args

def convert_and_write_band(band,xmax_cms,ofile):
    """
    Convert the X-coordinates of the band and write the resulting band
    to the outputfile 'ofile'.
    """
    length = len(band)
    if length == 0:
        return
    (Xmax,dummy) = band[length-1]
    for point in band:
        (X,E) = point
        X = xmax_cms*X/Xmax
        E = 27.2114*E
        ofile.write("%22.15f %22.15f\n" % (X,E))

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

def execute_with_arguments(args):
    """
    Convert the Elk band structure to a Comsuite band structure.

    The conversion is rather simple. The Elk band structure contains a linear
    array of data points X,Y where X refers to a point a long the high symmetry path
    and Y is the energy of the band at the corresponding k-point. The differences
    between the Elk and the Comsuite band structure are
    - that with Elk X runs from 0 to Xmax, whereas with Comsuite X run from 0 to 1;
    - Elk band structures are reported in units of Hartree whereas Comsuite uses 
      energies in the unit of eV.
    So we need to read the data corresponding to a band, convert
    the X coordinates as X_cms = X_i/Xmax, and then write the data out to the 
    Comsuite band structure file. Then rinse and repeat until all bands have been
    transformed.
    Also we need to convert the Energy from Hartree to eV which means scaling
    by a factor 27.2114.
    """
    inputFilename = args.input
    outputFilename = args.output
    if args.xmax_cms:
        xmax_cms = float(args.xmax_cms)
    elif args.match:
        refFilename = args.match
        reffile = open(refFilename,'rb')
        xmax_cms = extract_xmax(reffile)
    else:
        error = "No information about the new Xmax.\nTry --xmax or --match"
        print(error)
        return
    if inputFilename == "-":
        inputfile = sys.stdin
    else:
        inputfile = open(inputFilename,'rb')
    if outputFilename == "-":
        outputfile = sys.stdout
    else:
        outputfile = open(outputFilename,'w')
    band = []
    for line in inputfile:
        words = line.split()
        if len(words) == 0:
            convert_and_write_band(band,xmax_cms,outputfile)
            outputfile.write("\n")
            band = []
        else:
            band.append((float(words[0]),float(words[1])))
    convert_and_write_band(band,xmax_cms,outputfile)
    

def main():
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
