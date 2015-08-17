#!/usr/bin/env python
"""
  Routines and functions for :
  Read in .xvg files as created by Gromacs
  (contain #comments and @Grace-headers)
"""

import sys, os # math # os.path, string, re
import numpy as np

def read_xvg(infilename):
   """Reads in .xvg file as created by Gromacs without #comments and @Grace-headers 
      returns list of strings (lines)
   """
   #check filename length
   if len(infilename)==0 :
      print 'Provide filename (none given) \n' ;  sys.exit(1)

   try: infile = open(infilename, 'r') # open file for reading 
   except: print '%s:  reading file "%s" failed; does it exist?' % (sys.argv[0],infilename) ; sys.exit(1)

   #list of lines without # or @ comments (just values - floats)
   return  [ n for n in infile.readlines() if not (n.startswith("#") or n.startswith("@")) ]
   infile.close()   # close the inputfile


def lines_to_list_of_lists(lines):
   """Converts lines (strings) into a list of lists containing numbers-floats.
      Suitable for converting into numpy arrays
   """
   #create and fill list of lists with float numbers
   return [ [float(strnum) for strnum in line.split()] for line in lines ]

def PN_angle(angles):
   """Provides the mean value of P-N angle with std. as a dictionary
      claculated from the input np.array angles.
      The angle is defined as between 0-90deg, shifting is made accordingly
   """
   angles_mean = angles.mean()
   if angles_mean>90.0: return { 'mean': 180.0-angles_mean, 'std': angles.std() }
   else:                return { 'mean':       angles_mean, 'std': angles.std() }


def main():
   from optparse import OptionParser
   # help message is automatically provided
   # type=string, action=store is default
   parser = OptionParser()
   parser.add_option('-f', '--infile', dest='infilename', help='input filename of *.xvg file', default='')
   parser.add_option('-a', '--pnang' ,   action='store_true', dest='doPNang',   help='calculate the average value of all columns (useful for P-N angle, hence the nameing)', default=False)
   parser.add_option('-n', '--nomeans' , action='store_true', dest='dontmeans', help='do not calculate the average value of individual columns', default=False)
   options, args = parser.parse_args(sys.argv[1:])  #argument 'sys.argv..' not necessary, it's the default

   # a bit of abstractions
   doPNang   = options.doPNang
   dontmeans = options.dontmeans

   #create a numpy array with the contents of .xvg file;  nested: file -> lines -> list-of-lists(floats) -> np.array
   arrvals = np.array( lines_to_list_of_lists(read_xvg(options.infilename)) )

   if not dontmeans:
      print 'Mean values: ', arrvals.mean(axis=0)
      print 'Stdevs     : ', arrvals.std(axis=0)
      print

   if doPNang:
      t, angles = np.hsplit(arrvals,[1])
      PNang = PN_angle(angles)
      print 'P-N angle = %4.2f +- %4.2f deg' % (PNang['mean'],PNang['std'])

   #import IPython
   #IPython.embed()

   sys.exit(0)

if __name__ == '__main__':
    main()
