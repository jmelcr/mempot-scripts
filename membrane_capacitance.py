#!/usr/bin/python
#  -*- mode: python; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*-
#

__doc__ = """
This module contains classes and routines for analysis charge profiles and the resulting field and potential

------------------------------------------------------------
 Made by Joe,  Last edit 2015/08/17
------------------------------------------------------------
------------------------------------------------------------
          """

#import sys
import numpy as np
import read_xvg_calc_mean as rxvg   # Joe's library containing routines for reading .xvg files
import membrane_potential as mempot # Joe's library containing routines for integrating membrane potential of double-membrane systems
import IPython
#import cPickle
from scipy.signal import butter, filtfilt

# "constants" are capitalized
NDIM = 3
GRO2uFCM = 16.02 # conversion factor for Gromacs units e/V.nm2 --> uF/V.cm2


class Densities():
   ''' class for storage of density profiles from molecular simulation
       Originally designed to be used along with double membrane systems (e.g. potential integrator assumes that).
       Optional arguments during init:
          Contains an instance of q_profile as self.qp (stores charge, field, potential and capacitance profiles and so).
             As an option, it's possible to provide q_profile filename instead of q_profile instance and call q_profile init
          Stores arbitrary number of density profiles contained in one density file (Gromacs g_density output) as one big numpy array.
       TODO: Header, Gromacs (self.comments) and Grace (self.grace) comments, are stored in respective lists of lines.
   '''

   # class variables
   Xref = 14.8 # average z-box-length [nm] from simulations 

   def __init__(self, fname=None, qprofile=None, qp_fname=None):
       self.fname = fname       
       # Load densities from fname to self.dens
       if not fname==None:
          if isinstance(fname, str):
             #load the densities
             self.dens = self.load_dens()
             if self.dens.shape[1]==10:
                print "Splitting the loaded density into 10 individual densities assuming a particular order."
                (self.x, self.oil, self.headgr, self.wat, self.na, self.cl, self.phos, self.chol, self.carbonyl, self.ester) = np.hsplit(self.dens, self.dens.shape[1])
                self.scale_x()
          else:
             raise RuntimeError, "Parameter filename %s provided, but it is not a string. Densities not loaded. " % (fname,)
       else:
          print "Filename not provided, no densities loaded."

       if qp_fname!=None and qprofile!=None:
          print "Both q_profile filename and instance provided -- can't load both -- choose one!"
          self.qp = None
       else:
          # load q_profile instance from qprofile optional argument if provided
          if isinstance(qprofile, q_profile):
             self.qp = qprofile
          else:
             if qprofile!=None:
                print "Provided charge profile is not an intance of q_profile. Loading nothing. "
                self.qp = None
          # create a new q_profile instance based on qp_filename if provided
          if isinstance(qp_fname, str):
             self.qp = q_profile(qp_fname)
          else:
             if qp_fname!=None: raise RuntimeError, "Parameter filename %s provided, but it is not a string. Densities not loaded. " % (fname,)

   def load_dens(self):
       """ Returns multi-dim numpy array containing all density profiles from file self.fname """
       try: open(self.fname,"r")
       except: 
          raise IOError, "File %s could not be opened." % (self.fname,)
       return np.array( rxvg.lines_to_list_of_lists(rxvg.read_xvg(self.fname)) )

   def scale_x(self):
       """ Rescales the spacial coordinates self.x so that it spans the average box-z-length as set in Densities.Xref """
       self.x *= Densities.Xref/self.x.max()  # scale to the average dimension Xref (class variable)



class q_profile():
   ''' class for storage of charge, capacitance and other profiles of membrane (and possibly other) simulations
       the first instance to be initialized also sets the reference (class variables)
       some degree of consistency is assumed throughout all instances (like no. slices and so)
   '''

   # class variables
   q_ref_fname = None
   q_ref       = None # e/nm3
   x_ref       = None # nm
   slices      = 0
   Xref        = 14.8 # average z-box-length [nm] from simulations 
   

   def __init__(self, fname):
       if not isinstance(fname, str):
          print "Init argument must be a string -- path to file with xvg data."
       else:
          self.fname   = fname
          self.x       = None 
          self.q_noisy = None
          self.q       = None
          self.load_x_q_as_nparray() 
          self.scale_x()                     # scale to the average dimension Xref (class variable)
          self.q = noisefilt(self.q_noisy, critfreq=0.05)   # get smoother q-profile
          self.qcorr()                       # correct total charge to 0

          # get the field and potential profiles and the Voltage-drop
          (self.f, self.p, self.V) = mempot.int_q_twice(self.q_noisy, self.x, fcorr_type="cent_plat_zero")  # V/nm, V, Volts

          if q_profile.slices==0: 
             print """This, file %s, is the first instance of class q_profile. 
                      Setting it as a reference (class variable).""" % (self.fname,)
             self.set_q_ref()
          else:
             # old way of obtainind V: self.V  = mempot.get_potential(fname, q_profile.Xref, fcorr_type="cent_plat_zero")  # Volts
             self.dq = None
             self.c  = None
             self.calc_dq_c()  # get dq and capacity

   def scale_x(self):
       """ Rescales the spacial coordinates self.x so that it spans the average box-z-length as set in q_profile.Xref """
       self.x *= q_profile.Xref/self.x.max()  # scale to the average dimension Xref (class variable)

   def qcorr(self):
       """ Corrects the q-profile so that the total charge is 0. """
       self.q -= self.q.mean()

   def calc_dq_c(self):
       """ (re)calculates dq and capacitance profiles """
       if len(self.q)!=q_profile.slices:
          print """The number of slices %d of the data from file %s do not agree with 
                   the number of slices %d of the reference %s! Processing stopped.""" %\
                   (len(self.x), self.fname, q_profile.slices, q_profile.q_ref_fname)
       else:
          self.dq = self.q - q_profile.q_ref
          self.c  = GRO2uFCM * self.dq / self.V # uF/V.cm2

   def symmetrize_q(self, symmtype):
       """ function that symmetrizes the charge profile according to the specified type of symmetry: mirror/transl 
           Does not change the q-profile in-place.
       """
       symmtypes = ('mirror', 'transl')
       if not symmtype in symmtypes:
          raise RuntimeError, 'specified symmetry type %s is misspelled or not supported' % (symmtype,)
       else:
          l = len(self.q)
          if symmtype==symmtypes[0]:
             return ( self.q[:l/2] + np.flipud(self.q[l/2:]) )/2
          elif symmtype==symmtypes[1]:
             return ( self.q[:l/2] +           self.q[l/2:]  )/2
          else:
             raise RuntimeError, "Something unexpected happened!!"

   def set_q_ref(self):
       """ (re)sets a reference q-profile (e.g. zero voltage) for this class """
       q_profile.q_ref_fname = self.fname
       q_profile.x_ref       = self.x
       q_profile.q_ref       = self.q
       q_profile.slices      = len(q_profile.x_ref)

   def load_x_q_as_nparray(self):
       """ loads x and q profile from a xvg file self.fname into self.x,q in-place """
       try: open(self.fname,"r")
       except: 
          print "File %s could not be opened." % (self.fname,)
          raise IOError
       try: self.x, self.q_noisy = np.split(np.array( rxvg.lines_to_list_of_lists(rxvg.read_xvg(self.fname)) ).ravel(order="F") , 2 )
       except: raise RuntimeError, "File %s could not be loaded, possibly because of wrong formatting." % (self.fname,)



def noisefilt(data, order=8, critfreq=0.022):
    """ function that smoothens given data by applying butterworth's noise filter """
    b, a = butter(order, critfreq)  # beautiful smoothening filter from scipy.signal; its init
    return filtfilt(b, a, data)  





if __name__ == '__main__':
   # say "hello"
   print __doc__

   # just a shorthand for a long class name
   print "Aliasing classname q_profile to qp"
   qp = q_profile

   IPython.embed()

   #print "Pickling reslist ..."
   #with open("reslistwtraj.cpickle", "w") as file: cPickle.dump(reslist,file)
   # depickle with reslist = cPickle.load(file)

