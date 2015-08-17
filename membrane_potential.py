#!/usr/bin/env python

__doc__= """
  Calculate elstat. potential profile for double-membrane systems
  special symmetry assumed (2 membs with their COM centered)
  trapez integration scheme from numpy taken

  Stepan Timr and Josef Joe Melcr 12/4/2015
         """

import sys
import numpy as np
import math as mt
import read_xvg_calc_mean as rxvg  # Joe's library containing routines for reading .xvg files
#import IPython  # good for debugging, invoke IPython shell any time by calling IPython.embed()
import numpy as np

#Global Constants
EPS0 = 8.854187817e-12
Q = 1.602176565E-019
KON = Q/EPS0


def find_threshold(Array,start,threshold,reverse=1):
   """ Finds the first occurence of a value larger then threshold 
       in Array to the right (larger index) from start. (or left if reverse=-1)
       Returns the integer index of the cell last before the threshold
   """
   fast_step = 50
   i = start
   while abs(Array[i])<threshold: i+=fast_step*reverse
   while abs(Array[i])<threshold: i+=1*reverse
   return i


def integrate_from_centre(ArrayX,ArrayY):
   """ returns numpy array containing integral of function ArrayY of ArrayX
       uses trapezoid rule for integration
   """
   ArrLen = len(ArrayX)
   centre = int(ArrLen/2)
   Integral = np.zeros(ArrLen)
   #Integrate to the left
   for i in np.arange(0, centre):
       Integral[i] = np.trapz(np.flipud(ArrayY[i:centre+1]), x=np.flipud(ArrayX[i:centre+1]))
   #Integrate to the right
   for i in np.arange(centre, ArrLen):
       Integral[i] = np.trapz(ArrayY[centre:i+1], x=ArrayX[centre:i+1])

   return Integral


def int_q_twice(Charge, X, fcorr_type="cent_plat_zero", field_intensity=0.0, qcorr=True):
   """ Function int_q (integrate charge) returns the field and potential profile for simulations with two membranes.
       It uses trapezoid integration scheme,
       Starts the integration from the centre,
       applies different corrections for field ("no_work", "cent_plat_zero", "all_plat_zero") and charge (q_tot->0),
       allows the field_intensity (V/nm) to be set, so that the corrections can be used also for the constant-Efield simulations.
       It returns field and potential profiles (numpy arrays) and the average voltage drop (all as a tuple) in standard units (V/nm, V, V).
   """
   #Corrections
   #print "Corrections:"
   #Charge correction
   if qcorr:
      ChargeSum = np.mean(Charge)
      #print "Charge Correction: %6.4e e/nm^3" % (ChargeSum,)
      Charge -= ChargeSum

   #Integration from the centre of the box
   #FIELD [V/nm]
   #print "Integrating charge ..."
   Field = integrate_from_centre(X,1e9*KON*Charge)

   #get estimates for furhter potential-plateau search
   wat_thick_min = 1.0 #nm, minimal assumed thickness of water slab
   #find corresponding indices in X
   bounds_int_min = ( X.searchsorted(-wat_thick_min/2, side="right") , X.searchsorted( wat_thick_min/2, side="left") )
   #Calc fluctuation in the region - estimate for furhter potential-plateau search
   FieldStd = Field[bounds_int_min[0]:bounds_int_min[1]].std()
   FieldFluctThreshold=0.06  # just an estimate of an uncommonly large fluctuation
   if FieldStd >= FieldFluctThreshold: print "Beware! Large fluctuations in Field in central plateau, std = %6.4f V/nm" % (FieldStd,)

   #find the central plateau of potential (where water is and the potential is ~constant)
   #water thickness in no. cells of an array
   centre_int = int(len(Field)/2)
   FieldFluct_thres_factor = 5.0
   confidence_factor=0.75
   wat_rbound = find_threshold(Field, centre_int, FieldFluct_thres_factor*FieldStd)
   #wat_bounds = [ find_threshold(Field, centre_int, FieldFluct_thres_factor*FieldStd, reverse=i) for i in (1,-1) ] #for non-symmetric profile
   #print "water bounds indices: ", wat_bounds
   wat_thickness_int = int( 2*abs(wat_rbound-centre_int) * confidence_factor )  #assuming symmetry
   #wat_thickness_int = int( abs(wat_bounds[0]-wat_bounds[1]) * confidence_factor )
   #print "Central region after scaling: ", X[centre_int-int(wat_thickness_int/2):centre_int+int(wat_thickness_int/2)]
   
   #field_corrections = { "no_work"        : no_work_corr,
   #                      "cent_plat_zero" : cent_plat_zero_corr,
   #                      "all_plat_zero"  : all_plat_zero_corr   }

   # branched if for a CASE equivalent (simpler then using dictionary and functions)
   if fcorr_type=="no_work":
      FieldMean = Field.mean()
      #print "Field Correction: %6.4e V/nm" % (FieldMean,)
      Field -= FieldMean   
   elif fcorr_type=="cent_plat_zero":
      #print "Correcting field so that the central plateau is zero intensity on average"
      fieldCorrection_centre = Field[centre_int-int(wat_thickness_int/2):centre_int+int(wat_thickness_int/2)].mean()
      #print "Field correction: %6.4e V/nm" % (fieldCorrection_centre,)
      Field -= fieldCorrection_centre
   elif fcorr_type=="all_plat_zero":
      #print "Correcting field so that all plateaus (water regions - not true!!!!) have zero intensity on average"
      selection = np.select( [abs(Field)<FieldFluctThreshold], [1])
      lowfield = np.extract(selection,Field)
      #print "Field correction: %6.4e V/nm" % (lowfield.mean(),)
      Field -= lowfield.mean()


   # add up the specified homog. field intensity if non-zero (works well with fcorr_type=="no_work")
   if abs(field_intensity)>0.0 :
      print "Correcting field with homogenous intensity %6.4f V/nm" % (field_intensity,)
      Field += field_intensity   # shift the corrected field by the value of the external homog. field


   #POTENTIAL FROM UNBIASED FIELD [V]
   #print "Integrating field ..."
   Potential = integrate_from_centre(X,-Field)

   #calculate mean values for the plateaus in the double-bilayer setting
   Pot_centre = Potential[centre_int-int(wat_thickness_int/2):centre_int+int(wat_thickness_int/2)].mean()
   # start=25 for cropping edges that tend to have large overshoots
   Pot_left  =           Potential[ 20:int(find_threshold(Field,            25, FieldFluct_thres_factor*FieldStd)*confidence_factor)].mean()
   Pot_right = np.flipud(Potential)[20:int(find_threshold(np.flipud(Field), 25, FieldFluct_thres_factor*FieldStd)*confidence_factor)].mean()
   print
   print "V(left -centre) = %5.1f mV" % ((Pot_centre-Pot_left)*1000.0,)
   print "V(right-centre) = %5.1f mV" % ((Pot_centre-Pot_right)*1000.0,)
   if not abs(Pot_left)-abs(Pot_right)*1000.0<=10.0 : print "Beware: potential drop to the right and to the left is too different (>10mV)." # in mV
   #Potential_drop_aver = (2*Pot_centre-Pot_right-Pot_left)*1000.0/2.0
   Potential_drop_aver = (abs(Pot_centre-Pot_right) + abs(Pot_centre-Pot_left)) /2.0 # avg of the two drops (positive value)
   print "V( average    ) = %5.1f mV" % (Potential_drop_aver*1000.0,)
   # Statistical dispersion turns out to be negligible
   
   return (Field, Potential, Potential_drop_aver)




def get_potential(charge_file, Xref=14.8, crop=0, fcorr_type="cent_plat_zero", field_intensity=0.0, qcorr=True, plot=False, **kwargs):
   """ Function get_potential returns the voltage drop (in V) across the membrane from simulations with two membranes.
       It uses function int_q for integration, which starts the integration from the centre,
       Scales the box length to the reference value Xref (approximate average from the simulations),
       allows cropping from left/right sides,
       It returns the absolute (positive) average value of the voltage drop (in V) taken from drops centre-left and centre-right.
   """
   # Xref is in nm, value taken from the average of the longest simulations


   #create a numpy array with the contents of .xvg file;  nested: file -> lines -> list-of-lists(floats) -> np.array -> flatten(ravel) -> split in 2
   X, Charge = np.split(np.array( rxvg.lines_to_list_of_lists(rxvg.read_xvg(charge_file)) ).ravel(order="F") , 2 )
   #                            V---  list of lists     <-----------    list of lines
   #                      numpy array    ---------------------------------------------->    flattened  ...
   #     .... and split     ....       .......        in           .............                 ..........   two!

   #Scale potential to the average Z-length
   box_fluct_thres=0.2
   if abs(Xref-X.max())>box_fluct_thres:
      print "Scaling box length from %4.2f to %4.2f nm (average length from the longest simulation)" % (X.max(),Xref)
      print "This is quite a large difference. \n"
   X *= Xref/X.max()

   #Crop box ends
   if crop>=0 :
      X = X[crop:len(X)-crop]
      Charge = Charge[crop:len(Charge)-crop]

   #shift Box-Z centre to 0.0
   #assumes that Box starts at z=0 and ends at some z=Length
   X -= X[int(len(X)/2)]

   # total charge and dipole of the system
   #Qtot = np.trapz(Charge, x=X)
   #print "Net charge: %5.4f q/nm^2" % (Qtot,)
   #Ptot = np.trapz(np.multiply(Charge,X), x=X)
   #print "Net dipole: %5.4f D/nm^2 \n" % (1/0.393456*1.889725989e1*Ptot,) # Magic constants by Stepan

   (Field, Potential, V) = int_q_twice(Charge, X, fcorr_type, field_intensity, qcorr)



   if plot:
      import matplotlib.pyplot as plt
      #Read-in densities
      densities_file = kwargs.get("densfile")
      Xdens, phos, tail, wat = np.genfromtxt(densities_file, comments="#", skip_header=24, unpack=True)
      #Scaling to the same average Z-length (see above)
      Xdens *= Xref/Xdens.max()
      Numdens = len(Xdens)
      cendens = int(Numdens/2)
      Xdens_shift = Xdens - Xdens[cendens]

      #Plots
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      ax1.plot(Xdens_shift, phos/256, 'b', label="P")
      ax1.plot(Xdens_shift, tail/512, 'm', label="Terminal C")
      ax1.plot(Xdens_shift, wat/10128, 'c', label="Water")
      ax1.set_xlabel('Z [nm]')
      ax1.set_ylabel(r'Number density normalized to unity [nm$^{-3}$]')
      plt.legend(loc="upper left")


      ax2 = ax1.twinx()
      ax2.plot(X, Potential, 'k', label="Electrostatic potential", lw=2)
      ax2.plot([X[0],X[len(X)-1]], [0,0], '--', label="Zero", lw=2)
      ax2.set_ylabel('Potential [V]')
      plt.legend(loc="upper center")
      plt.show()


      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      ax1.plot(Xdens_shift, phos/256, 'b', label="P")
      ax1.plot(Xdens_shift, tail/512, 'm', label="Terminal C")
      ax1.plot(Xdens_shift, wat/10128, 'c', label="Water")
      ax1.set_xlabel('Z [nm]')
      ax1.set_ylabel(r'Number density normalized to unity [nm$^{-3}$]')
      plt.legend(loc="upper left")

      ax2 = ax1.twinx()
      ax2.plot(X, Field, 'k', label="Intensity of electric field", lw=2)
      ax2.set_ylabel('Field [V/nm]')
      plt.legend(loc="upper center")
      plt.show()

   return V



#--------------------------
if __name__ == '__main__':
   """ Calculates the electrostatic potential profile for a given charge density profile.
       Wokrs the best with double bilayers (integration starts from the centre). 
       This is the abstracted main programme that contains the parser of CMDline arguments. 
       The actual calculation is provided by the function get_potential, which can be loaded as a module to other Python codes.
   """
   import argparse

   parser = argparse.ArgumentParser()
   parser.add_argument("charge", help="charge.xvg from get_pot_inputs.sh")
   parser.add_argument("densities", help="density_profiles.xvg from get_pot_inputs.sh")
   parser.add_argument("-crop", type=int, help="Crop this number of slices from the beginning and end of the box", default=0)
   parser.add_argument("-fcorr",       action='store_true', help="field correction to zero on average")
   parser.add_argument("-fcorrcent",   action='store_true', help="field correction to zero on average in the central plateau")
   parser.add_argument("-fcorrplat",   action='store_true', help="field correction to zero on average in all plateaus (water regions)")
   parser.add_argument("-field", type=float, help="Set the intensity of external homog. electric field for field correction.", default=0.0)
   parser.add_argument("-zref",  type=float, help="Set the reference Z-box dimension.", default=14.8)
   parser.add_argument("-noqcorr", action='store_true', help="charge correction to zero on average")
   parser.add_argument("-plot", action='store_true', help="show plots of potential and field intensity")
   args = parser.parse_args()

   charge_file    = args.charge
   densities_file = args.densities
   z_reference    = args.zref

   if   args.fcorr     : fcorr_type ="no_work"
   elif args.fcorrcent : fcorr_type ="cent_plat_zero"
   elif args.fcorrplat : fcorr_type ="all_plat_zero"
   else :                fcorr_type =None  #no correction
   field_intensity = args.field
   
   if args.noqcorr: qcorr = False
   else           : qcorr = True
   crop = args.crop
   plot = args.plot

   # call with all arguments specified
   print """                                           ||-<--0V-->-||
  Membrane potential integrator for double-||-membrane || systems
<--------------------------------------V1--||          ||--V2--------->
  working with file: %s
         """ % (charge_file,)
   if plot: get_potential(charge_file, z_reference, crop, fcorr_type, field_intensity, qcorr, plot, densfile=densities_file)
   else:    get_potential(charge_file, z_reference, crop, fcorr_type, field_intensity, qcorr)

