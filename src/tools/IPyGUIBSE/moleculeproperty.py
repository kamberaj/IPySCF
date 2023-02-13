__doc__="""
# Id: moleculeproperty.py,v 1.0
#
#                This source code is part of
#
#   Interactive Quantum Mechanics Simulations
#
# Written by Hiqmet Kamberaj.
# Copyright (C) 2021 Hiqmet Kamberaj.
# Check out h.kamberaj@gmail.com for more information.
#
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software Foundation; 
# GPL-3.0
#
"""
import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import messagebox

import numpy as np
from .helpdoc   import *

from . utils   import *
import basis_set_exchange as bse
from . sto_ng import *
from . x_ng import *

from . gaussian_primitives import VOLUME, orbital_viz_dim

import os
from pyscf import gto, scf

class molecule_hfscf_energy(object):
   """
      Calculate Molecule Hartree-Fock SCF Energy class
   """
   def __init__ (self, master):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)
      self.symbol = "C"
      self.number = 7
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="Molecular Electronic Quantum Calaculations at HF Level", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)
        
      self.basisCheckVal = BooleanVar(master)
      self.basisCheckVal.set(FALSE)
      Label(spacerFrame, text="Basis Names:").place(x=0, y=20) 
      self.basis = StringVar(self.master) 
      self.basisChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.basis) 
      self.basisChosen['values'] = tuple(bse.api.get_all_basis_names())
      self.basisChosen.current(0) 
      self.basisChosen.place(x=110, y=20)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.basisCheckVal, command=self.updatedbasis).place(x=250, y=20)


      self.moleculeCheckVal = BooleanVar(master)
      self.moleculeCheckVal.set(FALSE)
      Checkbutton(spacerFrame, text="Atoms of Molecule:", indicatoron=TRUE, variable=self.moleculeCheckVal, command=self.updatedatoms).place(x=0, y=70)
      self.atomsTextVal = StringVar(self.master)
      self.atomsTextVal.set("H 0 0 0; H 0 0 1.2")
      Entry(spacerFrame, textvariable=self.atomsTextVal).place(x=150, y=70)      

      Label(spacerFrame, text="Total HF energy:").place(x=0, y=120) 
      self.energyTextVal = StringVar(self.master)
      self.energyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.energyTextVal).place(x=150, y=120)      

      Label(spacerFrame, text="Molecular Orbital energies:").place(x=0, y=170) 
      self.moenergyTextVal = StringVar(self.master)
      self.moenergyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.moenergyTextVal).place(x=170, y=170)      

      Label(spacerFrame, text="Molecular Orbital occupancy:").place(x=0, y=220) 
      self.mooccupTextVal = StringVar(self.master)
      self.mooccupTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.mooccupTextVal).place(x=170, y=220)      
 
      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=330)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=80, y=330)

      # Create the Run Button
      Button(spacerFrame, text="Plot", command=self.plot).place(x=150, y=330)

      # Create the Reference Button
      Button(spacerFrame, text="run", command=self.hfscf_energy).place(x=220, y=330)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

   def hfscf_energy(self):
       self.mol = gto.M(atom=self.updatedatoms(), basis=self.updatedbasis())
       e = self.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       myhf = self.mol.HF()
       myhf.kernel()
       # Orbital energies, Mulliken population etc.
       myhf.analyze()

   def plot(self):
      coords  = self.mol.atom_coords()
      charges = self.mol.atom_charges()
      self.moleculeName = self.atomsTextVal.get()
      nuc_coord_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
      
      n=30
      X,Y,Z = VOLUME(nuc_coord_center, -10.0,10.0, -10.0,10.0, -10.0,10.0, n)
      coords = []
      for i in range(n):
          for j in range(n):
              for k in range(n):
                  coords.append([X[i][j][k], Y[i][j][k], Z[i][j][k]])

          
      ao_value = self.mol.eval_gto("GTOval_sph", coords)
      mo_value = ao_value.dot( scf.RHF(self.mol).run().mo_coeff )
      Vol = np.sum(mo_value, axis = 1)
      
      orbital_viz_dim(X, Y, Z, Vol, -10, 10, self.moleculeName)
      
      mo_value = np.transpose( mo_value )
      for i in range(mo_value.shape[0]):
         orbital_viz_dim(X, Y, Z, mo_value[i], -10, 10, self.moleculeName + ": Orbital-"+str(i))
     
      return None

   def updatedbasis(self):
       return self.basis.get()

   def updatedatoms(self):
       return self.atomsTextVal.get()

   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()
