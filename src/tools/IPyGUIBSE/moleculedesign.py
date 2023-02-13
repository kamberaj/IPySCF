__doc__="""
# Id: moleculedesign.py,v 1.0
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
from . gaussian_primitives import VOLUME, orbital_viz_dim

import os
import basis_set_exchange as bse

from pyscf import gto, scf, mcscf, cc
from . periodictable_en import periodictable_en

from pyscf.geomopt import geometric_solver, as_pyscf_method
##from pyscf.geomopt import berny_solver
from pyscf.cc import ccsd_t_lambda_slow as ccsd_t_lambda
from pyscf.grad import ccsd_t as ccsd_t_grad

conv_params = {
 'convergence_energy': 1e-4,  # Eh
 'convergence_grms': 3e-3,    # Eh/Bohr
 'convergence_gmax': 4.5e-3,  # Eh/Bohr
 'convergence_drms': 1.2e-2,  # Angstrom
 'convergence_dmax': 1.8e-2,  # Angstrom
}
Charge   = 0
Spin     = None
Symmetry = True
Verbose  = 2
Coords_unit = 'Angstrom'

##mol.charge = 0
##mol.spin = 0 # 2j == nelec_alpha - nelec_beta
##mol.symmetry = 1  # Allow the program to apply point group symmetry if possible
### .unit can be 'bohr', 'ang' to indicate the coordinates unit of the input mol.atom
### If a number is assigned to unit, this number will be used as the length of
### 1 Bohr (in Angstrom).  Eg you can double the bond length of a system by
### setting mol.unit = 0.529*.5.
##mol.unit = 'Ang'    # (New in version 1.1)
##
### Output
### ------
### To write output on disk, assign a filename to Mole.output
##mol.output = 'path/to/my_out.txt'
### if Mole.output is not given, the default output would be stdout
##
### Print level
### -----------
### Mole.verbose is used to control print level.  The print level can be 0 (quite,
### no output) to 9 (very noise).  The default level is 1, which only outputs the
### error message, it works almost the same as level 0.  Level 4 (info), or 5 (debug)
### are recommended value if some calculation detials are needed.
##mol.verbose = 4
### level 4 hides some details such as CPU timings, the orbital energies during
### the SCF iterations.
##
### max memory to use
### -----------------
##mol.max_memory = 1000 # in MB
### or use evnrionment  PYSCF_MAX_MEMORY  to control the memory usage
### (New in PySCF-1.3) eg
###    export PYSCF_MAX_MEMORY=10000 # 10 GB
###    python 00-input_mole.py
##
### Whether to use Cartesian GTOs (New since version 1.5)
### -----------------------------------------------------
### default: False
##mol.cart = True

####################################################################################################
class molecule_small_design_level_5(object):
   """
      Geometry Optimization of a Molecule Using Hartree-Fock SCF Energy class
   """
   def __init__ (self, master, mol, atIndex, molBasis, atomBasis, moleculeName):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)

      # Copy previous level data
      self.mol = mol.copy()
      self.atIndex  = atIndex
      self.molBasis = molBasis
      self.atomBasis = atomBasis
      self.moleculeName = moleculeName
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="Geometry Optimization of Molecule - MOLECULAR GAME LEVEL 5", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      # create win log 
      self.WinLog = Text(spacerFrame, bd=0, bg="aliceblue", height="8", width="50", font="Arial",)
      self.WinLog.config(state=DISABLED)

      scrollbar = Scrollbar(spacerFrame, command=self.WinLog.yview, cursor="heart")
      self.WinLog['yscrollcommand'] = scrollbar.set

      ###
      Label(spacerFrame, text="Atoms to Molecule:").place(x=0, y=10)
      self.atomTextVal = StringVar(self.master)
      self.atomTextVal.set(self.mol.atom)
      Entry(spacerFrame, textvariable=self.atomTextVal).place(x=150, y=10)      

      Label(spacerFrame,text="Molecule:").place(x=0, y=60)
      self.moleculeTextVal = StringVar(self.master)
      self.moleculeTextVal.set(self.moleculeName)
      Entry(spacerFrame, textvariable=self.moleculeTextVal).place(x=150, y=60)      
       
      self.WinLog.config(state=DISABLED)
      self.WinLog.yview(END)

      scrollbar.place(x=920,  y=25, height=400)
      self.WinLog.place(x=500,y=10, height=400, width=400)
      
      Label(spacerFrame, text="Total energy:").place(x=0, y=120) 
      self.energyTextVal = StringVar(self.master)
      self.energyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.energyTextVal).place(x=190, y=120)
      

      Label(spacerFrame, text="Molecular Orbital energies:").place(x=0, y=170) 
      self.moenergyTextVal = StringVar(self.master)
      self.moenergyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.moenergyTextVal).place(x=190, y=170)      

      Label(spacerFrame, text="Molecular Orbital occupancy:").place(x=0, y=220) 
      self.mooccupTextVal = StringVar(self.master)
      self.mooccupTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.mooccupTextVal).place(x=190, y=220)      

      self.optimizerCheckVal = BooleanVar(master)
      self.optimizerCheckVal.set(FALSE)
      Label(spacerFrame, text="Chose Optimizer:").place(x=0, y=270) 
      self.optimizer = StringVar(self.master) 
      optimizerChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.optimizer) 
      optimizerChosen['values'] = ('BERNYSOLVER', 'GEOMETRICSOLVER', 'GRADIENTSCLASS')
      optimizerChosen.current(0) 
      optimizerChosen.place(x=190, y=270)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.optimizerCheckVal,
                  command=self.updatedOptimizer).place(x=410, y=270)

      self.qmCheckVal = BooleanVar(master)
      self.qmCheckVal.set(FALSE)
      Label(spacerFrame, text="QM Method:").place(x=0, y=320) 
      self.qm = StringVar(self.master) 
      qmChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.qm) 
      qmChosen['values'] = ('HF', 'CASSCF', 'PYSCF', 'CCSD')
      qmChosen.current(0) 
      qmChosen.place(x=190, y=320)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.qmCheckVal,
                  command=self.updatedQM).place(x=410, y=320)
 
      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=380)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=50, y=380)
      
      # Create the Reference Button
      Button(spacerFrame, text="run (low cost)", command=self.hfscf_energy).place(x=110, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="run (high cost)", command=self.optimize_at_high_cost).place(x=220, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="Molecule Orbitals", command=self.plot_mol_orb).place(x=340, y=380)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)


   def plot_mol_orb(self):
      coords  = self.mol.atom_coords()
      charges = self.mol.atom_charges()

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
   
   def updatedOptimizer(self):
       return self.optimizer.get()

   def updatedQM(self):
       return self.qm.get()

   def optimize_at_high_cost(self):
       self.mol.build()
       if self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          sol = geometric_solver.GeometryOptimizer(mf).set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          sol = mf.Gradients().optimizer(solver='geomeTRIC').set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = geometric_solver.optimize(mc, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = mc.Gradients().optimizer(solver='geomeTRIC').kernel(conv_params)
          sol.max_cycle = 100
          sol.run()
          
       elif self.qm.get().upper() == 'PYSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          f_method = as_pyscf_method(self.mol, self.f_pyscf)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CCSD' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          cc_scan = cc.CCSD(mf).as_scanner() 
          f_method = as_pyscf_method(self.mol, self.f_ccsd)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()
          
       else:
          messagebox.showinfo("Warning", "The Methods Combination is not supported.")
          sys.exit(0)
          
       e = sol.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords  = sol.mol.atom_coords(unit=Coords_unit)
       charges = sol.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       sol.mol.set_common_orig_(nuc_charge_center)
       dip_ints = sol.mol.intor('cint1e_r_sph', comp=3)
      
       print( "Before Optimization:"+ '\n', scf.RHF(sol.mol).kernel() )
        
       self.mol = sol.mol.copy()

       print( "After Optimization: "+ '\n', scf.RHF(self.mol).kernel() )

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'High Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return None

   def f_ccsd(self, mol):

       # This calculates CCSD(t) energy
       mf = scf.RHF(mol).run()
       mycc = cc.CCSD(mf).run()
       et_correction = mycc.ccsd_t()
       e = mycc.e_tot + et_correction

       # This calculates CCSD(t) gradients
       g = ccsd_t_grad.Gradients(mycc).kernel()
       
       return e, g

   def f_pyscf(self, mol):
       k = 0.1
       grad_scan = scf.RHF(mol).nuc_grad_method().as_scanner()
       e, g = grad_scan(mol)
       xyz = mol.atom_coords()
       e_penalty = ( np.linalg.norm( xyz[0] - xyz[1] )**2 ) * k
       e += e_penalty
       g[0] += (xyz[0] - xyz[1]) * 2.0 * k
       g[1] -= (xyz[0] - xyz[1]) * 2.0 * k
       return e, g
      
   def hfscf_energy(self):
       self.moleculeTextVal.set( self.moleculeName )
       self.mol.build()
       self.mol_opt = scf.RHF(self.mol)
       
       e = self.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords = self.mol.atom_coords()
       charges = self.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       self.mol.set_common_orig_(nuc_charge_center)
       dip_ints = self.mol.intor('cint1e_r_sph', comp=3)

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))

       return None

   def readBasis(self):
       f = open('src/tools/IPyGUIBSE/docs/scratch/data.doc', 'r')
       lines = f.readlines()
       f.close()
       logoutText = ''
       element= ' '
       i = 0
       for line in lines:
           if i == 0:
              tokens = line.split('#')
              name =  tokens[1].strip()
              
           tokens = line.split(' ')
           if tokens[0] != ' ' and len(tokens[0]) == 1:
              element = tokens[0]
           logoutText += line
           i += 1
           
       return logoutText, element, name
      
   def NewAtoms(self):
       periodictable_en()       
       return None

   def updatedatom(self):
       self.atomBasis, self.atomSymbol, self.atomBasisName = self.readBasis()
       self.molBasis += self.atomBasis
       self.moleculeName+=self.atomSymbol
       self.atomTextVal.set(self.atomSymbol)
       if self.atIndex == 1:
          self.mol.atom = self.atomSymbol + ' 0 0 0'
       else:   
          self.mol.atom += ' '.join(['; ' + self.atomSymbol + ' '+('%f  ' %y)*3 for y in (0.5*np.random.rand(1)+1.0)])

       if self.atomSymbol.upper() == 'H':
          augment_diffuse = 0
       else:
          augment_diffuse = 1
          
       self.mol.basis[self.atomSymbol] = gto.load(bse.api.get_basis(self.atomBasisName, \
                                                                    elements=self.atomSymbol, \
                                                                    augment_diffuse = augment_diffuse, \
                                                                    remove_free_primitives=False, \
                                                                    fmt='nwchem'), self.atomSymbol)
       
       self.atIndex += 1

       print("Mol basis:  ", self.mol.basis)
       print("Mol atoms:  ", self.mol.atom)
      
       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, self.atomBasis + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return self.atomTextVal.get()
      
   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()

####################################################################################################
class molecule_small_design_level_4(object):
   """
      Geometry Optimization of a Molecule Using Different Costly Methods class
   """
   def __init__ (self, master, mol, atIndex, molBasis, atomBasis, moleculeName):
      # Info needs to know the master 
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)

      # setup a random numbers sequence
      np.random.seed(1999)

      # Copy previous level data
      self.mol = mol.copy()
      self.atIndex  = atIndex
      self.molBasis = molBasis
      self.atomBasis = atomBasis
      self.moleculeName = moleculeName
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="Geometry Optimization of Molecule - MOLECULAR GAME LEVEL 4", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      # create win log 
      self.WinLog = Text(spacerFrame, bd=0, bg="aliceblue", height="8", width="50", font="Arial",)
      self.WinLog.config(state=DISABLED)

      scrollbar = Scrollbar(spacerFrame, command=self.WinLog.yview, cursor="heart")
      self.WinLog['yscrollcommand'] = scrollbar.set

      ###
      Label(spacerFrame, text="Add Atoms to Molecule:").place(x=0, y=10)
      Button(spacerFrame, text="New Selection", command=self.NewAtoms).place(x=190, y=10)

      Label(spacerFrame,text="Atom:").place(x=50, y=40)
      self.atomTextVal = StringVar(self.master)
      self.atomTextVal.set(self.mol.atom)
      Entry(spacerFrame, textvariable=self.atomTextVal).place(x=190, y=40)      
      Button(spacerFrame, text="Update", command=self.updatedatom).place(x=410, y=40)

      Label(spacerFrame,text="Molecule:").place(x=50, y=70)
      self.moleculeTextVal = StringVar(self.master)
      self.moleculeTextVal.set(self.moleculeName)
      Entry(spacerFrame, textvariable=self.moleculeTextVal).place(x=190, y=70)      
       
      self.WinLog.config(state=DISABLED)
      self.WinLog.yview(END)

      scrollbar.place(x=920,  y=25, height=400)
      self.WinLog.place(x=500,y=10, height=400, width=400)
      
      Label(spacerFrame, text="Total energy:").place(x=0, y=120) 
      self.energyTextVal = StringVar(self.master)
      self.energyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.energyTextVal).place(x=190, y=120)
      

      Label(spacerFrame, text="Molecular Orbital energies:").place(x=0, y=170) 
      self.moenergyTextVal = StringVar(self.master)
      self.moenergyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.moenergyTextVal).place(x=190, y=170)      

      Label(spacerFrame, text="Molecular Orbital occupancy:").place(x=0, y=220) 
      self.mooccupTextVal = StringVar(self.master)
      self.mooccupTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.mooccupTextVal).place(x=190, y=220)      

      self.optimizerCheckVal = BooleanVar(master)
      self.optimizerCheckVal.set(FALSE)
      Label(spacerFrame, text="Chose Optimizer:").place(x=0, y=270) 
      self.optimizer = StringVar(self.master) 
      optimizerChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.optimizer) 
      optimizerChosen['values'] = ('BERNYSOLVER', 'GEOMETRICSOLVER', 'GRADIENTSCLASS')
      optimizerChosen.current(0) 
      optimizerChosen.place(x=190, y=270)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.optimizerCheckVal,
                  command=self.updatedOptimizer).place(x=410, y=270)

      self.qmCheckVal = BooleanVar(master)
      self.qmCheckVal.set(FALSE)
      Label(spacerFrame, text="QM Method:").place(x=0, y=320) 
      self.qm = StringVar(self.master) 
      qmChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.qm) 
      qmChosen['values'] = ('HF', 'CASSCF', 'PYSCF', 'CCSD')
      qmChosen.current(0) 
      qmChosen.place(x=190, y=320)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.qmCheckVal,
                  command=self.updatedQM).place(x=410, y=320)
 
      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=380)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=50, y=380)
      
      # Create the Reference Button
      Button(spacerFrame, text="run (low cost)", command=self.hfscf_energy).place(x=110, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="run (high cost)", command=self.optimize_at_high_cost).place(x=230, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="Next Level ...", command=self.next_level).place(x=360, y=380)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

   def updatedOptimizer(self):
       return self.optimizer.get()

   def updatedQM(self):
       return self.qm.get()

   def next_level(self):
       fe = Tk()   
       feFrame = molecule_small_design_level_5(fe, self.mol, self.atIndex, self.molBasis, self.atomBasis, self.moleculeName)
       fe.title("Interactive Quantum Mechanics Molecular Design GAMING")
       fe.geometry("950x440")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None

   def optimize_at_high_cost(self):
       self.mol.build()
       if self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          sol = geometric_solver.GeometryOptimizer(mf).set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          sol = mf.Gradients().optimizer(solver='geomeTRIC').set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = geometric_solver.optimize(mc, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = mc.Gradients().optimizer(solver='geomeTRIC').kernel(conv_params)
          sol.max_cycle = 100
          sol.run()
          
       elif self.qm.get().upper() == 'PYSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          f_method = as_pyscf_method(self.mol, self.f_pyscf)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CCSD' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          cc_scan = cc.CCSD(mf).as_scanner() 
          f_method = as_pyscf_method(self.mol, self.f_ccsd)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()
          
       else:
          messagebox.showinfo("Warning", "The Methods Combination is not supported.")
          sys.exit(0)
          
       e = sol.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords  = sol.mol.atom_coords(unit=Coords_unit)
       charges = sol.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       sol.mol.set_common_orig_(nuc_charge_center)
       dip_ints = sol.mol.intor('cint1e_r_sph', comp=3)
      
       print( "Before Optimization:"+ '\n', scf.RHF(sol.mol).kernel() )
        
       self.mol = sol.mol.copy()

       print( "After Optimization: "+ '\n', scf.RHF(self.mol).kernel() )

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'High Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return None

   def f_ccsd(self, mol):

       # This calculates CCSD(t) energy
       mf = scf.RHF(mol).run()
       mycc = cc.CCSD(mf).run()
       et_correction = mycc.ccsd_t()
       e = mycc.e_tot + et_correction

       # This calculates CCSD(t) gradients
       g = ccsd_t_grad.Gradients(mycc).kernel()
       
       return e, g

   def f_pyscf(self, mol):
       k = 0.1
       grad_scan = scf.RHF(mol).nuc_grad_method().as_scanner()
       e, g = grad_scan(mol)
       xyz = mol.atom_coords()
       e_penalty = ( np.linalg.norm( xyz[0] - xyz[1] )**2 ) * k
       e += e_penalty
       g[0] += (xyz[0] - xyz[1]) * 2.0 * k
       g[1] -= (xyz[0] - xyz[1]) * 2.0 * k
       return e, g
      
   def hfscf_energy(self):
       self.moleculeTextVal.set( self.moleculeName )
       self.mol.build()
       self.mol_opt = scf.RHF(self.mol)
       
       e = self.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords = self.mol.atom_coords()
       charges = self.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       self.mol.set_common_orig_(nuc_charge_center)
       dip_ints = self.mol.intor('cint1e_r_sph', comp=3)

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))

       return None

   def readBasis(self):
       f = open('src/tools/IPyGUIBSE/docs/scratch/data.doc', 'r')
       lines = f.readlines()
       f.close()
       logoutText = ''
       element= ' '
       i = 0
       for line in lines:
           if i == 0:
              tokens = line.split('#')
              name =  tokens[1].strip()
              
           tokens = line.split(' ')
           if tokens[0] != ' ' and len(tokens[0]) == 1:
              element = tokens[0]
           logoutText += line
           i += 1
           
       return logoutText, element, name
      
   def NewAtoms(self):
       periodictable_en()       
       return None

   def updatedatom(self):
       self.atomBasis, self.atomSymbol, self.atomBasisName = self.readBasis()
       self.molBasis += self.atomBasis
       self.moleculeName+=self.atomSymbol
       self.atomTextVal.set(self.atomSymbol)
       if self.atIndex == 1:
          self.mol.atom = self.atomSymbol + ' 0 0 0'
       else:   
          self.mol.atom += ' '.join(['; ' + self.atomSymbol + ' '+('%f  ' %y)*3 for y in (0.5*np.random.rand(1)+1.0)])

       if self.atomSymbol.upper() == 'H':
          augment_diffuse = 0
       else:
          augment_diffuse = 1
          
       self.mol.basis[self.atomSymbol] = gto.load(bse.api.get_basis(self.atomBasisName, \
                                                                    elements=self.atomSymbol, \
                                                                    augment_diffuse = augment_diffuse, \
                                                                    remove_free_primitives=False, \
                                                                    fmt='nwchem'), self.atomSymbol)
       
       self.atIndex += 1

       print("Mol basis:  ", self.mol.basis)
       print("Mol atoms:  ", self.mol.atom)
      
       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, self.atomBasis + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return self.atomTextVal.get()
      
   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()

####################################################################################################
class molecule_small_design_level_3(object):
   """
      Geometry Optimization of a Molecule Using Different Costly Methods class
   """
   def __init__ (self, master, mol, atIndex, molBasis, atomBasis, moleculeName):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)

      # setup a random numbers sequence
      np.random.seed(9999)

      # Copy previous level data
      self.mol = mol.copy()
      self.atIndex  = atIndex
      self.molBasis = molBasis
      self.atomBasis = atomBasis
      self.moleculeName = moleculeName
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="Geometry Optimization of Molecule - MOLECULAR GAME LEVEL 3", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      # create win log 
      self.WinLog = Text(spacerFrame, bd=0, bg="aliceblue", height="8", width="50", font="Arial",)
      self.WinLog.config(state=DISABLED)

      scrollbar = Scrollbar(spacerFrame, command=self.WinLog.yview, cursor="heart")
      self.WinLog['yscrollcommand'] = scrollbar.set

      ###
      Label(spacerFrame, text="Add Atoms to Molecule:").place(x=0, y=10)
      Button(spacerFrame, text="New Selection", command=self.NewAtoms).place(x=190, y=10)

      Label(spacerFrame,text="Atom:").place(x=50, y=40)
      self.atomTextVal = StringVar(self.master)
      self.atomTextVal.set(self.mol.atom)
      Entry(spacerFrame, textvariable=self.atomTextVal).place(x=190, y=40)      
      Button(spacerFrame, text="Update", command=self.updatedatom).place(x=410, y=40)

      Label(spacerFrame,text="Molecule:").place(x=50, y=70)
      self.moleculeTextVal = StringVar(self.master)
      self.moleculeTextVal.set(self.moleculeName)
      Entry(spacerFrame, textvariable=self.moleculeTextVal).place(x=190, y=70)      
       
      self.WinLog.config(state=DISABLED)
      self.WinLog.yview(END)

      scrollbar.place(x=920,  y=25, height=400)
      self.WinLog.place(x=500,y=10, height=400, width=400)
      
      Label(spacerFrame, text="Total energy:").place(x=0, y=120) 
      self.energyTextVal = StringVar(self.master)
      self.energyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.energyTextVal).place(x=190, y=120)
      

      Label(spacerFrame, text="Molecular Orbital energies:").place(x=0, y=170) 
      self.moenergyTextVal = StringVar(self.master)
      self.moenergyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.moenergyTextVal).place(x=190, y=170)      

      Label(spacerFrame, text="Molecular Orbital occupancy:").place(x=0, y=220) 
      self.mooccupTextVal = StringVar(self.master)
      self.mooccupTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.mooccupTextVal).place(x=190, y=220)      

      self.optimizerCheckVal = BooleanVar(master)
      self.optimizerCheckVal.set(FALSE)
      Label(spacerFrame, text="Chose Optimizer:").place(x=0, y=270) 
      self.optimizer = StringVar(self.master) 
      optimizerChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.optimizer) 
      optimizerChosen['values'] = ('BERNYSOLVER', 'GEOMETRICSOLVER', 'GRADIENTSCLASS')
      optimizerChosen.current(0) 
      optimizerChosen.place(x=190, y=270)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.optimizerCheckVal,
                  command=self.updatedOptimizer).place(x=410, y=270)

      self.qmCheckVal = BooleanVar(master)
      self.qmCheckVal.set(FALSE)
      Label(spacerFrame, text="QM Method:").place(x=0, y=320) 
      self.qm = StringVar(self.master) 
      qmChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.qm) 
      qmChosen['values'] = ('HF', 'CASSCF', 'PYSCF', 'CCSD')
      qmChosen.current(0) 
      qmChosen.place(x=190, y=320)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.qmCheckVal,
                  command=self.updatedQM).place(x=410, y=320)
 
      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=380)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=50, y=380)
      
      # Create the Reference Button
      Button(spacerFrame, text="run (low cost)", command=self.hfscf_energy).place(x=110, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="run (high cost)", command=self.optimize_at_high_cost).place(x=230, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="Next Level ...", command=self.next_level).place(x=360, y=380)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

   def updatedOptimizer(self):
       return self.optimizer.get()

   def updatedQM(self):
       return self.qm.get()

   def next_level(self):
       fe = Tk()   
       feFrame = molecule_small_design_level_4(fe, self.mol, self.atIndex, self.molBasis, self.atomBasis, self.moleculeName)
       fe.title("Interactive Quantum Mechanics Molecular Design GAMING")
       fe.geometry("950x440")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None

   def optimize_at_high_cost(self):
       self.mol.build()
       if self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          sol = geometric_solver.GeometryOptimizer(mf).set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          sol = mf.Gradients().optimizer(solver='geomeTRIC').set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = geometric_solver.optimize(mc, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = mc.Gradients().optimizer(solver='geomeTRIC').kernel(conv_params)
          sol.max_cycle = 100
          sol.run()
          
       elif self.qm.get().upper() == 'PYSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          f_method = as_pyscf_method(self.mol, self.f_pyscf)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CCSD' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          cc_scan = cc.CCSD(mf).as_scanner() 
          f_method = as_pyscf_method(self.mol, self.f_ccsd)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()
          
       else:
          messagebox.showinfo("Warning", "The Methods Combination is not supported.")
          sys.exit(0)
          
       e = sol.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords  = sol.mol.atom_coords(unit=Coords_unit)
       charges = sol.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       sol.mol.set_common_orig_(nuc_charge_center)
       dip_ints = sol.mol.intor('cint1e_r_sph', comp=3)
      
       print( "Before Optimization:"+ '\n', scf.RHF(sol.mol).kernel() )
        
       self.mol = sol.mol.copy()

       print( "After Optimization: "+ '\n', scf.RHF(self.mol).kernel() )

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'High Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return None

   def f_ccsd(self, mol):

       # This calculates CCSD(t) energy
       mf = scf.RHF(mol).run()
       mycc = cc.CCSD(mf).run()
       et_correction = mycc.ccsd_t()
       e = mycc.e_tot + et_correction

       # This calculates CCSD(t) gradients
       g = ccsd_t_grad.Gradients(mycc).kernel()
       
       return e, g

   def f_pyscf(self, mol):
       k = 0.1
       grad_scan = scf.RHF(mol).nuc_grad_method().as_scanner()
       e, g = grad_scan(mol)
       xyz = mol.atom_coords()
       e_penalty = ( np.linalg.norm( xyz[0] - xyz[1] )**2 ) * k
       e += e_penalty
       g[0] += (xyz[0] - xyz[1]) * 2.0 * k
       g[1] -= (xyz[0] - xyz[1]) * 2.0 * k
       return e, g
      
   def hfscf_energy(self):
       self.moleculeTextVal.set( self.moleculeName )
       self.mol.build()
       self.mol_opt = scf.RHF(self.mol)
       
       e = self.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords = self.mol.atom_coords()
       charges = self.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       self.mol.set_common_orig_(nuc_charge_center)
       dip_ints = self.mol.intor('cint1e_r_sph', comp=3)

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))

       return None

   def readBasis(self):
       f = open('src/tools/IPyGUIBSE/docs/scratch/data.doc', 'r')
       lines = f.readlines()
       f.close()
       logoutText = ''
       element= ' '
       i = 0
       for line in lines:
           if i == 0:
              tokens = line.split('#')
              name =  tokens[1].strip()
              
           tokens = line.split(' ')
           if tokens[0] != ' ' and len(tokens[0]) == 1:
              element = tokens[0]
           logoutText += line
           i += 1
           
       return logoutText, element, name
      
   def NewAtoms(self):
       periodictable_en()       
       return None

   def updatedatom(self):
       self.atomBasis, self.atomSymbol, self.atomBasisName = self.readBasis()
       self.molBasis += self.atomBasis
       self.moleculeName+=self.atomSymbol
       self.atomTextVal.set(self.atomSymbol)
       if self.atIndex == 1:
          self.mol.atom = self.atomSymbol + ' 0 0 0'
       else:   
          self.mol.atom += ' '.join(['; ' + self.atomSymbol + ' '+('%f  ' %y)*3 for y in (0.5*np.random.rand(1)+1.0)])

       if self.atomSymbol.upper() == 'H':
          augment_diffuse = 0
       else:
          augment_diffuse = 1
          
       self.mol.basis[self.atomSymbol] = gto.load(bse.api.get_basis(self.atomBasisName, \
                                                                    elements=self.atomSymbol, \
                                                                    augment_diffuse = augment_diffuse, \
                                                                    remove_free_primitives=False, \
                                                                    fmt='nwchem'), self.atomSymbol)
       
       self.atIndex += 1

       print("Mol basis:  ", self.mol.basis)
       print("Mol atoms:  ", self.mol.atom)
      
       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, self.atomBasis + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return self.atomTextVal.get()
      
   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()

####################################################################################################
class molecule_small_design_level_2(object):
   """
      Geometry Optimization of a Molecule Using Different Costly Methods class
   """
   def __init__ (self, master, mol, atIndex, molBasis, atomBasis, moleculeName):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)

      # setup a random numbers sequence
      np.random.seed(2999)

      # Copy previous level data
      self.mol = mol.copy()
      self.atIndex  = atIndex
      self.molBasis = molBasis
      self.atomBasis = atomBasis
      self.moleculeName = moleculeName
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="Geometry Optimization of Molecule - MOLECULAR GAME LEVEL 2", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      # create win log 
      self.WinLog = Text(spacerFrame, bd=0, bg="aliceblue", height="8", width="50", font="Arial",)
      self.WinLog.config(state=DISABLED)

      scrollbar = Scrollbar(spacerFrame, command=self.WinLog.yview, cursor="heart")
      self.WinLog['yscrollcommand'] = scrollbar.set

      ###
      Label(spacerFrame, text="Add Atoms to Molecule:").place(x=0, y=10)
      Button(spacerFrame, text="New Selection", command=self.NewAtoms).place(x=190, y=10)

      Label(spacerFrame,text="Atom:").place(x=50, y=40)
      self.atomTextVal = StringVar(self.master)
      self.atomTextVal.set(self.mol.atom)
      Entry(spacerFrame, textvariable=self.atomTextVal).place(x=190, y=40)      
      Button(spacerFrame, text="Update", command=self.updatedatom).place(x=410, y=40)

      Label(spacerFrame,text="Molecule:").place(x=50, y=70)
      self.moleculeTextVal = StringVar(self.master)
      self.moleculeTextVal.set(self.moleculeName)
      Entry(spacerFrame, textvariable=self.moleculeTextVal).place(x=190, y=70)      
       
      self.WinLog.config(state=DISABLED)
      self.WinLog.yview(END)

      scrollbar.place(x=920,  y=25, height=400)
      self.WinLog.place(x=500,y=10, height=400, width=400)
      
      Label(spacerFrame, text="Total energy:").place(x=0, y=120) 
      self.energyTextVal = StringVar(self.master)
      self.energyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.energyTextVal).place(x=190, y=120)
      

      Label(spacerFrame, text="Molecular Orbital energies:").place(x=0, y=170) 
      self.moenergyTextVal = StringVar(self.master)
      self.moenergyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.moenergyTextVal).place(x=190, y=170)      

      Label(spacerFrame, text="Molecular Orbital occupancy:").place(x=0, y=220) 
      self.mooccupTextVal = StringVar(self.master)
      self.mooccupTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.mooccupTextVal).place(x=190, y=220)      

      self.optimizerCheckVal = BooleanVar(master)
      self.optimizerCheckVal.set(FALSE)
      Label(spacerFrame, text="Chose Optimizer:").place(x=0, y=270) 
      self.optimizer = StringVar(self.master) 
      optimizerChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.optimizer) 
      optimizerChosen['values'] = ('BERNYSOLVER', 'GEOMETRICSOLVER', 'GRADIENTSCLASS')
      optimizerChosen.current(0) 
      optimizerChosen.place(x=190, y=270)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.optimizerCheckVal,
                  command=self.updatedOptimizer).place(x=410, y=270)

      self.qmCheckVal = BooleanVar(master)
      self.qmCheckVal.set(FALSE)
      Label(spacerFrame, text="QM Method:").place(x=0, y=320) 
      self.qm = StringVar(self.master) 
      qmChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.qm) 
      qmChosen['values'] = ('HF', 'CASSCF', 'PYSCF', 'CCSD')
      qmChosen.current(0) 
      qmChosen.place(x=190, y=320)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.qmCheckVal,
                  command=self.updatedQM).place(x=410, y=320)
 
      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=380)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=50, y=380)
      
      # Create the Reference Button
      Button(spacerFrame, text="run (low cost)", command=self.hfscf_energy).place(x=110, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="run (high cost)", command=self.optimize_at_high_cost).place(x=230, y=380)

      # Create the Reference Button
      Button(spacerFrame, text="Next Level ...", command=self.next_level).place(x=360, y=380)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

   def updatedOptimizer(self):
       return self.optimizer.get()

   def updatedQM(self):
       return self.qm.get()

   def next_level(self):
       fe = Tk()   
       feFrame = molecule_small_design_level_3(fe, self.mol, self.atIndex, self.molBasis, self.atomBasis, self.moleculeName)
       fe.title("Interactive Quantum Mechanics Molecular Design GAMING")
       fe.geometry("950x440")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None

   def optimize_at_high_cost(self):
       self.mol.build()
       if self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          sol = geometric_solver.GeometryOptimizer(mf).set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'HF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          sol = mf.Gradients().optimizer(solver='geomeTRIC').set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = geometric_solver.optimize(mc, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = mc.Gradients().optimizer(solver='geomeTRIC').kernel(conv_params)
          sol.max_cycle = 100
          sol.run()
          
       elif self.qm.get().upper() == 'PYSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          f_method = as_pyscf_method(self.mol, self.f_pyscf)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qm.get().upper() == 'CCSD' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          cc_scan = cc.CCSD(mf).as_scanner() 
          f_method = as_pyscf_method(self.mol, self.f_ccsd)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()
          
       else:
          messagebox.showinfo("Warning", "The Methods Combination is not supported.")
          sys.exit(0)
          
       e = sol.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords  = sol.mol.atom_coords(unit=Coords_unit)
       charges = sol.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       sol.mol.set_common_orig_(nuc_charge_center)
       dip_ints = sol.mol.intor('cint1e_r_sph', comp=3)
      
       print( "Before Optimization:"+ '\n', scf.RHF(sol.mol).kernel() )
        
       self.mol = sol.mol.copy()

       print( "After Optimization: "+ '\n', scf.RHF(self.mol).kernel() )

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'High Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'High Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return None

   def f_ccsd(self, mol):

       # This calculates CCSD(t) energy
       mf = scf.RHF(mol).run()
       mycc = cc.CCSD(mf).run()
       et_correction = mycc.ccsd_t()
       e = mycc.e_tot + et_correction

       # This calculates CCSD(t) gradients
       g = ccsd_t_grad.Gradients(mycc).kernel()
       
       return e, g

   def f_pyscf(self, mol):
       k = 0.1
       grad_scan = scf.RHF(mol).nuc_grad_method().as_scanner()
       e, g = grad_scan(mol)
       xyz = mol.atom_coords()
       e_penalty = ( np.linalg.norm( xyz[0] - xyz[1] )**2 ) * k
       e += e_penalty
       g[0] += (xyz[0] - xyz[1]) * 2.0 * k
       g[1] -= (xyz[0] - xyz[1]) * 2.0 * k
       return e, g
      
   def hfscf_energy(self):
       self.moleculeTextVal.set( self.moleculeName )
       self.mol.build()
       self.mol_opt = scf.RHF(self.mol)
       
       e = self.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords = self.mol.atom_coords()
       charges = self.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       self.mol.set_common_orig_(nuc_charge_center)
       dip_ints = self.mol.intor('cint1e_r_sph', comp=3)

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))

       return None

   def readBasis(self):
       f = open('src/tools/IPyGUIBSE/docs/scratch/data.doc', 'r')
       lines = f.readlines()
       f.close()
       logoutText = ''
       element= ' '
       i = 0
       for line in lines:
           if i == 0:
              tokens = line.split('#')
              name =  tokens[1].strip()
              
           tokens = line.split(' ')
           if tokens[0] != ' ' and len(tokens[0]) == 1:
              element = tokens[0]
           logoutText += line
           i += 1
           
       return logoutText, element, name
      
   def NewAtoms(self):
       periodictable_en()       
       return None

   def updatedatom(self):
       self.atomBasis, self.atomSymbol, self.atomBasisName = self.readBasis()
       self.molBasis += self.atomBasis
       self.moleculeName+=self.atomSymbol
       self.atomTextVal.set(self.atomSymbol)
       if self.atIndex == 1:
          self.mol.atom = self.atomSymbol + ' 0 0 0'
       else:   
          self.mol.atom += ' '.join(['; ' + self.atomSymbol + ' '+('%f  ' %y)*3 for y in (0.5*np.random.rand(1)+1.0)])

       if self.atomSymbol.upper() == 'H':
          augment_diffuse = 0
       else:
          augment_diffuse = 1
          
       self.mol.basis[self.atomSymbol] = gto.load(bse.api.get_basis(self.atomBasisName, \
                                                                    elements=self.atomSymbol, \
                                                                    augment_diffuse = augment_diffuse, \
                                                                    remove_free_primitives=False, \
                                                                    fmt='nwchem'), self.atomSymbol)
       
       self.atIndex += 1

       print("Mol basis:  ", self.mol.basis)
       print("Mol atoms:  ", self.mol.atom)
      
       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, self.atomBasis + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return self.atomTextVal.get()
      
   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()

####################################################################################################       
class molecule_small_design_level_1(object):
   """
      Design a Small Molecule and Calculate Hartree-Fock SCF Energy class
   """
   def __init__ (self, master):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)

      # setup a random numbers sequence
      np.random.seed(999)
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="Small Molecule Design and Quantum Calculations at HF Level - MOLECULAR GAME LEVEL 1", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      # create win log 
      self.WinLog = Text(spacerFrame, bd=0, bg="aliceblue", height="8", width="50", font="Arial",)
      self.WinLog.config(state=DISABLED)

      scrollbar = Scrollbar(spacerFrame, command=self.WinLog.yview, cursor="heart")
      self.WinLog['yscrollcommand'] = scrollbar.set

      ###
      Label(spacerFrame, text="Atoms of Molecule:").place(x=0, y=10)
      Button(spacerFrame, text="New Selection", command=self.NewAtoms).place(x=150, y=10)

      Label(spacerFrame,text="Atom:").place(x=50, y=40)
      self.atomTextVal = StringVar(self.master)
      self.atomTextVal.set("H 0 0 0")
      Entry(spacerFrame, textvariable=self.atomTextVal).place(x=190, y=40)      
      Button(spacerFrame, text="Update", command=self.updatedatom).place(x=410, y=40)

      Label(spacerFrame,text="Molecule:").place(x=50, y=70)
      self.moleculeTextVal = StringVar(self.master)
      self.moleculeTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.moleculeTextVal).place(x=190, y=70)      
       
      self.WinLog.config(state=DISABLED)
      self.WinLog.yview(END)

      scrollbar.place(x=920,  y=25, height=400)
      self.WinLog.place(x=500,y=10, height=400, width=400)
      
      Label(spacerFrame, text="Total HF energy:").place(x=0, y=120) 
      self.energyTextVal = StringVar(self.master)
      self.energyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.energyTextVal).place(x=190, y=120)
      

      Label(spacerFrame, text="Molecular Orbital energies:").place(x=0, y=170) 
      self.moenergyTextVal = StringVar(self.master)
      self.moenergyTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.moenergyTextVal).place(x=190, y=170)      

      Label(spacerFrame, text="Molecular Orbital occupancy:").place(x=0, y=220) 
      self.mooccupTextVal = StringVar(self.master)
      self.mooccupTextVal.set(" ")
      Entry(spacerFrame, textvariable=self.mooccupTextVal).place(x=190, y=220)      
 
      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=330)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=80, y=330)
      
      # Create the Reference Button
      Button(spacerFrame, text="run", command=self.hfscf_energy).place(x=180, y=330)

      # Create the Reference Button
      Button(spacerFrame, text="Next Level ...", command=self.next_level).place(x=280, y=330)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

   def next_level(self):
       fe = Tk()   
       feFrame = molecule_small_design_level_2(fe, self.mol, self.atIndex, self.molBasis, self.atomBasis, self.moleculeName)
       fe.title("Interactive Quantum Mechanics Molecular Design GAMING")
       fe.geometry("950x440")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None
      
   def hfscf_energy(self):
       self.moleculeTextVal.set( self.moleculeName )
       self.mol.build()
       self.mol_opt = scf.RHF(self.mol)
       
       e = self.mol.apply(scf.RHF).run()
       self.energyTextVal.set(e.e_tot)
       self.moenergyTextVal.set(e.mo_energy)
       self.mooccupTextVal.set(e.mo_occ)

       coords = self.mol.atom_coords()
       charges = self.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       self.mol.set_common_orig_(nuc_charge_center)
       dip_ints = self.mol.intor('cint1e_r_sph', comp=3)

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'Before Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'Before Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'Before Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))

       return None

   def readBasis(self):
       f = open('src/tools/IPyGUIBSE/docs/scratch/data.doc', 'r')
       lines = f.readlines()
       f.close()
       logoutText = ''
       element= ' '
       i = 0
       for line in lines:
           if i == 0:
              tokens = line.split('#')
              name =  tokens[1].strip()
              
           tokens = line.split(' ')
           if tokens[0] != ' ' and (len(tokens[0]) == 1 or len(tokens[0]) == 2):
              element = tokens[0]
           logoutText += line
           i += 1

       return logoutText, element, name
      
   def NewAtoms(self):
       self.atIndex  = 1
       self.molBasis = " "
       self.mol = gto.Mole(charge=Charge, spin=Spin, symmetry=Symmetry, verbose=Verbose)
       self.mol.basis = {}
       self.moleculeName = ''
       periodictable_en()       
       return None

   def updatedatom(self):
       self.atomBasis, self.atomSymbol, self.atomBasisName = self.readBasis()
       self.molBasis += self.atomBasis
       self.moleculeName+=self.atomSymbol
       self.atomTextVal.set(self.atomSymbol)
       if self.atIndex == 1:
          self.mol.atom = self.atomSymbol + ' 0 0 0'
       else:   
          self.mol.atom += ' '.join(['; ' + self.atomSymbol + ' '+('%f  ' %y)*3 for y in (0.5*np.random.rand(1)+1.0)])

       if self.atomSymbol.upper() == 'H':
          augment_diffuse = 0
       else:
          augment_diffuse = 1
          
       self.mol.basis[self.atomSymbol] = gto.load(bse.api.get_basis(self.atomBasisName, \
                                                                    elements=self.atomSymbol, \
                                                                    augment_diffuse = augment_diffuse, \
                                                                    remove_free_primitives=False, \
                                                                    fmt='nwchem'), self.atomSymbol)
       
       self.atIndex += 1

       print("Mol atoms:  ", self.mol.atom)
      
       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, self.atomBasis + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return self.atomTextVal.get()
      
   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()
