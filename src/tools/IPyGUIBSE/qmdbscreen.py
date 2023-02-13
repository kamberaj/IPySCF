#!/usr/bin/env python
__doc__="""
# Id: qmdbscreen.py,v 1.0
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
from tkinter import Tk
from tkinter import messagebox

from tkinter import ttk

import os
import basis_set_exchange as bse

from pyscf import gto, scf, mcscf, cc

from pyscf.geomopt import geometric_solver, as_pyscf_method
from pyscf.cc import ccsd_t_lambda_slow as ccsd_t_lambda
from pyscf.grad import ccsd_t as ccsd_t_grad

from pyscf.hessian import thermo

from . gaussian_primitives import VOLUME, orbital_viz_dim

from . compound import *

from . utils   import *

conv_params = {
 'convergence_energy': 1e-4,  # Eh
 'convergence_grms': 3e-3,    # Eh/Bohr
 'convergence_gmax': 4.5e-3,  # Eh/Bohr
 'convergence_drms': 1.2e-2,  # Angstrom
 'convergence_dmax': 1.8e-2,  # Angstrom
}
Charge   = 0
Spin     = None
Symmetry = False
Verbose  = 0
Coords_unit = 'Angstrom'

############################################################################     
class ioForm(object):
   """
      ioForm for Molecular Processing 
   """
   def __init__ (self, master):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="Input/Output Configuration Options", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      # Create another frame to hold the center part of the form
      centerFrame = Frame(spacerFrame)
      leftColumn  = Frame(centerFrame, relief=GROOVE, borderwidth=2)
      rightColumn = Frame(centerFrame, relief=GROOVE, borderwidth=2)
      
      ######################################## Input #####################################################
       # Button for configuration file
      Label(leftColumn,  text="Input Configuration File:",).place(x=10, y=10)
      Button(leftColumn, text="Browse...", command=self.getConfigFileName).place(x=190, y=10)

      ########################## Output ##################################################################
      # Button for configuration file
      Checkbutton(rightColumn, text="Configuration File:", indicatoron=TRUE, command=self.setConfigFileName).place(x=5, y=20)
      self.setConfigFileNameText = StringVar(self.master)
      self.setConfigFileNameText.set("struct.mol")
      self.setConfigFileNameTextEntry = Entry(rightColumn, textvariable=self.setConfigFileNameText)
      self.setConfigFileNameTextEntry.place(x=150, y=20)
	   
      # Button for Energy file
      Checkbutton(rightColumn, text="Energy File:", indicatoron=TRUE, command=self.setEnergyFileName).place(x=5, y=60)
      self.setEnergyFileNameText = StringVar(self.master)
      self.setEnergyFileNameText.set("struct.ene")
      self.setEnergyFileNameTextEntry = Entry(rightColumn, textvariable=self.setEnergyFileNameText)
      self.setEnergyFileNameTextEntry.place(x=150, y=60)
	       
      # Make the frames visible
      leftColumn.pack(side=LEFT,  expand=YES, fill=BOTH)
      rightColumn.pack(side=LEFT, expand=YES, fill=BOTH)
      centerFrame.pack(side=TOP,  expand=YES, fill=BOTH)
      
      # Create two labels
      Label(spacerFrame, text="Inputs").place(x=2,    y=-13)
      Label(spacerFrame, text="Outputs").place(x=390, y=-13)

      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=530)

      # Create the Quit Button
      Button(spacerFrame, text="CANCEL", command=self.Cancel).place(x=80, y=530)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

############################################ Inputs ######################################################################
   def getConfigFileName(self):
       from tkinter import filedialog as fd
       from os import path
       self.inputConfigFile = []
       self.inputFileFmt = []
       fDir = fd.askdirectory()
       for root, dirs, files in os.walk(fDir):
           for file in files:
              tokens = file.split('.')
              if tokens[1] == 'xyz' or tokens[1] == 'mol' or tokens[1] == 'mol2':
                 token = file.split('.')
                 self.inputFileFmt.append( token[1].lower() )
                 self.inputConfigFile.append( fDir+'/'+file )
                 
       print( self.inputConfigFile )
       
       return None
       
###################################### outputs #########################################################################
   def setConfigFileName(self):
       from tkinter import filedialog as fd
       from os import path
       fDir = path.dirname(__file__)
       fName = fDir + '/scratch/geomopt/'+self.setConfigFileNameText.get()
       self.outputConfigFile = fName
       return None

   def setEnergyFileName(self):
       from tkinter import filedialog as fd
       from os import path
       fDir = path.dirname(__file__)
       fName = fDir + '/scratch/geomopt/'+self.setEnergyFileNameText.get()
       self.outputEnergyFile = fName
       return None
      
   def Okay(self):
       self.master.destroy()
       
   def Cancel(self):
       self.master.destroy()
 
############################################################################     
class qmdbscreenDialog(object):
   """
      QM Geometry Optimization Screening of Database class
   """
   def __init__ (self, master):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="QM Geometry Optimization of a Database Options", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      # Create another frame to hold the center part of the form
      centerFrame = Frame(spacerFrame)
      leftColumn  = Frame(centerFrame, relief=GROOVE, borderwidth=2)
      rightColumn = Frame(centerFrame, relief=GROOVE, borderwidth=2)

      # create win log 
      self.WinLog = Text(rightColumn, bd=0, bg="aliceblue", height="8", width="50", font="Arial",)
      self.WinLog.config(state=DISABLED)

      scrollbar = Scrollbar(rightColumn, command=self.WinLog.yview, cursor="heart")
      self.WinLog['yscrollcommand'] = scrollbar.set

      self.WinLog.config(state=DISABLED)
      self.WinLog.yview(END)

      scrollbar.place(x=420,  y=25, height=600)
      self.WinLog.place(x=0,  y=10, height=600, width=410)

      # Create widgets
      ypos = 10
      dy   = 60
 

      # Create some Radiobuttons
      ypos = 10
      dy   = 100
      self.qmmodelCheckVal = BooleanVar(master)
      self.qmmodelCheckVal.set(FALSE)
      Label(leftColumn, text="QM Model:").place(x=0, y=ypos) 
      self.qmmodel = StringVar(self.master) 
      qmmodelChosen = ttk.Combobox(leftColumn, width=12, textvariable=self.qmmodel) 
      qmmodelChosen['values'] = ('HF', 'CASSCF', 'PYSCF', 'CCSD')
      qmmodelChosen.current(0) 
      qmmodelChosen.place(x=150, y=ypos)
      Checkbutton(leftColumn, text="Selected", indicatoron=TRUE, variable=self.qmmodelCheckVal, \
                  command=self.updatedQmmodel).place(x=295, y=ypos)
      Button(leftColumn, text='(I)', command=self.infoQmmodel).place(x=395, y=ypos)
      ypos += dy
      
      self.optimizerCheckVal = BooleanVar(master)
      self.optimizerCheckVal.set(FALSE)
      Label(leftColumn, text="Chose Optimizer:").place(x=0, y=ypos) 
      self.optimizer = StringVar(self.master) 
      optimizerChosen = ttk.Combobox(leftColumn, width=12, textvariable=self.optimizer) 
      optimizerChosen['values'] = ('BERNYSOLVER', 'GEOMETRICSOLVER', 'GRADIENTSCLASS')
      optimizerChosen.current(0) 
      optimizerChosen.place(x=150, y=ypos)
      Checkbutton(leftColumn, text="Selected", indicatoron=TRUE, variable=self.optimizerCheckVal,\
                  command=self.updatedOptimizer).place(x=295, y=ypos)
      ypos += dy

      self.basisCheckVal = BooleanVar(master)
      self.basisCheckVal.set(FALSE)
      Label(leftColumn, text="Basis Names:").place(x=0, y=ypos) 
      self.basis = StringVar(self.master) 
      self.basisChosen = ttk.Combobox(leftColumn, width=12, textvariable=self.basis) 
      self.basisChosen['values'] = tuple(bse.api.get_all_basis_names())
      self.basisChosen.current(0) 
      self.basisChosen.place(x=150, y=ypos)
      Checkbutton(leftColumn, text="Selected", indicatoron=TRUE, variable=self.basisCheckVal,\
                  command=self.updatedbasis).place(x=295, y=ypos)

      ypos += dy
      
      self.basisfmtCheckVal = BooleanVar(master)
      self.basisfmtCheckVal.set(FALSE)
      Label(leftColumn, text="Basis Formats:").place(x=0, y=ypos) 
      self.basisfmt = StringVar(self.master) 
      self.basisfmtChosen = ttk.Combobox(leftColumn, width=12, textvariable=self.basisfmt) 
      self.basisfmtChosen['values'] = tuple(bse.api.get_formats())
      self.basisfmtChosen.current(0) 
      self.basisfmtChosen.place(x=150, y=ypos)
      Checkbutton(leftColumn, text="Selected", indicatoron=TRUE, variable=self.basisfmtCheckVal,\
                  command=self.updatedbasisfmt).place(x=295, y=ypos)

      # Button Callback
      ypos += dy
      Label(spacerFrame,  text="Input/Output Configuration Files:",).place(x=5, y=ypos)
      Button(spacerFrame, text="Choose ...", command=self.ioFiles).place(x=240, y=ypos)

      ypos += dy
      Checkbutton(leftColumn, text="Number of CPUs:", indicatoron=TRUE, command=self.updateCPUs).place(x=5, y=ypos)
      self.CPUs = StringVar(self.master)
      self.CPUs.set('1')
      self.setCPUsEntry = Entry(leftColumn, textvariable=self.CPUs)
      self.setCPUsEntry.place(x=150, y=ypos)
     
      # Make the frames visible
      leftColumn.pack(side=LEFT,  expand=YES, fill=BOTH)
      rightColumn.pack(side=LEFT, expand=YES, fill=BOTH)
      centerFrame.pack(side=TOP,  expand=YES, fill=BOTH)

      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=600)

      # Create the Upload DataBase of Molecules Button
      Button(spacerFrame, text="Upload", command=self.upload).place(x=45, y=600)

      # Create the High Cost Geometry Optimization Button
      Button(spacerFrame, text="High Cost", command=self.optimize_at_high_cost).place(x=105, y=600)

      # Create the Low Cost Geometry Optimization Button
      Button(spacerFrame, text="Low Cost", command=self.optimize_at_low_cost).place(x=185, y=600)

      # Create the Thermodynamics Property Button
      Button(spacerFrame, text="Properties", command=self.get_thermo_prop).place(x=260, y=600)

      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=350, y=600)

      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

############################################################################     
   def ioFiles(self):
       ioDialog = Tk()
       self.ioFrame = ioForm(ioDialog)       
       ioDialog.title("Interactive Molecular Structure Input/Output Files")
       ioDialog.geometry("800x600")
       ioDialog.resizable(width=TRUE, height=TRUE)       
       ioDialog.mainloop()
       return None
      
   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()

   def updatedbasis(self):
       return self.basis.get()

   def updatedbasisfmt(self):
       return self.basisfmt.get()

   def updateCPUs(self):
       return self.CPUs.get()

   def _get_thermo_prop(self):
      # Frequency analysis
      self.mf = self.mol.apply(scf.RHF).run()
      hessian = self.mf.Hessian().kernel()
      freq_info = thermo.harmonic_analysis(self.mf.mol, hessian)
      
      # Thermochemistry analysis at 298.15 K and 1 atmospheric pressure
      print('Do thermochemistry analysis ....')
      thermo_info = thermo.thermo(self.mf, freq_info['freq_au'], 298.15, 101325)

      print('Rotation constant')
      print(thermo_info['rot_const'])

      print('Zero-point energy')
      print(thermo_info['ZPE'])

      print('Internal energy at 0 K')
      print(thermo_info['E_0K'])

      print('Internal energy at 298.15 K')
      print(thermo_info['E_tot'])

      print('Enthalpy energy at 298.15 K')
      print(thermo_info['H_tot'])

      print('Gibbs free energy at 298.15 K')
      print(thermo_info['G_tot'])

      print('Heat capacity at 298.15 K')
      print(thermo_info['Cv_tot'])
      
      return_thermo_info = {}
      for key, value in thermo_info.items():
          return_thermo_info[key+' in '+ value[1]] = value[0]
          
      return return_thermo_info 
   
   def upload(self):
       self.db_mol  = []
       self.opt = False
       for inputConfigFile, inputFileFmt in zip(self.ioFrame.inputConfigFile, self.ioFrame.inputFileFmt):
          self.molBasis = " "
          self.mol = gto.Mole(charge=Charge, spin=Spin, symmetry=Symmetry, verbose=Verbose)
          self.mol.basis = {}
          self.moleculeName = ''

          print( inputConfigFile, inputFileFmt )
          
          if inputFileFmt == 'xyz':
             self.molecule = Compound(xyz=inputConfigFile)
          elif inputFileFmt == 'mol':  
             self.molecule = Compound(mol=inputConfigFile)
          elif inputFileFmt == 'mol2':  
             self.molecule = Compound(mol2=inputConfigFile)
             
          for i in range( self.molecule.natoms ):
              self.updatedatom(i, self.molecule.atomtypes[i], self.molecule.coordinates[i])

          self.mol.build()
          self.db_mol.append( self.mol )
   
       return None

   def get_thermo_prop(self):
       db_thermo_info = []
       if self.opt :
          for self.mol in self.db_opt_mol:
              db_thermo_info.append( self._get_thermo_prop() )
       else:
          for self.mol in self.db_mol:
              db_thermo_info.append( self._get_thermo_prop() )
          
       writeCSV( 'src/data.csv', db_thermo_info )
       tablePlot = plotApp()
       tablePlot.mainloop()
   
       return None

   def optimize_at_low_cost(self):
       self.db_opt_mol  = []
       self.opt = True
       for self.mol in self.db_mol:
           self._optimize_at_low_cost()
           self.db_opt_mol.append( self.mol )
         
       return None

   def optimize_at_high_cost(self):
       self.db_opt_mol  = []
       self.opt = True
       for self.mol in self.db_mol:
           self._optimize_at_high_cost()
           self.db_opt_mol.append( self.mol )
          
       return None
      
   def _optimize_at_high_cost(self):
       if self.qmmodel.get().upper() == 'HF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          sol = geometric_solver.GeometryOptimizer(mf).set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qmmodel.get().upper() == 'HF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          sol = mf.Gradients().optimizer(solver='geomeTRIC').set(params=conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qmmodel.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = geometric_solver.optimize(mc, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qmmodel.get().upper() == 'CASSCF' and self.optimizer.get().upper() == 'GRADIENTSCLASS':        
          mf = scf.RHF(self.mol)
          mc = mcscf.CASSCF(mf, 4, 4)
          sol = mc.Gradients().optimizer(solver='geomeTRIC').kernel(conv_params)
          sol.max_cycle = 100
          sol.run()
          
       elif self.qmmodel.get().upper() == 'PYSCF' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
          mf = scf.RHF(self.mol)
          f_method = as_pyscf_method(self.mol, self.f_pyscf)
          sol = geometric_solver.optimize(f_method, **conv_params)
          sol.max_cycle = 100
          sol.run()

       elif self.qmmodel.get().upper() == 'CCSD' and self.optimizer.get().upper() == 'GEOMETRICSOLVER':
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
 
       coords  = sol.mol.atom_coords(unit=Coords_unit)
       charges = sol.mol.atom_charges()

       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       sol.mol.set_common_orig_(nuc_charge_center)
       dip_ints = sol.mol.intor('cint1e_r_sph', comp=3)
      
       print( "Before Optimization:"+ '\n', scf.RHF(sol.mol).kernel() )
        
       self.mol = sol.mol.copy()

       print( "After Optimization: "+ '\n', scf.RHF(self.mol).kernel() )

       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'High Cost Total Energy (in hartree):' + '\n\n')
       self.WinLog.insert(END, str(e.e_tot) + '\n\n')
       self.WinLog.insert(END, 'High Cost Molecular Orbital Energies (in hartree):' + '\n\n')
       self.WinLog.insert(END, str(e.mo_energy) + '\n\n')
       self.WinLog.insert(END, 'High Cost Orbital Occupancies:' + '\n\n')
       self.WinLog.insert(END, str(e.mo_occ) + '\n\n')
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
      
   def _optimize_at_low_cost(self):
       self.mol_opt = scf.RHF(self.mol)    
       self.mf = self.mol.apply(scf.RHF).run()
       coords  = self.mol.atom_coords()
       charges = self.mol.atom_charges()
       nuc_charge_center = np.einsum('z,zx->x', charges, coords) / charges.sum()
       self.mol.set_common_orig_(nuc_charge_center)
       dip_ints = self.mol.intor('cint1e_r_sph', comp=3)
       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, 'Low Cost Total Energy (in hartree):' + '\n\n')
       self.WinLog.insert(END, str(self.mf.e_tot) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Molecular Orbital Energies (in hartree):' + '\n\n')
       self.WinLog.insert(END, str(self.mf.mo_energy) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Orbital Occupancies:' + '\n\n')
       self.WinLog.insert(END, str(self.mf.mo_occ) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Coordinates:' + '\n\n')
       self.WinLog.insert(END, str(coords) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Charges:' + '\n\n')
       self.WinLog.insert(END, str(charges) + '\n\n')
       self.WinLog.insert(END, 'Low Cost Geometry Optimization Dipole:' + '\n\n')
       self.WinLog.insert(END, str(dip_ints) + '\n\n')
       self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))

       return None
      
   def updatedatom(self, i, atomSymbol, coords):        
       x = []
       for c in coords: x.append(str(c))
       self.atomBasis     = self.basis.get()
       self.atomSymbol    = atomSymbol
       self.atomBasisName = self.basis.get()
       
       self.molBasis += self.atomBasis
       self.moleculeName+=self.atomSymbol
      
       if i == 0:
          self.mol.atom = self.atomSymbol + '  ' + x[0] + '  ' + x[1] + '  ' + x[2]
       else:
          self.mol.atom += ' '.join(['; ' + self.atomSymbol + '  ' + x[0] + '  ' + x[1] + '  ' + x[2]])

       if self.atomSymbol.upper() == 'H':
          augment_diffuse = 0
       else:
          augment_diffuse = 1
          
       self.mol.basis[self.atomSymbol] = gto.load(bse.api.get_basis(self.atomBasisName, \
                                                                    elements=self.atomSymbol, \
                                                                    augment_diffuse = augment_diffuse, \
                                                                    remove_free_primitives=False, \
                                                                    fmt='nwchem'), self.atomSymbol)

       print("Mol atoms:  ", self.mol.atom)

       if i == 0:
          self.WinLog.config(state=NORMAL)
          self.WinLog.insert(END, "Molecule Basis: " + self.atomBasis + '\n\n')
          self.WinLog.config(foreground="#000000", font=("Verdana", 10 ))
       
       return None

   def updatedOptimizer(self):
       return self.optimizer.get()


   def updatedQmmodel(self):
       return self.qmmodel.get()
      
        
   def infoQmmodel(self):
       messagebox.showinfo("Information", "Semiempirical QM model.")

   def infoOptimizer(self):
       messagebox.showinfo("Information", "Dynamics: Langevin Dynamics (LANG); BERENDSEN provides NVT and NPT dynamics for controlling temperature and pressure. NOSE dynamics are used only in NVT ensemble. NOSEHOOVER dynamics are used for both NVT and NPT ensembles. LEAP integrator provides increased accuracy.  It allows: Langevin dynamics (LANG) and Constant Temperature and Pressure (CPT) (based on Berendsen's method). VVER integrator also provides increase accuracy. It allows: Constant Temperature (NOSE) (Nose-Hoover method) Multiple Time Step (MTS). VV2 is a velocity-Verlet integrator based on the operator-splitting technique.  It works in conjunction with the TPCONTROL command.  It allows temperature and pressure control, and provides an efficient integration algorithm for polarizable force fields based on Drude oscillators.")
