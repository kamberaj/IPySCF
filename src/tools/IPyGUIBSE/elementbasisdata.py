__doc__="""
# Id: elementbasisdata.py,v 1.0
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

import numpy as np
from .helpdoc   import *

from .utils   import *
import basis_set_exchange as bse
from .sto_ng import *
from .x_ng import *

import os

class elementbasisDataDialog(object):
   """
      BSE basis Data Information for periodic table of elements class
   """
   from tkinter import ttk
   def __init__ (self, master, symbol, number):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)
      self.symbol = symbol
      self.number = number
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="BSE basis Data Information for periodic table of elements", relief=RIDGE)
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


      self.basisfmtCheckVal = BooleanVar(master)
      self.basisfmtCheckVal.set(FALSE)
      Label(spacerFrame, text="Basis Formats:").place(x=400, y=20) 
      self.basisfmt = StringVar(self.master) 
      self.basisfmtChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.basisfmt) 
      self.basisfmtChosen['values'] = tuple(bse.api.get_formats())
      self.basisfmtChosen.current(0) 
      self.basisfmtChosen.place(x=500, y=20)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.basisfmtCheckVal, command=self.updatedbasisfmt).place(x=650, y=20)

      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=330)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=80, y=330)

      # Create the Run Button
      Button(spacerFrame, text="Plot", command=self.plot).place(x=150, y=330)

      # Create the Filter Button
      Button(spacerFrame, text="Filter", command=self.filter).place(x=220, y=330)

      # Create the Reference Button
      Button(spacerFrame, text="Reference", command=self.reference).place(x=290, y=330)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

   def reference(self):
       ref = bse.get_references(self.updatedbasis(), fmt='bib', elements=[1,9])
       if not os.path.exists("src/tools/IPyGUIBSE/docs/scratch"): os.system("mkdir src/tools/IPyGUIBSE/docs/scratch")
       fout = open('src/tools/IPyGUIBSE/docs/scratch/ref.doc', 'w')
       fout.write(ref)
       fout.close()
       
       helpMe = Tk()
       self.helpFrame = helpDocument(helpMe, 'scratch/ref.doc')
       helpMe.title("BSE Reference:  " + self.updatedbasis())
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   def filter(self):
       bse_basis = bse.get_basis(self.updatedbasis(), elements=self.symbol)
       print(bse_basis)
       if not os.path.exists("src/tools/IPyGUIBSE/docs/scratch"): os.system("mkdir src/tools/IPyGUIBSE/docs/scratch")
       bse.write_formatted_basis_file(basis_dict=bse_basis, outfile_path='src/tools/IPyGUIBSE/docs/scratch/data.doc', basis_fmt=self.updatedbasisfmt(),header=bse_basis['name'])
       helpMe = Tk()
       self.helpFrame = helpDocument(helpMe, 'scratch/data.doc')
       helpMe.title("BSE Information for the Chemical Element:  " + self.symbol)
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None
      
   def plot(self):
       base_name = self.updatedbasis().upper()
       print('Base-name:   ', base_name)
       bse_basis = bse.get_basis(base_name, elements=self.symbol)
       if base_name == "STO-2G" or base_name == "STO-3G"  or base_name == "STO-3G*"  or base_name == "STO-4G"  or base_name == "STO-5G"  or base_name == "STO-6G" :
          
          sgto = STO_NG(bse_basis['elements'][str(self.number)]['electron_shells'][0]['exponents'], \
                  bse_basis['elements'][str(self.number)]['electron_shells'][0]['coefficients'],\
                  orb_type="s")
          sgto._STO_NG()
       
          pgto = STO_NG(bse_basis['elements'][str(self.number)]['electron_shells'][1]['exponents'], \
                  bse_basis['elements'][str(self.number)]['electron_shells'][1]['coefficients'],\
                  orb_type="sp-inner")
          pgto._STO_NG()

       elif base_name == "2ZAPA-NR":
          i=0
          while len(bse_basis['elements'][str(self.number)]['electron_shells']) > i:
              exponen = bse_basis['elements'][str(self.number)]['electron_shells'][i]['exponents']
              coeffic = bse_basis['elements'][str(self.number)]['electron_shells'][i]['coefficients']
              print('exp+coef: ', exponen, coeffic, \
                    bse_basis['elements'][str(self.number)]['electron_shells'][i]['function_type'],\
                    bse_basis['elements'][str(self.number)]['electron_shells'][i]['angular_momentum']\
                   )
              if bse_basis['elements'][str(self.number)]['electron_shells'][i]['angular_momentum'][0]   == 0:
                 print("Plot 1s AO")
                 x_ng = X_NG(exponen, coeffic, orb_type="1s")
                 x_ng._X_NG()
              elif bse_basis['elements'][str(self.number)]['electron_shells'][i]['angular_momentum'][0] == 1:
                 print("Plot 2p AO")
                 x_ng = X_NG(exponen, coeffic, orb_type="2p")
                 x_ng._X_NG()

              elif bse_basis['elements'][str(self.number)]['electron_shells'][i]['angular_momentum'][0] == 2:
                 print("Plot 1d AO")
                 x_ng = X_NG(exponen, coeffic, orb_type="1d")
                 x_ng._X_NG()

              else:
                 print("Error: Atomic orbital type not implemented:   ", \
                       bse_basis['elements'][str(self.number)]['electron_shells'][i]['angular_momentum'][0])
                 sys.exit(1)
              i += 1
          
       else:
          exponen = bse_basis['elements'][str(self.number)]['electron_shells'][0]['exponents']
          coeffic = bse_basis['elements'][str(self.number)]['electron_shells'][0]['coefficients']
                    
          x_ng = X_NG(exponen, coeffic, orb_type="1s")
          x_ng._X_NG()

          if (len(bse_basis['elements'][str(self.number)]['electron_shells']) < 2):
             return None
          else:   
             exponen = bse_basis['elements'][str(self.number)]['electron_shells'][1]['exponents']
             coeffic = bse_basis['elements'][str(self.number)]['electron_shells'][1]['coefficients']

             x_ng = X_NG(exponen, coeffic, orb_type="sp-inner")
          
             if (len(bse_basis['elements'][str(self.number)]['electron_shells']) < 3):
                x_ng._X_NG()
                return None
             else:   
                exponen = bse_basis['elements'][str(self.number)]['electron_shells'][2]['exponents']
                coeffic = bse_basis['elements'][str(self.number)]['electron_shells'][2]['coefficients']

                if len(coeffic) > len(coeffic[0]):
                   x_ng.set_sp_outer1(exponen, coeffic)
                else:
                   x_ng.set_d(exponen, coeffic)
                   x_ng._X_NG()
                   return None

                if (len(bse_basis['elements'][str(self.number)]['electron_shells']) < 4):
                    x_ng._X_NG()
                    return None
                else:   
                    exponen = bse_basis['elements'][str(self.number)]['electron_shells'][3]['exponents']
                    coeffic = bse_basis['elements'][str(self.number)]['electron_shells'][3]['coefficients']

                    if len(coeffic) > len(coeffic[0]):
                       x_ng.set_sp_outer2(exponen, coeffic)
                    else:
                       x_ng.set_d(exponen, coeffic)
                       x_ng._X_NG()
                       return None

                    if (len(bse_basis['elements'][str(self.number)]['electron_shells']) < 5):
                        x_ng._X_NG()
                        return None
                    else:
                        exponen = bse_basis['elements'][str(self.number)]['electron_shells'][4]['exponents']
                        coeffic = bse_basis['elements'][str(self.number)]['electron_shells'][4]['coefficients']
          
                        if len(coeffic) > len(coeffic[0]):
                           x_ng.set_sp_outer3(exponen, coeffic)
                        else:
                           x_ng.set_d(exponen, coeffic)
                           x_ng._X_NG()
                           return None
                     
                        if (len(bse_basis['elements'][str(self.number)]['electron_shells']) < 6):
                            x_ng._X_NG()
                            return None
                        else:
                            exponen = bse_basis['elements'][str(self.number)]['electron_shells'][5]['exponents']
                            coeffic = bse_basis['elements'][str(self.number)]['electron_shells'][5]['coefficients']

                            x_ng.set_d(exponen, coeffic)    
                            x_ng._X_NG()
                            return None

   def updatedbasis(self):
       return self.basis.get()

   def updatedbasisfmt(self):
       return self.basisfmt.get()

   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()
