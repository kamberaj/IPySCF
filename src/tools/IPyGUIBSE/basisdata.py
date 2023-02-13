__doc__="""
# Id: basisdata.py,v 1.0
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

from . utils   import *
import basis_set_exchange as bse

class basisDataDialog(object):
   """
      BSE basis Data Information class
   """
   def __init__ (self, master):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)
      #self.logoutText = ''
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="BSE basis Data Information", relief=RIDGE)
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
      Button(spacerFrame, text="Run", command=self.run).place(x=150, y=330)

      # Create the Filter Button
      Button(spacerFrame, text="Filter", command=self.filter).place(x=220, y=330)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)


   def run(self):
       data = bse.get_basis(self.updatedbasis(),fmt=self.updatedbasisfmt(), header=True)
       print(data)
       return None

   def filter(self):
       data = bse.api.filter_basis_sets( self.updatedbasis() )
       writeCSV( 'src/data.csv', data )
       tablePlot = plotApp()
       tablePlot.mainloop()

   def updatedbasis(self):
       return self.basis.get()

   def updatedbasisfmt(self):
       return self.basisfmt.get()

   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()
