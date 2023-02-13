__doc__="""
# Id: familydata.py,v 1.0
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

from . utils   import *
import basis_set_exchange as bse

class familyDataDialog(object):
   """
      BSE Family Data Information class
   """
   from tkinter import ttk
   def __init__ (self, master):
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      self.frame = Frame(master)
      #self.logoutText = ''
      
      # Create a Label to display a header
      self.headLbl = Label(self.frame, text="BSE Family Data Information", relief=RIDGE)
      self.headLbl.pack(side=TOP, fill=X)

      # Create a border using a dummy frame with a large border width
      spacerFrame = Frame(self.frame, borderwidth=10)

      self.familyCheckVal = BooleanVar(master)
      self.familyCheckVal.set(FALSE)
      Label(spacerFrame, text="Family Names:").place(x=0, y=20) 
      self.family = StringVar(self.master) 
      self.familyChosen = ttk.Combobox(spacerFrame, width=12, textvariable=self.family) 
      self.familyChosen['values'] = tuple(bse.get_families())
      self.familyChosen.current(0) 
      self.familyChosen.place(x=110, y=20)
      Checkbutton(spacerFrame, text="Selected", indicatoron=TRUE, variable=self.familyCheckVal, command=self.updatedfamily).place(x=250, y=20)

      # Create the Quit Button
      Button(spacerFrame, text="OK", command=self.Okay).place(x=10, y=330)
  
      # Create the Cancel Button
      Button(spacerFrame, text="Cancel", command=self.Cancel).place(x=80, y=330)

      # Create the Cancel Button
      Button(spacerFrame, text="Run", command=self.run).place(x=150, y=330)
  
      spacerFrame.pack(side=TOP, expand=YES, fill=BOTH)
      self.frame.pack(expand=YES, fill=BOTH)

   def run(self):
       data =bse.api.filter_basis_sets(family=self.updatedfamily())
       writeCSV( 'src/data.csv', data )
       tablePlot = plotApp()
       tablePlot.mainloop()
       return None

   def updatedfamily(self):
       return self.family.get()
 
   def Okay(self):
       self.master.destroy()

   def Cancel(self):
       self.master.destroy()
