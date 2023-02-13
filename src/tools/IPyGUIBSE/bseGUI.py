__doc__="""
# Id: bseGUI.py,v 1.0
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

import webbrowser

from . helpdoc   import helpDocument

import basis_set_exchange as bse

import numpy as np

from . familydata import familyDataDialog
from . basisdata import basisDataDialog
from . periodictable import periodictable
from . periodictable_en import periodictable_en
from . moleculeproperty import molecule_hfscf_energy
from . moleculedesign import molecule_small_design_level_1

from . qm import qmDialog
from . qmdbscreen import qmdbscreenDialog


class WinFrame:
   """
       Main window frame
   """
   def __init__ (self, master):
       
      # showInfo needs to know the master 5
      self.master = master
      # Create a frame to hold the widgets
      frame = Frame(master)
      menuBar = Menu(master)
      
      # BSE content menu 
      bseMenu = Menu(menuBar, tearoff = 0)
      bseMenu.add_command(label='BSE Package Documentation', command=self.helpBSEDoc)
      bseMenu.add_separator()
      menuBar.add_cascade(label='BSE Package', menu=bseMenu)

      # BSE Package content menu 
      bse_contentMenu = Menu(menuBar, tearoff = 0)      
      bse_contentMenu.add_command(label='BSE API Documentation', command=self.helpBSEAPIDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE BUNDLE Documentation', command=self.helpBSEBUNDLEDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE COMPOSING Documentation', command=self.helpBSECOMPOSEDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE CONVERT Documentation', command=self.helpBSECONVERTDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE FILEIO Documentation', command=self.helpBSEFileIODoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE One-center integrals for Gaussian-type and Slater-type orbitals', command=self.helpBSEIntsDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE Functions for LUT', command=self.helpBSELUTDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE Common Basis Set Manipulations', command=self.helpBSEMANIPDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE Miscellaneous Helper Functions', command=self.helpBSEMISCDoc)
      bse_contentMenu.add_separator()

      bse_contentMenu.add_command(label='BSE Helpers for Printing Pieces of Basis Sets', command=self.helpBSEPrintingDoc)
      bse_contentMenu.add_separator()
    
      menuBar.add_cascade(label='BSE Package Content', menu=bse_contentMenu)

      # BSE Package File content menu 
      bse_filesMenu = Menu(menuBar, tearoff = 0)
      
      bse_bse_dataMenu = Menu(bse_filesMenu, tearoff = 0)
      bse_bse_dataMenu.add_command(label='BSE Get Family Data', command=self.helpBSEFamilyData)
      bse_bse_dataMenu.add_separator()
      bse_bse_dataMenu.add_command(label='BSE Get Basis Data', command=self.helpBSEBasisData)
      bse_bse_dataMenu.add_separator()
      bse_filesMenu.add_cascade(label='BSE Data Content', menu=bse_bse_dataMenu)
      menuBar.add_cascade(label='BSE Files/Data Content', menu=bse_filesMenu) 

      # Atomic Basis sets menu 
      bse_atomMenu = Menu(menuBar, tearoff = 0)
      
      bse_enge_atomMenu = Menu(bse_atomMenu, tearoff = 0)
      bse_enge_atomMenu.add_command(label='German Language', command=self.getBSEGEAtomData)
      bse_enge_atomMenu.add_separator()
      bse_enge_atomMenu.add_command(label='English Language', command=self.getBSEENAtomData)
      bse_enge_atomMenu.add_separator()
      bse_atomMenu.add_cascade(label='BSE Periodic Table',  menu=bse_enge_atomMenu)
      
      menuBar.add_cascade(label='Atomic Orbitals', menu=bse_atomMenu) 

      # Molecular Properties menu 
      bse_molMenu = Menu(menuBar, tearoff = 0)
      
      bse_hfscf_molMenu = Menu(bse_molMenu, tearoff = 0)
      bse_hfscf_molMenu.add_command(label='HF-SCF Energy', command=self.getBSEHFSCFEnergy)
      bse_hfscf_molMenu.add_separator()
      bse_molMenu.add_cascade(label='QM Molecular Properties',  menu=bse_hfscf_molMenu)
      bse_hfscf_molMenu.add_separator()

      bse_molMenu.add_separator()

      bse_molMenu.add_command(label='QM Molecule Design GAME', command=self.getSmallMoleculeDesign)
      bse_molMenu.add_separator()

      bse_geomopt_molMenu = Menu(bse_molMenu, tearoff = 0)
      bse_geomopt_molMenu.add_command(label='QM Molecule Design View (2D & 3D)', command=self.MoleculeDesignView)
      bse_geomopt_molMenu.add_separator()
      bse_geomopt_molMenu.add_command(label='QM Molecule Geometry Optimisation', command=self.MoleculeGeomOpt)
      bse_geomopt_molMenu.add_separator()
      bse_molMenu.add_cascade(label='QM Molecule Design',  menu=bse_geomopt_molMenu)
      bse_molMenu.add_separator()

      bse_molMenu.add_command(label='QM Molecular Database Optimization Screening', command=self.qmMolDatabaseScreening)
      bse_molMenu.add_separator()

      bse_molMenu.add_separator()
      bse_molMenu.add_command(label='PySCF Documentation', command=self.helpPySCFDoc)
      
      menuBar.add_cascade(label='Molecular Properties', menu=bse_molMenu) 

      self.master.config(menu=menuBar)

   @classmethod
   def helpBSEDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'help.doc')
       helpMe.title("BSE Package Documentation")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSEAPIDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'api/help.doc')
       helpMe.title("BSE API Documentation")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSEBUNDLEDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'bundle/help.doc')
       helpMe.title("BSE BUNDLE Documentation")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSECOMPOSEDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'compose/help.doc')
       helpMe.title("BSE COMPOSING Documentation")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSECONVERTDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'convert/help.doc')
       helpMe.title("BSE CONVERT Documentation")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSEFileIODoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'fileio/help.doc')
       helpMe.title("BSE FileIO Documentation")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSEIntsDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'ints/help.doc')
       helpMe.title("BSE One-center integrals for Gaussian-type and Slater-type orbitals")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSELUTDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'lut/help.doc')
       helpMe.title("BSE Functions for LUT")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSEMANIPDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'manip/help.doc')
       helpMe.title("BSE Common Basis Set Manipulations")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSEMISCDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'misc/help.doc')
       helpMe.title("BSE Miscellaneous Felper Functions")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

   @classmethod
   def helpBSEPrintingDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'printing/help.doc')
       helpMe.title("BSE Helpers for Printing Pieces of Basis Sets")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None


   @classmethod
   def helpBSEFamilyData(cls):
       fe = Tk()   
       feFrame = familyDataDialog(fe)
       fe.title("Interactive BSE Family Data Dialog")
       fe.geometry("600x450")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None

   @classmethod   
   def helpBSEBasisData(cls):
       fe = Tk()   
       feFrame = basisDataDialog(fe)
       fe.title("Interactive BSE Basis Data Dialog")
       fe.geometry("900x450")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None

   @classmethod   
   def getBSEGEAtomData(cls):
       periodictable()
       return None

   @classmethod   
   def getBSEENAtomData(cls):
       periodictable_en()
       return None

   @classmethod
   def getBSEHFSCFEnergy(cls):
       fe = Tk()   
       feFrame = molecule_hfscf_energy(fe)
       fe.title("Interactive Electronic Molecular Properties Dialog")
       fe.geometry("900x450")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None

   @classmethod
   def getSmallMoleculeDesign(cls):
       fe = Tk()   
       feFrame = molecule_small_design_level_1(fe)
       fe.title("Interactive Quantum Mechanics Molecular Design GAMING")
       fe.geometry("950x440")
       fe.resizable(width=TRUE, height=TRUE)
       fe.mainloop()
       return None

   @classmethod
   def MoleculeDesignView(cls):
       webbrowser.open("http://localhost//icons/bsann/molview-master/index.php")
       return None

   @classmethod
   def MoleculeDesignView(cls):
       webbrowser.open("http://localhost//icons/bsann/molview-master/index.php")
       return None

   @classmethod 
   def MoleculeGeomOpt(self):
       qm = Tk()
       self.qmFrame = qmDialog(qm)
       qm.title("Interactive QM Geometry Optimization Dialog")
       qm.geometry("900x670")
       qm.resizable(width=TRUE, height=TRUE)
       qm.mainloop()
       return None


   @classmethod 
   def qmMolDatabaseScreening(self):
       qm = Tk()
       self.qmFrame = qmdbscreenDialog(qm)
       qm.title("Interactive QM Molecular Database Optimization Screening Dialog")
       qm.geometry("900x670")
       qm.resizable(width=TRUE, height=TRUE)
       qm.mainloop()
       return None


   @classmethod
   def helpPySCFDoc(cls):
       helpMe = Tk()
       helpFrame = helpDocument(helpMe, 'pyscf/help.doc')
       helpMe.title("Pyscf Documentation")
       helpMe.geometry("700x600")
       helpMe.resizable(width=TRUE, height=TRUE)
       helpMe.mainloop()
       return None

import tkinter 
from tkinter import *
from PIL import Image, ImageTk  # this is pillow
import os, sys
import tkinter as tk
     
def call_pyscf():
      """
       Main Window Frame
      """
      base_pyscf = Tk()

      myframe = WinFrame(base_pyscf)

      base_pyscf.title("Interactive Molecular Quantum Mechanics Simulations")
      base_pyscf.geometry("900x900")
      base_pyscf.resizable(width=TRUE, height=TRUE)

      if not os.path.exists("src/images"):
         print("Folder does not exist:   " + "images")
         sys.exit(1)
         
      photo = PhotoImage(file = "src/images/qmlogo.png")
      photoimage = photo.subsample(1,1)
      Button(base_pyscf, image = photoimage, ).pack(expand=NO,fill=BOTH,ipadx=120,ipady=120,anchor=CENTER)
      
      base_pyscf.mainloop()
      
