#!/usr/bin/env python
__doc__="""
# Id: program.py,v 1.0
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

import tkinter
from tkinter import *
from tkinter import Tk

from tkinter import messagebox


from . tools.IPyGUIBSE.bseGUI import call_pyscf

def program():
   """
       Main Window Frame
   """
   call_pyscf()
