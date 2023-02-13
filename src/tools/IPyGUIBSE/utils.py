__doc__="""
# Id: utils.py,v 1.0
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
import numpy as np
from tkinter import *
from pandastable import Table, TableModel
import pandas as pd

import sys

import tkinter
from tkinter import *
from tkinter import Tk
from tkinter import messagebox
_s = re.compile('\s+')

class dataInput(object):
    def __init__(self, name):
          self.data = name 		
    def loadData(self): #import csv file path
        data = self.data
        # import CSV
        strflag = 0
        if type(data) == None:
            return
        if isinstance(data,str) and data != 'None':
            self.CSV_PATH = data
            try:
                data = pd.read_csv(data)
            except FileNotFoundError as e:
                raise e
            except:
                print("Unexpected error:",sys.exc_info()[0])
                raise
            finally:
                strflag = 1

            self.isPandas = True
            
        self.csv = data
        self.csv["Row_INDEX_"] = np.linspace(1,len(self.csv),len(self.csv))
        self.data = data
        self.data["Row_INDEX_"] = np.linspace(1,len(self.data),len(self.data))
        return self.csv 
		
class plotApp(Frame):
        """Basic test frame for the table"""
        def __init__(self, parent=None):
            self.parent = parent
            Frame.__init__(self)
            self.main = Toplevel()
            self.main.geometry('600x400+200+100')
            #self.main.title('Table Data Plots GUI')
            f = Frame(self.main)
            f.pack(fill=BOTH,expand=1)
            inp = dataInput("src/data.csv")
            df = inp.loadData()
            pt = Table(f, dataframe=df, showtoolbar=True, showstatusbar=True)
            pt.show()
       
            return None

def writeCSV(filename, names):
    """
    """

    df = pd.DataFrame(names)
    df.to_csv(filename)
    return None


