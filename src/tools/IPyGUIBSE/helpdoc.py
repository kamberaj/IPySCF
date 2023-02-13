__doc__="""
# Id: helpdoc.py,v 1.0
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

from sys import *

class helpDocument(object):
   """
         Help Documentation of BSE program
   """
   def __init__ (self, master, file):
       # showInfo needs to know the master
       self.master = master
       # Create a frame to hold the widgets
       self.frame = Frame(master)
       self.WinLog = Text(self.master, bd=0, bg="white", height="8", width="50", font="Arial",)
       self.WinLog.config(state=DISABLED)

       #Bind scrollbar to Chat window
       scrollbar = Scrollbar(self.master, command=self.WinLog.yview, cursor="heart")
       self.WinLog['yscrollcommand'] = scrollbar.set

       #Create Button to send message
       btn = Button(self.master, font=("Verdana",10,'bold'), text="Leave", 
                       width="12", height=5,
                       bd=0, bg="gray", activebackground="red",fg="blue",
                       command=self.quit)

       self.ReadDoc('src/tools/IPyGUIBSE/docs/'+file)   
	   
       #Place all components on the screen
       scrollbar.place(x=685,y=15, height=585)
       self.WinLog.place(x=15,y=5, height=585, width=685)
       btn.place(x=550, y=550, height=40)
   

   def ReadDoc(self, filename):
       """
          Read Mindoc file from documentation
       """
       f = open(filename, 'r')
       lines = f.readlines()
       f.close()
       logoutText = ''
       for line in lines:
           logoutText += line
       	   
       self.WinLog.config(state=NORMAL)
       self.WinLog.insert(END, "Logo: " + logoutText + '\n\n')
       self.WinLog.config(foreground="#442265", font=("Verdana", 10 ))
 
       self.WinLog.config(state=DISABLED)
       self.WinLog.yview(END)
                        
       return None
       
   def quit(self):
       self.master.destroy()

  
