__doc__="""
# Id: compound.py,v 1.0
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
import collections

from . data import NUCLEAR_CHARGE

class Compound(object):
    """ The ``Compound`` class is used to store data from  

        :param xyz: Option to initialize the ``Compound`` with data from an XYZ file.
        :type xyz: string
    """

    def __init__(self, xyz = None, mol = None, mol2 = None):

        empty_array = np.asarray([], dtype = float)

        self.molid = float("nan")
        self.name = None

        # Information about the compound
        self.natoms = float("nan")
        self.natypes = {}
        self.atomtypes = []
        self.atomtype_indices = collections.defaultdict(list)
        self.nuclear_charges = empty_array
        self.coordinates = empty_array
        self.active_atoms = empty_array
        self.unit_cell = None

        # Container for misc properties
        self.energy = float("nan")
        self.properties = empty_array
        self.properties2 = empty_array

        # Representations:
        self.representation = empty_array

        if xyz is not None:
            self.read_xyz(xyz)
        elif mol is not None:
            self.ReadMol(mol)
        elif mol2 is not None:
            self.ReadMol2(mol2)

               
    def read_xyz(self, filename):
        """(Re-)initializes the Compound-object with data from an xyz-file.

    :param filename: Input xyz-filename.
    :type filename: string
    """
        try:
            f = open(filename, "r")
            lines = f.readlines()
            f.close()

            self.natoms = int(lines[0])
            self.atomtypes = []
            self.nuclear_charges = np.empty(self.natoms, dtype=int)
            self.coordinates = np.empty((self.natoms, 3), dtype=float)

            self.name = filename

            for i, line in enumerate(lines[2:self.natoms+2]):
                tokens = line.split()

                if len(tokens) < 4:
                    break

                self.atomtypes.append(tokens[0])
                self.atomtype_indices[tokens[0]].append(i)
                self.nuclear_charges[i] = NUCLEAR_CHARGE[tokens[0]]
        
                self.coordinates[i] = np.asarray(tokens[1:4], dtype=float)
       
            self.natypes = dict([(key, len(value)) for key,value in self.atomtype_indices.items()])

        except IOError as e:
             errno, strerror = e.args
             print( "(Readxyz) I/O error(%s):%s"%(errno,strerror) )
        except ValueError:
             print( "(Readxyz) Could not convert data to an integer." )
        except:
             print( "(Readxyz) Unexpected error:", sys.exc_info()[0] )
             raise

        
    def ReadMol(self, filename):
        """
          Read molecule structure from MOL format
        """
        try:
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            tokens = lines[3].split()
            self.natoms = int(tokens[0])
            self.atomtypes = []
            self.nuclear_charges = np.empty(self.natoms, dtype=int)
            self.coordinates = np.empty((self.natoms, 3), dtype=float)

            self.name = filename 
            for i, line in enumerate(lines[4:self.natoms+4]):
                tokens = line.split()
                if len(tokens) < 16: 
                   break
                else:
                   self.coordinates[i] = np.asarray(tokens[0:3], dtype=float)
                   self.atomtypes.append(tokens[3])
                   self.atomtype_indices[tokens[3]].append(i)
                   self.nuclear_charges[i] = NUCLEAR_CHARGE[tokens[3]]
                  
            self.natypes = dict([(key, len(value)) for key,value in self.atomtype_indices.items()])
            
        except IOError as e:
             errno, strerror = e.args
             print( "(ReadMol) I/O error(%s):%s"%(errno,strerror) )
        except ValueError:
             print( "(ReadMol) Could not convert data to an integer." )
        except:
             print( "(ReadMol) Unexpected error:", sys.exc_info()[0] )
             raise

    def ReadMol2(self, filename):
        """
          Read molecule structure from MOL format
        """
        try:
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            tokens = lines[2].split()
            self.natoms = int(tokens[0])
            self.atomtypes = []
            self.nuclear_charges = np.empty(self.natoms, dtype=int)
            self.coordinates = np.empty((self.natoms, 3), dtype=float)

            self.name = filename
            for i, line in enumerate(lines[7:self.natoms+7]):
                tokens = line.split()
                if len(tokens) < 8: 
                   break
                else:
                   self.coordinates[i] = np.asarray(tokens[2:5], dtype=float)
                   atype = tokens[5].split(".")
                   self.atomtypes.append(atype[0])
                   self.atomtype_indices[atype[0]].append(i)
                   self.nuclear_charges[i] = NUCLEAR_CHARGE[atype[0]]
                  
            self.natypes = dict([(key, len(value)) for key,value in self.atomtype_indices.items()])
            
        except IOError as e:
             errno, strerror = e.args
             print( "(ReadMol2) I/O error(%s):%s"%(errno,strerror) )
        except ValueError:
             print( "(ReadMol2) Could not convert data to an integer." )
        except:
             print( "(ReadMol2) Unexpected error:", sys.exc_info()[0] )
             raise


