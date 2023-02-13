__doc__="""
# Id: sto_ng.py,v 1.0
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
import sys
from .gaussian_primitives import *

class STO_NG(object):
      def __init__ (self, exponen, coeff, orb_type="s"):
          self.orb_type = orb_type
          if orb_type == "s": 
             self.exponents, self.coefficients_1s = [], []
             for i in range( len(exponen) ):
                 self.exponents.append( float(exponen[i]) )
                 self.coefficients_1s.append( float(coeff[0][i])   )
          elif orb_type == "sp":       
             self.exponents, self.coefficients_2s, self.coefficients_2p = [], [], []
             for i in range( len(exponen) ):
                  self.exponents.append( float(exponen[i]) )
                  self.coefficients_2s.append( float(coeff[0][i])   )
                  self.coefficients_2p.append( float(coeff[1][i])   )
          else:
                print("Orbital type not implemented yet ...")
                sys.exit(1)
                
      def _STO_NG(self):
            """ Plot GTO atomic orbital
            """
            X,Y,Z,R = self._VOLUME()
            if self.orb_type == "s":
               title = '1s-type Atomic Orbital'
               Vol = self.coefficients_1s[0] * s_gaussian( self.exponents[0], R )
               for i in range(1, len(self.exponents) ):
                   Vol += self.coefficients_1s[i] * s_gaussian( self.exponents[i], R )

               orbital_viz(X,Y,Z,Vol,title)
               
            elif self.orb_type == "sp":
               title = '2s-type Atomic Orbital'
               Vol = self.coefficients_2s[0] * s_gaussian( self.exponents[0], R )
               for i in range(1, len(self.exponents) ):
                   Vol += self.coefficients_2s[i] * s_gaussian( self.exponents[i], R )

               orbital_viz(X,Y,Z,Vol,title)
               
               title = '2p-type Atomic Orbital'    
               Vol = self.coefficients_2p[0] * p_gaussian( self.exponents[0], X, R )
               for i in range(1, len(self.exponents) ):
                   Vol += self.coefficients_2p[i] * p_gaussian( self.exponents[i], X, R )

               orbital_viz(X,Y,Z,Vol,title)
            

      def _VOLUME(self):
            X_c, Y_c, Z_c = 0.0, 0.0, 0.0
            X, Y, Z    = np.meshgrid(np.linspace(-5.0, 5.0, 40),
                                     np.linspace(-5.0, 5.0, 40),
                                     np.linspace(-5.0, 5.0, 40))
            Xc = X - X_c
            Yc = Y - Y_c
            Zc = Z - Z_c

            R = Xc**2 + Yc**2 + Zc**2
            return Xc,Yc,Zc,R
      
      
      
