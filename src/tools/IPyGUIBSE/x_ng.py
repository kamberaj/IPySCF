__doc__="""
# Id: x_ng.py,v 1.0
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

class X_NG(object):
      def __init__ (self, exponen, coeff, orb_type="1s"):
          if orb_type == "1s": 
             self.orb_1s_type        = False
             self.orb_2p_type        = False
             self.orb_sp_inner_type  = False
             self.orb_sp_outer1_type = False
             self.orb_sp_outer2_type = False
             self.orb_sp_outer3_type = False
             self.orb_1d_type        = False
             self.set_s(exponen, coeff)
             
          elif orb_type == "sp-inner":       
             self.set_sp_inner(exponen, coeff)

          elif orb_type == "sp-outer1":       
             self.set_sp_outer1(exponen, coeff)

          elif orb_type == "sp-outer2":       
             self.set_sp_outer2(exponen, coeff)

          elif orb_type == "sp-outer3":       
             self.set_sp_outer3(exponen, coeff)

          elif orb_type == "1d":       
             self.set_d(exponen, coeff)
             
          elif orb_type == "2p":       
             self.set_p(exponen, coeff)
 
          else:
                print("Orbital type not implemented yet ...")
                sys.exit(1)

      def set_d(self, exponen, coeff):
          self.orb_1d_type        = True
          self.orb_1s_type        = False
          self.orb_2p_type        = False
          self.orb_sp_inner_type  = False
          self.exponents_1d, self.coefficients_1d = [], []
          for i in range( len(exponen) ):
              self.exponents_1d.append( float(exponen[i]) )
              self.coefficients_1d.append( float(coeff[0][i])   )
          return None

      def set_s(self, exponen, coeff):
          self.orb_1s_type = True
          self.exponents, self.coefficients_1s = [], []
          for i in range( len(exponen) ):
              self.exponents.append( float(exponen[i]) )
              self.coefficients_1s.append( float(coeff[0][i])   )
          return None

      def set_p(self, exponen, coeff):
          self.orb_1d_type        = False
          self.orb_1s_type        = False
          self.orb_sp_inner_type  = False
          self.orb_2p_type        = True
          self.exponents_2p, self.coefficients_2p = [], []
          for i in range( len(exponen) ):
              self.exponents_2p.append( float(exponen[i]) )
              self.coefficients_2p.append( float(coeff[0][i])   )
          return None

      def set_sp_inner(self, exponen, coeff):
          self.orb_1d_type        = False
          self.orb_1s_type        = False
          self.orb_2p_type        = False
          self.orb_sp_outer1_type = False
          self.orb_sp_outer2_type = False
          self.orb_sp_outer3_type = False
          self.orb_sp_inner_type  = True
          self.exponents_inner, self.coefficients_2s_inner, self.coefficients_2p_inner = [], [], []
          for i in range( len(exponen) ):
              self.exponents_inner.append( float(exponen[i]) )
              self.coefficients_2s_inner.append( float(coeff[0][i])   )
              self.coefficients_2p_inner.append( float(coeff[1][i])   )
          return None

      def set_sp_outer1(self, exponen, coeff):
          self.orb_1d_type        = False
          self.orb_1s_type        = False
          self.orb_2p_type        = False
          self.orb_sp_outer2_type = False
          self.orb_sp_outer3_type = False
          self.orb_sp_outer1_type = True
          self.exponents_outer1, self.coefficients_2s_outer1, self.coefficients_2p_outer1 = [], [], []
          for i in range( len(exponen) ):
              self.exponents_outer1.append( float(exponen[i]) )
              self.coefficients_2s_outer1.append( float(coeff[0][i])   )
              self.coefficients_2p_outer1.append( float(coeff[1][i])   )
          return None

      def set_sp_outer2(self, exponen, coeff):
          self.orb_1d_type        = False
          self.orb_1s_type        = False
          self.orb_2p_type        = False
          self.orb_sp_outer1_type = False
          self.orb_sp_outer3_type = False
          self.orb_sp_outer2_type = True
          self.exponents_outer2, self.coefficients_2s_outer2, self.coefficients_2p_outer2 = [], [], []
          for i in range( len(exponen) ):
              self.exponents_outer2.append( float(exponen[i]) )
              self.coefficients_2s_outer2.append( float(coeff[0][i])   )
              self.coefficients_2p_outer2.append( float(coeff[1][i])   )
          return None

      def set_sp_outer3(self, exponen, coeff):
          self.orb_1d_type        = False
          self.orb_1s_type        = False
          self.orb_2p_type        = False
          self.orb_sp_outer1_type = False
          self.orb_sp_outer2_type = False
          self.orb_sp_outer3_type = True
          self.exponents_outer3, self.coefficients_2s_outer3, self.coefficients_2p_outer3 = [], [], []
          for i in range( len(exponen) ):
              self.exponents_outer3.append( float(exponen[i]) )
              self.coefficients_2s_outer3.append( float(coeff[0][i])   )
              self.coefficients_2p_outer3.append( float(coeff[1][i])   )
          return None
      
      def _X_NG(self):
            """ Plot X-NG atomic orbital
            """
            X,Y,Z,R = self._VOLUME()
            if self.orb_1s_type:
               title = '1s-type Atomic Orbital'
               Vol = self.coefficients_1s[0] * s_gaussian( self.exponents[0], R )
               for i in range(1, len(self.exponents) ):
                   Vol += self.coefficients_1s[i] * s_gaussian( self.exponents[i], R )

               orbital_viz(X,Y,Z,Vol,title)
               
            if self.orb_sp_inner_type:
                  title = '2s-type Atomic Orbital'
                  Vol = self.coefficients_2s_inner[0] * s_gaussian( self.exponents_inner[0], R )
                  for i in range(1, len(self.exponents_inner) ):
                      Vol += self.coefficients_2s_inner[i] * s_gaussian( self.exponents_inner[i], R )

                  if self.orb_sp_outer1_type:
                     for i in range( len(self.exponents_outer1) ):
                         Vol += self.coefficients_2s_outer1[i] * s_gaussian( self.exponents_outer1[i], R )
                  if self.orb_sp_outer2_type:
                     for i in range( len(self.exponents_outer2) ):
                         Vol += self.coefficients_2s_outer2[i] * s_gaussian( self.exponents_outer2[i], R )
                  if self.orb_sp_outer3_type:
                     for i in range( len(self.exponents_outer3) ):
                         Vol += self.coefficients_2s_outer3[i] * s_gaussian( self.exponents_outer3[i], R )

                  orbital_viz(X,Y,Z,Vol,title)
               
                  title = '2p-type Atomic Orbital'    
                  Vol = self.coefficients_2p_inner[0] * p_gaussian( self.exponents_inner[0], X, R )
                  for i in range(1, len(self.exponents_inner) ):
                      Vol += self.coefficients_2p_inner[i] * p_gaussian( self.exponents_inner[i], X, R )

                  if self.orb_sp_outer1_type:
                     for i in range( len(self.exponents_outer1) ):
                         Vol += self.coefficients_2p_outer1[i] * p_gaussian( self.exponents_outer1[i], X, R )
                  if self.orb_sp_outer2_type:
                     for i in range( len(self.exponents_outer2) ):
                         Vol += self.coefficients_2p_outer2[i] * p_gaussian( self.exponents_outer2[i], X, R )
                  if self.orb_sp_outer3_type:
                     for i in range( len(self.exponents_outer3) ):
                         Vol += self.coefficients_2p_outer3[i] * p_gaussian( self.exponents_outer3[i], X, R )

                  orbital_viz(X,Y,Z,Vol,title)
      
            if self.orb_1d_type:
                  title = '1d-type Atomic Orbital'    
                  Vol = self.coefficients_1d[0] * d_gaussian( self.exponents_1d[0], X, R )
                  for i in range(1, len(self.exponents_1d) ):
                      Vol += self.coefficients_1d[i] * d_gaussian( self.exponents_1d[i], X, R )

                  orbital_viz(X,Y,Z,Vol,title)

            if self.orb_2p_type:
                  title = '2p-type Atomic Orbital'    
                  Vol = self.coefficients_2p[0] * p_gaussian( self.exponents_2p[0], X, R )
                  for i in range(1, len(self.exponents_2p) ):
                      Vol += self.coefficients_2p[i] * p_gaussian( self.exponents_2p[i], X, R )

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
      
      
      
