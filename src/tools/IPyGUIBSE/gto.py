__doc__="""
# Id: gto.py,v 1.0
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
import plotly.graph_objects as g
import sys

class GTO(object):
      def __init__ (self, exponen, coeff, orb_type="s"):
          self.orb_type = orb_type
          self.exponents, self.coefficients = [], []
          print(exponen, coeff)
          if orb_type == "s": 
             for i in range( len(exponen) ):
                 self.exponents.append( float(exponen[i]) )
                 self.coefficients.append( 0.5*float(coeff[i])   )
          elif orb_type == "p":       
             for i in range( len(exponen) ):
                 for j in range( len(exponen) ):
                     self.exponents.append( float(exponen[i]) )
                     self.coefficients.append( float(coeff[i][j])   )
          else:
                print("Orbital type not implemented yet ...")
                sys.exit(1)
                
      def _GTO(self):
            """ Plot GTO atomic orbital
            """
            X,Y,Z,R = self._VOLUME()
            if self.orb_type == "s":
               c = 2.0 / np.pi   
               title = 's-type Atomic Orbital'
               Vol = self.coefficients[0] * ( (c * self.exponents[0])**0.75 ) * np.exp( -self.exponents[0] * R )
               for i in range(1, len(self.exponents) ):
                   Vol += self.coefficients[i] * ( (c * self.exponents[i])**0.75 ) * np.exp( -self.exponents[i] * R )
            elif self.orb_type == "p":
               c = 128.0 / (np.pi**3.0)
               title = 'p-type Atomic Orbital'
               Vol = self.coefficients[0] * ((c * self.exponents[0]**5)**0.25) * X * Y * np.exp( -self.exponents[0] * R )
               for i in range(1, len(self.exponents) ):
                   Vol += self.coefficients[i] * ((c * self.exponents[i]**5)**0.25) * X * Y * np.exp( -self.exponents[i] * R )
                   
            fig = g.Figure(data=g.Volume(\
                x=X.flatten(),
                y=Y.flatten(),
                z=Z.flatten(),
                value=Vol.flatten(),
                isomin= -1.0,
                isomax=  1.0,
                opacity=0.3,
                opacityscale=[[0,1], [0.5, 0.2], [1, 1]],
                surface_count=100,
                colorscale=[[0, 'rgb(0,0,255)'], [1,'rgb(255,0,0)']],
                ),
                layout={'xaxis': {'color': 'green','range': [0, 1]}}
                )

            fig.update_layout(
                                height = 750,
                                width  = 750,
                                title=dict(
                                      text= title,
                                      font=dict(
                                            family="Arial",
                                            color='#000000',
                                            size=25
                                       )),
                                font=dict(
                                     family="Arial",
                                     color='#000000',
                                     size=25),                       
            ##                          xaxis_tickmode='linear',
            ##                          xaxis_tick0 = -4.0,
            ##                          xaxis_dtick = 2,
                                scene = dict(                            
                                    xaxis_title = 'X-axis',
                                    xaxis = dict(
                                        nticks=4,
                                        range=[-5,5],
                                        tickfont=dict(
                                                 family="Arial",
                                                 color='#000000',
                                                 size=25),                       
                                    ),
                                    yaxis_title = 'Y-axis',
                                    yaxis = dict(
                                        nticks=4,
                                        range=[-5,5],
                                        tickfont=dict(
                                                 family="Arial",
                                                 color='#000000',
                                                 size=25),                       
                                       
                                    ),
                                    zaxis_title = 'Z-axis',
                                    zaxis = dict(
                                        nticks=4,
                                        range=[-5,5],
                                        tickfont=dict(
                                                 family="Arial",
                                                 color='#000000',
                                                 size=25),                       
                                    ),
                                    
                                )
            )
            fig.show()
            return None

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
      
      
      
