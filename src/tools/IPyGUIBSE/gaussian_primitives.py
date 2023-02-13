__doc__="""
# Id: gaussian_primitives.py,v 1.0
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
import plotly.graph_objects as g


def VOLUME(c, _x, x, _y, y, _z, z, n):
    X, Y, Z = np.meshgrid(np.linspace(_x, x, n),
                             np.linspace(_y, y, n),
                             np.linspace(_z, z, n))
    
    return X+c[0], Y+c[1], Z+c[2]

def s_gaussian(a, r2):
      """
         s-type gaussian function primitive:
         g = N * (x^m * y^n * z^l) * exp(-a*r^2)
         r^2 = x^2 + y^2 + z^2
         m+n+l = 0
         N: normalization constant (spherical coords)
      """
      N  = ( 2.0*a/np.pi )**0.75
      gs = N * np.exp( -a*r2 ) 
      return gs

def p_gaussian(a, x, r2):
      """
         p-type gaussian function primitive:
         g = N * (x^m * y^n * z^l) * exp(-a*r^2)
         m+n+l = 1
         N: normalization constant (spherical coords)
      """
      a5  = a**5.0
      pi3 = np.pi**3.0      
      N = ( 128.0 * a5 / pi3 )**0.25
      gp = N * x * np.exp( -a*r2 ) 
      return gp

def d_gaussian(a, x, r2):
      """
         d-type gaussian function primitive:
         g = N * (x^m * y^n * z^l) * exp(-a*r^2)
         m+n+l = 2
         N: normalization constant
      """
      a7  = a**7.0
      pi3 = np.pi**3.0     
      N = ( 2048.0 * a7 / pi3 )**0.25
      gp = N * x * x * np.exp( -a*r2 ) 
      return gp

def orbital_viz_dim(X,Y,Z,Vol,_m, m, title):
    """
        Volumetric plot for visualisation of orbitals
    """
    fig = g.Figure(data=g.Volume(\
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=Vol.flatten(),
        isomin= Vol.min()-0.3,
        isomax= Vol.max()+0.3,
        opacity=0.5,
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
                        scene = dict(                            
                            xaxis_title = 'X-axis',
                            xaxis = dict(
                                nticks=4,
                                range=[_m,m],
                                tickfont=dict(
                                         family="Arial",
                                         color='#000000',
                                         size=25),                       
                            ),
                            yaxis_title = 'Y-axis',
                            yaxis = dict(
                                nticks=4,
                                range=[_m,m],
                                tickfont=dict(
                                         family="Arial",
                                         color='#000000',
                                         size=25),                       
                               
                            ),
                            zaxis_title = 'Z-axis',
                            zaxis = dict(
                                nticks=4,
                                range=[_m,m],
                                tickfont=dict(
                                         family="Arial",
                                         color='#000000',
                                         size=25),                       
                            ),
                            
                        )
    )
    fig.show()
    return None


def orbital_viz(X,Y,Z,Vol,title):
    """
        Volumetric plot for visualisation of orbitals
    """
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
