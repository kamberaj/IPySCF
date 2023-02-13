import logging
import pydoc

from . import program

__all__ = (
    'program',
)

__version__ = '0.1a0'
__author__ = 'The IPySCF Project Developers'
__copyright__ = 'Copyright (C)  The IPySCF Project Developers'
__license__ = 'GPLv3+'

LICENSE = __doc__="""
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


def print_license():
    """Print the license."""
    pydoc.pager(LICENSE)
