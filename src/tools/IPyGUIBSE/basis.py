__doc__="""
# Id: basis.py,v 1.0
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

import os
from typing import Dict, List, Optional, Union

import basis_set_exchange as bse


class Basis:
    """Specify the AO basis."""

    def __init__(self,
                 basis: Union[Dict[str, str], str],
                 ri: Optional[Union[Dict[str, str], str]] = None,
                 admm: Optional[Union[Dict[str, str], str]] = None) -> None:
        """Initialize Basis instance.

        Args:
            basis: Basis set for molecule. Basis can be a string, if the same basis set
                should be applied to all atoms. Basis can be a dict of the type
                {atom_name: basis} if different basis sets are to be used for
                different atoms.
            ri: Auxiliary basis set for RI (optional). Follows the same format as basis.
            admm: Auxiliary basis set for ADMM (optional). Follows the same format as basis.
        """
        validate_basis(basis)
        self.basis = basis
        if ri is not None:
            validate_basis(ri)
        self.ri = ri
        if admm is not None:
            validate_basis(admm)
        self.admm = admm

    def write(self, basis_format: str = 'dalton', data_dir: str = 'data') -> None:
        """Write basis set to file.

        Args:
            basis_format: Format of the basis set file.
        """
        for basis_set in (self.basis, self.ri, self.admm):
            # ri and admm sets can be empty
            if not basis_set:
                continue
            # basis can be either a str or dict
            # so we make a dummy dict for when it is str
            if isinstance(basis_set, str):
                basis_dict = {'dummy': basis_set}
            else:
                basis_dict = basis_set
            for basis in basis_dict.values():
                if os.path.isfile(basis):
                    continue
                bse_basis = bse.get_basis(name=basis,fmt=basis_format, data_dir=data_dir)
                bse.write_formatted_basis_file(basis_dict=bse_basis, outfile_path=basis, basis_fmt=basis_format)
 
                
    def read(self, basis_format: str = 'nwchem', data_dir: str = 'data') -> None:
        """read basis set.

        Args:
            basis_format: Format of the basis set file.
        """
        for basis_set in (self.basis, self.ri, self.admm):
            # ri and admm sets can be empty
            if not basis_set:
                continue
            # basis can be either a str or dict
            # so we make a dummy dict for when it is str
            if isinstance(basis_set, str):
                basis_dict = {'dummy': basis_set}
            else:
                basis_dict = basis_set
            for basis in basis_dict.values():
                if os.path.isfile(basis):
                    continue
                bse_basis = bse.get_basis(name=basis,fmt=basis_format, data_dir=data_dir)



def validate_basis(basis: Union[Dict[str, str], str]) -> None:
    """Validate basis set input.

    Args:
        basis: Basis set.
    """
    if not isinstance(basis, str) and not isinstance(basis, dict):
        raise TypeError(f'Unsupported basis set type: {type(basis)}.')
    elif isinstance(basis, dict):
        for label, basis_type in basis.items():
            if not isinstance(label, str):
                raise TypeError(f'Invalid atom label type: {type(label)}.')
            if not isinstance(basis_type, str):
                raise TypeError(f'Unknown basis type: {type(basis_type)}.')


def get_atom_basis(basis: Union[str, Dict[str, str]], num_atoms: int, labels: List[str]) -> List[str]:
    """Process basis set input.

    Args:
        basis: Basis set.
        num_atoms: Number of atoms.
        labels: Atom labels.

    Returns:
        List containing the basis set of each individual atoms.
    """
    if isinstance(basis, str):
        atom_basis = [basis] * num_atoms
    elif isinstance(basis, dict):
        atom_basis = []
        for label in labels:
            try:
                atom_basis.append(basis[label])
            except KeyError:
                raise KeyError(f'no basis for {label} was provided.')
    else:
        raise TypeError(f'basis must be of type str or dict and not {type(basis)}.')
    return atom_basis


