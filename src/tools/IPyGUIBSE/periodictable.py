__doc__="""
# Id: periodictable.py,v 1.0
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
import json

import tkinter
from tkinter import *
import tkinter as tk

import webbrowser
from locale import setlocale, LC_ALL, format_string
from typing import Optional, Union, Iterable, Iterator, NamedTuple

import basis_set_exchange as bse

from . elementbasisdata import elementbasisDataDialog


from tkinter import messagebox

import numpy as np

from .utils   import *

class Element(NamedTuple):
    symbol: str
    name: str
    number: int
    category: str
    group: Union[str, int]
    period: int
    block: str
    mass: float
    phase: Optional[str]
    density: Optional[float]
    electronegativity: Optional[float]


class PlacedElement(NamedTuple):
    row: int
    column: int
    element: Element


def format_float(s: Optional[float]) -> str:
    if s is None:
        return 'n.A'
    return format_string('%.2f', s)


def place_elements(elements: Iterable[Element]) -> Iterator[PlacedElement]:
    OFFSET = 2
    la_offset = 2
    ac_offset = 2

    for element in elements:
        period, group_name = element.period, element.group

        if group_name == 'La':
            group = la_offset + OFFSET
            la_offset += 1
            period += OFFSET
        elif group_name == 'Ac':
            group = ac_offset + OFFSET
            ac_offset += 1
            period += OFFSET
        else:
            group = group_name

        yield PlacedElement(row=period - 1, column=group - 1, element=element)


def load_json(filename: str = 'src/tools/IPyGUIBSE/elements.json') -> Iterator[Element]:
    with open(filename, encoding='utf-8') as f:
        for element_dict in json.load(f):
            yield Element(**element_dict)


class ElementButton:
    BORDER = 3

    CATEGORY_COLORS = {
        'Alkalimetalle': '#fe6f61',
        'Erdalkalimetalle': '#6791a7',
        'Übergangsmetalle': '#83b8d0',
        'Metalle': '#cae2ed',
        'Halbmetalle': '#a7d6bc',
        'Nichtmetalle': '#ffde66',
        'Halogene': '#e9aa63',
        'Edelgase': '#e29136',
        'Unbekannt': '#cec0bf',
        'Lanthanoide': '#696071',
        'Actinoide': '#5b4c68',
    }

    PHASE_COLORS = {
        'fest': 'black',
        'flüssig': 'blue',
        'gasförmig': 'red',
        None:  'grey',
    }

    def __init__(self, parent: tk.Widget, placed_element: PlacedElement) -> None:
        self.element = placed_element.element
        self.background = self.CATEGORY_COLORS[self.element.category]
        self.frame = frame = tk.Frame(
            parent, relief=tk.RAISED,
            name=f'frame_{self.element.symbol}',
            background=self.background,
            border=self.BORDER,
        )
        self.frame.grid_columnconfigure(1, weight=2)
        self.frame.grid(row=placed_element.row, column=placed_element.column, sticky=tk.EW)

        self.populate()

        frame.bind('<ButtonPress-1>', self.press)
        frame.bind('<ButtonRelease-1>', self.release)
        for child in frame.winfo_children():
            child.bindtags((frame,))

    def populate(self) -> None:
        prefix = f'label_{self.element.symbol}_'

        tk.Label(
            self.frame, name=prefix + 'number',
            text=self.element.number, background=self.background,
        ).grid(row=0, column=0, sticky=tk.NW)

        tk.Label(
            self.frame, name=prefix + 'mass',
            text=format_float(self.element.mass), background=self.background,
        ).grid(row=0, column=2, sticky=tk.NE)

        tk.Label(
            self.frame, name=prefix + 'name',
            text=self.element.name, background=self.background,
        ).grid(row=1, column=0, sticky=tk.EW, columnspan=3)

        tk.Label(
            self.frame, name=prefix + 'symbol',
            text=self.element.symbol, font='bold', background=self.background,
            foreground=self.PHASE_COLORS[self.element.phase],
        ).grid(row=2, column=0, sticky=tk.EW, columnspan=3)

        tk.Label(
            self.frame, name=prefix + 'electronegativity',
            text=format_float(self.element.electronegativity), background=self.background,
        ).grid(row=3, column=0, sticky=tk.SW)

        tk.Label(
            self.frame, name=prefix + 'density',
            text=format_float(self.element.density), background=self.background,
        ).grid(row=3, column=2, sticky=tk.SE)

    def press(self, event: tk.Event) -> None:
        self.frame.configure(relief='sunken')
        
    def release(self, event: tk.Event) -> None:
        self.frame.configure(relief='raised')
        print(self.element.name, " , ", self.element.symbol)
        self.BSEElementBasisData()
        webbrowser.open(
            url=f'https://de.wikipedia.org/wiki/{self.element.name}',
            new=2,
        )

    def BSEElementBasisData(self):
        fe = Tk()   
        feFrame = elementbasisDataDialog(fe, self.element.symbol, self.element.number)
        fe.title("Interactive BSE Basis Data Dialog")
        fe.geometry("900x400")
        fe.resizable(width=TRUE, height=TRUE)
        fe.mainloop()
        return None

def periodictable() -> None:
    setlocale(LC_ALL, 'en_US.UTF-8')

    root = tk.Tk()
    root.title('Periodic system of Elements')

    frame = tk.Frame(root, name='grid_container')
    frame.pack_configure(fill=tk.BOTH)

    elements = tuple(place_elements(load_json()))
    for element in elements:
        ElementButton(frame, element)

    columns = {elm.column for elm in elements}
    for x in columns:
        frame.grid_columnconfigure(index=x, weight=1)

    root.mainloop()


if __name__ == '__main__':
    periodictable()
