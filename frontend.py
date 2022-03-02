"""library of widgets exposed to user via Compute_Spectra_Workflow notebook"""
#data tools
import pandas as pd
import numpy as np
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure

#simtool loading and interface
from simtool import findInstalledSimToolNotebooks,searchForSimTool
from simtool import getSimToolInputs,getSimToolOutputs,Run

#user interface utilities 
import os, stat
import ipywidgets as widgets
from IPython.display import display
from IPython.display import clear_output

#visualization tools
import qgrid
import plotly
#import matplotlib.pyplot as plt

# modification of Hublib's HideCodeButton
class HideCodeButton(object):

    def __init__(self, **kwargs):

        label = kwargs.get('label', None)
        icon = kwargs.get('icon', '')
        tooltip = kwargs.get('tooltip', '')
        style = kwargs.get('style', '')
        bcb = kwargs.get('cb', None)
        runall = kwargs.get('RunAll', False)

        self.func = 'hide'

        if label is None:
            self.hlabel = "Hide Code Cells"
            self.slabel = "Show Code Cells"
        elif isinstance(label, list):
            self.hlabel, self.slabel = label
        elif isinstance(label, str):
            self.hlabel = self.slabel = label
        else:
            raise ValueError('label should be a string or list containing two strings')

        if style is None:
            self.hstyle = ''
            self.sstyle = ''
        elif isinstance(style, list):
            self.hstyle, self.sstyle = style
        elif isinstance(style, str):
            self.hstyle = self.sstyle = style
        else:
            raise ValueError('style should be a string or list containing two strings')

        self.w = widgets.Button(
            description=self.hlabel,
            icon=icon,
            tooltip=tooltip,
            button_style=self.hstyle
            )

        def button_cb(ignore):
            js = Javascript("$('div.input').%s()" % self.func)
            if self.func == 'hide':
                self.func = 'show'
                self.w.description = self.slabel
                self.w.button_style = self.sstyle
            else:
                self.func = 'hide'
                self.w.description = self.hlabel
                self.w.button_style = self.hstyle

            if bcb is not None:
                bcb()
            display(js)

        self.w.on_click(button_cb)

    def _ipython_display_(self):
        self.w._ipython_display_()
