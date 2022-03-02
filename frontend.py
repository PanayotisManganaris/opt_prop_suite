"""library of widgets exposed to user via Compute_Spectra_Workflow notebook"""
#data tools
import pandas as pd
import numpy as np
import qgrid
from pymatgen.ext.matproj import MPRester

#simtool loading and interface
from simtool import findInstalledSimToolNotebooks,searchForSimTool
from simtool import getSimToolInputs,getSimToolOutputs,Run

#user interface utilities 
import os, stat
import plotly.express as px
import ipywidgets as widgets
from IPython.display import display
from IPython.display import clear_output

class Authenticate:
    """
    Container class for user authentication functions

    eventually provide an interface to authenticate into multiple data providers
    - materials project
    - citrine informatics
    - nomad?

    for now, just materials project is provided.
    (current remote query apps work only through MP anyway)
    (look into OPTIMATE?)
    """
    @staticmethod
    def mpkey():
        """
        write user's materials project key to a file which can be
        accessed by later instances of MPRester. Delete key from
        kernel vars dict.
        """
        try:
            mpkey = str(input('Paste MP API key: '))
            clear_output()
            if not mpkey.isalnum():
                raise TypeError('Wrong Key')
            if mpkey is None:
                raise TypeError('Empty')
            with open(os.path.expanduser('~/.mpkey.txt'), 'w') as keyfile:
                keyfile.write(mpkey)
            del mpkey
            os.chmod(os.path.expanduser('~/.mpkey.txt'), stat.S_IREAD | stat.S_IWRITE)
            print("Success")
        except:
            print("Something seems wrong with your key")


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
