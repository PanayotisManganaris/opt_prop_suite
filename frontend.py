"""
library of widgets facilitating material structure input intended for
nanohub simtools to be published using the Anaconda2020.11 jupyter
tool default kernel.
"""
#data tools
import pandas as pd
import numpy as np
import qgrid
from pymatgen.ext.matproj import MPRester

#data handling
from pymatgen.core import Structure

#simtool loading and interface
from simtool import findInstalledSimToolNotebooks,searchForSimTool
from simtool import getSimToolInputs,getSimToolOutputs,Run

#user interface utilities 
import os, stat, io
import plotly.express as px
import ipywidgets as widgets
from IPython.display import display, Javascript, clear_output

from functools import partial
import re
import codecs

class Authenticate():
    """
    Nanohub friendly user authentication widgets.
    If users already have a key, it is not then requested by default.

    Provides interfaces to authenticate into multiple data providers
    - materials project
    - citrine informatics?
    - nomad?
    """
    def __init__(self):
        #display managers
        self.mpout = widgets.Output()
        #AutoAuth Checks
        self.force = False
        self.checkmp = os.path.isfile(os.path.expanduser("~/.mpkey.txt"))
        #Control and Base Widgets
        self.authbox = widgets.Text(value="", description='Token', disabled=False)
        self.good = widgets.IntProgress(value=10, min=0, max=10,
                                        description="You're In!", bar_style="success")
        self.renew_button = widgets.Button(description='Renew Tokens',
                                           disabled=False,
                                           button_style='warning')
        self.renew_button.on_click(self.forceAuth)
        #Auth widgets
        self.mpkey = partial(self._token, out=self.mpout, authmethod=self._mpkey, check=self.checkmp)

    def forceAuth(self, event):
        self.force = True
        self.mpkey()
        self.force = False        

    def _token(self, out, authmethod, check):
        """
        generic conditional authentication
        """
        with out:
            clear_output(wait=True)
            if check and not self.force:
                display(self.good)
            else:
                self.authbox.on_submit(authmethod)
                display(self.authbox)
        return out
    
    def _mpkey(self, text):
        """
        mpkey authmethod -- write user's materials project key to a
        file which can be accessed by later instances of
        MPRester. Delete key from kernel vars dict.
        """
        mpkey = text.value
        text.value = ""
        with self.mpout:
            clear_output(wait=True)
            assert mpkey, 'MP keys cannot be empty strings, entry not written'
            assert mpkey.isalnum(), 'MP keys should consist of numbers and letters, entry not written'
            with open(os.path.expanduser('~/.mpkey.txt'), 'w') as keyfile:
                keyfile.write(mpkey)
            del mpkey
            os.chmod(os.path.expanduser('~/.mpkey.txt'), stat.S_IREAD | stat.S_IWRITE)
            display(self.good)
            print("Entry written. If your MP queries return authentication errors, try reentering your key")

class MPQuery():
    """
    The Materials Project Mutable Query Object. Used to get and apply MP data.
    Note: structure retrieval requires "query" return list include the "task_id" key
    """
    def __init__(self, content=[], **kwargs):
        """current stopgap:
        if OPTIMATE is hard to implement, just make case-by-case query
        object with the following API

        authentication is pulled from the user's filesystem to offload security onto them...
        frame attribute returns a dataframe of requested values
        """
        self.authfile = kwargs.get('authfile', "~/.mpkey.txt")
        self.query = kwargs.get('query', ({ "crystal_system": "cubic"},
                                          ["task_id","formula", "elements", "e_above_hull",
                                           "spacegroup.number", "band_gap", "crystal_system"]))
        self.dbcode_key = "task_id"
        self.rester = self._make_rester()
        self.frame = pd.DataFrame(self.get_content(content)) #mutable

    def _make_rester(self):
        with open(os.path.expanduser(self.authfile), "r+") as file:
            apikey = file.readline()
        return MPRester(apikey) #investigate other apis, plug and play parent? OPTIMATE?

    def get_content(self, content):
        """mutate list of MP results"""
        if not content: #should persist through every instantiation of QueryObject
            content += self.rester.query(*self.query)
        return content

    def get_structure(self, mpid):
        return self.rester.get_structure_by_material_id(mpid,
                                                        final=False,
                                                        conventional_unit_cell=False)

class FakeQuery():
    """ A minimal Query Object serving both default and debugging purposes """
    def __init__(self, content=[], Debug=False, **kwargs):
        self.dbcode_key = "task_id"
        self.frame = pd.DataFrame(["You must authenticate into a remote data source and run a query to use this interface."],
                                  index=["alert!"], columns=["attention!"])
        if Debug:
            self.frame = pd.read_pickle("./cubics.pkl")
        monitor = widgets.Output()

    def get_structure(self, mpid):
        raise NotImplementedError("No rest api needed to see this output!")

def struct_plot(struct):
    """ 3-D Interactive plot function for inspecting the unit cell of a material """
    POSCAR_str = struct.to(fmt = "poscar")
    lines = POSCAR_str.split('\n')
    # get the lattice information
    lattice = lines[1]
    cell_vectors = np.array([lines[2].split() , lines[3].split() , lines[4].split()]).astype(float)
    # get the list of sites
    sites = []
    for line in lines[8:]:
        if not line:
            break
        sites.append([line.split()[0],line.split()[1],line.split()[2]])
    # convert from fractional to xyz
    sites = np.array(sites).astype(float)
    xyz = np.matmul(sites,cell_vectors).transpose()
    # get the coordinates of the box
    corners = np.array([[0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,0],
                        [0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,1],
                        [0,0,0,0,0,1,1,1,1,1,1,1,0,0,1,1,0]]).T
    cell = np.matmul(corners,cell_vectors).transpose()
    box = pd.DataFrame({'x':cell[0],'y':cell[1],'z':cell[2]})
    # get a color dictionary
    elements = lines[5].split()
    atoms = lines[6].split()
    hues = ['tab:blue','tab:orange','tab:green','tab:red']
    colors = []
    for i, atom in enumerate(atoms):
        colors.extend([elements[i]]*int(atoms[i]))

    zip_iterator = zip(elements,hues)
    color_dict = dict(zip_iterator)

    size = 72
    POSCAR_df = pd.DataFrame({'x':xyz[0],'y':xyz[1],'z':xyz[2],'element':colors,'size':size})    
    fig = px.scatter_3d(POSCAR_df, x='x', y='y', z='z',
              color='element',labels={'x':'','y':'','z':''},size='size',size_max=size)

    box = px.line_3d(box,x='x', y='y', z='z')
    fig.add_traces(list(box.select_traces()))

    fig.update_layout(scene = dict(
                    xaxis = dict(
                        nticks=0,showbackground=False,showticklabels=False,),
                    yaxis = dict(
                        nticks=0,showbackground=False,showticklabels=False,),
                    zaxis = dict(
                        nticks=0,showbackground=False,showticklabels=False,),),
                    width=700,
                    margin=dict(r=10, l=10, b=10, t=10)
                  )
    fig.show()

class QueryPanel():
    """
    widget class for performing queries. Currently necessary to use
    even if a tool does not require remote access. InputSuite must be
    instantiated using QueryPanel
    
    1. interactively instantiate and swap between query objects using the Toggles widget
    2. Collect and display query status with the progressout widget
    wishlist: customize query scope using a panel of interactive sliders
    feed the result to the respective query objects.
    probably a post-OPTIMATE thing
    """
    def __init__(self):
        self.progressout = widgets.Output(layout={'border': "5px solid black"})
        self.toggles = widgets.ToggleButtons(
            options=["None",
                     "Debug", #comment for release
                     "Materials Project"],
            description="Remote Data: ",
            disabled=False,
            button_style="",
            tooltips=["load a saved query for quick debugging",
                      "Contact the Materials Project REST API for semiconductor information"]
        )
        self.Q = FakeQuery() #default
        self.toggles.observe(self._assign_data_on_pick, names='value')
        self.InputSuite = None

    def _assign_data_on_pick(self, change):
        """ observes update to the traitlets bunch """
        with self.progressout: #stdout captured by default
            if self.toggles.value == "None":
                self.Q = FakeQuery()
            elif self.toggles.value == "Materials Project":
                self.Q = MPQuery()
            elif self.toggles.value == "Debug":
                self.Q = FakeQuery(Debug=True)
        self.InputSuite.remoteout #Must instantiate QueryPanel then InputSuite

class ThisString():
    """
    contains logic for determining a crystallographic file format string from the string content
    """
    def __init__(self, string):
        self.string = string

    def is_cif(self):
        return False

    def is_poscar(self):
        return True

    def format_is(self):
        if self.is_cif():
            return 'cif'
        if self.is_poscar():
            return 'poscar'
        if self.is_cssr():
            return 'cssr'
        if self.is_json():
            return 'json'
        if self.is_yaml():
            return 'yaml'
        if self.is_xsf():
            return 'xsf'
        if self.is_mcsqs():
            return 'mcsqs'

class InputSuite():
    """
    create and return widgets for selecting and reviewing a structure from various sources
    1. use instance plotout widget to view a selected structure
    2. use remote_menu to select a structure from an authorized and connected REST API
    3. use copybox to write or copy/paste a plain text structure acceptable by pymatgen
    4. use upload_button and file_menu to select from a batch of structure files of potentially mixed types
    """
    def __init__(self, QueryPanel):
        self._remoteout = widgets.Output()
        self._filesout = widgets.Output()
        self.plotout = widgets.Output(layout={'border': '5px solid black'})
        # remote interface
        self.QP = QueryPanel
        setattr(self.QP, "InputSuite", self)
        # write widgets
        self.copybox = widgets.Textarea(
            value="",
            placeholder="formula\nscale\nx y z\nel\na b c\netc... ",
            description="Structure: ",
            disabled=False,
            layout=widgets.Layout(height="300px", width="auto")
        )
        self.plotbutton = widgets.Button(
            description='Plot',
            disabled=False,
            button_style='info',
            tooltip='Generates Structure and Plots'
        )
        self.plotbutton.on_click(self._process_input)
        # file widgets
        self.upload_button = widgets.FileUpload(
            accept='',
            multiple=True,
        )
        self._get_uploads() #set default
        self.upload_button.observe(self._get_uploads, names='value') #also a callback

    #complex display handlers
    @property
    def remoteout(self):
        with self._remoteout:
            clear_output(wait=True)
            self._remote_menu = qgrid.show_grid(self.QP.Q.frame) #make protected grid widget with convenience method
            self._remote_menu.observe(self._process_pick_remote)
            display(self._remote_menu)
        return self._remoteout

    @property
    def filesout(self):
        with self._filesout:
            clear_output(wait=True)
            self._files_menu = qgrid.show_grid(self.files_frame)
            self._files_menu.observe(self._process_pick_file)
            display(self._files_menu)
        return self._filesout

    def _get_uploads(self, *args):
        if not self.upload_button.value:
            upload_dict = {'alert': 'upload a batch of files to select from here'}
        else:
            upload_dict = {k:codecs.decode(v['content']) for k,v in self.upload_button.value.items()}
        files = pd.Series(upload_dict)
        files.name = 'file content'
        self.files_frame = files.to_frame()
        self.filesout
        
    # string parsing utility
    def _induce_format(self, raw_string:str):
        """ pymatgen contains no logic for this, so here is a simple stuffer """
        sfmt = ThisString(raw_string).format_is()
        self.struct = Structure.from_str(raw_string,
                                         primitive=False, #only for cifs
                                         sort=False,
                                         fmt=sfmt) #currently hardcoded

    #plotting callbacks:
    def _process_input(self, event):
        """ callback for copybox """
        with self.plotout:
            clear_output(wait=True)
            self._induce_format(self.copybox.value)
            struct_plot(self.struct)

    def _process_pick_remote(self, event):
        """ callback for menu """
        grid = self._remote_menu
        with self.plotout:
            clear_output(wait=True)
            sid = grid.get_changed_df().iloc[grid.get_selected_rows()[0]].loc[self.QP.Q.dbcode_key]
            self.struct = self.QP.Q.get_structure(sid)
            struct_plot(self.struct)

    def _process_pick_file(self, event):
        """ callback for menu """
        grid = self._files_menu
        with self.plotout:
            clear_output(wait=True)
            raw_string = grid.get_changed_df().iloc[grid.get_selected_rows()[0], 0]
            self._induce_format(raw_string)
            struct_plot(self.struct)

class SimSettingSuite():
    def __init__(self): #users_pp_file_list
        installedSimToolNotebooks = findInstalledSimToolNotebooks(simToolName,returnString=True)
        print(installedSimToolNotebooks)

        inputs = getSimToolInputs(simToolLocation)
        inputs

        self.pp_list = []
        for filename in os.listdir("./simtool/pseudo/"):
            f = os.path.join("./simtool/pseudo/", filename)
            # get a list of all the PPs -- is this best instantiated here or globally?
            # if instanced here, the user could probably pass their own PPs to the constructor as well
            if os.path.isfile(f):
                self.pp_list.append(filename)

        # TODO: filter by selected compound compositions
        self.filtered_pp_list = self.pp_list

    def _create_widgets(self):
        self.log = widgets.Select(
            options=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            value='DEBUG',
            # rows=10,
            description='Log Level:',
            disabled=False
        )
        self.walltime = widgets.Text(
            value='01:00:00',
            placeholder='walltime',
            description='walltime:',
            disabled=False
        )
        self.numnodes = widgets.IntText(
            value=8,
            placeholder='nodes',
            description='nodes:',
            disabled=False
        ) 
        self.button = widgets.Button(
            description='run simtool',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            tooltip='run to submit qe simtool'
        )
        self.button.on_click(self._on_button_clicked)
        self.ecutwfc = widgets.BoundedFloatText(
            value=50,
            min=50,
            max=400,
            step=10,
            description='ecutwfc:',
            disabled=False
        )
        self.ecutrho = widgets.BoundedFloatText(
            value=200,
            min=200,
            max=1600,
            step=40,
            description='ecutrho:',
            disabled=False
        )
        self.smearing = widgets.Select(
            options=['smearing','fixed'],
            value='fixed',
            rows = 2,
            description='smearing:',
            disabled=False
        )
        self.pp_menu1 = widgets.Combobox(
            placeholder="choose a pseudopotential",
            options=filtered_pp_list,
            description='pseudopotential 1:',
            disabled=False
        )
        self.pp_menu2 = widgets.Combobox(
            placeholder="choose a pseudopotential",
            options=filtered_pp_list,
            description='pseudopotential 2:',
            disabled=False
        )
        self.output = widgets.Output()
        def _on_button_clicked(self, change):
            with output:
                print("submitting sim2l run with formula" , self.compound.value)
                runSim2l()
                r = Run(simToolLocation,inputs)

        # display(c, s, button, output)
        def runSim2l():
            inputs['loglevel'].value = log.value
            inputs['walltime'].value = walltime.value
            inputs['numnodes'].value = numnodes.value
            inputs['ecutwfc'].value = ecutwfc.value
            inputs['ecutrho'].value = ecutrho.value
            inputs['pps'].value = [pp_menu1.value, pp_menu2.value]
            inputs['smearing'].value = smearing.value
            inputs['struct_dict'].value = struct_dict
