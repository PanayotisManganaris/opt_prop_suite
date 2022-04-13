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
from IPython.display import display, Javascript, clear_output

class Authenticate():
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
    def __init__(self):
        self.checkmp = os.path.isfile(os.path.expanduser("~/.mpkey.txt"))

    def mpkey(self):
        """
        ask user for authentication if they do not already have mpkey file
        """
        if self.checkmp:
            pass
        else:
            self._mpkey()

    @staticmethod
    def _mpkey():
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

class MPQuery():
    """ might want to try implementing OPTIMATE """
    def __init__(self, content=[], **kwargs):
        """current stopgap:
        if OPTIMATE is hard to implement, just make case-by-case query
        object with the following API

        authentication is pulled from the user's filesystem to offload security onto them...
        frame returns a dataframe of requested values
        """
        self.authfile = kwargs.get('authfile', "~/.mpkey.txt")
        self.query = kwargs.get('query', ({ "crystal_system": "cubic"},
                                          ["task_id","pretty_formula","formula",
                                           "elements","e_above_hull", "spacegroup.number",
                                           "band_gap", "crystal_system"]))
        self.rester = self._make_rester()
        self.frame = pd.DataFrame(self.get_content(content)) #mutable

    def _make_rester(self):
        with open(os.path.expanduser(self.authfile), "r+") as file:
            apikey = file.readline()
        return MPRester(apikey) #investigate other apis, plug and play?... make a function probably. OPTIMATE?

    def get_content(self, content):
        """mutate list of MP results"""
        if not content: #should persist through every instantiation of QueryObject
            content += self.rester.query(*self.query)
        return content

    def get_structure(self, mpid):
        return self.rester.get_structure_by_material_id(mpid,
                                                        final=False,
                                                        conventional_unit_cell=False)

class QuerySuite():
    """widgets for making and reviewing queries"""
    def __init__(self, QueryObj):
        q = QueryObj
        self.menu = qgrid.show_grid(q.frame)
        
        self.plotbutton = widgets.Button(
            description='plot',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            tooltip='plot'
        )
        def click_to_plot_structure(b):
            selected_row = self.menu.get_changed_df().index[self.menu.get_selected_rows()][0]
            mpid_selected = q.frame.at[selected_row,'task_id']
            struct = q.get_structure(mpid_selected)
            self.mpid_plot(struct)

        self.plotbutton.on_click(click_to_plot_structure)
    
    def get_qgrid(self):
        return self.menu

    def get_plot(self):
        return self.plotbutton

    def _display_qgrid(self):
        """debugging qgrid"""
        display(self.menu)


    @staticmethod
    def mpid_plot(struct):
        # import POSCAR file
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

class InputSuite():
    def __init__(self, inputs, struct_dict): #users_pp_file_list

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
