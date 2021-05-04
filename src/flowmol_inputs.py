#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import time
import os
import shutil

import wx
import wx.propgrid as wxpg
import wx.html as html
from wx.adv import BitmapComboBox
import wx.lib.dialogs

import numpy as np
from itertools import product

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

from SetupInputs import SetupInputs

# Code to read input file
import sys
sys.path.append("/home/es205/codes/python/SimWrapPy/")
import simwraplib as swl

sys.path.insert(0, "/home/es205/codes/python/pyDataView/")
import postproclib as ppl
import postproclib.visualiser as pplv

class CanvasPanel(wx.Panel):
    def __init__(self, parent, tmpdir="temp"):
        wx.Panel.__init__(self, parent)
        self.parent = parent
        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.tmpdir = tmpdir
        self.ft = True
        self.resultsdir = self.tmpdir + "/results/"
        self.ThreeD = True
        if (self.ThreeD):
            self.axes = self.figure.add_subplot(111, projection='3d', proj_type = 'ortho')
        else:
            self.axes = self.figure.add_subplot(111)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def axisEqual3D(self, ax):
        extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        sz = extents[:,1] - extents[:,0]
        centers = np.mean(extents, axis=1)
        maxsize = max(abs(sz))
        r = maxsize/2
        for ctr, dim in zip(centers, 'xyz'):
            getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    def draw(self, plottype):

        if not self.ft:
            #Store current angle and delete axis
            azim = self.axes.azim
            elev = self.axes.elev
            self.axes.cla()
            self.draw_grid(self.gridcheck.IsChecked())
            #self.parent.gridcheck.SetValue(False)

            #Delete any existing colorbar
            try:
                self.cb.remove()
                del self.cb
            except AttributeError:
                pass

        try:
            header = ppl.MDHeaderData(self.resultsdir)
        except FileNotFoundError:
            return
        N = int(header.globalnp)
        r = np.fromfile(self.resultsdir + "/initial_dump_r",  dtype=np.float).reshape(N,3)
        #Limit plots to be fast
        skip = int(N/1000.)
        r = r[::skip,:]
        try:
            if plottype == "tags":
                tag = np.fromfile(self.resultsdir + "/initial_dump_tag",  dtype=np.int32)
                print("tags include ", np.unique(tag))
                if (self.ThreeD):
                    scatter = self.axes.scatter(r[:,0], r[:,1], r[:,2], ".", c=tag[::skip], cmap=cm.RdYlBu_r)
                else:
                    scatter = self.axes.scatter(r[:,0], r[:,1], c=tag[::skip], cmap=cm.RdYlBu_r)
                #Generate tag labels from data
                elems = scatter.legend_elements()
                #tagDict = {"free": 0, "fixed": 1, "fixed_slide": 2, "teth": 3, "thermo": 4, 
                #            "teth_thermo": 5, "teth_slide": 6, "teth_thermo_slide": 7}
                tagList = ["free", "fixed", "fixed_slide", "teth", "thermo", 
                            "teth_thermo", "teth_slide", "teth_thermo_slide"]
                tags = [int(t.replace('$\\mathdefault{','').replace('}$','')) for t in elems[1]]
                newelems = (elems[0], [tagList[t] for t in tags])
                legend = self.axes.legend(*newelems)
                self.axes.add_artist(legend)
            elif plottype == "moltype":
                moltype = np.fromfile(self.resultsdir + "/initial_dump_moltype",  dtype=np.int32)
                if (self.ThreeD):
                    self.axes.scatter(r[:,0], r[:,1], r[:,2], ".", c=moltype[::skip], cmap=cm.RdYlBu_r)
                else:
                    self.axes.scatter(r[:,0], r[:,1], c=moltype[::skip], cmap=cm.RdYlBu_r)
            elif plottype == "v":
                v = np.fromfile(self.resultsdir + "/initial_dump_v",  dtype=np.float).reshape(N,3)
                vmag = np.sqrt(v[::skip,0]**2 + v[::skip,1]**2 + v[::skip,2]**2)
                if (self.ThreeD):
                    cs = self.axes.scatter(r[:,0], r[:,1], r[:,2], ".", c=vmag, cmap=cm.RdYlBu_r)
                else:
                    cs = self.axes.scatter(r[:,0], r[:,1], c=vmag, cmap=cm.RdYlBu_r)
                self.cb = self.figure.colorbar(cs, ax=self.axes)
            elif plottype == "None":
                if (self.ThreeD):
                    self.axes.scatter(r[:,0], r[:,1], r[:,2], ".")
                else:
                    self.axes.scatter(r[:,0], r[:,1])
        except (FileNotFoundError, ValueError) as e:
            if (self.ThreeD):
                self.axes.scatter(r[:,0], r[:,1], r[:,2], ".")
            else:
                self.axes.scatter(r[:,0], r[:,1])

        if (self.ThreeD):
            #try:
            #    self.axes.set_box_aspect((np.ptp(r[:,0]), np.ptp(r[:,1]), np.ptp(r[:,2])))
            #except AttributeError:
            #    self.axisEqual3D(self.axes)

            if self.ft:
                self.axes.view_init(90, -90)
                self.ft=False
            else:
                self.axes.view_init(elev, azim)

        self.canvas.draw()


    def draw_grid(self, draw):

        if draw:
            try:
                header = ppl.MDHeaderData(self.resultsdir)
            except FileNotFoundError:
                return

            nx = int(header.gnbins1)
            ny = int(header.gnbins2)
            nz = int(header.gnbins3)

            Lx = float(header.globaldomain1)
            Ly = float(header.globaldomain2)
            Lz = float(header.globaldomain3)

            x = np.linspace(-Lx/2.,Lx/2., nx+1)
            y = np.linspace(-Ly/2, Ly/2., ny+1)
            z = np.linspace(-Lz/2, Lz/2., nz+1)
            #X, Y, Z = np.meshgrid(x,y,z)

            #Origin at zero
            ox = oy = oz = 0.

            x1, z1 = np.meshgrid(x, z)
            y11 = np.ones_like(x1)*(oy-Ly/2)
            y12 = np.ones_like(x1)*(oy+Ly/2)
            x2, y2 = np.meshgrid(x, y)
            z21 = np.ones_like(x2)*(oz-Lz/2)
            z22 = np.ones_like(x2)*(oz+Lz/2)
            y3, z3 = np.meshgrid(y, z)
            x31 = np.ones_like(y3)*(ox-Lx/2)
            x32 = np.ones_like(y3)*(ox+Lx/2)

            self.grid = []
            # outside surface
            self.grid.append(self.axes.plot_wireframe(x1, y11, z1, color='k', rstride=1, cstride=1, alpha=0.6))
            # inside surface
            self.grid.append(self.axes.plot_wireframe(x1, y12, z1, color='k', rstride=1, cstride=1, alpha=0.6))
            # bottom surface
            self.grid.append(self.axes.plot_wireframe(x2, y2, z21, color='k', rstride=1, cstride=1, alpha=0.6))
            # upper surface
            self.grid.append(self.axes.plot_wireframe(x2, y2, z22, color='k', rstride=1, cstride=1, alpha=0.6))
            # left surface
            self.grid.append(self.axes.plot_wireframe(x31, y3, z3, color='k', rstride=1, cstride=1, alpha=0.6))
            # right surface
            self.grid.append(self.axes.plot_wireframe(x32, y3, z3, color='k', rstride=1, cstride=1, alpha=0.6))

            self.canvas.draw()


        else:
            try:
                for g in self.grid:
                    g.remove()
                del self.grid
                self.canvas.draw()
            except AttributeError:
                pass
    

class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyFrame.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.width = 800
        self.height = 600
        self.SetSize((self.width, self.height))
        self.SetTitle("Flowmol Input")

        #Top menu
        self.InitUI()

        #Setup notebook pages
        self.notebook_1 = wx.Notebook(self, wx.ID_ANY)
        self.notebook_1_pane_1 = wx.Panel(self.notebook_1, wx.ID_ANY)
        self.notebook_1.AddPage(self.notebook_1_pane_1, "Setup")


        #Top sizer
        sizer_top = wx.BoxSizer(wx.HORIZONTAL)

        #Split display and properties window
        self.window_LR = wx.SplitterWindow(self.notebook_1_pane_1, wx.ID_ANY)
        self.window_LR.SetMinimumPaneSize(20)
        sizer_top.Add(self.window_LR, 1, wx.EXPAND, 0)

        #Left window
        self.create_input_panel()

        #Create matploltib panel
        self.create_matplot_panel()

        #Add all sizers to windows
        self.window_LR.SplitVertically(self.window_left, self.window_right)
        self.notebook_1_pane_1.SetSizer(sizer_top)
        self.Layout()

        #Add pyDataView to second tab
        self.notebook_1_pane_2 = pplv.MainPanel(self.notebook_1, "./", 
                                                catch_noresults=False)
        self.notebook_1.AddPage(self.notebook_1_pane_2, "Results")
    

        #Close handle
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def InitUI(self):    

        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        open = fileMenu.Append(wx.ID_OPEN, '&Open\tCtrl+O')
        save = fileMenu.Append(wx.ID_SAVE, '&Save\tCtrl+S')
        quit = fileMenu.Append(wx.ID_EXIT, '&Quit\tCtrl+Q')

        menubar.Append(fileMenu, '&File')
        self.SetMenuBar(menubar)
        
        self.Bind(wx.EVT_MENU, self.OnOpen, open)
        self.Bind(wx.EVT_MENU, self.OnSaveAs, save)
        self.Bind(wx.EVT_MENU, self.OnQuit, quit)


    def create_input_panel(self):

        self.window_left = wx.Panel(self.window_LR, wx.ID_ANY)
        grid_sizer_1 = wx.BoxSizer(wx.VERTICAL)

        search_run_sizer = wx.BoxSizer(wx.HORIZONTAL)
        grid_sizer_1.Add(search_run_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)

        #########################################
        # Search control for top level keywords #
        #########################################
        self.searchctrl = BitmapComboBox(self.window_left,
                              wx.ID_ANY,  
                              choices=[],
                              style= wx.TE_PROCESS_ENTER)# | wx.CB_SORT)
        self.searchctrl.SetMinSize((200, 30))
        search_run_sizer.Add(self.searchctrl, 1, 0, 0)

        self.searchctrl.Bind(wx.EVT_TEXT_ENTER, self.on_search)
        self.searchctrl.Bind(wx.EVT_COMBOBOX, self.on_search)
        self.searchctrl.Bind(wx.EVT_TEXT, self.on_search)
        self.searchctrl.Bind(wx.EVT_COMBOBOX_CLOSEUP, self.on_search)
        self.searchctrl.SetFocus()


        #########################################
        #             Run Button                #
        #########################################
        self.runbtn = wx.Button(self.window_left, wx.ID_ANY, "Run")
        self.runbtn.SetMinSize((30, 30))
        search_run_sizer.Add(self.runbtn, 1, 0, 0)
        self.runbtn.Bind(wx.EVT_BUTTON, self.run_btn)

        #########################################
        #           Setup Button                #
        #########################################
        self.button_2 = wx.Button(self.window_left, wx.ID_ANY, "Setup")
        self.button_2.SetMinSize((50, 30))
        search_run_sizer.Add(self.button_2, 1, 0, 0)
        self.button_2.Bind(wx.EVT_BUTTON, self.run_setup)

        #Setup adjustable window between properties and help
        self.window_help_props = wx.SplitterWindow(self.window_left, wx.ID_ANY)
        self.window_help_props.SetMinSize((-1, 600))
        self.window_help_props.SetMinimumPaneSize(20)
        grid_sizer_1.Add(self.window_help_props, 0, wx.EXPAND, 0)

        #########################################
        #             Help panel                #
        #########################################
        self.window_help = wx.Panel(self.window_help_props, wx.ID_ANY)
        sizer_help = wx.BoxSizer(wx.HORIZONTAL)
        self.helptxt = wx.TextCtrl(self.window_help, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(self.width/3, self.height/4))
        self.helptxt.SetValue("Help Text \n\n\n\n\n\n\n")
        sizer_help.Add(self.helptxt, 1, wx.EXPAND, 0)

        #########################################
        #           Property Grid               #
        #########################################
        self.window_props = wx.Panel(self.window_help_props, wx.ID_ANY)
        sizer_props = wx.BoxSizer(wx.HORIZONTAL)
        self.propgrid = wxpg.PropertyGridManager(self.window_props, wx.ID_ANY,
                                                 style=wxpg.PG_SPLITTER_AUTO_CENTER)
        sizer_props.Add(self.propgrid, 1, wx.EXPAND, 0)
        self.propgrid.Bind(wxpg.EVT_PG_CHANGED, self.change_propgrid)

        self.window_props.SetSizer(sizer_props)
        self.window_help.SetSizer(sizer_help)
        self.window_left.SetSizer(grid_sizer_1)
        self.window_help_props.SplitHorizontally(self.window_help, self.window_props)


        #Not sure if any or all of these are needed to create the panel 
        self.Layout()
        self.Refresh()
        self.Update()


    def create_matplot_panel(self):

        #Right window
        self.window_right = wx.Panel(self.window_LR, wx.ID_ANY)
        sizer_Matplot_panel = wx.BoxSizer(wx.VERTICAL)

        self.plotpanel = CanvasPanel(self.window_right)
        self.tmpdir = self.plotpanel.tmpdir
        sizer_Matplot_panel.Add(self.plotpanel, 1, wx.EXPAND, 0)

        sizer_plotcontrol_panel = wx.BoxSizer(wx.HORIZONTAL)


        #Create radiobox
        self.plottypes = ["v", "tags", "moltype", "None"]
        self.radio_box_1 = wx.RadioBox(self.window_right, wx.ID_ANY, "", 
                                       choices=self.plottypes, majorDimension=4, 
                                       style=wx.RA_SPECIFY_COLS)
        self.radio_box_1.SetSelection(0)
        self.plotype = "v"
        self.radio_box_1.Bind(wx.EVT_RADIOBOX, self.onRadioBox)
        sizer_plotcontrol_panel.Add(self.radio_box_1, 0, wx.EXPAND, 0)

        #Add grid tickbox
        self.gridcheck = wx.CheckBox(self.window_right, wx.ID_ANY, "Grid")
        sizer_plotcontrol_panel.Add(self.gridcheck, 0, wx.EXPAND, 0)
        self.gridcheck.Bind(wx.EVT_CHECKBOX , self.ongridtick)
        self.plotpanel.gridcheck = self.gridcheck

        #Add a loading gauge (FOR SOME REASON THIS DOESN'T WORK IN UBUNTU 18.04
        #                     JUST SITS IN THE BOTTOM LEFT)
        #self.gauge = wx.Gauge(self, range=100, style=wx.GA_HORIZONTAL)
        #self.gauge.SetValue(100)
        #sizer_Matplot_panel.Add(self.gauge)

        sizer_Matplot_panel.Add(sizer_plotcontrol_panel, 0, wx.EXPAND, 0)
        self.window_right.SetSizer(sizer_Matplot_panel)


    def populate_searchctrl(self):

        #Save choices
        self.choices = list(self.InputsDict.keys())

        #Images for optional and compulsory
        req = wx.Image('Required.png').ConvertToBitmap()
        ops = wx.Image('Optional.png').ConvertToBitmap()

        #Order by if they have values
        self.ChangeDict = {}
        self.set = []; self.unset = []
        for k in self.choices:
            try:
                set = self.InputFile.read_inputs(k)
            except AttributeError:
                set = []
            if not set:
                self.unset.append(k)
            else:
                self.set.append(k)
        for s in [self.set, self.unset]:
            for k in s:
                if (self.InputsDict[k]["Optional"]):
                    self.searchctrl.Append(k, req)
                else:
                    self.searchctrl.Append(k, ops)

        out = self.searchctrl.AutoComplete(list(self.InputsDict.keys()))
    

    def Create_InputsDict(self):

        """
            Setup input of the form

            FORMAT
            InputsDict = {"KEYWORD": 
                                   {"HELP":"Text to describe input", 
                                    "optional":True/False, 
                                    "vars":
                                           {"inputname":input value}
                                    }     
                          }
        """

        #Get src directory from local directory
        self.srcdir = os.path.dirname(os.path.abspath(__file__)) 
        fsetup = self.srcdir + "/setup_read_input.f90"
        if os.path.isfile(fsetup):
            #Code to get inputs from setup set parameters
            FlowmolInputs = SetupInputs(fsetup='./setup_read_input.f90')
            FlowmolInputDict = FlowmolInputs.get_items()
        else:
            raise FileNotFoundError("Error - setup_read_input.f90 not found " +
                                    "in current directory, edit location in code")

        InputsDict = {}
        for key, item in FlowmolInputDict.items():
            InputsDict[key] = {}
            #Get all variables
            helpstr = FlowmolInputs.get_helpstring(key)
            InputsDict[key]["HELP"] = helpstr
            optional = FlowmolInputs.get_optional(key)
            InputsDict[key]["Optional"] = optional
            EnumProperties = FlowmolInputs.variables_from_string(helpstr)

            #Get values from input file
            try:
                params = self.InputFile.read_inputs(key)
            except AttributeError:
                params = []

            if (len(EnumProperties) != 0):
                EnumDict = {}
                if (len(item) == len(params)):
                    EnumDict = {}
                    for i in range(len(item)):
                        try:
                            setval = int(params[i])
                        except ValueError:
                            try:
                                setval = float(params[i])
                            except ValueError:
                                setval = 0
                        EnumDict[item[i]] = {"names":[], "numbers":[], "set":setval}
                else:
                    for i in range(len(item)):
                        if (i < len(params)):
                            EnumDict[item[i]] = {"names":[], "numbers":[], "set":params[i]}
                        else:
                            EnumDict[item[i]] = {"names":[], "numbers":[], "set":0}
                for e in EnumProperties:
                    itemnum, val, name = e
                    try:
                        EnumDict[item[itemnum-1]]["names"].append(name)
                        EnumDict[item[itemnum-1]]["numbers"].append(val)
                    except IndexError:
                        print("IndexError", e)
                InputsDict[key]["vars"] = EnumDict 
            else:
                InputsDict[key]["vars"] = {}
                for i in range(len(item)):
                    if (i < len(params)):
                        InputsDict[key]["vars"][item[i]] = params[i]
                    else:
                        InputsDict[key]["vars"][item[i]] = "0"

        return InputsDict


    def on_search(self, event):

        # Get key and build property grid
        value = self.searchctrl.GetValue()

        #Look for value in Dict
        try:
            VarsDict = self.InputsDict[value]
        except KeyError:
            #print("Keyword ", value, "not found in ", self.InputsDict.keys())
            return 

        try:
            cleantext = VarsDict["HELP"].replace("#","").replace("!","").replace("--","").replace("\t","")
            self.helptxt.SetValue(cleantext)
            vars = VarsDict["vars"]
        except KeyError:
            print("InputDict entries should contain HELP and vars entries")
            raise

        self.propgrid.Clear()
        self.propgrid.AddPage( "Page" )

        if isinstance(vars, dict):
            for k, var in vars.items():
                try:
                    #Should be in form of nested dictonaries
                    if isinstance(var, dict):
                        #A list of numbers corresponding to flags
                        if isinstance(var["numbers"][0], int):
                            self.propgrid.Append( wxpg.EnumProperty(k, k, var["names"], var["numbers"], int(var["set"])) )
                        #Integers such as number of unit cells
                        elif ("int" in var["numbers"][0]):
                            self.propgrid.Append( wxpg.IntProperty(k, value=int(var["set"])) )
                        #Floats such as density of system
                        elif ("float" in var["numbers"][0]):
                            self.propgrid.Append( wxpg.FloatProperty(k, value=float(var["set"])) )
                        #Or a case with string based keywords 
                        #(so store a mapping to use enum list) 
                        elif isinstance(var["numbers"][0], str):
                            mapping = list(range(len(var["numbers"])))
                            print(k, var["numbers"], mapping, var["set"])

                            self.propgrid.Append( wxpg.EnumProperty(k, k, var["numbers"], mapping,   int(var["set"])) )
                        else:
                            raise KeyError("Dictonary format not known")
                        #print(n)
                    elif ".true." in var:
                        self.propgrid.Append( wxpg.BoolProperty(k, value=True) )
                    elif ".false." in var:
                        self.propgrid.Append( wxpg.BoolProperty(k, value=False) )
                    elif "." in var:
                        self.propgrid.Append( wxpg.FloatProperty(k, value=float(var)) )
                    else: 
                        self.propgrid.Append( wxpg.IntProperty(k, value=int(var)) )
                except ValueError:
                    self.propgrid.Append( wxpg.StringProperty(k,value=str(var)) )
                    raise

    def change_propgrid(self, event):
        key = self.searchctrl.GetValue()
        prop = event.GetPropertyName()
        val = event.GetPropertyValue()
        PropObj = event.GetProperty()
        #print(event.GetSelection(), event.GetProperty())
        if (isinstance(self.InputsDict[key]["vars"][prop], dict)):
            self.InputsDict[key]["vars"][prop]["set"] = val
            #print("Set value = ", val, PropObj.GetChoiceSelection(), PropObj.ValueToString(val))
        else:
            self.InputsDict[key]["vars"][prop] = str(val)

        keys = self.InputsDict[key]["vars"].keys() 
        try:
            if isinstance(self.ChangeDict[key], list):
                changes = self.ChangeDict[key] 
            else:
                changes = [None]*len(keys)
        except KeyError:
            self.ChangeDict[key] = [None]*len(keys)
            changes = [None]*len(keys)

        for i, k in enumerate(keys):
            if prop == k:
                changes[i] = val
                #print("change_propgrid vars", i, k, val)

        self.ChangeDict[key] = changes

        print("Changes = ", key, changes, self.ChangeDict)

        return

    def run_btn(self, event): 

        #btn = event.GetEventObject().GetLabel() 

        # otherwise ask the user what new file to open
        with wx.FileDialog(self, "Choose output directory",
                            defaultDir='./',
                            style=wx.FD_OPEN) as folderDiag:

            if folderDiag.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

            rundir = folderDiag.GetPath()
            #print("Label of button = ", btn, rundir)

        run = swl.MDRun(self.srcdir, self.srcdir, rundir,
                  "parallel_md.exe", 
                  self.inputfilename.split("/")[-1], "setup.out",
                  inputchanges=changes[0], finishargs = {},
                  dryrun=False)


        # Run the study
        runlist = [run]
        threadlist =[runlist]
        with wx.BusyInfo("Working, please wait", self):
            study = swl.Study(threadlist, ncpus)

    def run_setup(self, event): 

        btn = event.GetEventObject().GetLabel() 

        # Number of threads and runs per thread
        ncpus = 1
        maxlicenses = ncpus

        #Check input file has been loaded
        try:
            print("Running setup run", self.inputfilename)
        except AttributeError:
            msgbx = wx.MessageDialog(self, "Setup run with no input file specified",
                                    style=wx.OK|wx.ICON_ERROR)
            msgbx.ShowModal()
            msgbx.Destroy()
            return

        #Create list of changes
        try:
            #Concat dictonaries (incase item already in)
            changes = dict({'VMD_OUTFLAG': [5]}, **self.ChangeDict)
            #Changes needs to be a list of dictonaries
            #changes = [{k:v} for k,v in ChangeDict.items()]

        except IndexError:
            changes = {'VMD_OUTFLAG': [5]}


        print("Running setup run", self.inputfilename, " with changes", changes)

        #Remove and recreate temp directory
        try:
            shutil.rmtree(self.tmpdir)
        except FileNotFoundError:
            pass
        print("Making ", self.tmpdir)
        os.mkdir(self.tmpdir)
        os.mkdir(self.tmpdir+"/results")

        #Use current source code location
        #srcdir = "/home/es205/codes/flowmol/src/"
        inputfile = self.inputfilename
        if self.srcdir in inputfile:
            pathtoinput = inputfile.replace(self.srcdir,'')
            inputfile = inputfile.split("/")[-1]
            basedir = self.inputfilename.replace(inputfile,"")

        #print(inputfile)

        self.run = swl.MDRun(self.srcdir, basedir, self.srcdir + "/" + self.tmpdir,
                  "parallel_md.exe", 
                  inputfile, "setup.out",
                  inputchanges=changes, finishargs = {},
                  dryrun=False, minimalcopy=True)                

        self.run.setup()
        self.run.execute(print_output=False, out_to_file=False, blocking=False)
        self.Bind(wx.EVT_IDLE, self.OnIdle)


        #self.progress = wx.ProgressDialog("Running Setup", "please wait", 
        #                                   parent=self, 
        #                                   style=wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)

        # #Line by line step through
        # i = 0
        # for stdout_line in iter(run.proc.stdout.readline, ""):
            # lastline = stdout_line.replace("\n","")
            # errormsg = run.proc.stderr.read()
            # #time.sleep(0.1)
            # contnue, skp = self.progress.Update(i, lastline)
            # if contnue is False:
                # print("Cancel Pressed")
                # self.progress.Destroy()
                # self.proc.stdout.close()
                # run.proc.kill()
                # break

            # print(lastline, run.proc.returncode, run.proc.stderr.read())
            # if "Time taken" in lastline:
                # self.progress.Destroy()
                # run.proc.kill()
                # break

        # # except run.CalledProcessError:

            # if (run.proc.returncode or errormsg != ""):
                # #print remaining stdout
                # stdout = run.proc.stdout.read()
                # print(stdout)
                # self.progress.Destroy()
                # run.proc.kill()
                # #Clean stderr
                # print("stderr = ", errormsg)
                # errormsg_box = errormsg.split("\n")[0]
                # msgbx = wx.MessageDialog(self, errormsg_box 
                                         # +"\n Look at terminal for more error information",
                                         # style=wx.OK|wx.ICON_ERROR)
                # msgbx.ShowModal()
                # msgbx.Destroy()
                # return

        # #Redraw the figure with latest setup
        # self.plotpanel.draw(self.plotype)


    def OnIdle(self, event):

        run = self.run
        stdout = run.proc.stdout.read()
        errormsg = run.proc.stderr.read()
        if (stdout != ""):
            print(stdout)
        self.helptxt.SetValue(stdout)
        #self.gauge.SetValue(0)
        if (run.proc.returncode or errormsg != ""):
            run.proc.kill()
            #Clean stderr
            print("stderr = ", errormsg, "stdout =", stdout)
            errormsg_box = errormsg.split("\n")[0]
            msgbx = wx.MessageDialog(self, errormsg_box 
                                     +"\n Look at terminal for more error information",
                                     style=wx.OK|wx.ICON_ERROR)
            msgbx.ShowModal()
            msgbx.Destroy()
            self.Unbind(wx.EVT_IDLE)
            return

        if "Time taken" in stdout:
            run.proc.kill()
        else:
            event.RequestMore()
            return

        self.plotpanel.draw(self.plotype)
        self.Unbind(wx.EVT_IDLE)
        

    def OnOpen(self, event):

        # otherwise ask the user what new file to open
        with wx.FileDialog(self, "Open Flowmol input file", 
                           wildcard="in files (*.in)|*.in",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

            # Proceed loading the file chosen by the user
            fdir = fileDialog.GetPath()
            self.inputfilename = fdir
            try:
                with open(fdir, 'r') as file:
                    #Destroy current panel if existing
                    #try:
                    #    self.panel_1.destroy()
                    #except AttributeError:
                    #    print("Creating new panel")
                    #Load new input file
                    self.InputFile = swl.KeywordInputMod(fdir)
                    #Create dictonary of inputs from flowmol setup_read_inputs
                    self.InputsDict = self.Create_InputsDict()
                    #Create all details
                    self.populate_searchctrl()

            except IOError:
                wx.LogError("Cannot open file '%s'." % fdir)

    def OnSaveAs(self, event):

        with wx.FileDialog(self, "Save Flowmol input file", wildcard="in files (*.in)|*.in",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:

            # the user changed their mind
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return    

            # save the current contents in the file
            fdir = fileDialog.GetPath()
            try:
                with open(fdir, 'w') as file:
                    self.doSaveData(file)
            except IOError:
                wx.LogError("Cannot save current data in file '%s'." % fdir)

    def onRadioBox(self, e): 
        print(self.radio_box_1.GetStringSelection(),' is clicked from Radio Box')
        self.plotype = self.radio_box_1.GetStringSelection()
        self.plotpanel.draw(self.plotype)

    def ongridtick(self, event):
        cb = event.GetEventObject() 
        self.plotpanel.draw_grid(cb.GetValue())

    def OnClose(self, event):
        try:
            shutil.rmtree(self.tmpdir)
        except FileNotFoundError:
            pass        
        self.Destroy()

    def OnQuit(self, e):
        try:
            shutil.rmtree(self.tmpdir)
        except FileNotFoundError:
            pass        
        self.Close()


    def About(self, event):
        from platform import platform
        myos = platform()
        aboutInfo = wx.AboutDialogInfo()
        aboutInfo.SetName("Flowmol")
        aboutInfo.SetVersion("1.0")
        aboutInfo.SetDescription("Flowmol," \
            " Molecular Fluid Dynamics simulation toool: "+myos)
        aboutInfo.SetCopyright("(C) Edward Smith-2021")
        aboutInfo.SetLicense("https://www.gnu.org/licenses/gpl-3.0.html")
        aboutInfo.AddDeveloper("Edward Smith")
        aboutInfo.AddDocWriter("Edward Smith")
        aboutInfo.SetWebSite('https://www.edwardsmith.co.uk')
        wx.AboutBox(aboutInfo)


class MyApp(wx.App):
    def OnInit(self):
        self.frame = MyFrame(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True

#import wx.lib.inspection

if __name__ == "__main__":
    app = MyApp(0)
    #Use wxPython debugging tool
    #wx.lib.inspection.InspectionTool().Show()
    app.MainLoop()



    # def create_panel(self):
        # self.panel_1 = wx.Panel(self, wx.ID_ANY)
        # self.sizer_1 = wx.BoxSizer(wx.VERTICAL)

        # self.sizer_2 = wx.FlexGridSizer(1, 3, 0, 0)
        # self.sizer_1.Add(self.sizer_2, 1, wx.ALL | wx.EXPAND, 0)


        # #########################################
        # # Search control for top level keywords #
        # #########################################
        # self.choices = list(self.InputsDict.keys())
        # self.searchctrl = BitmapComboBox(self.panel_1,
                              # size=wx.DefaultSize,  
                              # choices=[],
                              # style= wx.TE_PROCESS_ENTER)# | wx.CB_SORT)

        # #Images for optional and compulsory
        # req = wx.Image('Required.png').ConvertToBitmap()
        # ops = wx.Image('Optional.png').ConvertToBitmap()

        # #Order by if they have values
        # self.ChangeDict = {}
        # self.set = []; self.unset = []
        # for k in self.choices:
            # try:
                # set = self.InputFile.read_inputs(k)
            # except AttributeError:
                # set = []
            # if not set:
                # self.unset.append(k)
            # else:
                # self.set.append(k)
        # for s in [self.set, self.unset]:
            # for k in s:
                # if (self.InputsDict[k]["Optional"]):
                    # self.searchctrl.Append(k, req)
                # else:
                    # self.searchctrl.Append(k, ops)

        # self.sizer_2.Add(self.searchctrl, 0, 0, 0)
        # out = self.searchctrl.AutoComplete(list(self.InputsDict.keys()))

        # self.searchctrl.Bind(wx.EVT_TEXT_ENTER, self.on_search)
        # self.searchctrl.Bind(wx.EVT_COMBOBOX, self.on_search)
        # self.searchctrl.Bind(wx.EVT_TEXT, self.on_search)
        # self.searchctrl.Bind(wx.EVT_COMBOBOX_CLOSEUP, self.on_search)

        # self.searchctrl.SetFocus()

        # #########################################
        # #             Run Button                #
        # #########################################
        # self.runbtn = wx.Button(self.panel_1, wx.ID_ANY, "Run")
        # self.sizer_2.Add(self.runbtn, 0, 0, 0)
        # self.runbtn.Bind(wx.EVT_BUTTON, self.run_btn)

        # #########################################
        # #             Help panel                #
        # #########################################
        # self.helptxt = wx.TextCtrl(self.panel_1, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(800, 400))
        # self.helptxt.SetValue("Help Text \n\n\n\n\n\n\n")
        # self.sizer_1.Add(self.helptxt, 0, wx.EXPAND, 0)

        # #########################################
        # #           Property Grid               #
        # #########################################
        # self.propgrid = wxpg.PropertyGridManager(self.panel_1, wx.ID_ANY,
                        # style=wxpg.PG_SPLITTER_AUTO_CENTER)
        # self.sizer_1.Add(self.propgrid, 1, wx.EXPAND, 0)

        # self.propgrid.Bind(wxpg.EVT_PG_CHANGED, self.change_propgrid)
        # self.panel_1.SetSizer(self.sizer_1)

        # self.Layout()
        # self.Refresh()
        # self.Update()
        # self.panel_1.Layout()
        # self.panel_1.Refresh()
        # self.panel_1.Update()
        # #self.Layout()
        # #self.Refresh()
        # #self.SetSize((self.width, self.height))
        # #self.Update()

