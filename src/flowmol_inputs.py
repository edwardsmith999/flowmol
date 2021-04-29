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
sys.path.append("/home/es205/codes/SimWrapPy/")
import simwraplib as swl

sys.path.append("/home/es205/codes/pyDataView/")
import postproclib as ppl


def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


class CanvasPanel(wx.Panel):
    def __init__(self, parent, tmpdir="temp"):
        wx.Panel.__init__(self, parent)
        self.parent = parent
        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.tmpdir = tmpdir

        self.axes = self.figure.add_subplot(111, projection='3d', proj_type = 'ortho')
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()


    def draw(self):

        self.axes.cla()
        resultsdir = self.tmpdir + "/results/"
        header = ppl.MDHeaderData(resultsdir)
        N = int(header.globalnp)
        r = np.fromfile(resultsdir + "/initial_dump_r",  dtype=np.float).reshape(N,3)
        try:
            tag = np.fromfile(resultsdir + "/initial_dump_tag",  dtype=np.int)
            c = cm.RdYlBu_r(tag/tag.max())
        except FileNotFoundError:
            c = "b"
        except ValueError:
            #print("Failed to load tags", tag.size, c.size)
            c = "k"

        self.axes.scatter(r[:,0], r[:,1], r[:,2], c=c)
        try:
            self.axes.set_box_aspect((np.ptp(r[:,0]), np.ptp(r[:,1]), np.ptp(r[:,2])))
        except AttributeError:
            axisEqual3D(self.axes)  
        self.axes.view_init(90, 90)
        #size = tuple(self.parent.GetClientSize())
        #self.figure.set_size_inches(float(size[0])/self.figure.get_dpi(),
        #                            float(size[1])/self.figure.get_dpi())
        #self.toolbar = NavigationToolbar2Wx(self.canvas)
        #self.toolbar.Realize()
        self.canvas.draw()

        #self.Fit()

class MyFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyFrame.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.width = 800
        self.height = 800
        self.SetSize((self.width, self.height))
        self.SetTitle("Flowmol Input")

        #Top menu
        self.InitUI()

        #Top sizer
        sizer_top = wx.BoxSizer(wx.HORIZONTAL)

        #Split display and properties window
        self.window_LR = wx.SplitterWindow(self, wx.ID_ANY)
        self.window_LR.SetMinimumPaneSize(20)
        sizer_top.Add(self.window_LR, 1, wx.EXPAND, 0)

        #Left window
        self.create_input_panel()

        #Right window
        self.window_right = wx.Panel(self.window_LR, wx.ID_ANY)
        sizer_2 = wx.BoxSizer(wx.VERTICAL)

        self.plotpanel = CanvasPanel(self.window_right)
        self.tmpdir = self.plotpanel.tmpdir
        #self.plotpanel.draw()
        #wx.Panel(self.window_right, wx.ID_ANY)
        sizer_2.Add(self.plotpanel, 1, wx.EXPAND, 0)

        self.radio_box_1 = wx.RadioBox(self.window_right, wx.ID_ANY, "", 
                                       choices=["tags", "moltype", "v"], majorDimension=1, style=wx.RA_SPECIFY_COLS)
        self.radio_box_1.SetSelection(0)
        sizer_2.Add(self.radio_box_1, 0, wx.EXPAND, 0)

        self.window_right.SetSizer(sizer_2)
        self.window_help_props.SplitHorizontally(self.window_help, self.window_props)
        self.window_LR.SplitVertically(self.window_left, self.window_right)
        self.SetSizer(sizer_top)
        self.Layout()

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

        #Not sure if any or all of these are needed to create the panel 
        self.Layout()
        self.Refresh()
        self.Update()

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
            try:
                InputsDict[key] = {}

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
                            print("param=item", i, item[i], params[i], setval)
                            EnumDict[item[i]] = {"names":[], "numbers":[], "set":setval}
                    else:
                        for i in range(len(item)):
                            if (i < len(params)):
                                EnumDict[item[i]] = {"names":[], "numbers":[], "set":params[i]}
                            else:
                                EnumDict[item[i]] = {"names":[], "numbers":[], "set":0}
                    # if (len(item) == len(params)):
                        # EnumDict = {}
                        # for i in range(len(item)):
                            # print(i, item[i], params[i])
                            # EnumDict[item[i]] = {"names":[], "numbers":[], "set":params[i]}
                    # else:
                        # EnumDict = {i:{"names":[], "numbers":[], "set":0} for i in item}
                    for e in EnumProperties:
                        itemnum, val, name = e
                        try:
                            EnumDict[item[itemnum-1]]["names"].append(name)
                            EnumDict[item[itemnum-1]]["numbers"].append(val)
                        except IndexError:
                            print("IndexError", e)
                        #EnumDict[item[itemnum-1]] = {"names":name, "numbers":val}
                    InputsDict[key]["vars"] = EnumDict 
                else:
                    InputsDict[key]["vars"] = {}
                    for i in range(len(item)):
                        if (i < len(params)):
                            #print(i, item[i], params[i])
                            InputsDict[key]["vars"][item[i]] = params[i]
                        else:
                            InputsDict[key]["vars"][item[i]] = "0"

                    #if (len(item) == len(params)):
                    #    InputsDict[key]["vars"] = {}
                    #    for i in range(len(item)):
                    #        InputsDict[key]["vars"][item[i]] = params[i]
                    #else:
                    #    InputsDict[key]["vars"] = {i:"0" for i in item}
            except KeyError:
                print("key ", key, " not found")
                raise

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
                changes = [[None]*len(keys)]
        except KeyError:
            self.ChangeDict[key] = [None]*len(keys)
            changes = [[None]*len(keys)]

        for i, k in enumerate(keys):
            if prop == k:
                changes[0][i] = val
                #print("change_propgrid vars", i, k, val)

        self.ChangeDict[key] = changes

        #print("Changes = ", self.ChangeDict)



        # try:
            # self.ChangeDict[key][prop] = val
        # except KeyError:
            # self.ChangeDict[key] = {prop:val}

        #print("after", key, prop, self.InputsDict[key]["vars"][prop], self.ChangeDict)
        return

    def run_btn(self, event): 

        btn = event.GetEventObject().GetLabel() 

        # otherwise ask the user what new file to open
        with wx.FileDialog(self, "Choose output directory",
                            defaultDir='./',
                            style=wx.FD_OPEN) as folderDiag:

            if folderDiag.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

            rundir = folderDiag.GetPath()
            print("Label of button = ", btn, rundir)

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
            changes = swl.InputDict({'VMD_OUTFLAG': [5]})+swl.InputDict(self.ChangeDict)
            print(self.ChangeDict, changes)

        except IndexError:
            changes = swl.InputDict({'VMD_OUTFLAG': [5]}).expand()


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

        run = swl.MDRun(self.srcdir, basedir, self.srcdir + self.tmpdir,
                  "parallel_md.exe", 
                  inputfile, "setup.out",
                  inputchanges=changes[0], finishargs = {},
                  dryrun=False, minimalcopy=True)                

        run.setup()
        run.execute(print_output=False, out_to_file=False, blocking=False)
        self.progress = wx.ProgressDialog("Running Setup", "please wait", 
                                           parent=self, 
                                           style=wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)

        #Line by line step through
        i = 0
        for stdout_line in iter(run.proc.stdout.readline, ""):
            lastline = stdout_line.replace("\n","")
            errormsg = run.proc.stderr.read()
            #time.sleep(0.1)
            contnue, skp = self.progress.Update(i, lastline)
            if contnue is False:
                print("Cancel Pressed")
                self.progress.Destroy()
                self.proc.stdout.close()
                run.proc.kill()
                break

            print(lastline, run.proc.returncode, run.proc.stderr.read())
            if "Time taken" in lastline:
                self.progress.Destroy()
                run.proc.kill()
                break

        # except run.CalledProcessError:

            if (run.proc.returncode or errormsg != ""):
                #print remaining stdout
                stdout = run.proc.stdout.read()
                print(stdout)
                self.progress.Destroy()
                run.proc.kill()
                #Clean stderr
                print("stderr = ", errormsg)
                errormsg_box = errormsg.split("\n")[0]
                msgbx = wx.MessageDialog(self, errormsg_box 
                                         +"\n Look at terminal for more error information",
                                         style=wx.OK|wx.ICON_ERROR)
                msgbx.ShowModal()
                msgbx.Destroy()
                return

        #Redraw the figure with latest setup
        self.plotpanel.draw()

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

    def OnClose(self, event):
        shutil.rmtree(self.tmpdir) 
        self.Destroy()



    def OnQuit(self, e):
        shutil.rmtree(self.tmpdir) 
        self.Close()


    def About(self, event):
        from platform import platform
        myos = platform()
        aboutInfo = wx.AboutDialogInfo()
        aboutInfo.SetName("My Application ")
        aboutInfo.SetVersion("1.0")
        aboutInfo.SetDescription("My Super App," \
            " That does amazing things\nRunning on: "+myos)
        aboutInfo.SetCopyright("(C) Joe Bloggs-2016")
        aboutInfo.SetLicense("https://www.gnu.org/licenses/gpl-2.0.html")
        aboutInfo.AddDeveloper("Joe Bloggs")
        aboutInfo.AddDocWriter("Joe Bloggs")
        aboutInfo.SetWebSite('https://www.JoeBlogs.com')
        wx.AboutBox(aboutInfo)


class MyApp(wx.App):
    def OnInit(self):
        self.frame = MyFrame(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True

# end of class MyApp

if __name__ == "__main__":
    app = MyApp(0)
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

