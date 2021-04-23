#!/usr/bin/env python3
# -*- coding: UTF-8 -*-




import wx
import wx.propgrid as wxpg
import wx.html as html
from wx.adv import BitmapComboBox

from SetupInputs import SetupInputs


# Code to read input file
import sys
sys.path.append("/home/es205/codes/SimWrapPy/")
import simwraplib as swl

fdir = "./MD_2.in"
InputFile = swl.KeywordInputMod(fdir)

#Code to 
FlowmolInputs = SetupInputs()
FlowmolInputDict = FlowmolInputs.get_items()


# InputsDict = {}
# for key, item in FlowmolInputDict.items():
    # try:
        # InputsDict[key] = {}
        # helpstr = FlowmolInputs.get_helpstring(key)
        # InputsDict[key]["HELP"] = helpstr
        # InputsDict[key]["vars"] = {i:"0" for i in item}
        # #InputsDict[key] = {"vars":{i:"0" for i in item}}
        # print(key, InputsDict[key])


    # except KeyError:
        # print("key ", key, " not found")



#
#InputsDict = {"INPUT":{"HELP":"THis is text to describe variable", "vars":{"name":"2"}}, 
#              "THING":{"HELP":"different help text", "vars":{"xcells":"1","ycells":"2","zcells":"3"}},
#"STR":{"HELP":"Example of a string with a really long help example to see if this fits in the box or needs to be wrapped", 
#                "vars":{"string":"Hello","logical":".true.","int":"2","float":"3.14159",
#                "List":{"names":['NVE', 'NVT', 'Tag Move system'], "numbers":[0,1,6]}}}}

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


class CanvasPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)

        self.axes = self.figure.add_subplot(111, projection='3d')
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def draw(self):
        N = 6831
        r = np.fromfile("results/initial_dump_r",  dtype=np.float).reshape(N,3)
        try:
            tag = np.fromfile("results/initial_dump_tag",  dtype=np.int)
            c = tag
        except FileNotFoundError:
            c = "b"
            
        self.axes.plot(r[:,0], r[:,1], r[:,2],'.', c=c)
        try:
            self.axes.set_box_aspect((np.ptp(r[:,0]), np.ptp(r[:,1]), np.ptp(r[:,2])))
        except AttributeError:
            axisEqual3D(self.axes)  
        self.axes.view_init(0, 90)



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

        #Create panel
        self.sizer_top = wx.BoxSizer(wx.HORIZONTAL)
        self.panel_1 = self.create_input_panel()
        self.sizer_top.Add(self.panel_1, 1, wx.EXPAND, 0)


        #Add display panel
        sizer_3 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_top.Add(sizer_3, 1, wx.EXPAND, 0)

        self.panel_2 = CanvasPanel(self)
        self.panel_2.draw()
        sizer_3.Add(self.panel_2, 1, wx.EXPAND, 0)

        choices = ["choice 1", "choice 2", "choice 3", "choice 4"]
        self.radio_box_1 = wx.RadioBox(self, wx.ID_ANY, "", 
                                       choices=choices, majorDimension=2, 
                                       style=wx.RA_SPECIFY_COLS)
        self.radio_box_1.SetSelection(0)
        sizer_3.Add(self.radio_box_1, 0, wx.EXPAND, 0)

        self.SetSizer(self.sizer_top)

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

        panel = wx.Panel(self, wx.ID_ANY)

        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2 = wx.FlexGridSizer(1, 3, 0, 0)
        self.sizer_1.Add(self.sizer_2, 1, wx.ALL | wx.EXPAND, 0)

        #########################################
        # Search control for top level keywords #
        #########################################
        self.searchctrl = BitmapComboBox(panel,
                              size=wx.DefaultSize,  
                              choices=[],
                              style= wx.TE_PROCESS_ENTER)# | wx.CB_SORT)

        self.sizer_2.Add(self.searchctrl, 0, 0, 0)

        self.searchctrl.Bind(wx.EVT_TEXT_ENTER, self.on_search)
        self.searchctrl.Bind(wx.EVT_COMBOBOX, self.on_search)
        self.searchctrl.Bind(wx.EVT_TEXT, self.on_search)
        self.searchctrl.Bind(wx.EVT_COMBOBOX_CLOSEUP, self.on_search)

        self.searchctrl.SetFocus()

        #########################################
        #             Run Button                #
        #########################################
        self.runbtn = wx.Button(panel, wx.ID_ANY, "Run")
        self.sizer_2.Add(self.runbtn, 0, 0, 0)
        self.runbtn.Bind(wx.EVT_BUTTON, self.run_btn)

        #########################################
        #             Help panel                #
        #########################################
        self.helptxt = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(self.width/3, self.height/4))
        self.helptxt.SetValue("Help Text \n\n\n\n\n\n\n")
        self.sizer_1.Add(self.helptxt, 0, wx.EXPAND, 0)

        #########################################
        #           Property Grid               #
        #########################################
        self.propgrid = wxpg.PropertyGridManager(panel, wx.ID_ANY,
                        style=wxpg.PG_SPLITTER_AUTO_CENTER)
        self.sizer_1.Add(self.propgrid, 1, wx.EXPAND, 0)

        self.propgrid.Bind(wxpg.EVT_PG_CHANGED, self.change_propgrid)
        panel.SetSizer(self.sizer_1)

        self.Layout()
        self.Refresh()
        self.Update()
        panel.Layout()
        panel.Refresh()
        panel.Update()
        #self.Layout()
        #self.Refresh()
        #self.SetSize((self.width, self.height))
        #self.Update()

        return panel

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

        try:
            self.ChangeDict[key][prop] = val
        except KeyError:
            self.ChangeDict[key] = {prop:val}

        print("after", key, prop, self.InputsDict[key]["vars"][prop], self.ChangeDict)
        return

    def run_btn(self, event): 

        btn = event.GetEventObject().GetLabel() 

        # otherwise ask the user what new file to open
        with wx.FileDialog(self, "Choose output directory",
                            defaultDir='./',
                            style=wx.FD_OPEN) as folderDiag:

            if folderDiag.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

            fdir = folderDiag.GetPath()
            print("Label of button = ", btn, fdir)

    def OnOpen(self, event):

        # otherwise ask the user what new file to open
        with wx.FileDialog(self, "Open Flowmol input file", 
                           wildcard="in files (*.in)|*.in",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

            # Proceed loading the file chosen by the user
            fdir = fileDialog.GetPath()
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
        self.Destroy()


    def OnQuit(self, e):
        self.Close()


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

