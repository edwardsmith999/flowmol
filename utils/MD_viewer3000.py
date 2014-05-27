#!/usr/bin/env python
# -*- coding: utf-8 -*-
# generated by wxGlade 0.6.4 on Sun May 25 18:19:27 2014
import numpy as np
import wx

import matplotlib
matplotlib.interactive(False)
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib.pyplot import gcf, setp
import matplotlib.pyplot as plt 

#from MDFields import *
#from HeaderData import *
from MD_PostProc import MD_PostProc

# begin wxGlade: extracode
# end wxGlade




# We are going to need 2 levels of frame
#   (o) Top level "FdirPanel" has the folder select dialogue. On folder select (or plotable files found in default directory) it instantiates the plotting frame
#         > Next level is the "PlotPanel" and includes all radio boxes etc. This is destroyed when new folder path is selected
# Questions -- How do we specify and fix the layout? Should these be frames or panels? 


class MyFrame(wx.Frame):

    def __init__(self, *args, **kwds):
        # begin wxGlade: MyFrame.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)

        self.xyz= 0
        self.normal = 0
        self.fdir = '../MD_dCSE/src_code/results/'
        self.fielddict = MD_PostProc(self.fdir)
        self.field = self.fielddict.plotlist.values()[0]
        maxpos = int(len(self.field.grid[self.normal])-1)
        maxrec = int(self.field.Raw.maxrec)
        self.pos = int(maxpos/2.)
        self.rec = int(maxrec/2.)
        self.get_orthogonal(self.normal,self.xyz,self.pos)
        self.redraw = self.redraw_plot
        self.update = self.update_plot
        print("Field name = ", self.field)

        #Textbox and buttons to select directory
        self.text_ctrl_1 = wx.TextCtrl(self, -1, self.fdir,style=wx.TE_PROCESS_ENTER)
        self.button_1 = wx.Button(self, -1, "...")

        self.radio_plotype_2 = wx.RadioBox(self, -1, "Plot_type", choices=["Profile", "Contour", "3D"], majorDimension=0, style=wx.RA_SPECIFY_COLS)
        self.radio_fieldtype_1 = wx.RadioBox(self, -1, "Field_type", choices=self.fielddict.plotlist.keys(), majorDimension=0, style=wx.RA_SPECIFY_ROWS)
        try:
            self.combo_xyz = wx.ComboBox(self, -1, choices=self.field.labels, style=wx.TE_PROCESS_ENTER)
        except AttributeError:
            self.combo_xyz = wx.ComboBox(self, -1, choices=[str(x) for x in range(self.field.nperbin)], style=wx.TE_PROCESS_ENTER)


        self.combo_norm = wx.ComboBox(self, -1, choices=["x", "y", "z"], style=wx.TE_PROCESS_ENTER)
        self.panel_2 = MatPlotPanel(self, -1)
        self.slider_pos_1  = wx.Slider(self, -1,  self.pos, 0, maxpos)
        self.slider_time_2 = wx.Slider(self, -1,  self.rec, 0, maxrec)

        self.Bind(wx.EVT_TEXT_ENTER, self.handle_changeFdir, self.text_ctrl_1)
        self.Bind(wx.EVT_BUTTON, self.OpenFileDialogue, self.button_1)
        self.Bind(wx.EVT_RADIOBOX, self.Handle_radio_plotype, self.radio_plotype_2)
        self.Bind(wx.EVT_RADIOBOX, self.Handle_radio_field, self.radio_fieldtype_1)
        self.Bind(wx.EVT_COMBOBOX, self.Handle_combo_xyz, self.combo_xyz)
        self.Bind(wx.EVT_COMBOBOX, self.handle_combo_normal, self.combo_norm)
        self.Bind(wx.EVT_COMMAND_SCROLL, self.handle_slider_pos, self.slider_pos_1)
        self.Bind(wx.EVT_COMMAND_SCROLL, self.handle_slider_time, self.slider_time_2)
        self.Bind(wx.EVT_LEFT_DCLICK,self.panel_2.debugtest, self.panel_2) 

        self.__set_properties()
        self.__do_layout()
        self.__setup_menubar()
        self.redraw()
        self.update()
        self.toggle_position_slider("Off") 
        # end wxGlade


    def __set_properties(self):
        # begin wxGlade: MyFrame.__set_properties
        self.SetTitle("MDViewer_3000")
        _icon = wx.EmptyIcon()
        _icon.CopyFromBitmap(wx.Bitmap("./moleculecell.png", wx.BITMAP_TYPE_ANY))
        self.SetIcon(_icon)

        self.radio_plotype_2.SetSelection(0)
        self.radio_fieldtype_1.SetSelection(0)
        self.combo_xyz.SetSelection(0)
        self.combo_norm.SetSelection(0)
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: MyFrame.__do_layout
        sizer_3 = wx.BoxSizer(wx.VERTICAL)
        sizer_4 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_5 = wx.BoxSizer(wx.VERTICAL)
        #sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7 = wx.BoxSizer(wx.VERTICAL)
        sizer_8 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_8.Add(self.text_ctrl_1, 8, 0, 0)
        sizer_8.Add(self.button_1, 1, wx.ALIGN_RIGHT, 0)
        sizer_7.Add(sizer_8, 1, wx.EXPAND, 0)
        sizer_3.Add(sizer_7, 1, wx.EXPAND, 0)
        self.sizer_5.Add(self.radio_plotype_2, 1, wx.EXPAND, 0)
        self.sizer_5.Add(self.radio_fieldtype_1, 5, wx.EXPAND, 0)
        self.sizer_5.Add(self.combo_xyz, 1, 0)
        self.sizer_5.Add(self.combo_norm, 1, 0)
        #sizer_5.Add(sizer_6, 1, 0, 3)
        sizer_4.Add(self.sizer_5, 1, 0, 0)
        sizer_4.Add(self.panel_2, 3, wx.EXPAND, 0)
        sizer_3.Add(sizer_4, 6, wx.EXPAND, 0)
        sizer_3.Add(self.slider_pos_1, 0, wx.EXPAND, 0)
        sizer_3.Add(self.slider_time_2, 0, wx.EXPAND, 0)
        self.SetSizer(sizer_3)
        sizer_3.Fit(self)
        self.Layout()
        # end wxGlade

    def __setup_menubar(self):
        # Menu Bar
        self.frame_2_menubar = wx.MenuBar()
        wxglade_tmp_menu = wx.Menu()
        self.frame_2_menubar.Append(wxglade_tmp_menu, "File")
        wxglade_tmp_menu = wx.Menu()
        self.frame_2_menubar.Append(wxglade_tmp_menu, "Edit")
        wxglade_tmp_menu = wx.Menu()
        self.frame_2_menubar.Append(wxglade_tmp_menu, "Help")
        self.SetMenuBar(self.frame_2_menubar)

    def Handle_radio_plotype(self, event):  # wxGlade: MyFrame.<event_handler>
        self.plottype = event.GetInt()
        if self.plottype == 0:
            self.redraw = self.redraw_plot
            self.update = self.update_plot
            self.toggle_position_slider("Off")
        elif self.plottype == 1:
            self.redraw = self.redraw_contour
            self.update = self.update_contour
            self.toggle_position_slider("On")
        elif self.plottype == 2:
            try:
                from mayavi import mlab
            except ImportError:
                self.showMessageDlg("3D plotting requires mayavi to be installed","Information", wx.OK|wx.ICON_INFORMATION)
        else:
            quit("Error in plotype specified")

        self.redraw()

    def Handle_radio_field(self, event):  # wxGlade: MyFrame.<event_handler>

        if (self.field == self.fielddict.plotlist[self.radio_fieldtype_1.GetStringSelection()]):
            pass
        else:
            self.field = self.fielddict.plotlist[self.radio_fieldtype_1.GetStringSelection()]
            print("field changed to ",self.radio_fieldtype_1.GetStringSelection(), "obj = ", self.field)
        self.combo_xyz.Clear()
        try:
            self.combo_xyz.AppendItems(self.field.labels)
        except AttributeError:
            self.combo_xyz.AppendItems([str(x) for x in range(self.field.nperbin)])
        #for label in self.field.labels:
        #    self.combo_xyz.append(label)
        self.combo_xyz.SetSelection(0)
        self.xyz = 0
        self.slider_time_2.SetMax(self.field.Raw.maxrec)
        if (self.field.Raw.maxrec < self.rec):
            print("Set max = " ,self.rec, self.field.Raw.maxrec,int(self.field.Raw.maxrec) < self.rec,self.field.Raw.maxrec > self.rec)
            self.rec = self.field.Raw.maxrec
            self.slider_time_2.SetValue(self.rec)

        self.redraw()

    def Handle_combo_xyz(self, event):  # wxGlade: MyFrame.<event_handler>

        self.xyz = event.GetInt()
        print("Component = ", self.xyz)
        self.get_orthogonal(self.normal,self.xyz,self.pos)
        self.redraw()

    def handle_combo_normal(self, event):  # wxGlade: MyFrame.<event_handler>
    
        self.normal = event.GetInt()
        self.get_orthogonal(self.normal,self.xyz,self.pos)
        self.slider_pos_1.SetMax(int(len(self.field.grid[self.normal])-1))
        print("Normal = ", self.normal, "Other components", self.naxes)
        self.redraw()

    def handle_slider_pos(self, event):  # wxGlade: MyFrame.<event_handler>

        self.pos = event.GetInt()
        print("Postition = ", self.pos)
        self.get_orthogonal(self.normal,self.xyz,self.pos)
        self.update()

    def handle_slider_time(self, event):  # wxGlade: MyFrame.<event_handler>

        self.rec = event.GetInt()
        print("Record = ", self.rec)
        self.update()

    def handle_changeFdir(self, event):  # wxGlade: MyFrame.<event_handler>
        self.fdir = self.text_ctrl_1.GetValue()
        self.fielddict = MD_PostProc(self.fdir)

        #Here we need to redefine the radio box
#        self.unbind(wx.EVT_RADIOBOX, self.Handle_radio_field, self.radio_fieldtype_1)
#        self.radio_fieldtype_1.Destroy()
#        self.radio_fieldtype_1 = wx.RadioBox(self, -1, "Field_type", choices=self.fielddict.plotlist.keys(), majorDimension=0, style=wx.RA_SPECIFY_ROWS)
#        self.sizer_5.Add(self.radio_fieldtype_1, 5, wx.EXPAND, 0)
#        self.Bind(wx.EVT_RADIOBOX, self.Handle_radio_field, self.radio_fieldtype_1)

        self.__set_properties()
        self.field = self.fielddict.plotlist.values()[0]
        self.slider_pos_1.SetMax(int(len(self.field.grid[0])-1))
        self.slider_time_2.SetMax(int(len(self.field.Raw.maxrec)))
        print("Directory = ", self.fdir)
        print(self.fielddict)

        self.redraw()
        self.Refresh()

    def OpenFileDialogue(self, event):  # wxGlade: MyFrame.<event_handler>

        fdir = ""  # Use  folder as a flag
        dlg = wx.DirDialog(self, message="Choose a file", defaultPath ="../MD_dCSE/src_code/results/")
 
        if dlg.ShowModal() == wx.ID_OK:
 
            # get the new folder from the dialog
            fdir = dlg.GetPath() + "/"
            dlg.SetPath(fdir)
        dlg.Destroy()  # best to do this sooner than later

        if fdir:
            # use the file name
            self.text_ctrl_1.SetValue(fdir)
            event = wx.PyCommandEvent(wx.EVT_TEXT_ENTER.typeId, self.text_ctrl_1.GetId())
            self.GetEventHandler().ProcessEvent(event)
            #self.handle_changeFdir(fdir)

        self.redraw()
        self.Refresh()

    def get_orthogonal(self,normal,xyz,pos):
        self.naxes = []
        for ixyz in [0,1,2]:
            if ixyz != normal:
                self.naxes.append(ixyz)
        self.binlimits = [None,None,None]
        self.binlimits[normal] = (pos,pos)

    def handle_forceredraw(self,event):
        self.redraw()
        self.Refresh()

    def update_plotdata(self):
        print("Plt Field name = ", self.field,self.normal,self.rec)
        self.ax1, self.data = self.field.profile(self.normal,startrec=self.rec,endrec=self.rec)

    def update_contourdata(self):
        print("Contour Field name = ", self.field,self.naxes,self.rec,self.binlimits)
        self.ax1, self.ax2, self.data = self.field.contour(axes=self.naxes,
                                                          startrec=self.rec,
                                                          endrec=self.rec,
                                                          binlimits=self.binlimits)

    def redraw_plot(self):
        self.panel_2.figure.clf(keep_observers=True)
        self.panel_2.subplot1 = self.panel_2.figure.add_subplot(111)
        self.update_plotdata()
        print(self.ax1.shape,self.data.shape)
        self.lines = self.panel_2.subplot1.plot(self.ax1, self.data[:,self.xyz], 'r-o', linewidth=2)
        self.panel_2.subplot1.set_xlim([self.ax1.min(), self.ax1.max()])
        #self.panel_2.subplot1.set_ylim([self.data.min(), self.data.max()])
        self.panel_2.subplot1.set_xlabel("x", fontsize = 12)
        self.panel_2.subplot1.set_ylabel("u", fontsize = 12)
        self.panel_2.canvas.draw()
        self.Refresh()

    def update_plot(self):
        self.update_plotdata()
        setp(self.lines, xdata=self.ax1, ydata=self.data[:,self.xyz])
        self.panel_2.canvas.draw()
        self.Refresh()

    def redraw_contour(self):

        self.panel_2.figure.clf(keep_observers=True)
        self.panel_2.subplot1 = self.panel_2.figure.add_subplot(111)
        self.update_contourdata()
        self.colormesh = self.panel_2.subplot1.pcolormesh(self.ax1, self.ax2, self.data[:,1:,self.xyz],cmap=self.panel_2.cmap)
        self.cbar = self.panel_2.figure.colorbar(self.colormesh)
        self.panel_2.subplot1.set_xlim([self.ax1.min(), self.ax1.max()])
        self.panel_2.subplot1.set_ylim([self.ax2.min(), self.ax2.max()])
        self.panel_2.subplot1.set_ylabel("z", fontsize = 12)
        self.panel_2.subplot1.set_xlabel("x", fontsize = 12)
        self.panel_2.figure.canvas.draw_idle()
        self.Refresh()

    def update_contour(self):

        """
            update canvas with pcolormesh and colorbar
        """

        #update the current frame
        self.update_contourdata()
        self.colormesh.set_array(self.data[:,1:,self.xyz].ravel())
        self.panel_2.figure.canvas.draw_idle()
        self.Refresh()


    def toggle_position_slider(self,switchon):
        if (switchon == "On"):
            self.slider_pos_1.Enable(True)
            self.slider_pos_1.SetBackgroundColour(self.slider_time_2.GetBackgroundColour())
            self.slider_pos_1.SetTransparent(255)
        elif(switchon == "Off"):
            self.slider_pos_1.SetValue(0)
            self.slider_pos_1.Enable(False)
            self.slider_pos_1.SetBackgroundColour("Black")
            self.slider_pos_1.SetTransparent(10)
            #self.slider_pos_1.SetForegroundColour("Blue")

        else:
            quit("Error - toggle_position_slider must be str On of Off")

    def showMessageDlg(self, msg, title, style):
        """"""
        dlg = wx.MessageDialog(parent=None, message=msg, 
                               caption=title, style=style)
        dlg.ShowModal()
        dlg.Destroy()


# end of class MyFrame


class MatPlotPanel(wx.Panel):

    """
        The panel where the matplotlib canvas is assigned and data plotted
    """

    def __init__(self, *args, **kwargs):

        wx.Panel.__init__(self, *args, **kwargs)

        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        rect = self.figure.patch
        rect.set_facecolor((0.9,0.9,0.9)) #(self.GetBackgroundColour().rgb())
        self.cmap = plt.cm.RdYlBu_r

    def debugtest(self,event):
        quit("CLICKY")






#Framework for custom event, from the answer at
# http://stackoverflow.com/questions/747781/wxpython-calling-an-event-manually

#myEVT_CUSTOM = wx.NewEventType()
#EVT_CUSTOM = wx.PyEventBinder(myEVT_CUSTOM, 1)

#event = MyEvent(myEVT_CUSTOM, self.GetId())
#event.SetMyVal('here is some custom data')
#self.GetEventHandler().ProcessEvent(event)

#self.Bind(EVT_CUSTOM, self.on_event)

#def on_event(self, e):
#    data = e.GetMyVal()
#    print 'custom data is: {0}'.format(data)

class MyEvent(wx.PyCommandEvent):
    def __init__(self, evtType, id):
        wx.PyCommandEvent.__init__(self, evtType, id)
        myVal = None

    def SetMyVal(self, val):
        self.myVal = val

    def GetMyVal(self):
        return self.myVal




if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    frame_1 = MyFrame(None, -1, "")
    app.SetTopWindow(frame_1)
    frame_1.Show()
    app.MainLoop()
