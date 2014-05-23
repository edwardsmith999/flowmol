import numpy as np
import wx

import matplotlib
matplotlib.interactive(False)
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
from matplotlib.pyplot import gcf, setp
import matplotlib.pyplot as plt

from MDFields import *
from HeaderData import *
from MD_PostProc import MD_PostProc

class Knob:
    """
    Knob - simple class with a "setKnob" method.  
    A Knob instance is attached to a Param instance, e.g., param.attach(knob)
    Base class is for documentation purposes.
    """
    def setKnob(self, value):
        pass


class Param:
    """
    The idea of the "Param" class is that some parameter in the GUI may have
    several knobs that both control it and reflect the parameter's state, e.g.
    a slider, text, and dragging can all change the value of the frequency in
    the waveform of this example.  
    The class allows a cleaner way to update/"feedback" to the other knobs when 
    one is being changed.  Also, this class handles min/max constraints for all
    the knobs.
    Idea - knob list - in "set" method, knob object is passed as well
      - the other knobs in the knob list have a "set" method which gets
        called for the others.
    """
    def __init__(self, initialValue=None, minimum=0, maximum=1):
        self.minimum = minimum
        self.maximum = maximum
        if initialValue != self.constrain(initialValue):
            raise ValueError('illegal initial value')
        self.value = initialValue
        self.knobs = []
        
    def attach(self, knob):
        self.knobs += [knob]
        
    def set(self, value, knob=None):
        self.value = value
        self.value = self.constrain(value)
        for feedbackKnob in self.knobs:
            if feedbackKnob != knob:
                feedbackKnob.setKnob(self.value)
        return self.value

    def constrain(self, value):
        if value <= self.minimum:
            value = self.minimum
        if value >= self.maximum:
            value = self.maximum
        return value


class SliderGroup(Knob):

    def __init__(self, parent, label, param):
        self.sliderLabel = wx.StaticText(parent, label=label)
        self.sliderText = wx.TextCtrl(parent, -1, style=wx.TE_PROCESS_ENTER)
        self.slider = wx.Slider(parent, -1)
        self.slider.SetMax(param.maximum)
        self.setKnob(param.value)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.sliderLabel, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.sliderText,  0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.slider, 1, wx.EXPAND)
        self.sizer = sizer

        self.slider.Bind(wx.EVT_SLIDER, self.sliderHandler)
        self.sliderText.Bind(wx.EVT_TEXT_ENTER, self.sliderTextHandler)

        self.param = param
        self.param.attach(self)

    def sliderHandler(self, evt):
        value = evt.GetInt()
        self.param.set(value)
        
    def sliderTextHandler(self, evt):
        value = float(self.sliderText.GetValue())
        self.param.set(value)
        
    def setKnob(self, value):
        self.sliderText.SetValue('%g'%value)
        self.slider.SetValue(value)


class ContourSliderFrame(wx.Frame):

    """
        The frame is the top level here and contains everything
        else (sizers and all that they contain). This is what
        would normally be called a window on a program
    """
    def __init__(self, fielddict, *args, **kwargs):

        #Take in field and call wx.Frame's constructor
        wx.Frame.__init__(self, *args, **kwargs)

        #Save field dictonary
        self.fielddict = fielddict

        # Setup frame layout with sizers
        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_1.Add(self.sizer_2, 1, wx.EXPAND, 0)

        #Add radiobox to select fieldtype
        self.radio_box_1 = wx.RadioBox(self, -1, "Field", choices=self.fielddict.plotlist.keys(), 
                                        majorDimension=15, style=wx.RA_SPECIFY_ROWS)
        self.Bind(wx.EVT_RADIOBOX, self.OnRadio,self.radio_box_1) 
        self.radio_box_1.SetSelection(1)
        self.field = self.fielddict.plotlist[self.radio_box_1.GetStringSelection()]

        self.sizer_2.Add(self.radio_box_1, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        self.setup_panel()

    def setup_panel(self):

        """
            Setup sliders and main matplotlib viewing window
        """

        #Setup main contour viewing window
        self.ContourSliderWindow = ContourSliderWindow(self)
        self.sizer_2.Add(self.ContourSliderWindow, 1, wx.EXPAND)

        #Add Sliders
        #Draw a component slider and add to sizer_1 if needed, otherwise don't
        if int(self.field.nperbin)-1 > 0:
            self.componentSliderGroup = SliderGroup(self,label='x,y,z:', \
                                                    param=self.ContourSliderWindow.comp)
            self.sizer_1.Add(self.componentSliderGroup.sizer, 0, \
                wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)
        else:
            self.componentSliderGroup = None

        #Draw a bin number slider and add to sizer_1
        self.positionSliderGroup = SliderGroup(self, label='Position: ', \
                                                param=self.ContourSliderWindow.pos)
        self.sizer_1.Add(self.positionSliderGroup.sizer, 0, \
            wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)

        #Draw a record number slider and add to sizer_1
        self.recordSliderGroup = SliderGroup(self,   label='Record:   ', \
                                                param=self.ContourSliderWindow.rec)
        self.sizer_1.Add(self.recordSliderGroup.sizer, 0, \
            wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)
        self.SetSizer(self.sizer_1)


    def reset_panel(self):

        # Setup frame layout with sizers
        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_1.Add(self.sizer_2, 1, wx.EXPAND, 0)

        #Add radiobox to select fieldtype
        self.radio_box_1 = wx.RadioBox(self, -1, "Field", choices=self.fielddict.plotlist.keys(), 
                                        majorDimension=15, style=wx.RA_SPECIFY_ROWS)
        self.Bind(wx.EVT_RADIOBOX, self.OnRadio,self.radio_box_1) 
        self.radio_box_1.SetSelection(0)
        self.sizer_2.Add(self.radio_box_1, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        
        self.setup_panel()
        
    def OnRadio(self, evt): 
        if (self.field == self.fielddict.plotlist[self.radio_box_1.GetStringSelection()]):
            pass
        else:
            #Ideally we should destroy current window and replace with a new one!!
            self.field = self.fielddict.plotlist[self.radio_box_1.GetStringSelection()]
            print("field changed to ",self.radio_box_1.GetStringSelection(), "obj = ", self.field)
            #self.ContourSliderWindow.Destroy()
            #self.reset_panel()

    def setup_menubar(self):
        # Menu Bar
        self.frame_2_menubar = wx.MenuBar()
        wxglade_tmp_menu = wx.Menu()
        self.frame_2_menubar.Append(wxglade_tmp_menu, "File")
        wxglade_tmp_menu = wx.Menu()
        self.frame_2_menubar.Append(wxglade_tmp_menu, "Edit")
        wxglade_tmp_menu = wx.Menu()
        self.frame_2_menubar.Append(wxglade_tmp_menu, "Help")
        self.SetMenuBar(self.frame_2_menubar)


class ContourSliderWindow(wx.Window, Knob):

    """
        The window is where the canvas is defined and data plotted
    """

    def __init__(self, *args, **kwargs):

        wx.Window.__init__(self, *args, **kwargs)
        self.frameimin = self.GetParent()
        self.field = self.frameimin.field 
        self.cmap = plt.cm.RdYlBu_r
        self.colormesh = []
        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)

        #Define behaviour of sliders
        print(self.frameimin.field)
        if int(self.field.nperbin)-1 > 0:
            self.comp = Param(0,        minimum=0, maximum=int(self.field.nperbin)-1)
        else:
            self.comp = Param(0,minimum=0, maximum=1)
            self.comp.value = 0
        self.pos  = Param(int(len(self.field.grid[1])/2.),  minimum=0, maximum=int(len(self.field.grid[1])-1))
        self.rec  = Param(int(self.field.maxrec/2.),        minimum=0, maximum=int(self.field.Raw.maxrec))

        #Bind the slider interface to sizehandler
        self.comp.attach(self)
        self.pos.attach(self)
        self.rec.attach(self)
        self.Bind(wx.EVT_SIZE, self.sizeHandler)

        #Draw initial contour
        self.draw_initial()
       
    def sizeHandler(self, *args, **kwargs):
        self.canvas.SetSize(self.GetSize())

    def draw_initial(self):

        """
            Draw intial canvas with pcolormesh and colorbar
        """

        if not hasattr(self, 'subplot1'):
            self.subplot1 = self.figure.add_subplot(111)
        component = int(self.comp.value)
        self.initpos = int(len(self.field.grid[1])/2.)
        self.initrec = int(self.field.maxrec/2.)
        self.ax1, self.ax2, data = self.field.contour(axes=naxes,startrec=self.initrec,endrec=self.initrec,
                                                      binlimits=[None,(self.initpos,self.initpos),None])
        self.colormesh = self.subplot1.pcolormesh(self.ax1, self.ax2, data[:,1:,component],cmap=self.cmap)

        #Set some plot attributes
        self.cbar = self.figure.colorbar(self.colormesh)
        self.subplot1.set_xlim([self.ax1.min(), self.ax1.max()])
        self.subplot1.set_ylim([self.ax2.min(), self.ax2.max()])
        self.subplot1.set_ylabel("z", fontsize = 12)
        self.subplot1.set_xlabel("x", fontsize = 12)
        rect = self.figure.patch
        rect.set_facecolor((0.9,0.9,0.9))#(self.GetBackgroundColour().rgb())

    def setKnob(self, value):

        """
            Update canvas with pcolormesh and colorbar
        """
        component = int(self.comp.value)
        rec = int(np.rint(self.rec.value))
        pos = int(np.rint(self.pos.value))
        if (self.field != self.frameimin.field):
            self.figure.clf(keep_observers=True)
            self.subplot1 = self.figure.add_subplot(111)
            self.field = self.frameimin.field
            ax1, ax2, data = self.field.contour(axes=naxes,startrec=rec,endrec=rec,binlimits=[None,(pos,pos),None])
            self.colormesh = self.subplot1.pcolormesh(self.ax1, self.ax2, data[:,1:,component],cmap=self.cmap)
            self.cbar = self.figure.colorbar(self.colormesh)
            self.subplot1.set_xlim([self.ax1.min(), self.ax1.max()])
            self.subplot1.set_ylim([self.ax2.min(), self.ax2.max()])
            self.subplot1.set_ylabel("z", fontsize = 12)
            self.subplot1.set_xlabel("x", fontsize = 12)
            rect = self.figure.patch
            rect.set_facecolor((0.9,0.9,0.9)) #(self.GetBackgroundColour().rgb())
        else:
            #Update the current frame
            ax1, ax2, data = self.field.contour(axes=naxes,startrec=rec,endrec=rec,binlimits=[None,(pos,pos),None])
            self.colormesh.set_array(data[:,1:,component].ravel())

        self.figure.canvas.draw_idle()

class App(wx.App):

    """
        The application -- calls the frame to initialise and run by call to Mainloop
    """

    def __init__(self,fielddict):

        self.fielddict = fielddict
        wx.App.__init__(self)

    def OnInit(self):
        self.frame1 = ContourSliderFrame(fielddict=self.fielddict ,parent=None, title="ContourSlider", size=(640, 480))
        self.frame1.Show()
        return True

component = 0

naxes = (0,2)

#Setup array of values to plot
#fdir = '../MD_dCSE/src_code/results/'
fdir = '/home/es205/scratch/Re400/iter1918000_to_2233899/'
fielddict = MD_PostProc(fdir)
print(fielddict)

#var = 'vel' #raw_input("Please enter choice of field: ")
#field = fielddict.plotlist[var]

app = App(fielddict=fielddict)
app.MainLoop()
