import wx
from PyplotPanel import PyplotPanel
from FieldChooser import FieldChooserPanel
from SliderPanels import RecordSliderPanel

from MD_PostProc import MD_PostProc  

class VisualiserPanel(wx.Panel):
 
    def __init__(self,parent,fdir,**kwargs):

        wx.Panel.__init__(self,parent,**kwargs)

        self.fdir = fdir
        self.MD_PP = MD_PostProc(self.fdir,cpol_bins=True)
        self.field = self.MD_PP.plotlist.values()[0]

        self.pyplotp = PyplotPanel(self)
        self.choosep = FieldChooserPanel(self)
        self.slidersp = RecordSliderPanel(self)
    
        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.choosep, 0, wx.EXPAND | wx.ALL)
        hbox.Add(self.pyplotp, 1, wx.EXPAND | wx.ALL)
        vbox.Add(hbox, 1, wx.EXPAND | wx.ALL)
        vbox.Add(self.slidersp, 0, wx.EXPAND | wx.ALL)
        self.SetSizer(vbox)

        self.set_bindings()
        self.set_defaults()

    def set_defaults(self):

        self.autoscale = False
        self.update_components()
        self.component = 0
        self.normal = 0
        self.redraw = self.redraw_plot
        self.update = self.update_plot

        self.maxbin = len(self.field.grid[self.normal])-1
        self.maxrec = self.field.Raw.maxrec
        self.slidersp.recslider.SetMax(self.maxrec)
        self.slidersp.binslider.SetMax(self.maxbin)

        self.bin = int(float(self.maxbin)/2.)
        self.rec = int(float(self.maxrec)/2.)
        self.recwidth = 0
        self.binwidth = 0
        self.redraw()
        self.SetBin(self.bin)
        self.SetRecord(self.rec)

    def set_bindings(self):

        self.Bind(wx.EVT_RADIOBOX, self.handle_plottype, 
                  self.choosep.plottype_p.fieldradiobox)
        self.Bind(wx.EVT_RADIOBOX, self.handle_fieldtype, 
                  self.choosep.fieldtype_p.fieldradiobox)
        self.Bind(wx.EVT_COMBOBOX, self.handle_component, 
                  self.choosep.component_p.componentcombobox)
        self.Bind(wx.EVT_COMBOBOX, self.handle_normal, 
                  self.choosep.component_p.normalcombobox)

        self.Bind(wx.EVT_COMMAND_SCROLL, self.handle_recslider, 
                  self.slidersp.recslider.slider)
        self.Bind(wx.EVT_COMMAND_SCROLL, self.handle_binslider, 
                  self.slidersp.binslider.slider)
        self.Bind(wx.EVT_TEXT_ENTER, self.handle_rectxt,
                  self.slidersp.recslider.slidertext)
        self.Bind(wx.EVT_TEXT_ENTER, self.handle_bintxt,
                  self.slidersp.binslider.slidertext)

        self.Bind(wx.EVT_SPINCTRL, self.handle_recspin,
                  self.slidersp.recslider.spin)
        self.Bind(wx.EVT_SPINCTRL, self.handle_binspin,
                  self.slidersp.binslider.spin)

        self.Bind(wx.EVT_CHECKBOX, self.handle_autoscale, 
                  self.choosep.autoscale_b) 


    def handle_plottype(self, event):
        print('Plot type changed!')  
        plottype = event.GetString()
        if plottype == 'Profile':
            self.redraw = self.redraw_plot
            self.update = self.update_plot
            #self.toggle_position_slider("Off")
        elif plottype == 'Contour':
            self.redraw = self.redraw_contour
            self.update = self.update_contour
            #self.toggle_position_slider("On")
        #else: 
            #try:
            #    from mayavi import mlab
            #except ImportError:
            #    self.showMessageDlg("3D plotting requires mayavi to be installed",
            #                        "Information", wx.OK|wx.ICON_INFORMATION)
        else:
            quit("Error in plotype specified")

        self.redraw()

    def handle_fieldtype(self, event):
        ftype = event.GetString()
        if (self.field == self.MD_PP.plotlist[ftype]):
            pass
        else:
            self.field = self.MD_PP.plotlist[ftype]
        self.update_components()

        self.maxbin = len(self.field.grid[self.normal]) - 1
        self.maxrec = self.field.Raw.maxrec
        self.slidersp.recslider.slider.SetMax(self.maxrec)
        self.slidersp.binslider.slider.SetMax(self.maxbin)
        if (self.rec > self.maxrec):
            self.SetRecord(self.maxrec)
        if (self.bin > self.maxbin):
            self.SetBin(self.maxbin)

        self.redraw()

    def update_components(self):
        self.choosep.component_p.componentcombobox.Clear()
        try:
            self.choosep.component_p.componentcombobox.AppendItems(self.field.labels)
        except AttributeError:
            self.choosep.component_p.componentcombobox.AppendItems(
                [str(x) for x in range(self.field.nperbin)])
        self.choosep.component_p.componentcombobox.SetSelection(0)
        self.component = 0

    def handle_component(self, event):
        self.component = event.GetInt()
        self.redraw()
    def handle_normal(self, event):
        self.normal = event.GetInt()
        self.slidersp.binslider.SetMax(len(self.field.grid[self.normal])-1) 
        self.redraw()
    def handle_autoscale(self,event):
        self.autoscale = event.GetInt()
        self.redraw()


    def handle_recslider(self, event):
        self.SetRecord(event.GetInt())
    def handle_rectxt(self, event):
        rec = int(event.GetString())
        if (rec > self.maxrec):
            self.SetRecord(self.maxrec)
        else:
            self.SetRecord(rec)
    def handle_recspin(self, event):
        width = event.GetInt()
        self.SetRecordWidth(width)


    def handle_binslider(self, event):
        self.SetBin(event.GetInt()) 
    def handle_bintxt(self, event):
        bin = int(event.GetString())
        if (bin > self.maxbin):
            self.SetBin(self.maxbin)
        else:
            self.SetBin(bin)
    def handle_binspin(self, event):
        width = event.GetInt()
        self.SetBinWidth(width)

    def SetRecord(self, rec):
        self.rec = rec
        if (self.rec + self.recwidth > self.maxrec):
            self.rec = self.maxrec-self.recwidth
        elif (self.rec - self.recwidth < 0):
            self.rec = self.recwidth
        self.slidersp.recslider.SetValue(self.rec)
        if self.autoscale:
            self.redraw()
        else:
            self.update()
    def SetRecordWidth(self, width):
        self.recwidth = width
        if (self.recwidth > self.maxrec/2):
            self.recwidth = self.maxrec/2
            self.SetRecord(self.maxrec/2)
        self.SetRecord(self.rec)

    def SetBin(self, bin):
        self.bin = bin
        if (self.bin + self.binwidth > self.maxbin):
            self.bin = self.maxbin - self.binwidth
        elif (self.bin - self.binwidth < 0):
            self.bin = self.binwidth
        self.slidersp.binslider.SetValue(self.bin)
        if self.autoscale:
            self.redraw()
        else:
            self.update()
    def SetBinWidth(self, width):
        self.binwidth = width
        if (self.binwidth > self.maxbin/2):
            self.binwidth = self.maxbin/2
            self.SetBin(self.maxbin/2)
        self.SetBin(self.bin)

    # Below this point are the data updating/plotting routines

    def get_contour_data(self):
        naxes = [0,1,2]
        naxes.remove(self.normal)
        binlimits = [None]*3
        binlimits[self.normal] = (self.bin-self.binwidth, self.bin+self.binwidth)
        ax1, ax2, data = self.field.contour(axes=naxes, 
                                            startrec=self.rec-self.recwidth,
                                            endrec=self.rec+self.recwidth,
                                            binlimits=binlimits,
                                            quit_on_error=False)
        return ax1, ax2, data

    def redraw_plot(self):
        ax, data = self.field.profile(self.normal, 
                                      startrec=self.rec-self.recwidth, 
                                      endrec=self.rec+self.recwidth, 
                                      quit_on_error=False)
        self.pyplotp.redraw_plot(ax, data[:,self.component])
        self.Refresh()
    def update_plot(self):
        ax, data = self.field.profile(self.normal, 
                                      startrec=self.rec-self.recwidth, 
                                      endrec=self.rec+self.recwidth, 
                                      quit_on_error=False)
        self.pyplotp.update_plot(ax, data[:,self.component])
        self.Refresh()

    def redraw_contour(self):
        ax1, ax2, data = self.get_contour_data()
        self.pyplotp.redraw_contour(ax1, ax2, data[:,:,self.component])
        self.Refresh()

    def update_contour(self):
        ax1, ax2, data = self.get_contour_data()
        self.pyplotp.update_contour(data[:,:,self.component])
        self.Refresh()

    def toggle_position_slider(self,switchon):
        pass
#        if (switchon == "On"):
#            self.slider_pos_1.Enable(True)
#            self.slider_pos_1.SetBackgroundColour(self.slidersp.recslider.slider.GetBackgroundColour())
#            self.slider_pos_1.SetTransparent(255)
#        elif(switchon == "Off"):
#            self.slider_pos_1.SetValue(0)
#            self.slider_pos_1.Enable(False)
#            self.slider_pos_1.SetBackgroundColour("Black")
#            self.slider_pos_1.SetTransparent(10)
#            #self.slider_pos_1.SetForegroundColour("Blue")#
#
#        else:
#            quit("Error - toggle_position_slider must be str On of Off")

