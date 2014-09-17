import wx
import wx.lib.scrolledpanel as scrolled

class PlotTypePanel(wx.Panel):

    def __init__(self,parent,**kwargs):
        wx.Panel.__init__(self,parent,**kwargs)
        choices = ['Profile','Contour','CPL']
        self.fieldradiobox = wx.RadioBox(self,label='Plot Type',    
                                    style=wx.RA_SPECIFY_COLS,
                                    choices=choices)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.fieldradiobox, 0, wx.EXPAND|wx.ALL, 10)
        self.SetSizer(vbox)

class FieldTypePanel(scrolled.ScrolledPanel):

    def __init__(self,parent,**kwargs):
        scrolled.ScrolledPanel.__init__(self, parent,**kwargs)
        choices = sorted(parent.parent.PP.plotlist.keys())
        self.fieldradiobox = wx.RadioBox(self,label='Field',    
                                    style=wx.RA_SPECIFY_ROWS,
                                    choices=choices)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.fieldradiobox, 0, wx.EXPAND|wx.ALL, 10)
        self.SetSizer(vbox)
        self.SetAutoLayout(1)
        self.SetupScrolling()

class SaveFigurePanel(wx.Panel):

    def __init__(self,parent,**kwargs):
        wx.Panel.__init__(self,parent,**kwargs)
        self.savebutton = wx.Button(self,-1,"Save figure...")

class FieldComponentPanel(wx.Panel):
    
    def __init__(self,parent,**kwargs):
        wx.Panel.__init__(self,parent,**kwargs)

        self.componenttitle = wx.StaticText(self,-1,label='Component',size=(100,-1))
        self.componentcombobox = wx.ComboBox(self, size=(10,-1), value='0')

        choices = ['0','1','2']
        self.normaltitle = wx.StaticText(self,-1,label='Normal', size=(100,-1))
        self.normalcombobox = wx.ComboBox(self, choices=choices, size=(10,-1), 
                                          value='0')

        grid = wx.GridBagSizer(hgap=3)
        grid.Add(self.componenttitle,    (0,0), flag=wx.EXPAND)
        grid.Add(self.componentcombobox, (1,0), flag=wx.EXPAND)
        grid.Add(self.normaltitle,       (0,1), flag=wx.EXPAND)
        grid.Add(self.normalcombobox,    (1,1), flag=wx.EXPAND)
        grid.AddGrowableCol(0)
        grid.AddGrowableCol(1)
        self.SetSizer(grid)

class FieldChooserPanel(wx.Panel):
    
    def __init__(self,parent,**kwargs):
        wx.Panel.__init__(self,parent,**kwargs)
        self.parent = parent 
        # Plot type chooser box
        self.plottype_p = PlotTypePanel(self)    
        # Field type chooser box
        self.fieldtype_p = FieldTypePanel(self, size = (-1, 400))
        # Component chooser combo box
        self.component_p = FieldComponentPanel(self)
        # Autoscale button
        self.autoscale_b = wx.CheckBox(self,-1,label='Autoscale')
        # Save button
        self.save_b = SaveFigurePanel(self)

        # Sizer
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.plottype_p, 0,wx.EXPAND, 0)
        vbox.Add(self.fieldtype_p,0,wx.EXPAND, 0)
        vbox.Add(self.component_p,0,wx.EXPAND, 0)
        vbox.Add(self.autoscale_b,0,wx.EXPAND, 0)
        vbox.Add(self.save_b,     0,wx.EXPAND, 0)
        self.SetSizer(vbox) 
