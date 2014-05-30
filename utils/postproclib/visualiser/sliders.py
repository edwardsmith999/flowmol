import wx

class SliderPlusWidth(wx.Panel):
    
    def __init__(self,parent,slidername,minval=0,maxval=100,**kwargs):

        wx.Panel.__init__(self,parent,**kwargs)
        sliderlabel = wx.StaticText(self,-1,label=slidername+':',size=(50,-1))
        self.slidertext = wx.TextCtrl(self,-1,style=wx.TE_PROCESS_ENTER,size=(50,-1))
        self.slider = wx.Slider(self, minValue=minval,maxValue=maxval)
        spintext = wx.StaticText(self,-1,label=u"\u00B1",size=(10,-1),)
        self.spin = wx.SpinCtrl(self,value='0',initial=0,size=(50,-1))

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(sliderlabel,0,wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 10)
        hbox.Add(self.slidertext,0,wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 10)
        hbox.Add(self.slider,1,wx.EXPAND,0)
        hbox.Add(spintext,0,wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 10)
        hbox.Add(self.spin,0,wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(hbox,1,wx.EXPAND,0)
        self.SetSizer(vbox) 

        self.SetValue(0)

    def SetValue(self,pos):
        self.slidertext.SetValue(str(pos))
        self.slider.SetValue(pos)
   
    def GetValue(self):
        return self.slider.GetValue() 

    def SetMax(self,maximum):
        self.slider.SetMax(maximum)

class RecordSliderPanel(wx.Panel):

    def __init__(self,parent,**kwargs):

        wx.Panel.__init__(self,parent,**kwargs)

        self.binslider = SliderPlusWidth(self, 'Bin')
        self.recslider = SliderPlusWidth(self, 'Rec')

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.binslider,1,wx.EXPAND,0)
        vbox.Add(self.recslider,1,wx.EXPAND,0)
        self.SetSizer(vbox)
