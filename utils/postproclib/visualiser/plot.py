import wx
import numpy as np
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.backends.backend_wxagg as wxaggb
import matplotlib.backends.backend_wx as wxb
import matplotlib.pyplot as plt

class PyplotPanel(wx.Panel):

    def __init__(self,parent,**kwargs):
        wx.Panel.__init__(self,parent,**kwargs)
        self.figure = matplotlib.figure.Figure()
        self.canvas = wxaggb.FigureCanvasWxAgg(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.Bind(wx.EVT_SIZE, self.sizeHandler)

        #self.toolbar = wxb.NavigationToolbar2Wx(self.canvas) 
        #self.toolbar.Realize()
        #self.toolbar.update()

        self.cmap = matplotlib.cm.RdYlBu_r

    def sizeHandler(self, event):
        self.canvas.SetSize(self.GetSize())

    def redraw_plot(self, ax, data, xlabel=None, ylabel=None):
        self.figure.clf(keep_observers=True)
        self.ax = self.figure.add_subplot(111)
        self.lines = self.ax.plot(ax, data, 'r-o', linewidth=2)
        self.ax.set_xlim(ax.min(), ax.max())
        if (xlabel): self.ax.set_xlabel(xlabel)
        if (ylabel): self.ax.set_ylabel(ylabel)
        self.canvas.draw()

    def redraw_plot_many(self, axs, datas, styles=None, xlabel=None, ylabel=None):
        self.figure.clf(keep_observers=True)
        self.ax = self.figure.add_subplot(111)

        if (styles == None):
            styles = [{}]*len(axs)

        self.lines = [] 
        for ax,data,style in zip(axs,datas,styles):
            line = self.ax.plot(ax, data, **style)
            self.lines.append(line)

        # Maximum and minimum grid values
        maxval = np.max([np.max(ax) for ax in axs])
        minval = np.min([np.min(ax) for ax in axs])
        self.ax.set_xlim(minval,maxval)

        if (xlabel): self.ax.set_xlabel(xlabel)
        if (ylabel): self.ax.set_ylabel(ylabel)
        self.canvas.draw()
   
    def update_plot_many(self, axs, datas):
        for line, ax, data in zip(self.lines, axs, datas):
            plt.setp(line, xdata=ax, ydata=data)
        self.canvas.draw() 
    
    def update_plot(self, ax, data):
        plt.setp(self.lines, xdata=ax, ydata=data)
        self.canvas.draw()

    def redraw_contour(self, ax1, ax2, data, xlabel=None, ylabel=None):
        self.figure.clf(keep_observers=True)
        self.ax = self.figure.add_subplot(111)
        self.colormesh = self.ax.pcolormesh(ax1, ax2, data[:,1:], cmap=self.cmap)
        self.cbar = self.figure.colorbar(self.colormesh)
        self.ax.set_xlim(ax1.min(), ax1.max())
        self.ax.set_ylim(ax2.min(), ax2.max())
        if (xlabel): self.ax.set_xlabel(xlabel)
        if (ylabel): self.ax.set_ylabel(ylabel)
        self.canvas.draw()

    def update_contour(self, data):
        self.colormesh.set_array(data[:,1:].ravel())          
        self.canvas.draw()

    def savefigure(self,fpath):
        self.figure.savefig(str(fpath),dpi=300)
