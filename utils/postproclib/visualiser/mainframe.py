# /usr/bin/env python
import wx
import sys
from postproclib.visualiser import __path__ as pplvpath
from postproclib.postproc import NoResultsInDir

from visuals import VisualiserPanel
from directory import DirectoryChooserPanel

def showMessageDlg(msg, title='Information', style=wx.OK|wx.ICON_INFORMATION):
    """"""
    dlg = wx.MessageDialog(parent=None, message=msg, 
                           caption=title, style=style)
    dlg.ShowModal()
    dlg.Destroy()

class MainFrame(wx.Frame):

    def __init__(self, parent=None, fdir='./', title='DataViewer 5000', 
                 size=(800,600), **kwargs):

        wx.Frame.__init__(self,parent,title=title,size=size)
        try:
            _icon = wx.EmptyIcon()
            _icon.CopyFromBitmap(
                # postproclib.visualiser.__path__ (pplvpath) is a list
                wx.Bitmap(pplvpath[0]+"/logo.gif", wx.BITMAP_TYPE_ANY)
            )
            self.SetIcon(_icon)
        except IOError:
            print('Couldn\'t load icon')

        self.dirchooser = DirectoryChooserPanel(self, fdir)

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.dirchooser, 0, wx.EXPAND, 0)
        self.SetSizer(self.vbox)

        self.visualiserpanel = None
        self.new_visualiserpanel(fdir) 

        self.set_bindings()

    def set_bindings(self):

        self.Bind(wx.EVT_TEXT_ENTER, self.handle_chdir, 
                  self.dirchooser.textctrl)
        self.Bind(wx.EVT_BUTTON, self.fdir_dialogue, 
                  self.dirchooser.changebutton)

    def destroy_visualiserpanel(self):

        if (self.visualiserpanel != None):
            self.visualiserpanel.Destroy()
            self.visualiserpanel = None
        
    def new_visualiserpanel(self, fdir):

        self.destroy_visualiserpanel()
        try:
            newvp = VisualiserPanel(self, fdir)
        except IOError:
            raise
            #showMessageDlg('Invalid directory.')
        except NoResultsInDir:
            showMessageDlg('No results in this directory.')
        else:
            self.visualiserpanel = newvp 
            self.vbox.Add(self.visualiserpanel, 1, wx.EXPAND, 0)
            self.SetSizer(self.vbox)
            self.Layout()
            print('New visualiser file directory: ' + fdir)
    
    def handle_chdir(self, event):

        fdir = self.dirchooser.textctrl.GetValue()
        self.new_visualiserpanel(fdir)

    def fdir_dialogue(self, event):

        fdir = ""  # Use  folder as a flag
        dlg = wx.DirDialog(self, defaultPath ="../")
        if dlg.ShowModal() == wx.ID_OK:
            fdir = dlg.GetPath() + "/"
            dlg.SetPath(fdir)
        dlg.Destroy()  # best to do this sooner than later

        if fdir:
            self.dirchooser.textctrl.SetValue(fdir)
            event = wx.PyCommandEvent(wx.EVT_TEXT_ENTER.typeId, 
                                      self.dirchooser.textctrl.GetId())
            self.GetEventHandler().ProcessEvent(event)
