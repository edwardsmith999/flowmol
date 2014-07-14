import wx
import sys
import postproclib.visualiser as pplv

if __name__ == "__main__":

    if (len(sys.argv) > 1):
        direc = sys.argv[1]
    else:
        direc = '../MD_dCSE/src_code/results/'

    app = wx.App()
    fr = pplv.MainFrame(None, fdir=direc)
    fr.Show()
    app.MainLoop()
