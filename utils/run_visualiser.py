import wx
import sys
import postproclib.visualiser as pplv
import argparse

def run_visualiser(parent_parser=argparse.ArgumentParser(add_help=False)):

    #Keyword arguments
    parser = argparse.ArgumentParser(description="""Runs visualiser MD XXXX where XXXX is an 
                                                  increasingly more futuristic and exciting number""",
                                     parents=[parent_parser])
    parser.add_argument('-f', '--fdir', dest='fdir', help='Directory containing results', default='../MD_dCSE/src_code/results/')
    args = vars(parser.parse_args())

    app = wx.App()
    fr = pplv.MainFrame(None, fdir=args['fdir'])
    fr.Show()
    app.MainLoop()

if __name__ == "__main__":

    run_visualiser()
