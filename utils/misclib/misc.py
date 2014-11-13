#! /usr/bin/env python
# Routines which are not specific to MD/CFD or CPL code
import os
import numpy as np
from matplotlib.colors import colorConverter
import latex2utf
from math import log10, floor
import math as maths

class Chdir:          
    """
       Wrapper to move from current directory to new directory
       and return when using with

       Example usage:

       with Chdir('./../'):
           os.system('./a.out')
    """
    def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        self.newPath = newPath

    def __enter__( self ):
        os.chdir(self.newPath)

    def __exit__( self, etype, value, traceback):
        os.chdir( self.savedPath )


#Some simple functions to generate colours.
def pastel(colour, weight=2.4):
    """ Convert colour into a nice pastel shade"""
    rgb = np.asarray(colorConverter.to_rgb(colour))
    # scale colour
    maxc = max(rgb)
    if maxc < 1.0 and maxc > 0:
        # scale colour
        scale = 1.0 / maxc
        rgb = rgb * scale
    # now decrease saturation
    total = sum(rgb)
    slack = 0
    for x in rgb:
        slack += 1.0 - x

    # want to increase weight from total to weight
    # pick x s.t.  slack * x == weight - total
    # x = (weight - total) / slack
    x = (weight - total) / slack

    rgb = [c + (x * (1.0-c)) for c in rgb]

    return rgb

def update_line(hl, new_data):
    hl.set_xdata(np.append(hl.get_xdata(), new_data))
    hl.set_ydata(np.append(hl.get_ydata(), new_data))
    plt.draw()

def get_colours(n):
    """ Return n pastel colours. """
    base = np.asarray([[1,0,0], [0,1,0], [0,0,1]])

    if n <= 3:
        return base[0:n]

    # how many new colours to we need to insert between
    # red and green and between green and blue?
    needed = (((n - 3) + 1) / 2, (n - 3) / 2)

    colours = []
    for start in (0, 1):
        for x in np.linspace(0, 1, needed[start]+2):
            colours.append((base[start] * (1.0 - x)) +
                           (base[start+1] * x))

    return [pastel(c) for c in colours[0:n]]

def round_to_n(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(maths.log10(x))
    tens = maths.pow(10, e - p + 1)
    n = maths.floor(x/tens)

    if n < maths.pow(10, p - 1):
        e = e -1
        tens = maths.pow(10, e - p+1)
        n = maths.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= maths.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)


def latextounicode(strings):

    if type(strings) is unicode:
        string = strings.encode('utf8')
        strings = strings.replace('rho','\xcf\x81')
    if type(strings) is str:
        strings = strings.replace('rho','\xcf\x81')
    elif type(strings) is list:
        for i, string in enumerate(strings):
            strings[i] = string.replace('rho','\xcf\x81')
            #latex2utf.latex2utf(string)

    return strings

def unicodetolatex(strings):

    if type(strings) is unicode:
        string = strings.encode('utf8')
        strings = string.replace('\xcf\x81','rho')
    if type(strings) is str:
        strings = string.replace('\xcf\x81','rho')
    elif type(strings) is list:
        for i, string in enumerate(strings):
            string = string.encode('utf8')
            strings[i] = string.replace('\xcf\x81','rho')

    return strings

