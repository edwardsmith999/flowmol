import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.colors import colorConverter
from matplotlib.patches import Rectangle
from pylab import asarray

#Some simple functions to generate colours.
def pastel(colour, weight=2.4):
    """ Convert colour into a nice pastel shade"""
    rgb = asarray(colorConverter.to_rgb(colour))
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

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]


def pdf2eps(filename):
    #Strip file extention (if any)
    name = os.path.splitext(filename)[0]
    pdfilename = name + '.pdf'
    epsfilename = name + '.eps'
    pdfilenametemp = name + '_temp.pdf'
    cmd = 'pdfcrop ' + pdfilename + ' ' + pdfilenametemp
    os.system(cmd)
    cmd = 'pdftops -level3 -eps ' + pdfilenametemp +  ' ' +  epsfilename
    os.system(cmd)
    os.system('rm ' + pdfilenametemp)


def plot_CVconserveforce(outtype):

    assert(outtype in ['mass', 'momentum','energy'])

    a=np.genfromtxt('./' + outtype)
    advection = a[:,8]
    forcing = a[:,7]
    accumulation = a[:,9]
    value = a[:,11]

    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    delta_t = 0.005
    plotmin = 0
    plotmax = a.shape[0]
    x = np.linspace(plotmin*delta_t,plotmax*delta_t,(plotmax-plotmin))

    if (outtype.find('mass') == 0):
        gap = 5.0
    if (outtype.find('momentum') == 0):
        gap = 50.0
    elif(outtype.find('energy') == 0):
        gap = 50.0
    Accumulationshift = 0.0
    Forcingshift = Accumulationshift + gap
    Advectionshift = Forcingshift + gap
    AbsoluteValueshift = Accumulationshift - gap

    #ADVECTION
    plt.plot(x,advection[plotmin:plotmax]+Advectionshift,'-ro',alpha=0.5,label='Advection',markersize=4.5)
    plt.plot(x,np.ones(plotmax-plotmin)*Advectionshift, 'k',linewidth = 0.5) #ZERO LINE

    #FORCING
    plt.plot(x,forcing[plotmin:plotmax]+Forcingshift,'r',linewidth = 3.0,label='Forcing')
    plt.plot(x,np.ones(plotmax-plotmin)*Forcingshift, 'k',linewidth = 0.5) #ZERO LINE

    #ACCUMULATION
    plt.plot(x,accumulation[plotmin:plotmax]+Accumulationshift,'k',label='Accumulation')
    plt.fill_between(x,np.ones(plotmax-plotmin)*Accumulationshift,accumulation[plotmin:plotmax]+Accumulationshift,color='k',alpha=0.4)

    #MASS/MOMENTUM/ENERGY CHANGE
    if (outtype in ['mass']):
        factor = 0.5
        name = '(' + outtype.replace('_sin',' Sinusoid')+'-'+ str(np.around(np.mean(value[plotmin:plotmax]), decimals=3)) + r') $\times$ ' + str(factor)
        plt.plot(x,factor*(value[plotmin:plotmax]-np.mean(value[plotmin:plotmax]))+AbsoluteValueshift,':',color=[0.5,0.0,0.0],linewidth = 1.0,label='$\mathcal{CV}$ Mass')
    elif (outtype in ['momentum']):
        factor = 10
        plt.plot(x,factor*value[plotmin:plotmax]+AbsoluteValueshift,':',color=[0.5,0.0,0.0],linewidth = 1.0,label='$\mathcal{CV}$ Momentum')
        name = outtype.replace('_sin',' Sinusoid')+r' $\times$ ' + str(factor)
    elif (outtype in ['energy']):
        factor = 10
        plt.plot(x,factor*(value[plotmin:plotmax]-np.mean(value[plotmin:plotmax]))+AbsoluteValueshift,':',color=[0.5,0.0,0.0],linewidth = 1.0,label='$\mathcal{CV}$ Energy')
        name = '(' + outtype.replace('_sin',' Sinusoid')+'-'+ str(np.around(np.mean(value[plotmin:plotmax]), decimals=2)) + r') $\times$ ' + str(factor)

    fsize = 16
    plt.text((plotmin+2.0)*delta_t,Advectionshift+gap/2.0,r'Advection',fontsize=fsize)
    plt.text((plotmin+2.0)*delta_t,Forcingshift+1.5,r'Forcing',fontsize=fsize)
    txt = plt.text((plotmin+2.0)*delta_t,Accumulationshift+1.8,r'Accumulation',fontsize=fsize)
    plt.text((plotmin+2.0)*delta_t,AbsoluteValueshift+2.0,name,fontsize=fsize)

    plt.setp(ax, yticks=np.linspace(-1.5*gap,4*gap,12), yticklabels=[ str(-gap/2.),'0.0','$\pm$'+str(gap/2.), '0', '$\pm$'+str(gap/2.),
                                                                 '0', '$\pm$'+str(gap/2.),
                                                                '0', '$\pm$'+str(gap/2.),'0', '$\pm$'+str(gap/2.),str(gap)])

    plt.ylim(-1.5*gap,3*gap)
    plt.xlim(plotmin*delta_t,plotmax*delta_t)

    plt.rc('font', size=fsize)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rc('text', usetex=True)

    plt.xlabel('Time',fontsize=fsize)
    plt.ylabel(outtype.replace('_sin',' Sinusoid').capitalize(),fontsize=fsize)

    filename = '../CV' + outtype + '.pdf'
    plt.savefig(filename)
    plt.clf()
    pdf2eps(filename)


if len(sys.argv) == 1:
    types = ['mass','momentum','energy']
else:
    outtype = sys.argv[1]
    types = [outtype]
    print('Using supplied type ' + outtype)

filenames=[]
for outtype in types:
    plot_CVconserveforce(outtype)
    filenames.append('../CV' + outtype + '.pdf')

os.system('acroread ' + ' '.join(filenames))



