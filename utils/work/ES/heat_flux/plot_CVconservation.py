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

    assert(outtype in ['momentum','momentum_sin','energy','energy_sin'])

    a=np.genfromtxt('./' + outtype)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    delta_t = 0.005
    plotmin = 0
    plotmax = 400
    x = np.linspace(plotmin*delta_t,plotmax*delta_t,(plotmax-plotmin))

    if (outtype.find('momentum') == 0):
        gap = 6.0
    elif(outtype.find('energy') == 0):
        gap = 6.0
    Accumulationshift = 0.0
    Constraintshift = Accumulationshift+gap
    Forcingshift = Constraintshift + gap
    Advectionshift = Forcingshift + gap
    AbsoluteValueshift = Accumulationshift - gap

    #ADVECTION
    plt.plot(x,a[plotmin:plotmax,8]+Advectionshift,'-ro',alpha=0.5,label='Advection',markersize=4.5)
    plt.plot(x,np.ones(plotmax-plotmin)*Advectionshift, 'k',linewidth = 0.5) #ZERO LINE

    #FORCING
    plt.plot(x,a[plotmin:plotmax,7]+Forcingshift,'r',linewidth = 3.0,label='Forcing')
    plt.plot(x,np.ones(plotmax-plotmin)*Forcingshift, 'k',linewidth = 0.5) #ZERO LINE

    #ACCUMULATION
    plt.plot(x,a[plotmin:plotmax,9]+Accumulationshift,'k',label='Accumulation')
    plt.fill_between(x,np.ones(plotmax-plotmin)*Accumulationshift,a[plotmin:plotmax,9]+Accumulationshift,color='k',alpha=0.4)

    #CV CONSTRAINT
    plt.plot(x,a[plotmin:plotmax,10]+Constraintshift,'b--',linewidth = 1.0,label='$\mathcal{CV}$ Constraint')
    plt.fill_between(x,Constraintshift,a[plotmin:plotmax,10]+Constraintshift,alpha=0.2)
    plt.plot(x,np.ones(plotmax-plotmin)*Constraintshift, 'k',linewidth = 0.5) #ZERO LINE

    #MOMENTUM/ENERGY CHANGE
    if (outtype in ['momentum','momentum_sin']):
        factor = 10
        plt.plot(x,factor*a[plotmin:plotmax,11]+AbsoluteValueshift,':',color=[0.5,0.0,0.0],linewidth = 1.0,label='$\mathcal{CV}$ Momentum')
        name = outtype.replace('_sin',' Sinusoid')+r' $\times$ ' + str(factor)
    elif (outtype in ['energy','energy_sin']):
        factor = 10
        plt.plot(x,factor*(a[plotmin:plotmax,11]-np.mean(a[plotmin:plotmax,11]))+AbsoluteValueshift,':',color=[0.5,0.0,0.0],linewidth = 1.0,label='$\mathcal{CV}$ Energy')
        name = '(' + outtype.replace('_sin',' Sinusoid')+'-'+ str(np.around(np.mean(a[plotmin:plotmax,11]), decimals=2)) + r') $\times$ ' + str(factor)
    #plt.text((plotmin+(plotmax-plotmin)/2.0)*delta_t,Advectionshift-gap/2.0,r'$+$',fontsize=24)
    #plt.text((plotmin+(plotmax-plotmin)/2.0)*delta_t,Forcingshift-gap/2.0,r'$-$',fontsize=24)
    #plt.text((plotmin+(plotmax-plotmin)/2.0)*delta_t,Constraintshift-gap/2.0,r'$=$',fontsize=24)
    #lg = plt.legend(loc=2, borderaxespad=0.0, ncol=3)
    #lg.draw_frame(False)

    fsize = 16
    plt.text((plotmin+2.0)*delta_t,Advectionshift+gap/2.0,r'Advection',fontsize=fsize)
    plt.text((plotmin+2.0)*delta_t,Forcingshift+1.5,r'Forcing',fontsize=fsize)
    plt.text((plotmin+2.0)*delta_t,Constraintshift+0.5,r'$\mathcal{CV}$ Constraint',fontsize=fsize)
    txt = plt.text((plotmin+2.0)*delta_t,Accumulationshift+1.8,r'Accumulation',fontsize=fsize)
    #box = txt.get_clip_box(); Rectangle(box.extents,box.size[0],box.size[1], angle=0.0, fc ='w')
    #txt = plt.text((plotmin+2.0)*delta_t,Accumulationshift+1.8,r'Accumulation',fontsize=fsize)
    plt.text((plotmin+2.0)*delta_t,AbsoluteValueshift+2.0,name,fontsize=fsize)



    plt.text((plotmin+plotmax)*delta_t/2.-0.095,AbsoluteValueshift-2.5,r'$\underbrace{\phantom{a + b \; }}_{Sigmoid}$',fontsize=fsize)

    #PLOT ALL ON SAME LINE
    #plt.plot(x,zero_to_nan(a[plotmin:plotmax,8]),'ro',alpha=0.5,label='Advection',markersize=4.5)
    #plt.plot(x,a[plotmin:plotmax,7],'r',linewidth = 3.0,label='Forcing')
    #plt.plot(x,a[plotmin:plotmax,9],'k',label='Accumulation')
    #plt.fill_between(x,np.zeros(plotmax-plotmin),a[plotmin:plotmax,9],color='k',alpha=0.4)
    #plt.plot(x,a[plotmin:plotmax,10],'b--',linewidth = 1.0,label='Constraint')
    #plt.fill_between(x,0.0,a[plotmin:plotmax,10],alpha=0.2)
    #plt.plot(x,20*a[plotmin:plotmax,11]-5.0,':',color=[0.5,0.0,0.0],linewidth = 1.0,label='CV Momentum')


    #plt.setp(ax, yticks=np.linspace(-4.0,16.0,11))

    plt.setp(ax, yticks=np.linspace(-1.5*gap,4*gap,12), yticklabels=[ str(-gap/2.),'0.0','$\pm$'+str(gap/2.), '0', '$\pm$'+str(gap/2.),
                                                                 '0', '$\pm$'+str(gap/2.),
                                                                '0', '$\pm$'+str(gap/2.),'0', '$\pm$'+str(gap/2.),str(gap)])

    plt.ylim(-1.5*gap,4*gap)
    #plt.xlim(473.0,550.0)
    plt.xlim(plotmin*delta_t,plotmax*delta_t)
    #plt.xlim(300.0,700.0)

    plt.rc('font', size=fsize)
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    plt.rc('text', usetex=True)

    plt.xlabel('Time',fontsize=fsize)
    plt.ylabel(outtype.replace('_sin',' Sinusoid').capitalize(),fontsize=fsize)

    #ax.set_yscale('log')
    #filename = '../CV' + outtype + '.png'
    #plt.savefig(filename)
    #os.system('eog ' + filename)

    filename = '../CV' + outtype + '.pdf'
    plt.savefig(filename)
    plt.clf()
    pdf2eps(filename)
    #os.system('evince ' + filename)



if len(sys.argv) == 1:
    types = ['momentum','momentum_sin','energy','energy_sin']
else:
    outtype = sys.argv[1]
    types = [outtype]
    print('Using supplied type ' + outtype)

filenames=[]
for outtype in types:
    plot_CVconserveforce(outtype)
    filenames.append('../CV' + outtype + '.pdf')

os.system('acroread ' + ' '.join(filenames))



