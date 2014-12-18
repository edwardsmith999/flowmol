import numpy as np
import matplotlib.pyplot as plt

def zero_one_test(x,Nrandnos=100,alpha=2.5):

    """
        zero_one_test(X) applies the Gottwald-Melbourne 0-1 test
        to time history vector X. Result should be 0.0 for 
        non-chaotic data and 1.0 for chaotic data.
        
    """

    #Check it is a numpy array
    if type(x).__module__ == np.__name__:
        print('Numpy array in z1test')
    else:
        print('Warning, not numpy array')

    N=len(x)
    #From Detection of low-dimensional chaos in quasi-periodic
    # time series: The 0-1 test by Francois Ascani, Paul Keeler,
    #Ruben Kubiak, Samuel Laney, Srideep Musuvathy, David Papo
    #It is expected that n << N so not all data points are used in
    # the calculation. In practice, M is computed over roughly the 
    #first 10% of the time series
    j = np.arange(0,N); n= np.round(N/10)
    t = np.arange(0,n); M = np.zeros(n)
    D = np.zeros(n)
    kcorr = np.empty(Nrandnos)
    chist = np.empty(Nrandnos)

    # To address noise in data, it is best to 
    #Loop over a range of random c values in range [pi/5,4pi/5]
    for its in range(0,Nrandnos):
        c=np.pi/5.+np.random.rand()*3.*np.pi/5.
        chist[its] = c
        # Get p and q, the cumulative sum of x times sin/cosine  
        p=np.cumsum(x*np.cos(j*c))
        q=np.cumsum(x*np.sin(j*c))
        # Calculate mean square displacement (MSD) of p and q
        # In a chaotic system, the MSD grows faster than the number 
        # of data points
        Expectedx = -np.power(np.mean(x),2.)
        for n_ in range(0,n): 
            #MSD of p and q
            M[n_] = ( np.mean( np.power((p[n_:N]-p[0:N-n_]),2) 
                              +np.power((q[n_:N]-q[0:N-n_]),2) ))
            #Damps the expression reducing susceptibility to noise 
            D[n_] = M[n_]-alpha*Expectedx*np.sin(np.sqrt(2.)*n_)
            #D[n_]=M[n_]-alpha*Expectedx*(1.-np.cos(n_*c))/(1.-np.cos(c)))

        kcorr[its]=np.corrcoef(t,D)[0,1]


    fig, ax = plt.subplots(2,sharex=False)
    ax[0].plot(chist,kcorr,'x')
    ax[1].plot(t,M)
    plt.show()

    return np.median(kcorr)


#temp = np.arange(0,1000)
#print(zero_one_test(temp))

#temp = np.random.rand(1000)
#print(zero_one_test(temp))
