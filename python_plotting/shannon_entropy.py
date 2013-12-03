import numpy as np
import sys

def shannon_entropy(pmf_samp):
    """
    Shannon entropy of a probability mass function.  Needs pmf_samp to be normalized such that Sum(pmf_samp)=1.0
    """

    #Bad error handling
    if np.sum(pmf_samp)!=1.0:
        sys.exit('pmf_samp not normalized')

    with np.errstate(divide='ignore',invalid='ignore'):
        #Uses masked array (ma), which eventually gives NaN, which is ignored in nansum below
        log_contrib = np.ma.log2(pmf_samp)
        ent_contrib = - pmf_samp * np.log2(pmf_samp)

    return np.nansum(ent_contrib)

#Program
if __name__=='__main__':
    #Cosmology data
    #data=np.loadtxt('nsralpha.txt')
    #
    #ns=data[:,0]
    #r=data[:,1]
    #alpha=1000*data[:,2]

    #Test data
    numb=10000
    ns = np.random.randn(numb)
    r = np.random.randn(numb)
    alpha = np.random.randn(numb)



    #Create normalized PMF
    bins = 20
    norm=False

    nsr, xedges1, yedges1 = np.histogram2d(ns, r, bins=bins, normed=norm)
    pmf_samp = nsr/len(ns)

    #pmf_samp = np.array([[1,2],[3,4]])/10.0

    print shannon_entropy(pmf_samp)
