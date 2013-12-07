import numpy as np
import sys
import scipy.stats as stats

def shannon_entropy(p_samp,q_samp=None,base=2):
    """
    Shannon entropy of a probability mass function.  If q_samp is supplied, then calculates the Kullback-Leibler divergence.  Needs p_samp to be normalized such that Sum(p_samp)=1.0
    """

    def normed(samp, tol=1e-12):
        return np.abs(1.0-np.sum(samp)) < tol

    p_samp.flatten()

    if base != 2:
        base_convert = 1.0/np.log2(base)
    else:
        base_convert = 1.0

    if not normed(p_samp) or (q_samp != None and not normed(q_samp)):
        #Replace by raising valueerror
        print 'Sample not normalized.'
        print 'Sum(P) =', np.sum(p_samp)
        if q_samp != None: print 'Sum(Q) =', np.sum(q_samp)
        sys.exit()

    if q_samp!= None and len(q_samp) != len(p_samp):
        print 'Q sample and P sample not of same length.'
        print 'len(P) =', len(p_samp)
        print 'len(Q) =', len(q_samp)
        sys.exit()

    with np.errstate(divide='ignore',invalid='ignore'):
        #Uses masked array (ma), which eventually gives NaN for log(0)
        #which is ignored in nansum below

        if q_samp == None:
            #Shannon entropy: H = - \Sum P log P
            log_contrib = np.ma.log2(p_samp)
        else:
            #Kullback-Leibler: D = \Sum P log (P/Q)
            #NB: Picks up extra minus
            log_contrib = - np.ma.log2(p_samp/q_samp)

        ent_contrib = - p_samp * log_contrib
        ent_contrib *= base_convert

    return np.nansum(ent_contrib)

def histogram_probability(data,bins):

    numb = float(len(data))

    hist, a = np.histogramdd(data, bins=bins)

    return hist/numb


def resample_with_repl(array,numb):

    return array[np.random.randint(array.shape[0],size=numb),:]

#def scott_binning(samp):


nats_to_bits = 1.0/np.log(2.0)


#Program
if __name__=='__main__':

    import matplotlib.pyplot as plt

    #Cosmology data
    data=np.loadtxt('nsralpha.txt')

    ns=data[:,0]
    r=data[:,1]
    alpha=1000*data[:,2]



    nbins = 10


    #Test data

    #ns = np.random.randn(numb)
    #r = np.random.randn(numb)
    #alpha = np.random.randn(numb)
    #mu, sigma = 100, 15
    #x = mu + sigma * np.random.randn(numb)


    #entropy= [ 'Sample size', 'Bin #', 'Entropy 1D ns', 'Entropy 1D r', 'Entropy 2D ns r']
    entropy= []

    for sampsize in np.arange(10000,len(data),1000):

        for bin_numb in np.arange(10,1010,100):


            print 'sampsize', sampsize, 'bin_numb', bin_numb

            ns=data[:sampsize,0]
            r=data[:sampsize,1]
            alpha=1000*data[:sampsize,2]

            hist_1dr = histogram_probability(r,bin_numb)
            hist_1dns = histogram_probability(ns,bin_numb)
            hist_1dalpha = histogram_probability(alpha,bin_numb)
            hist_nsr = histogram_probability(np.column_stack((ns,r)),bin_numb)
            hist_nsalpha = histogram_probability(np.column_stack((ns,alpha)),bin_numb)

            rand = 100*np.random.rand(len(ns))
            hist_rand = histogram_probability(rand,bin_numb)


            #print 'Mine 1D: ', shannon_entropy(hist_1dns,hist_1dr,base=np.e)
            #print 'Mine 2D: ', shannon_entropy(hist_nsalpha)

            entropy.append([
                sampsize,
                bin_numb,
                shannon_entropy(hist_1dns),
                shannon_entropy(hist_1dr),
                shannon_entropy(hist_1dalpha),
                shannon_entropy(hist_nsr),
                shannon_entropy(hist_nsalpha)
                ])


    np.array(entropy)
    np.savetxt('entropy.txt',entropy)
    print 'Entropy printed.'
