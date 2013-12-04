import numpy as np
import sys

def shannon_entropy(pmf_samp):
    """
    Shannon entropy of a probability mass function.  Needs pmf_samp to be normalized such that Sum(pmf_samp)=1.0
    """

    pmf_samp.flatten()

    if np.abs(1.0-np.sum(pmf_samp)) > 1e-12:
        #Bad error handling
        print 'pmf_samp not normalized. Sum =', np.sum(pmf_samp)
        sys.exit()

    with np.errstate(divide='ignore',invalid='ignore'):
        #Uses masked array (ma), which eventually gives NaN for log(0)
        #which is ignored in nansum below
        log_contrib = np.ma.log2(pmf_samp)
        ent_contrib = - pmf_samp * np.log2(pmf_samp)

    return np.nansum(ent_contrib)

def get_normed_hist(data,bins):

    numb = float(len(data))

    if len(data.shape)==1 or data.shape[1]==1:

        hist, bins = np.histogram(r, bins=bins)

        return hist/numb


    elif data.shape[1]==2:

        hist, xedges1, yedges1 = np.histogram2d(data[:,0], data[:,1],
                bins=bins, normed=False)

        return hist/numb


    elif data.shape[1]>2:
        print 'Can\'t get a normed histogram for more \n than 2 dimensions and d =',data.shape[1]
        sys.exit()



#Program
if __name__=='__main__':

    import matplotlib.pyplot as plt

    #Cosmology data
    data=np.loadtxt('nsralpha_100k.txt')

    print 'Done loading txt'

    numb=50000

    ns=data[0:numb,0]
    r=data[0:numb,1]
    alpha=1000*data[0:numb,2]


    print 'Numb of points: ', numb

    nbins = 200


    #Test data

    #ns = np.random.randn(numb)
    #r = np.random.randn(numb)
    #alpha = np.random.randn(numb)
    #mu, sigma = 100, 15
    #x = mu + sigma * np.random.randn(numb)

    #print np.arange(10,500,10)

    entropy_r = {}
    entropy_nsr = {}
    for size in np.arange(10,1010,100):

        hist1d = get_normed_hist(r,size)
        hist2d = get_normed_hist(np.column_stack((ns,r)),size)

        entropy_r[size]=shannon_entropy(hist1d)
        entropy_nsr[size]=shannon_entropy(hist2d)

        #print 'H 1D r: ', shannon_entropy(hist1d)
        #print 'H 2D ns-r: ', shannon_entropy(hist2d)

    print 'entropy_r', entropy_r
    print 'entropy_nsr', entropy_nsr
