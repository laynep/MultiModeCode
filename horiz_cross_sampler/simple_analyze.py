#!/usr/bin/python

"""Program to sample the horizon crossing surface."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import scipy.special as sp
import myplotting as my
import sys

def build_big_dataframe(fileroot):

    files = np.arange(0,100,1)

    data=[]
    for numb in files:
        data_file = pd.load(fileroot+'output'+str(numb)+'.pkl')
        #data_file['n_t/r'] = pd.Series(data_file['n_t']/(data_file['p'][0]*4.0/55.0))
        data_file['n_t/r'] = pd.Series(data_file['n_t']/data_file['r'])
        data.append(data_file)
    data = pd.concat(data)
    data.save(fileroot+'full_output.pkl')
    print data


def load_data(fileroot):

    data = pd.load(fileroot+'full_output.pkl')

    #print data['nfields']
    #data.sort('nfields')

    #print data['nfields']

    return data

def make_plot():

    my.seaborn_style('dark')

    fig = my.make_fig("PRL")

    ax1 = fig.add_subplot(111)

    return fig, ax1

def make_scatter(data):

    plt.scatter(np.array(data['nfields']),
            np.array(data['n_t/r']))

    plt.ylim(-0.05,-0.025)
    plt.xlim(10,1000)

    #plt.show()


def make_hist(data):

    nbins = [40, 100]

    #datarange = [[0,1000], [-0.05, -0.025] ]
    datarange = [[0,1000], [-0.45, -0.2] ]


    counts, xedges1, yedges1 = np.histogram2d(
            np.array(data['nfields']),
            np.array(data['n_t/r']),
            range=datarange,
            bins=nbins,
            normed=True)


    extent= [xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]]

    plt.imshow(counts.T, origin='lower', extent=extent, aspect='auto', interpolation='nearest')

    #plt.show()


def make_kde(data):

    new_data = np.vstack([np.array(data['nfields']),np.array(data['n_t/r'])])

    kernel = stats.gaussian_kde(new_data)

    datarange = [[0,1100], [-0.45, -0.2] ]
    #datarange = [[0,1100], [-5, 5] ]
    #datarange = [[0,1100], [-2, -0.2] ]

    x, y = np.mgrid[
            datarange[0][0]:datarange[0][1]:20j,
            datarange[1][0]:datarange[1][1]:20j]


    positions = np.vstack([x.ravel(), y.ravel()])

    f = np.reshape(kernel(positions).T, x.shape)

    extent = [datarange[0][0], datarange[0][1], datarange[1][0], datarange[1][1]]

    #color_map = plt.get_cmap('RdYlBu_r')
    #color_map = plt.get_cmap('Blues_r')
    #color_map = plt.get_cmap('hot_r')
    color_map = plt.get_cmap('coolwarm')
    #color_map = plt.get_cmap('Blues')

    plt.imshow(f.T, origin='lower', extent=extent, aspect='auto',
            cmap=color_map,
            alpha=1.0)

    #plt.show()


def get_ranges_log(pval,low,high,nfields):

    #Intermediate var
    zed = np.sqrt(2.0)*sp.erfinv(2.0*pval-1.0)

    #Correlation
    #gamma = -0.978941
    gamma = -0.95

    #Moments of log dist
    def mu_n(n):
        mu = (1.0/np.log(10.0)/(high-low)/nfields)*(10**(nfields*high)-10**(nfields*low))
        return mu

    mu_1 = mu_n(1.0)
    mu_2 = mu_n(2.0)
    mu_4 = mu_n(4.0)

    mean = mu_1
    sig1 = np.sqrt(mu_2 - mu_1**2)
    sig2 = np.sqrt(mu_4 - mu_2**2)

    mean_num = -(1.0/8.0)*nfields**2*(mu_2)
    std_num = (1.0/8.0)*nfields**(3.0/2.0)*sig2
    mean_denom = nfields**2*mean**2
    std_denom = 2.0*nfields**(3.0/2.0)*mean*sig1

    A = zed**2*std_denom**2 - mean_denom**2
    B = 2.0*mean_num*mean_denom - 2.0*zed**2*std_num*std_denom*gamma
    C = zed**2*std_num**2 - mean_num**2

    discriminant = B**2 - 4.0*A*C
    #discriminant = 0.0
    if discriminant < 0.0:
        if all(abs(discriminant*np.array([1.0/A,1.0/B,1.0/C]))) < 1e-6 :
            discriminant = 0.0
        else:
            #DEBUG
            ans = [-1.0/8.0]
            return ans
            raise Exception("Error: discriminant = %s" %discriminant)

    ans = [(-B + np.sqrt(discriminant))/(2.0*A),(-B - np.sqrt(discriminant))/(2.0*A)]
    #ans = (-B + np.sqrt(discriminant))/(2.0*A)

    return ans

def make_pred():

    #Unif pred --- mean (many) and value (single)

    #high = 1.0e-13
    #low = 1.0e-14

    high = -12
    low = -14

    x = np.linspace(2,1500,200)
    #y= (-1.0/8.0)*np.ones(200)
    y2= (-1.0/16.0)*np.log(10.0)*(high-low)*((10.0**high + 10.0**low)/\
            (10.0**high-10.0**low))*np.ones(200)

    pval = 0.25
    #pval = (1.0 - 0.68)/2.0
    y4 = [min(get_ranges_log(pval,low,high,me)) for me in x]
    pval = 0.75
    #pval = 1.0 -(1.0-0.68)/2.0
    y3 = [max(get_ranges_log(pval,low,high,me)) for me in x]

    #pval = 0.25 - 0.25/2.0
    pval = 0.0027/2.0
    #pval = 0.0455/2.0
    y6 = [min(get_ranges_log(pval,low,high,me)) for me in x]
    #pval = 0.75 + 0.25/2.0
    pval = 1.0 - 0.0027/2.0
    #pval = 1.0 - 0.0455/2.0
    y5 = [max(get_ranges_log(pval,low,high,me)) for me in x]

    #siglow_pred, = plt.plot(x,y3, 'k--', linewidth=1.0, alpha=0.75)
    #sighigh_pred, = plt.plot(x,y4,'k--', linewidth=1.0, alpha=0.75)


    #sf_pred, = plt.plot(x,y, linestyle='--', color='Chocolate', linewidth=1.0, alpha=0.70,zorder=1)
    mf_pred, = plt.plot(x,y2, 'k-', linewidth=2.0, alpha=0.50,zorder=1)

    #ax1.fill_between(x,y4,y3,alpha=0.20,edgecolor='b',lw=0.0)
    #ax1.fill_between(x,y6,y5,alpha=0.10,edgecolor='b',lw=0.0)



def master():

    #fileroot='./data/pandas/log_p1.5/'
    #fileroot='./data/pandas/log_p0.66/'
    #fileroot='./data/pandas/log_p2/'
    #fileroot='./data/pandas/log_p1.0/'
    #fileroot='./data/pandas/log_p4.0/'
    #fileroot='./data/pandas/log_p8.0/'

    fileroot='./data/pandas/unif_p0.66/'

    build_big_dataframe(fileroot)
    data = load_data(fileroot)
    fig, ax1 = make_plot()

    make_kde(data)
    #make_hist(data)
    #make_pred()

    plt.xlim(0,1100)
    plt.show()





if __name__=="__main__":
    master()
