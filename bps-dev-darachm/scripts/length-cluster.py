#!/usr/bin/env python3

import argparse
import sklearn as sk
import os
import pathlib
import textwrap
from Bio import SeqIO
import numpy as np
from sklearn.mixture import GaussianMixture
from matplotlib import pyplot as plt
import multiprocessing
import scipy 
from scipy.optimize import minimize

parser = argparse.ArgumentParser('''
    thing to split fastq files into clusters by read length by GMM for now
    ''')
parser.add_argument('input_dir')
parser.add_argument('--output-dir',default='./')
parser.add_argument('--plots-dir',default='./')
parser.add_argument('--delim',default='---',type=str)
parser.add_argument('--min-reads',default=5,type=int)
parser.add_argument('--max-clusters',default=4,type=int)
parser.add_argument('--keep-tail',action="store_true")
parser.add_argument('--threads',default=1,type=int)
args = parser.parse_args()

# could re-write this as a minibatch, maybe in TF for speeeeeed

if args.max_clusters > 5:
    print("hey can't handle that many clusters, shorten it up")

class readLengthUnmixer():
    
    def __init__(self,components=1,scaling_factor=1e4,maxiter=5e2):
        self.components = components
        self.clusters = list(range(1,components+1))
        self.scaling_factor = scaling_factor
        self.maxiter = maxiter
        # [fractions], [sizes], [sds]
        self.bounds = ( 
            [(1e-6, 5e-1)] + [(1e-9, 1)] * components +
            [ (5e2 / scaling_factor, 1e5 / scaling_factor) ] * (components+1) + 
            [(1 / scaling_factor, 2e2 / scaling_factor )] * (components)
            )
        
    def fit(self,x):
        
        self.x = x
        self.x_as_zeros = np.zeros(x.shape)
        
        def calc_pdfs(fractions,sizes,sds,fit_x):
            #print(fractions,sizes,sds)
            distz = [
                    scipy.stats.uniform.pdf( fit_x, 
                        #loc=np.min(fit_x)-1, scale=sizes[0] )
                        loc=1e3, scale=sizes[0] )
                    ] + [
                    scipy.stats.norm.pdf( fit_x,
                        loc=sizes[i+1], scale=sds[i] )
                    for i in range(self.components) 
                    ]
            prob_together = np.array(
                    [ distz[i] * fractions[i] 
                        for i in range(len(distz)) ]
                    )
            return np.apply_along_axis(lambda y: np.choose( np.isfinite(y),
                        [np.zeros(len(y)), y ]) , 0, prob_together)
        self.calc_pdfs = calc_pdfs

        def calc_ll(fractions,sizes,sds,fit_x):
            pdfs = calc_pdfs(fractions,sizes,sds,fit_x)
            log_pdfs = np.log( np.apply_along_axis(lambda y: np.sum( y ), 
                                    0, pdfs) )
            return np.sum(
                    np.choose( np.isfinite(log_pdfs),
                        [np.zeros(len(log_pdfs)), log_pdfs ]
                        ) )

        def fitting_func(params,fit_x=x):
            # [fractions], [sizes], [sds]
            fractions = params[0:(self.components+1)]
            sizes = params[(self.components+1):(2*(self.components+1))] * self.scaling_factor
            sds = params[(2*(self.components+1)):] * self.scaling_factor
            return -calc_ll(fractions,sizes,sds,fit_x)
        self.fitting_func = fitting_func

        # stub of idea to take the biggest peaks, but the noise...
        #lenz_hist, lenz_bins = np.histogram(lenz_array,20,density=True)
        #print( sorted( zip(lenz_hist, lenz_bins) , key=lambda x: x[0], reverse=True ) )

        fracs0 = np.linspace(0.5,1,self.components+1)
        x_mode = scipy.stats.mode(x,keepdims=False)[0]
        x0 = (  [ i for i in fracs0 / np.sum(fracs0) ] + 
                list( x_mode / self.scaling_factor * 0.95 ) + 
                list( x_mode / self.scaling_factor * 
                    [   np.array([1]) ,
                        np.array([1,2]) ,
                        np.array([0.75,1,2]) ,
                        np.array([0.75,1,1.5,2]) ,
                        np.array([0.75,1,1.5,2,3]) ,
                        np.array([0.75,1,1.5,2,3,4]) 
                        ][self.components-1]
                    ) + 
                [5e1 / self.scaling_factor] * self.components
                #list( np.linspace(1e1,1e2,self.components) / 
                #            self.scaling_factor )
                )
        #print(x0)

        results = minimize(fitting_func, x0, method='trust-constr',
                bounds=self.bounds,
                options={'maxiter': self.maxiter, 'verbose': 0} ,
                constraints=[
                        {'type': 'eq', 'fun': 
                            lambda z: np.sum( z[0:(1+self.components)] ) - 1 } 
                        ] )
        self.fit_params = results.x
        self.scaled_params = np.concatenate( [
                results.x[0:(1+self.components)] ,
                results.x[(self.components+1):(2*(self.components+1))] * 
                    self.scaling_factor ,
                results.x[(2*(self.components+1)):] * self.scaling_factor
                ] )
        self.scaled_params = np.choose( np.less(0,self.scaled_params), 
                                [np.zeros(len(self.scaled_params)),self.scaled_params] )

        def prob_func(x):
            fractions = self.scaled_params[0:(self.components+1)]
            sizes = self.scaled_params[(self.components+1):(2*(self.components+1))]
            sds = self.scaled_params[(2*(self.components+1)):] 
            together = calc_pdfs(fractions,sizes,sds,fit_x=x)
            return together
        self.prob_func = prob_func

    def prob_sum(self,x):
        return np.apply_along_axis(lambda y: np.sum(y), 0, self.prob_func(x))
        
    def predict(self,x):
        probs = self.prob_func(x)
        if not args.keep_tail:
            probs[0,:] = 0
        return np.apply_along_axis(
                lambda y: list(y).index(np.max(y)),
                0, probs )
    
    def aic(self,x):
        ll = -self.fitting_func(self.fit_params,fit_x=x)
        return 2*len(self.fit_params) - 2 * ll

    def bic(self,x):
        ll = -self.fitting_func(self.fit_params,fit_x=x)
        return len(self.fit_params)*np.log(len(x)) - 2 * ll




def proc_the_file(filename):

    with open(filename) as f:

        basename = filename.stem
        
        recordz = list(SeqIO.parse(f,'fastq'))
        lenz = [ len(j.seq) for j in recordz ]
        # store lengths of all records, and reshape as numpy array ready for fit
        lenz_array = np.array([lenz]).reshape(-1,1)

        models = []
        for n in range(1,args.max_clusters+1):

            if len(lenz) / args.min_reads < n :
                # not enough reads, return
                return

            mixer = readLengthUnmixer(n,scaling_factor=1e0,maxiter=1e2)
            mixer.fit(lenz_array)
            the_aic = mixer.aic(lenz_array)

            print(basename)
            [ print(z) for z in zip( 
                    ["tail frac"] + ["component frac"] * n + 
                        ["tail size"] + ["component size"] * n + 
                        ["component sd"] * n
                    , mixer.scaled_params) ]
            print(the_aic)
            print("")

            plot_x = np.linspace(min(lenz), max(lenz),1000)
            plot_y = mixer.prob_sum(plot_x)
            lenz_hist, lenz_bins = np.histogram(lenz_array,50,density=True)

            # pretty much copied from: 
            # https://www.astroml.org/book_figures/chapter4/fig_GMM_1D.html
            fig = plt.figure(figsize=(9,5),layout='constrained')
            
            ax = fig.add_subplot(121)
            ax.plot(lenz_bins[1:],lenz_hist)
            ax.plot(plot_x, plot_y, '--k')
            plt.ylim(0,max(lenz_hist))
            ax.set_title(basename+'\n'+str(n)+' components/clusters, '+
                    'AIC: '+','.join(textwrap.wrap(str(int(the_aic))[::-1],width=3))[::-1]
                    )
            ax.set_xlabel('read length')
            ax.set_ylabel('density observed (blue lin),\nor sum of fits (dashed lines)')
            
            ax = fig.add_subplot(122)
            ax.plot(lenz_bins[1:],lenz_hist)
            ax.plot(plot_x, plot_y, '--k')
            plt.ylim(1e-6,max(lenz_hist))
            ax.set_yscale('log')
            ax.set_xlabel('read length')
            ax.set_ylabel('density observed (blue lin),\nor sum of fits (dashed lines)')
            
            plt.savefig(args.plots_dir+'/'+basename+'_'+str(n)+'clusters.jpeg')
            plt.close()

            models.append((n,the_aic,mixer,fig))
            
        
        try:

            this_min_aic = min([ aic for (n,aic,mod,plot) in models if n > 1])
            # pull out the minimum aic, then select the model with it
            n_and_mod = [ (n,mod) for (n,aic,mod,plot) in models if aic == this_min_aic ][0]
            # use it to label all the records
            labelz = n_and_mod[1].predict(lenz_array)

            output_bins = [list()]*(n_and_mod[0])
            # for each cluster, store the records/seqs that are in that cluster
            for j in range(n_and_mod[0]):
                output_bins[j] = [ z[1] for z in zip(labelz,recordz) if z[0] == j]

            # for each cluster, write those out as new file
            for j in range(n_and_mod[0]):
                SeqIO.write(output_bins[j], 
                    args.output_dir+'/'+basename+args.delim+'cluster'+
                        str(j)+'.fastq', 'fastq'
                    )

            return [ (n,aic,mod.scaled_params,plot) for (n,aic,mod,plot) in models if aic == this_min_aic ][0]

        except:
            return 'fail'


with multiprocessing.Pool(args.threads) as p:
    p.map(proc_the_file,pathlib.Path(args.input_dir).iterdir())
#proc_the_file(list(pathlib.Path(args.input_dir).iterdir())[0])
#this_mod = proc_the_file(pathlib.Path('examples/bc_bc_5---CTATAAAACAATTTAAGAT---donor_barcode_148.fastq'))
#this_mod = proc_the_file(pathlib.Path('togo/bc_bc_10---ACATGAAACCTGTTGGCGA---donor_barcode_47.fastq'))

