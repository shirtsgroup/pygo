#!/usr/bin/python2.4

import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import cPickle
import optparse
import plot_dG_solution
import pdb

def parse_args():
    parser = optparse.OptionParser(description='Plots dG, dS, dH of folding and adsorption')
    parser.add_option('--direc', dest='direc', help = 'directory of bootstraps')
    parser.add_option('-n', type='int', dest='n', help = 'number of bootstraps')
    parser.add_option('--all', action="store_true", default=False, help='show all bootstraps in plot')
    parser.add_option('--show', action="store_true", default=False, help='show plots')
    (options,args) = parser.parse_args()
    return options

def plot_avg_bootstrap(boot_direc, a, b):
    dG = []
    dS = []
    dH = []
    for i,direc in enumerate(boot_direc):
        f = open(direc,'rb')
        temp = cPickle.load(f)
        bin = cPickle.load(f)
        G = cPickle.load(f)
        dG_i = G[:,a] - G[:,b]
        dS_i = (dG_i[2::] - dG_i[0:-2])/-10
        dH_i = dG_i[1:-1] + temp[1:-1]*dS_i
        dG.append(dG_i)
        dS.append(dS_i)
        dH.append(dH_i)
    dG = numpy.array(dG)
    dS = numpy.array(dS)
    dH = numpy.array(dH)
    
    dG_boot = numpy.average(dG,axis=0)
    ddG_boot = numpy.std(dG,axis=0)#/float(len(boot_direc)**.5)
    dS_boot = numpy.average(dS,axis=0)
    ddS_boot = numpy.std(dS,axis=0)#/float((len(boot_direc)-2)**.5)
    dH_boot = numpy.average(dH,axis=0)
    ddH_boot = numpy.std(dH,axis=0)#/float((len(boot_direc)-2)**.5)
    
    plt.errorbar(temp,dG_boot,ddG_boot)
    
    return temp, dG_boot, ddG_boot, dS_boot, ddS_boot, dH_boot, ddH_boot

def rel_ads_bootstrap(boot_direc):
    dG = []
    dS = []
    dH = []
    for i,direc in enumerate(boot_direc):
        f = open(direc,'rb')
        temp = cPickle.load(f)
        bin = cPickle.load(f)
        G = cPickle.load(f)
        dG_i = (G[:,1] - G[:,3])-(G[:,0]-G[:,2])
        dS_i = (dG_i[2::] - dG_i[0:-2])/-10
        dH_i = dG_i[1:-1] + temp[1:-1]*dS_i
        dG.append(dG_i)
        dS.append(dS_i)
        dH.append(dH_i)
    dG = numpy.array(dG)
    dS = numpy.array(dS)
    dH = numpy.array(dH)
    
    dG_boot = numpy.average(dG,axis=0)
    ddG_boot = numpy.std(dG,axis=0)#/float(len(boot_direc)**.5)
    dS_boot = numpy.average(dS,axis=0)
    ddS_boot = numpy.std(dS,axis=0)#/float((len(boot_direc)-2)**.5)
    dH_boot = numpy.average(dH,axis=0)
    ddH_boot = numpy.std(dH,axis=0)#/float((len(boot_direc)-2)**.5)
    
    plt.errorbar(temp,dG_boot,ddG_boot)
    
    return temp, dG_boot, ddG_boot, dS_boot, ddS_boot, dH_boot, ddH_boot


def plot_all_bootstrap(boot_direc, a, b):
    colors = cm.cool(numpy.linspace(0,1,len(boot_direc)))
    for i,direc in enumerate(boot_direc):
        f = open(direc,'rb')
        temp = cPickle.load(f)
        bin = cPickle.load(f)
        dG = cPickle.load(f)
        ddG = cPickle.load(f)
        dGf = dG[:,a] - dG[:,b]
        ddGf = [numpy.sqrt(x[a,b]) for x in ddG]
        plt.errorbar(temp,dGf,ddGf,color=colors[i])

def plot_orginal(direc, a, b):
    f = open(direc,'rb')
    temp = cPickle.load(f)
    bin = cPickle.load(f)
    dG = cPickle.load(f)
    ddG = cPickle.load(f)
    dGf = dG[:,a] - dG[:,b]
    ddGf = [numpy.sqrt(x[a,b]) for x in ddG]
    plt.errorbar(temp,dGf,ddGf,color='k')

def save(filename,temp,X,dX):
    f = open(filename,'wb')
    cPickle.dump(temp,f)
    cPickle.dump(X,f)
    cPickle.dump(dX,f)
    f.close()

def main():
    args = parse_args()
    
    org_file = '%s/dG_raw_noint_2.pkl' % args.direc
    nboots = args.n
    boot_files = ['%s/bootstrap/dG_raw_boot%i.pkl' %(args.direc,x) for x in range(nboots)]

        ### thermodynamics of folding ###
    plt.figure(1)
    a = 1
    b = 0
    temp, dGf_boot, ddGf_boot, dSf_boot, ddSf_boot, dHf_boot, ddHf_boot = plot_avg_bootstrap(boot_files,a,b)
    if args.all:
        plot_all_bootstrap(boot_files,a,b)
    plot_orginal(org_file,a,b)

    save('%s/bootstrap/dGf_boot.pkl' % args.direc, temp, dGf_boot, ddGf_boot)
    save('%s/bootstrap/dSf_boot.pkl' % args.direc, temp[1:-1], dSf_boot, ddSf_boot)
    save('%s/bootstrap/dHf_boot.pkl' % args.direc, temp[1:-1], dHf_boot, ddHf_boot)

    ### thermodynamics of adsorption, folded ###
    plt.figure(2)
    a = 1
    b = 3
    temp, dGadsf_boot, ddGadsf_boot, dSadsf_boot, ddSadsf_boot, dHadsf_boot, ddHadsf_boot = plot_avg_bootstrap(boot_files,a,b)
    if args.all:
        plot_all_bootstrap(boot_files,a,b)
    plot_orginal(org_file,a,b)

    save('%s/bootstrap/dGadsf_boot.pkl' % args.direc, temp, dGadsf_boot, ddGadsf_boot)
    save('%s/bootstrap/dSadsf_boot.pkl' % args.direc, temp[1:-1], dSadsf_boot, ddSadsf_boot)
    save('%s/bootstrap/dHadsf_boot.pkl' % args.direc, temp[1:-1], dHadsf_boot, ddHadsf_boot)

    ### thermodynamics of adsorption, unfolded ###
    plt.figure(3)
    a = 0
    b = 2
    temp, dGadsu_boot, ddGadsu_boot, dSadsu_boot, ddSadsu_boot, dHadsu_boot, ddHadsu_boot = plot_avg_bootstrap(boot_files,a,b)
    if args.all:
        plot_all_bootstrap(boot_files,a,b)
    plot_orginal(org_file,a,b)

    save('%s/bootstrap/dGadsu_boot.pkl' % args.direc, temp, dGadsu_boot, ddGadsu_boot)
    save('%s/bootstrap/dSadsu_boot.pkl' % args.direc, temp[1:-1], dSadsu_boot, ddSadsu_boot)
    save('%s/bootstrap/dHadsu_boot.pkl' % args.direc, temp[1:-1], dHadsu_boot, ddHadsu_boot)

    plt.figure(4)
    temp, dGads_boot, ddGads_boot, dSads_boot, ddSads_boot, dHads_boot, ddHads_boot = rel_ads_bootstrap(boot_files)
    save('%s/bootstrap/ddGads_boot.pkl' % args.direc, temp, dGads_boot, ddGads_boot)
    save('%s/bootstrap/ddSads_boot.pkl' % args.direc, temp[1:-1], dSads_boot, ddSads_boot)
    save('%s/bootstrap/ddHads_boot.pkl' % args.direc, temp[1:-1], dHads_boot, ddHads_boot)
    
    if args.show: 
        plt.show()

if __name__ == '__main__':
    main()
