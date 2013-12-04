import numpy
import matplotlib
import matplotlib.pyplot as plt

def main():
    plt.rc('text',usetex=True)
    matplotlib.rc('font', family = 'serif')
    font = {'family' : 'serif',
            'size'   : 'larger'}

    adsorbed = numpy.load('adsorbed_rate.out.npy')
    desorbed = numpy.load('desorbed_rate.out.npy')
    N = len(adsorbed)
    width = 0.35
    ind = numpy.arange(N)
    fig,ax = plt.subplots()
    rects1 = ax.bar(ind,adsorbed[:,0],width,yerr=adsorbed[:,1]/N**.5)
    rects2 = ax.bar(ind+width,desorbed[:,0],width,color='g',yerr=desorbed[:,1]/N**.5)
    
    lam = numpy.arange(0.1,0.65,.05)
    ax.set_ylabel('Number of folding and unfolding transitions')
    ax.set_xlabel(r'$\lambda$')
    ax.set_xticks(ind+width)
    ax.set_xticklabels([str(x) for x in lam])
    
    ax.legend((rects1[0],rects2[0]),('adsorbed states','desorbed states'))
    plt.xlim((-.3,N))
    
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Tab2.pdf')
    plt.savefig('/home/edz3fz/proteinmontecarlo/manuscripts/figures/Tab2.eps')
    plt.show()

if __name__ == '__main__':
    main()
