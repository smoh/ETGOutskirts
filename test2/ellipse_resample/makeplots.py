from glob import glob
from pylab import *
from astropy.table import Table, Column

def read_ellipse_output(filename):
    return Table.read(
        filename, data_start=0,
        format='ascii.commented_header', header_start=1, fill_values=['INDEF', '-99'])

def build_data_array(dirname):
    colnames = ['PA', 'ELLIP', 'STOP']
    buffsize = 100
    d = Table()
    d.add_column(Column(name='SMA', data=arange(0, buffsize*2, 2)))
    for i, fn in enumerate(glob(dirname+'/*.dat')):
        t = read_ellipse_output(fn)
        if i is 0:
            d.add_columns([Column(name=name, data=zeros([buffsize, 100])-99) for name in colnames])
        assert array_equiv(d['SMA'][:len(t)].data, t['SMA'].data), "SMA array do not match"
        for col in colnames:
            d[col][:len(t),i] = t[col]
    return d


rc('lines', linewidth=1)
# d = build_data_array('J1003')

def plotAllData(dirname):
    fig, axs = subplots(2, 1, figsize=(6, 6), sharex=True)
    fig.suptitle(dirname, fontsize=15)
    fig.subplots_adjust(hspace=0.05, bottom=0.15, top=0.9, left=0.15, right=0.95)
    axs[0].set_ylabel('PA')
    axs[1].set_ylabel('ELLIP')
    axs[1].set_xlabel('SMA pixels')
    for fn in glob(dirname+'/*.dat'):
        t = read_ellipse_output(fn)
        bad = t['STOP'] != 0
        axs[0].plot(t['SMA'], t['PA'], 'k-')
        axs[0].plot(t['SMA'][bad], t['PA'][bad], 'ro', mfc='None', mec='r')
        axs[1].plot(t['SMA'], t['ELLIP'], 'k-')
        axs[1].plot(t['SMA'][bad], t['ELLIP'][bad], 'ro', mfc='None', mec='r')
    fig.savefig('%s_raw.png' % (dirname), dpi=80)
    close()

def plotPercentiles(dirname):
    fig, axs = subplots(3, 1, figsize=(6, 10), sharex=True)
    fig.suptitle(dirname, fontsize=15)
    fig.subplots_adjust(hspace=0.08, bottom=0.15, top=0.9, left=0.15, right=0.95)
    axs[0].set_ylabel('PA')
    axs[1].set_ylabel('ELLIP')
    axs[2].set_ylabel('N records')
    axs[2].set_xlabel('SMA pixels')
    d = build_data_array(dirname)

    for ax, col in zip(axs, ['PA', 'ELLIP']):
        sma = []
        prc = []
        nrecords = []
        for smatmp, drow, stoprow in zip(d['SMA'], d[col], d['STOP']):
            filtered = drow[(drow!=-99)&(stoprow==0)]
            if count_nonzero(filtered) is 0: 
                continue
            nrecords.append(len(filtered))
            prc.append(percentile(filtered, [0, 30, 50, 70, 100]))
            sma.append(smatmp)
        sma = atleast_1d(sma)
        nrecords = atleast_1d(nrecords)
        prc = atleast_2d(prc)
        ax.plot(sma, prc, 'k-')
        ax.fill_between(sma, prc[:,1], prc[:,3], alpha=0.5)
        ax.set_xlim(sma.min(), sma.max())
        axs[2].plot(sma, nrecords, label=col)
    fig.savefig('%s_prc.png' % (dirname), dpi=80)
    close()


if __name__ == '__main__':
    for dirname in glob('J[0-9]*/')[:]:
        print dirname
        try:
            plotAllData(dirname.strip('/'))
            plotPercentiles(dirname.strip('/'))
        except:
            print dirname, "failed"
