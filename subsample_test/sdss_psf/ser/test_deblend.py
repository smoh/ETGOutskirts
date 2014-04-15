
from pylab import *
from astropy.table import Table

import mpltools
mpltools.set_pub_single()
rc('savefig', dpi=50)
nsa = Table.read('../../SampleZMprobaEllSub_visual.fits')
traw = Table.read('raw/RAWFIT00000.00048.fits')
tmasked = Table.read('masked/RAWFIT00000.00048.fits')
tdeblended = Table.read('deblended/RAWFIT00000.00048.fits')

praw, pmasked, pdeblended = {}, {}, {}
perr_raw, perr_masked, perr_deblended = {}, {}, {}
for i, name in enumerate(["Ie", "Reff", "N", "Q_BA", "C", "XC", "YC", "PHI"]):
    praw[name] = traw['FIT_SER'][:,i].data
    pmasked[name] = tmasked['FIT_SER'][:,i].data
    pdeblended[name] = tdeblended['FIT_SER'][:,i].data
    perr_raw[name] = traw['PERR_SER'][:,i].data
    perr_masked[name] = tmasked['PERR_SER'][:,i].data
    perr_deblended[name] = tdeblended['PERR_SER'][:,i].data


def effective_radius():
    errorbar(pdeblended['Reff'], praw['Reff'],
                xerr=perr_deblended['Reff'],
                yerr=perr_raw['Reff'],
                fmt='bo', label='raw')
    errorbar(pdeblended['Reff'], pmasked['Reff'],
                xerr=perr_deblended['Reff'],
                yerr=perr_masked['Reff'],
                fmt='ro', label='masked')
    errorbar(pdeblended['Reff'], nsa['SERSIC_TH50']/0.396,
                xerr=perr_deblended['Reff'],
                fmt='ko', label='NSA')
    xscale('log')
    yscale('log')
    xlim(0., 120)
    ylim(0., 200)
    v = linspace(0., 120., 100)
    plot(v, v, c='k', ls='--')
    xlabel("Reff deblended")
    ylabel("Reff other")
    legend(loc='lower right')
    savefig('reff.png')
    close()


def Ie():
    errorbar(log10(pdeblended['Ie']), log10(praw['Ie']),
                    xerr=perr_deblended['Ie']/pdeblended['Ie'],
                    yerr=perr_raw['Ie']/praw['Ie'],
                    fmt='bo', label='raw')
    errorbar(log10(pdeblended['Ie']), log10(pmasked['Ie']),
                    xerr=perr_deblended['Ie']/pdeblended['Ie'],
                    yerr=perr_masked['Ie']/pmasked['Ie'],
                    fmt='ro', label='masked')

    xlim(-2.0, 0.5)
    ylim(-2.5, 0.5)
    v = linspace(-2.5, 0.5, 100)
    plot(v, v, c='k', ls='--')
    xlabel("log(Ie) deblended")
    ylabel("log(Ie) other")
    legend(loc='lower right')
    savefig('Ie.png')
    close()


def center():
    errorbar(praw['XC'] - pdeblended['XC'], praw['YC'] - pdeblended['YC'],
                    xerr=sqrt(perr_raw['XC']**2 + perr_deblended['XC']**2),
                    yerr=sqrt(perr_raw['YC']**2 + perr_deblended['YC']**2),
                    fmt='bo', label='raw')
    errorbar(pmasked['XC'] - pdeblended['XC'], pmasked['YC'] - pdeblended['YC'],
                    xerr=sqrt(perr_masked['XC']**2 + perr_deblended['XC']**2),
                    yerr=sqrt(perr_masked['YC']**2 + perr_deblended['YC']**2),
                    fmt='ro', label='masked')
    axvline(0, c='k', ls='--')
    axhline(0, c='k', ls='--')
    xlabel("XC other - deblended")
    ylabel("YC other - deblended")
    legend(loc='lower right')
    savefig('center.png')
    close()


def sersic_N():
    errorbar(pdeblended['N'], praw['N'],
                xerr=perr_deblended['N'],
                yerr=perr_raw['N'],
                fmt='bo', label='raw')
    errorbar(pdeblended['N'], pmasked['N'],
                xerr=perr_deblended['N'],
                yerr=perr_masked['N'],
                fmt='ro', label='masked')
    errorbar(pdeblended['N'], nsa['SERSIC_N'],
                xerr=perr_deblended['N'],
                fmt='ko', label='NSA')

    xlim(2, 9.5)
    ylim(2, 9.5)
    v = linspace(2, 9.5, 100)
    plot(v, v, c='k', ls='--')
    xlabel("N deblended")
    ylabel("N other")
    legend(loc='lower right')
    savefig('sersicN.png')
    close()


def q_BA():
    errorbar(pdeblended['Q_BA'], praw['Q_BA'],
                xerr=perr_deblended['Q_BA'],
                yerr=perr_raw['Q_BA'],
                fmt='bo', label='raw')
    errorbar(pdeblended['Q_BA'], pmasked['Q_BA'],
                xerr=perr_deblended['Q_BA'],
                yerr=perr_masked['Q_BA'],
                fmt='ro', label='masked')
    errorbar(pdeblended['Q_BA'], nsa['SERSIC_BA'],
                xerr=perr_deblended['Q_BA'],
                fmt='ko', label='NSA')
    xlim(0.65, 1)
    ylim(0.65, 1)
    v = linspace(0.65, 1., 1000)
    plot(v, v, c='k', ls='--')
    xlabel("Q_BA deblended")
    ylabel("Q_BA other")
    legend(loc='lower right')
    savefig('q_BA.png')
    close()


def posang():
    phi_deblended = pdeblended['PHI'] % pi*0.5
    phi_raw = praw['PHI'] % pi*0.5
    phi_masked = pmasked['PHI'] % pi*0.5
    errorbar(phi_deblended, phi_raw,
                xerr=perr_deblended['PHI'],
                yerr=perr_raw['PHI'],
                fmt='bo', label='raw')
    errorbar(phi_deblended, phi_masked,
                xerr=perr_deblended['PHI'],
                yerr=perr_masked['PHI'],
                fmt='ro', label='masked')
    # xlim(0.65, 1)
    # ylim(0.65, 1)
    # v = linspace(0.65, 1., 1000)
    # plot(v, v, c='k', ls='--')
    xlabel("PHI deblended")
    ylabel("PHI other")
    legend(loc='upper left')
    savefig('posang.png')
    close()


def makeImageTable():
    import jinja2
    template = jinja2.Template(
        """
        <!DOCTYPE html>
        <html>
        <body>
        <h2>masked/deblended</h2>
        <table border="1">
            <thead>
                <td>NAME</td>
                <td>Ie</td>
                <td>N</td>
                <td>Reff</td>
                <td>Image</td>
            </thead>
            {% for row in table %}
            <tr height="20%">
                <td>{{  row['NAME'] }}</td>
                {% for p in ['Ie', 'N', 'Reff'] %}
                <td><font color="{{ paint(row[p]) }}">{{ '%4.3f' % (row[p]) }}</font></td>
                {% endfor %}
                <td>
                    <img width="80%" src="{{ '%s' % (plotdir) }}/{{ row['NAME'] }}.png">
                </td>
            </tr>
            {% endfor %}
        </table>
        </body>
        </html>
        """)
    tt = Table()
    tt['NAME'] = tdeblended['NAME'].copy()
    tt['Ie'] = tmasked['FIT_SER'][:,0].data / tdeblended['FIT_SER'][:,0].data
    tt['Reff'] = tmasked['FIT_SER'][:,1].data / tdeblended['FIT_SER'][:,1].data
    tt['N'] = tmasked['FIT_SER'][:,2].data / tdeblended['FIT_SER'][:,2].data

    # paint cells
    def paint(x):
        if x > 2 or x < 0.5:
            return "blue"
        elif x > 3 or x < 0.33:
            return "red"
    templatevars = {'table':tt, "plotdir":"deblended/plots", "paint":paint}
    with open('imagetable.html', 'w') as f:
        f.write(template.render(templatevars))

if __name__ == '__main__':
    makeImageTable()
    # effective_radius()
    # Ie()
    # center()
    # sersic_N()
    # q_BA()
    # posang()