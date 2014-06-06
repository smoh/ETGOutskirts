import tempfile
import subprocess
import os
from astropy.table import Table


def kinemetry_phot(imagefile, x0, y0, initial_pa, initial_q,
                   verbose=False, ntrm=10):
    """
    Do ellipse fitting to an image using IDL KINEMETRY_PHOP

    imagefile : file name of the image
    x0, y0 : position of the center (fixed for all ellipses)
    initial_pa : initial value of PA (-90 to 90 degrees, 0 along the x-axis)
    initial_q ; initial value of flattening

    Returns Table instance containing these columns
        RAD : 1-d array of radius
        PA : 1-d array of PA at each radius
        Q : 1-d array of Q at each radius
        CF : 2-d array containing coefficients of the Fourier expansion at
            each radius. (NRAD, NCOEFF)
        ER_PA : 1-sigma error of PA
        ER_Q : 1-sigma error of Q
        ER_CF : 1-sigma error of CF
    """
    # make temporary output file
    dum, tempname = tempfile.mkstemp(suffix='.fits')
    outname = tempname

    cmd = "KINEMETRY_PHOT, '{:s}', {:f}, {:f}, {:f}, {:f}, '{:s}' ".format(
        imagefile, x0, y0, initial_pa, initial_q, outname)
    cmd += ', NTRM=%i' % (ntrm)
    if verbose:
        cmd += ', /VERBOSE'

    proc = subprocess.Popen(["idl", "-e", cmd],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in iter(proc.stdout.readline, b''):
        print line
    out = proc.communicate()
    if proc.returncode:
        print proc.stderr

    out = Table.read(outname, format='fits')
    os.remove(outname)
    return out
