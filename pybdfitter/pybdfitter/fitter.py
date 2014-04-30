import subprocess
import datetime

def mkdir_p(path):
    import os, errno
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def build_command_str(input, profile, start, end, outdir, datadir='data', imgdir=None,
        filter='r', residual=True, debug=True, savepsf=False):
    """
    Build idl command string
    """
    cmd = "fit_sample, '{:s}', {:d}, {:d}, '{:s}', '{:s}'".format(
            input, start, end+1, outdir+'/', datadir+'/')
    cmd = cmd + ", profiles={" + profile + "}"
    if filter:
        cmd += ', filter=%s' % (filter)
    if residual:
        cmd += ', /residuals'
    if debug:
        cmd += ', /debug'
    if savepsf:
        cmd += ', /savepsf'
    if imgdir:
        cmd += ", imgdir='{:s}'".format(imgdir + '/')
    return cmd

def fit_sample(input, profile, start, end, outdir, datadir='data', imgdir=None,
        filter='r', residual=True, debug=True, savepsf=False,
        stdout=None, stderr=None):
    
    # prepare directories
    mkdir_p(outdir)
    if residual:
        mkdir_p(outdir + '/models')
    if savepsf:
        mkdir_p(datadir + '/sdss_psf')

    cmd = build_command_str(input, profile, start, end, outdir,
            datadir=datadir, imgdir=imgdir, filter=filter, residual=residual,
            debug=debug, savepsf=savepsf)
    with open(stdout, 'w') as stdout:
        subprocess.call(["idl", "-e", cmd], stdout=stdout)

    outTable = outdir + '/RAWFIT%05i.%05i.fits' % (start, end)
    return outTable

