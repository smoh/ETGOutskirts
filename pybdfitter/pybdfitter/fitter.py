import subprocess
import datetime
import os

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

def check_IDL():
    """ return 1 if no IDL command found """
    proc = subprocess.Popen(
        ["which", "idl"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.communicate()
    return proc.returncode


def fit_sample(input, profile, start, end, outdir,
        datadir='data', imgdir=None, output=None,
        filter='r', residual=True, debug=True, savepsf=False,
        stdout=None, stderr=None):
    
    if check_IDL():
        print "NO IDL found. Nothing else done"
        return 0
    
    # prepare directories
    mkdir_p(outdir)
    if residual:
        mkdir_p(outdir + '/models')
    if savepsf:
        mkdir_p(datadir + '/sdss_psf')

    cmd = build_command_str(input, profile, start, end, outdir,
            datadir=datadir, imgdir=imgdir, filter=filter, residual=residual,
            debug=debug, savepsf=savepsf)
    proc = subprocess.Popen(["idl", "-e", cmd],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = proc.communicate()
    if out[0]:
        print out[0]
    if proc.returncode:
        print out[1]
        return 0

    out_expected = 'RAWFIT%05i.%05i.fits' % (start, end)
    if output:
        os.rename(out_expected, output)
    outTable = outdir +'/'+ [out_expected if not output else output][0]
    return outTable

def fit_sample_stream(input, profile, start, end, outdir,
        datadir='data', imgdir=None, output=None,
        filter='r', residual=True, debug=True, savepsf=False,
        stdout=None, stderr=None):
    
    if check_IDL():
        print "NO IDL found. Nothing else done"
        return 0
    
    # prepare directories
    mkdir_p(outdir)
    if residual:
        mkdir_p(outdir + '/models')
    if savepsf:
        mkdir_p(datadir + '/sdss_psf')

    cmd = build_command_str(input, profile, start, end, outdir,
            datadir=datadir, imgdir=imgdir, filter=filter, residual=residual,
            debug=debug, savepsf=savepsf)
    proc = subprocess.Popen(["idl", "-e", cmd],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in iter(proc.stdout.readline, b''):
        print line,
    out = proc.communicate() # close p.stdout, wait for the subprocess to exit
    print out[1]

    out_expected = 'RAWFIT%05i.%05i.fits' % (start, end)
    if output:
        os.rename(out_expected, output)
    outTable = outdir +'/'+ [out_expected if not output else output][0]
    return outTable
