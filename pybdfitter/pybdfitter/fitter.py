import subprocess
import os, sys
import errno

__all__ = ['fit_sample']


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def build_command_str(input, profile, start, end, outdir, datadir='data',
                      imgdir=None, filter='r', residual=True, debug=True):
    """
    Build idl command string for FIT_SAMPLE
    """
    cmd = "fit_sample, '{:s}', {:d}, {:d}, '{:s}', '{:s}'".format(
        input, start, end + 1, outdir + '/', datadir + '/')
    cmd = cmd + ", profiles={" + profile + "}"
    if filter:
        cmd += ', filter=%s' % (filter)
    if residual:
        cmd += ', /residuals'
    if debug:
        cmd += ', /debug'
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
                      filter='r', residual=True, debug=True, 
                      stdout=None, stderr=None, quiet=False):
    """
    Run bdfitter

    input : str
        input file name
    profile : str
        'profile_name:number_of_params' e.g., 'DVC:8'
    start : int
        start row
    end : int
        end row (inclusive)
    outdir : str
        output directory
    datadir : str, default='data'
        data root directory
    imgdir : str, optional
        set in case images are not in [datadir]/images
    output : str, optional
        output table name. If not given, it will be name RAWFITXXX..
    residual : bool, default=True
        save model images?
    """
    if check_IDL():
        print "NO IDL found. Nothing else done"
        return 0
    # prepare directories
    mkdir_p(outdir)
    if residual:
        mkdir_p(outdir + '/models')

    cmd = build_command_str(
        input, profile, start, end, outdir,
        datadir=datadir, imgdir=imgdir, filter=filter, residual=residual,
        debug=debug)
    if quiet:
        proc = subprocess.Popen(["idl", "-e", cmd, "-quiet"],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        proc = subprocess.Popen(["idl", "-e", cmd],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in iter(proc.stdout.readline, b''):
            # print ' STDOUT '.center(80, '=')
            print line,
            sys.stdout.flush()
            sys.stderr.flush()
    out = proc.communicate()  # close p.stdout, wait for the subprocess to exit
    # if proc.returncode:  # if fails
    print proc.stderr

    out_expected = 'RAWFIT%05i.%05i.fits' % (start, end)
    if output:
        os.rename(out_expected, output)
    outTable = outdir + '/' + [out_expected if not output else output][0]
    return outTable


def test_fit_sample_stream():
    from pybdfitter.make_input import single_DVC
    input = 'test/sample_test.fits'
    profile = 'DVC:8'
    start = 0
    end = 1
    outdir = 'test/out'
    datadir = '../../subsampletest/data'
    imgdir = '../../subsampletest/data/deblended'

    single_DVC(input, 'test/sample_test_DVC.fits')
    fit_sample_stream('test/sample_test_DVC.fits', profile, start, end, outdir,
                      datadir=datadir, imgdir=imgdir)


if __name__ == '__main__':
    test_fit_sample_stream()
