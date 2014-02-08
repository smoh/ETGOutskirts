""" Script to download NSA deblend, ivar, and psf images """

import urllib
import os, errno
from urllib2 import urlopen, URLError, HTTPError


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def downloadfile(url, savename=None, prefix=None):
    """
    Download a file from url

    savename : target filename (default=basename of url)
    """
    if not savename:
        savename = os.path.basename(url)
    if not prefix:
        prefix = './'

    mkdir_p(prefix)
    try:
        f = urlopen(url)
        print "downloading " + url

        # Open our local file for writing
        with open(prefix + savename, "wb") as local_file:
            local_file.write(f.read())
            local_file.close()

    #handle errors
    except HTTPError, e:
        print "HTTP Error:", e.code, url
    except URLError, e:
        print "URL Error:", e.reason, url


def downloadNSA(subdir, iauname, pid, aid, prefix, savename):
    """
    Download NSA child, parent, and psf images
    """
    #TODO: add PSF images

    url = 'http://sdss.physics.nyu.edu/mblanton/v0/detect/v0_1/%s/atlases/%s/' % (subdir, pid)
    child = '%s-%s-atlas-%s.fits.gz' % (iauname, pid, aid)
    parent = '%s-parent-%s.fits.gz' % (iauname, pid)

    for f, sub in zip([child, parent], ['image/', 'ivar/']):
        downloadfile(url+f, savename=savename, prefix=prefix + '/' + sub)



#TODO: write argument parser if needed


if __name__ == '__main__':
    nsaid = '124654.fits'
    subdir = '16h/p48/J161618.57+495542.2'
    iauname = 'J161618.57+495542.2'
    pid = '37'
    aid = '3'
    downloadNSA(subdir, iauname, pid, aid, 'test', nsaid)

