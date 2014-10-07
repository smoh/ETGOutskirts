import numpy as np
from scipy.stats import binned_statistic

__all__ = ["nanomaggie2mag",
           "mag2nanomaggie",
           "radialprofile"]


def nanomaggie2mag(nnmg):
    """ nnmg : flux in nanomaggie """
    return 22.5 - 2.5*np.log10(nnmg)


def mag2nanomaggie(mag):
    return 10**((22.5 - mag)/2.5)


def radialprofile(image, center=None, binsize=2):
    """Calculate radial profile from an image

    image : 2d array
        input image array
    center : (x, y)
        center coordinate in pixels
        if None, use the center of the image
    binsize : float
        size of radius bin

    Returns
    bin_centers, values

    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    bins = np.arange(0, r.max(), binsize)
    bin_centers = (bins[1:] + bins[:-1]) * 0.5

    p = binned_statistic(r.ravel(), image.ravel(), bins=bins)[0]

    return bin_centers, p


