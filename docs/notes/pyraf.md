
Mon Oct 27 21:46:49 EDT 2014

- pyraf programmer's guide](http://stsdas.stsci.edu/stsci_python_sphinxdocs_2.13/docs/pyraf_guide.pdf)
- From pyraf tutorial, you can start pyraf from anywhere unlike iraf. It will
  first look for login.cl in the current directory, and then in ~/iraf.
- pyraf 2.1.7 does not set parameters correctly for some tasks including
  ellipse. I got a response to look at this [ticket](https://aeon.stsci.edu/ssb/trac/pyraf/ticket/217) that says, for now, set parameters directly from the main task as in `iraf.ellipse.x0 = 250` instead of `iraf.geompar.x0 = 250`.

