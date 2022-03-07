# SALT_RSS_Extract_SNR
Extract regions (outward from the object's centre) for a FITS file processed by the SALT RSS (Southern African Large Telescope Robert Stobie Spectrograph) pipeline.

What this program does:
========================================

The SALT RSS data reduction procedure is well documented.

See the data reduction procedure here:
http://mips.as.arizona.edu/~khainline/salt_redux.html

Also see the data reduction FAQ:
https://astronomers.salt.ac.za/data/data-reduction-faq/

This program takes a FITS image processed by the SALT-RSS pipeline and extracts regions from the image above a desired
minimum signal-to-noise ratio.

Prerequisites for using this program:
========================================

This program requires a programming environment with Python 3.6 installed.
With the following external libraries:
- PyRAF, 
- Astropy
- pathlib
- see: https://faculty1.coloradocollege.edu/~sburns/courses/18-19/pc362/Anaconda_IRAF_install.html for installation instructions.

