'''
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SALT RSS Aperture Extraction (With Minimum Desired Single to Noise Ratio)

Extract regions (outward from the object's centre) for a FITS file
processed by the SALT RSS pipeline.

The target FITS file is typically a 2D flux spectrum of a galaxy.

The SALT RSS Data Reduction procedure is described in:
http://mips.as.arizona.edu/~khainline/salt_redux.html


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
'''

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Import Libraries
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

import os  # For bash commands
import numpy as np # For array handling
from pyraf import iraf # For IRAF commands
import astropy.io.fits as fits # For FITS file handling
from pathlib import Path # To extract filenames

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Load IRAF Libraries
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

iraf.images()
iraf.images.imutil()

'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Logarithmic Rebinning of Galaxy Spectrum 
(Taken from Cappellari & Emsellem (2004).
Please cite and credit the original paper
from which the log_rebin command comes from.

See: https://iopscience.iop.org/article/10.1086/381875/meta
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''


def log_rebin(lamRange, spec, oversample=1, velscale=None, flux=False):
    """
    Logarithmically rebin a spectrum, while rigorously conserving the flux.
    Basically the photons in the spectrum are simply redistributed according
    to a new grid of pixels, with non-uniform size in the spectral direction.

    When the `flux` keyword is set, this program performs an exact integration
    of the original spectrum, assumed to be a step function within the
    linearly-spaced pixels, onto the new logarithmically-spaced pixels.
    The output was tested to agree with the analytic solution.

    :param lamRange: two elements vector containing the central wavelength
        of the first and last pixels in the spectrum, which is assumed
        to have constant wavelength scale! E.g. from the values in the
        standard FITS keywords: LAMRANGE = CRVAL1 + [0, CDELT1*(NAXIS1 - 1)].
        It must be LAMRANGE[0] < LAMRANGE[1].
    :param spec: input spectrum.
    :param oversample: can be used, not to loose spectral resolution,
        especally for extended wavelength ranges and to avoid aliasing.
        Default: OVERSAMPLE=1 ==> Same number of output pixels as input.
    :param velscale: velocity scale in km/s per pixels. If this variable is
        not defined, then it will contain in output the velocity scale.
        If this variable is defined by the user it will be used
        to set the output number of pixels and wavelength scale.
    :param flux: (boolean) True to preserve total flux. In this case the
        log rebinning changes the pixels flux in proportion to their
        dLam so the following command will show large differences
        beween the spectral shape before and after LOG_REBIN:

           plt.plot(exp(logLam), specNew)  # Plot log-rebinned spectrum
           plt.plot(np.linspace(lamRange[0], lamRange[1], spec.size), spec)

        By defaul, when this is False, the above two lines produce
        two spectra that almost perfectly overlap each other.
    :return: [specNew, logLam, velscale] where logLam is the natural
        logarithm of the wavelength and velscale is in km/s.

    """
    lamRange = np.asarray(lamRange)
    assert len(lamRange) == 2, 'lamRange must contain two elements'
    assert lamRange[0] < lamRange[1], 'It must be lamRange[0] < lamRange[1]'
    assert spec.ndim == 1, 'input spectrum must be a vector'
    n = spec.shape[0]
    m = int(n * oversample)

    dLam = np.diff(lamRange) / (n - 1.)  # Assume constant dLam
    lim = lamRange / dLam + [-0.5, 0.5]  # All in units of dLam
    borders = np.linspace(*lim, num=n + 1)  # Linearly
    logLim = np.log(lim)

    c = 299792.458  # Speed of light in km/s
    if velscale is None:  # Velocity scale is set by user
        velscale = np.diff(logLim) / m * c  # Only for output
    else:
        logScale = velscale / c
        m = int(np.diff(logLim) / logScale)  # Number of output pixels
        logLim[1] = logLim[0] + m * logScale

    newBorders = np.exp(np.linspace(*logLim, num=m + 1))  # Logarithmically
    k = (newBorders - lim[0]).clip(0, n - 1).astype(int)

    specNew = np.add.reduceat(spec, k)[:-1]  # Do analytic integral
    specNew *= np.diff(k) > 0  # fix for design flaw of reduceat()
    specNew += np.diff((newBorders - borders[k]) * spec[k])

    if not flux:
        specNew /= np.diff(newBorders)

    # Output log(wavelength): log of geometric mean
    logLam = np.log(np.sqrt(newBorders[1:] * newBorders[:-1]) * dLam)

    return specNew, logLam, velscale


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Aperture Extraction

The following methods are quite similar,
with minor differences.

extract(): Extracts regions from the centre symmetrically
and then combines them.

extract_1kpc(): Extract 0-1 kiloparsec regions from the centre.

extract_left(): Extract regions left (-) from the centre of the object,
such that these regions have a signal-to-noise ratio (SNR) above the desired
minimum.

extract_right(): Extract regions right (+) from the centre of the object,
such that these regions have a signal-to-noise ratio (SNR) above the desired
minimum.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''


def extract(path, file, desired_SNR, centre, pix_per_kpc, Re_kpc):
    '''

    Extract all regions for a galaxy spectra within the effective radius, such that
    these regions have a signal-to-noise ratio above the desired minimum.

    :param centre [Float]: Centre pixel of the galaxy.
    :param path [String]: File path of the galaxy's FITS file.
    :param desired_SNR [Float]: Desired minimum signal to noise ratio.
    :param pix_per_kpc [Float]: Pixels per kiloparsec ratio for the galaxy observation.
    :param Re_kpc [Float]: Effective radius of the galaxy in kiloparsec.

    :return: None.
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Loop Condition
    flag = False

    # Set size of aperture (in pixels)
    ap_size = 1

    # Set starting pixel
    start_pix_R = centre
    start_pix_L = centre

    # Set run number
    i = 1

    # Get name
    name = Path(file).stem
    print('name = {}'.format(name))

    while flag is False:

        try:

            '''
            Extract Right-Hand Side
            >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            '''
            # Perform block average calculation to get **SUM** total flux
            # Save this as a temporary FITS file
            ymin_R = int(start_pix_R)
            ymax_R = int(start_pix_R + ap_size)

            # Set name of FITS file
            kpc = (((ymax_R + ymin_R) / 2.0) - centre) / pix_per_kpc
            kpc = np.round(kpc, 2)
            if kpc > Re_kpc:
                print('Effective Radius Reached')
                break;

            ap_R = '{}_{}kpc_R.fits'.format(name, kpc)
            ap_L = '{}_{}kpc_L.fits'.format(name, kpc)
            ap = '{}_{}kpc.fits'.format(name, kpc)

            print('ymin_R = {}'.format(ymin_R))
            print('ymax_R = {}'.format(ymax_R))
            print('apR = {}'.format(ap_R))
            iraf.images.imgeom.blkavg(input='{}[*,{}:{}]'.format(file, ymin_R, ymax_R),
                                      output=ap_R, option='sum',
                                      b1=1, b2=ymax_R - ymin_R + 1, b3=1, b4=1, b5=1, b6=1, b7=1, mode='ql')
            '''
            Extract Left-Hand Side
            >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            '''
            ymin_L = int(start_pix_L - ap_size)
            ymax_L = int(start_pix_L)
            print('ymin_L = {}'.format(ymin_L))
            print('ymax_L = {}'.format(ymax_L))
            print('apL = {}'.format(ap_L))
            iraf.images.imgeom.blkavg(input='{}[*,{}:{}]'.format(file, ymin_L, ymax_L),
                                      output=ap_L, option='sum',
                                      b1=1, b2=ymax_L - ymin_L + 1, b3=1, b4=1, b5=1, b6=1, b7=1, mode='ql')
            '''
            Combine Left-Hand Side and Right-Hand Side
            >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            '''
            iraf.imutil.imarith(operand1=ap_L, op='+', operand2=ap_R, result=ap,
                                verbose='no',
                                noact='no')
            iraf.imutil.imdelete(images=ap_L, verify='no', default_action='no', go_ahead='no', mode='ql')
            iraf.imutil.imdelete(images=ap_R, verify='no', default_action='no', go_ahead='no', mode='ql')

        except:
            print('--------------------------')
            print('End of FITS file reached.')
            print('--------------------------')
            # Break the loop when the effective radius is reached and stop extracting regions.
            try:
                os.remove(ap_R)
            except:
                flag = True
                break
            try:
                os.remove(ap_L)
            except:
                flag = True
                break
            try:
                os.remove(ap)
            except:
                flag = True
                break
            # End loop
            flag = True
            break

        '''
        Open Temporary FITS File
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        hdu = fits.open(ap)
        gal_lin = hdu[0].data
        try:
            if len(gal_lin[0] > 1):
                gal_lin = gal_lin[0]
        except:
            pass
        hdr = hdu[0].header
        lamRange1 = hdr['CRVAL1'] + np.array([0., hdr['CDELT1'] * (hdr['NAXIS1'] - 1)])
        galaxy, logLam1, velscale = log_rebin(lamRange1, gal_lin)
        galaxy = galaxy / np.median(galaxy)
        lam = np.exp(logLam1)
        print('lam = {}'.format(lam))
        print('galaxy = {}'.format(galaxy))

        '''
        Calculate Signal-To-Noise Ratio
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        signal = np.median(galaxy)
        n = len(galaxy)
        noise = 0.6052697 * np.median(np.abs(2.0 * galaxy[2:n - 2] - galaxy[0:n - 4] - galaxy[4:n]))
        SNR = signal / noise
        print('Singal-To-Noise Ratio = {}'.format(SNR))

        # Delete aperture if Signal-To-Noise Ratio is too low
        # Then further increase aperture size
        # Else, keep FITS file and adjust starting pixel
        if SNR < desired_SNR:
            iraf.images.imutil.imdelete(images=ap)
            ap_size = ap_size + 1
        else:
            # Add info to header of new FITS file
            iraf.imutil.hedit(images=ap, fields='SNR', value=SNR)
            iraf.imutil.hedit(images=ap, fields='ymin_pixel', value=ymin_R)
            iraf.imutil.hedit(images=ap, fields='ymax_pixel', value=ymax_R)
            # Set new starting pixel
            start_pix_R = start_pix_R + 1
            start_pix_L = start_pix_L - 1
            # Reset aperture size
            ap_size = 1
            # Increment run number
            i = i + 1

    return None


def extract_1kpc(path, file, centre, pix_per_kpc):
    '''
    Extract the 0kpc - 1kpc regions for a galaxy spectra within the effective radius.

    :param centre [Float]: Centre pixel of the galaxy.
    :param path [String]: File path of the galaxy's FITS file.
    :param pix_per_kpc [Float]: Pixels per kiloparsec ratio for the galaxy observation.
    :param Re_kpc [Float]: Effective radius of the galaxy in kiloparsec.

    :return: None.
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Get name
    name = Path(file).stem
    print('name = {}'.format(name))

    # Set starting pixel
    start_pix_R = centre
    start_pix_L = centre

    ap_size = pix_per_kpc

    '''
    Extract Right-Hand Side
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    # Perform block average calculation to get **SUM** total flux
    # Save this as a temporary FITS file
    ymin_R = int(start_pix_R)
    ymax_R = int(start_pix_R + ap_size)

    # Set name of FITS file
    kpc = 0.5
    ap_R = '{}_{}kpc_R.fits'.format(name, kpc)
    ap_L = '{}_{}kpc_L.fits'.format(name, kpc)
    ap = '{}_{}kpc.fits'.format(name, kpc)

    print('ymin_R = {}'.format(ymin_R))
    print('ymax_R = {}'.format(ymax_R))
    print('apR = {}'.format(ap_R))
    iraf.images.imgeom.blkavg(input='{}[*,{}:{}]'.format(file, ymin_R, ymax_R),
                              output=ap_R, option='sum',
                              b1=1, b2=ymax_R - ymin_R + 1, b3=1, b4=1, b5=1, b6=1, b7=1, mode='ql')
    '''
    Extract Left-Hand Side
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    ymin_L = int(start_pix_L - ap_size)
    ymax_L = int(start_pix_L)
    print('ymin_L = {}'.format(ymin_L))
    print('ymax_L = {}'.format(ymax_L))
    print('apL = {}'.format(ap_L))
    iraf.images.imgeom.blkavg(input='{}[*,{}:{}]'.format(file, ymin_L, ymax_L),
                              output=ap_L, option='sum',
                              b1=1, b2=ymax_L - ymin_L + 1, b3=1, b4=1, b5=1, b6=1, b7=1, mode='ql')
    '''
    Combine Left-Hand Side and Right-Hand Side
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    iraf.imutil.imarith(operand1=ap_L, op='+', operand2=ap_R, result=ap,
                        verbose='no',
                        noact='no')
    iraf.imutil.imdelete(images=ap_L, verify='no', default_action='no', go_ahead='no', mode='ql')
    iraf.imutil.imdelete(images=ap_R, verify='no', default_action='no', go_ahead='no', mode='ql')

    return None


def extract_left(path, file, desired_SNR, centre, pix_per_kpc, Re_kpc):
    '''

    Extract all regions for a galaxy spectra **on its left-hand-side** within the effective radius, such that
    these regions have a signal-to-noise ratio above the desired minimum.

    :param centre [Float]: Centre pixel of the galaxy.
    :param path [String]: File path of the galaxy's FITS file.
    :param desired_SNR [Float]: Desired minimum signal to noise ratio.
    :param pix_per_kpc [Float]: Pixels per kiloparsec ratio for the galaxy observation.
    :param Re_kpc [Float]: Effective radius of the galaxy in kiloparsec.

    :return: None.
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Loop Condition
    flag = False

    # Set size of aperture (in pixels)
    ap_size = 1

    # Set starting pixel
    start_pix_L = centre

    # Set run number
    i = 1

    # Get name
    name = Path(file).stem
    print('name = {}'.format(name))

    while flag is False:

        try:

            '''
            Extract Left-Hand Side (ONLY)
            >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            '''
            # Perform block average calculation to get **SUM** total flux
            # Save this as a temporary FITS file
            ymin_L = int(start_pix_L - ap_size)
            ymax_L = int(start_pix_L)

            # Set name of FITS file
            kpc = (((ymax_L + ymin_L) / 2.0) - centre) / pix_per_kpc
            kpc = np.round(kpc, 2)
            if np.abs(kpc) > Re_kpc:
                print('Effective Radius Reached')
                break;

            ap_L = '{}_{}kpc_L.fits'.format(name, kpc)

            print('ymin_L = {}'.format(ymin_L))
            print('ymax_L = {}'.format(ymax_L))
            print('apL = {}'.format(ap_L))
            iraf.images.imgeom.blkavg(input='{}[*,{}:{}]'.format(file, ymin_L, ymax_L),
                                      output=ap_L, option='sum',
                                      b1=1, b2=ymax_L - ymin_L + 1, b3=1, b4=1, b5=1, b6=1, b7=1, mode='ql')

        except:
            print('--------------------------')
            print('End of FITS file reached.')
            print('--------------------------')
            # Delete last FITS from failed blkavg run
            try:
                os.remove(ap_L)
            except:
                flag = True
                break
            flag = True
            break

        '''
        Open Temporary FITS File
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        hdu = fits.open(ap_L)
        gal_lin = hdu[0].data
        try:
            if len(gal_lin[0] > 1):
                gal_lin = gal_lin[0]
        except:
            pass
        hdr = hdu[0].header
        lamRange1 = hdr['CRVAL1'] + np.array([0., hdr['CDELT1'] * (hdr['NAXIS1'] - 1)])
        galaxy, logLam1, velscale = log_rebin(lamRange1, gal_lin)
        galaxy = galaxy / np.median(galaxy)
        lam = np.exp(logLam1)
        print('lam = {}'.format(lam))
        print('galaxy = {}'.format(galaxy))

        '''
        Calculate Signal-To-Noise Ratio
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        signal = np.median(galaxy)
        n = len(galaxy)
        noise = 0.6052697 * np.median(np.abs(2.0 * galaxy[2:n - 2] - galaxy[0:n - 4] - galaxy[4:n]))
        SNR = signal / noise
        print('Singal-To-Noise Ratio = {}'.format(SNR))

        # Delete aperture if Signal-To-Noise Ratio is too low
        # Then further increase aperture size
        # Else, keep FITS file and adjust starting pixel
        if SNR < desired_SNR:
            iraf.images.imutil.imdelete(images=ap_L)
            ap_size = ap_size + 1
        else:
            # Add info to header of new FITS file
            iraf.imutil.hedit(images=ap_L, fields='SNR', value=SNR)
            iraf.imutil.hedit(images=ap_L, fields='ymin_pixel', value=ymin_L)
            iraf.imutil.hedit(images=ap_L, fields='ymax_pixel', value=ymax_L)
            # Set new starting pixel
            start_pix_L = start_pix_L - 1
            # Reset aperture size
            ap_size = 1
            # Increment run number
            i = i + 1

    return None


def extract_right(path, file, desired_SNR, centre, pix_per_kpc, Re_kpc):
    '''

    Extract all regions for a galaxy spectra **on its right-hand-side** within the effective radius, such that
    these regions have a signal-to-noise ratio above the desired minimum.

    :param centre [Float]: Centre pixel of the galaxy.
    :param path [String]: File path of the galaxy's FITS file.
    :param desired_SNR [Float]: Desired minimum signal to noise ratio.
    :param pix_per_kpc [Float]: Pixels per kiloparsec ratio for the galaxy observation.
    :param Re_kpc [Float]: Effective radius of the galaxy in kiloparsec.

    :return: None.
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Loop Condition
    flag = False

    # Set size of aperture (in pixels)
    ap_size = 1

    # Set starting pixel
    start_pix_R = centre

    # Set run number
    i = 1

    # Get name
    name = Path(file).stem
    print('name = {}'.format(name))

    while flag is False:

        try:

            '''
            Extract Right-Hand Side (ONLY)
            >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            '''
            # Perform block average calculation to get **SUM** total flux
            # Save this as a temporary FITS file
            ymin_R = int(start_pix_R)
            ymax_R = int(start_pix_R + ap_size)

            # Set name of FITS file
            kpc = (((ymax_R + ymin_R) / 2.0) - centre) / pix_per_kpc
            kpc = np.round(kpc, 2)
            if kpc > Re_kpc:
                print('Effective Radius Reached')
                break;

            ap_R = '{}_{}kpc_R.fits'.format(name, kpc)

            print('ymin_R = {}'.format(ymin_R))
            print('ymax_R = {}'.format(ymax_R))
            print('apR = {}'.format(ap_R))
            iraf.images.imgeom.blkavg(input='{}[*,{}:{}]'.format(file, ymin_R, ymax_R),
                                      output=ap_R, option='sum',
                                      b1=1, b2=ymax_R - ymin_R + 1, b3=1, b4=1, b5=1, b6=1, b7=1, mode='ql')

        except:
            print('--------------------------')
            print('End of FITS file reached.')
            print('--------------------------')
            # Break the loop when the effective radius is reached and stop extracting regions.
            try:
                os.remove(ap_R)
            except:
                flag = True
                break
            flag = True
            break

        '''
        Open Temporary FITS File
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        hdu = fits.open(ap_R)
        gal_lin = hdu[0].data
        try:
            if len(gal_lin[0] > 1):
                gal_lin = gal_lin[0]
        except:
            pass
        hdr = hdu[0].header
        lamRange1 = hdr['CRVAL1'] + np.array([0., hdr['CDELT1'] * (hdr['NAXIS1'] - 1)])
        galaxy, logLam1, velscale = log_rebin(lamRange1, gal_lin)
        galaxy = galaxy / np.median(galaxy)
        lam = np.exp(logLam1)
        print('lam = {}'.format(lam))
        print('galaxy = {}'.format(galaxy))

        '''
        Calculate Signal-To-Noise Ratio
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        signal = np.median(galaxy)
        n = len(galaxy)
        noise = 0.6052697 * np.median(np.abs(2.0 * galaxy[2:n - 2] - galaxy[0:n - 4] - galaxy[4:n]))
        SNR = signal / noise
        print('Singal-To-Noise Ratio = {}'.format(SNR))

        # Delete aperture if Signal-To-Noise Ratio is too low
        # Then further increase aperture size
        # Else, keep FITS file and adjust starting pixel
        if SNR < desired_SNR:
            iraf.images.imutil.imdelete(images=ap_R)
            ap_size = ap_size + 1
        else:
            # Add info to header of new FITS file
            iraf.imutil.hedit(images=ap_R, fields='SNR', value=SNR)
            iraf.imutil.hedit(images=ap_R, fields='ymin_pixel', value=ymin_R)
            iraf.imutil.hedit(images=ap_R, fields='ymax_pixel', value=ymax_R)
            # Set new starting pixel
            start_pix_R = start_pix_R + 1
            # Reset aperture size
            ap_size = 1
            # Increment run number
            i = i + 1

    return None


'''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Run Program
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

# Set extraction mode
mode = str(raw_input("Please enter extraction mode. Enter 'A' for All, 'L' for Left and 'R' for Right:  "))
# Set filepath of data
path = raw_input('Please enter file path of the galaxy\'s FITS file:  ')
# Set filename
file = raw_input('Please enter file name of the galaxy\'s FITS file (with complete extension):  ')
# Set desired signal-to-noise ratio
desired_SNR = float(input('Please enter the desired Signal-To-Noise Ratio:  '))
# Set centre pixel number
centre = float(input('Please enter the center pixel of the galaxy:  '))
# Set pixel per kiloparsec ratio
pix_per_kpc = float(input('Please enter the pixel per kpc value for this galaxy:  '))
# Set effective radius
Re_kpc = float(input('Please enter the effective radius of this galaxy (in kpc):  '))

# Extract regions
if mode == 'A':
    extract(path, file, desired_SNR, centre, pix_per_kpc, Re_kpc)
    extract_1kpc(path, file, centre, pix_per_kpc)
elif mode == 'L':
    extract_left(path, file, desired_SNR, centre, pix_per_kpc, Re_kpc)
elif mode == 'R':
    extract_right(path, file, desired_SNR, centre, pix_per_kpc, Re_kpc)
