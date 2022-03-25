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
Aperture Extraction

The following methods are quite similar,
with minor differences.

extract_1kpc(): Extract 0-1 kiloparsec regions from the centre.

extract_left(): Extract regions left (-) from the centre of the object,
such that these regions have a signal-to-noise ratio (SNR) above the desired
minimum.

extract_right(): Extract regions right (+) from the centre of the object,
such that these regions have a signal-to-noise ratio (SNR) above the desired
minimum.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''

def extract_1kpc_left(path, file, centre, pix_per_kpc):

    '''
    :param path [String]: File path containing the galaxy's FITS file.
    :param file [String]: The galaxy's FITS file.
    :param centre [float]: Central pixel of the galaxy's FITS file.
    :param pix_per_kpc [float]: Pixel per kiloparsec ratio for the galaxy.
    :return:
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Set pixel ranges
    pix_end = int(centre)
    pix_start = pix_end - 1

    # Set operands
    op1 = '{}[*,{}]'.format(file, pix_start)
    op2 = '{}[*,{}]'.format(file, pix_end)

    # Set file names
    res = 'result.fits'
    temp = 'temp.fits'

    '''
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Extract Apertures
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    flag = False;
    while flag is False:

        # Set kpc
        kpc = (((pix_start + pix_end) / 2.0) - centre) / pix_per_kpc

        print('pix_start = {}'.format(pix_start))
        print('pix_end = {}'.format(pix_end))
        iraf.images.imutil.imarith(operand1=op1, op='+', operand2=op2, result=res, verbose='No', mode='ql')

        if kpc > -0.5:
            iraf.images.imutil.imcopy(input=res, output=temp)
            iraf.images.imutil.imdelete(images=res)
            pix_start = pix_start - 1
            op1 = '{}[*,{}]'.format(file, pix_start)
            op2 = '{}[0]'.format(temp)
        else:
            out_file = '{}_-0.5kpcL.fits'.format(file)
            iraf.images.imutil.imcopy(input=res, output=out_file)
            iraf.images.imutil.imdelete(images=res)
            # Update Header Information
            iraf.imutil.hedit(images=out_file, fields='ymin_pixel', value=pix_start)
            iraf.imutil.hedit(images=out_file, fields='ymax_pixel', value=pix_end)
            pix_end = pix_start - 1
            pix_start = pix_end - 1
            flag = True
            iraf.images.imutil.imdelete(images=temp)
            break

    return None

def extract_1kpc_right(path, file, centre, pix_per_kpc):

    '''
    :param path [String]: File path containing the galaxy's FITS file.
    :param file [String]: The galaxy's FITS file.
    :param centre [float]: Central pixel of the galaxy's FITS file.
    :param pix_per_kpc [float]: Pixel per kiloparsec ratio for the galaxy.
    :return:
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Set pixel ranges
    pix_start = int(centre)
    pix_end = int(centre + 1)

    # Set operands
    op1 = '{}[*,{}]'.format(file, pix_start)
    op2 = '{}[*,{}]'.format(file, pix_end)

    # Set file names
    res = 'result.fits'
    temp = 'temp.fits'

    '''
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Extract Apertures
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    flag = False;
    while flag is False:

        # Set kpc
        kpc = (((pix_start + pix_end) / 2.0) - centre) / pix_per_kpc

        print('pix_start = {}'.format(pix_start))
        print('pix_end = {}'.format(pix_end))
        iraf.images.imutil.imarith(operand1=op1, op='+', operand2=op2, result=res, verbose='No', mode='ql')

        if kpc < 0.5:
            iraf.images.imutil.imcopy(input = res, output = temp)
            iraf.images.imutil.imdelete(images = res)
            pix_end = pix_end + 1
            op1 = '{}[0]'.format(temp)
            op2 = '{}[*,{}]'.format(file, pix_end)
        else:
            out_file = '{}_0.5kpcR.fits'.format(file)
            iraf.images.imutil.imcopy(input=res, output=out_file)
            iraf.images.imutil.imdelete(images=res)
            # Update Header Information
            iraf.imutil.hedit(images=out_file, fields='ymin_pixel', value=pix_start)
            iraf.imutil.hedit(images=out_file, fields='ymax_pixel', value=pix_end)
            pix_start = pix_end + 1
            pix_end = pix_start + 1
            flag = True
            iraf.images.imutil.imdelete(images=temp)
            break

    return None

def extract_left(path, file, desired_SNR, centre, pix_per_kpc):

    '''
    :param path [String]: File path containing the galaxy's FITS file.
    :param file [String]: The galaxy's FITS file.
    :param desired_SNR [float]: Desired *minimum* signal-to-noise ratio for aperture extraction.
    :param centre [float]: Central pixel of the galaxy's FITS file.
    :param pix_per_kpc [float]: Pixel per kiloparsec ratio for the galaxy.
    :return:
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Set pixel ranges
    pix_end = int(centre)
    pix_start = pix_end - 1

    # Set operands
    op1 = '{}[*,{}]'.format(file, pix_start)
    op2 = '{}[*,{}]'.format(file, pix_end)

    # Set file names
    res = 'result.fits'
    temp = 'temp.fits'

    '''
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Extract Apertures
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    flag = False;
    while flag is False:

        # Set kpc
        kpc = (((pix_start + pix_end) / 2.0) - centre) / pix_per_kpc
        kpc = np.round(kpc, 2)

        if pix_end - pix_start > 50:
            print('Maximum aperture size of 50 exceeded.')
            try:
                iraf.images.imutil.imdelete(images=temp)
            except:
                pass
            try:
                iraf.images.imutil.imdelete(images=res)
            except:
                pass
            flag = True;
            break;

        print('pix_start = {}'.format(pix_start))
        print('pix_end = {}'.format(pix_end))
        iraf.images.imutil.imarith(operand1=op1, op='+', operand2=op2, result=res, verbose='No', mode='ql')

        '''
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        Open FITS File
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        hdu = fits.open(res)
        galaxy = hdu[0].data
        try:
            if len(galaxy[0] > 1):
                galaxy = galaxy[0]
        except:
            pass
        print('galaxy = {}'.format(galaxy))
        '''
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
            iraf.images.imutil.imcopy(input=res, output=temp)
            iraf.images.imutil.imdelete(images=res)
            pix_start = pix_start - 1
            op1 = '{}[*,{}]'.format(file, pix_start)
            op2 = '{}[0]'.format(temp)
        else:
            out_file = '{}_{}kpcL.fits'.format(file, kpc)
            iraf.images.imutil.imcopy(input=res, output=out_file)
            iraf.images.imutil.imdelete(images=res)
            # Update Header Information
            iraf.imutil.hedit(images=out_file, fields='SNR', value=SNR)
            iraf.imutil.hedit(images=out_file, fields='ymin_pixel', value=pix_start)
            iraf.imutil.hedit(images=out_file, fields='ymax_pixel', value=pix_end)
            pix_end = pix_start - 1
            pix_start = pix_end - 1
            op1 = '{}[*,{}]'.format(file, pix_start)
            op2 = '{}[*,{}]'.format(file, pix_end)

    return None

def extract_right(path, file, desired_SNR, centre, pix_per_kpc):

    '''
    :param path [String]: File path containing the galaxy's FITS file.
    :param file [String]: The galaxy's FITS file.
    :param desired_SNR [float]: Desired *minimum* signal-to-noise ratio for aperture extraction.
    :param centre [float]: Central pixel of the galaxy's FITS file.
    :param pix_per_kpc [float]: Pixel per kiloparsec ratio for the galaxy.
    :return:
    '''

    # Change to file path
    os.chdir(path)
    print('Changed to path: {}'.format(path))

    # Set pixel ranges
    pix_start = int(centre)
    pix_end = int(centre + 1)

    # Set operands
    op1 = '{}[*,{}]'.format(file, pix_start)
    op2 = '{}[*,{}]'.format(file, pix_end)

    # Set file names
    res = 'result.fits'
    temp = 'temp.fits'

    '''
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Extract Apertures
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    flag = False;
    while flag is False:

        # Set kpc
        kpc = (((pix_start + pix_end) / 2.0) - centre) / pix_per_kpc
        kpc = np.round(kpc, 2)

        if pix_end - pix_start > 50:
            print('Maximum aperture size of 50 exceeded.')
            try:
                iraf.images.imutil.imdelete(images = temp)
            except:
                pass
            try:
                iraf.images.imutil.imdelete(images = res)
            except:
                pass
            flag = True;
            break;

        print('pix_start = {}'.format(pix_start))
        print('pix_end = {}'.format(pix_end))
        iraf.images.imutil.imarith(operand1=op1, op='+', operand2=op2, result=res, verbose='No', mode='ql')

        '''
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        Open FITS File
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        '''
        hdu = fits.open(res)
        galaxy = hdu[0].data
        try:
            if len(galaxy[0] > 1):
                galaxy = galaxy[0]
        except:
            pass
        print('galaxy = {}'.format(galaxy))
        '''
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
            iraf.images.imutil.imcopy(input = res, output = temp)
            iraf.images.imutil.imdelete(images = res)
            pix_end = pix_end + 1
            op1 = '{}[0]'.format(temp)
            op2 = '{}[*,{}]'.format(file, pix_end)
        else:
            out_file = '{}_{}kpcR.fits'.format(file, kpc)
            iraf.images.imutil.imcopy(input=res, output=out_file)
            iraf.images.imutil.imdelete(images=res)
            # Update Header Information
            iraf.imutil.hedit(images=out_file, fields='SNR', value=SNR)
            iraf.imutil.hedit(images=out_file, fields='ymin_pixel', value=pix_start)
            iraf.imutil.hedit(images=out_file, fields='ymax_pixel', value=pix_end)
            pix_start = pix_end + 1
            pix_end = pix_start + 1
            op1 = '{}[*,{}]'.format(file, pix_start)
            op2 = '{}[*,{}]'.format(file, pix_end)

    return None

# Set extraction mode
mode = str(raw_input("Please enter extraction mode. Enter 'B' for Both, 'L' for Left and 'R' for Right:  "))
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

# Extract regions
if mode == 'B':
    extract_left(path, file, desired_SNR, centre, pix_per_kpc)
    extract_1kpc_left(path, file, centre, pix_per_kpc)
    extract_right(path, file, desired_SNR, centre, pix_per_kpc)
    extract_1kpc_right(path, file, centre, pix_per_kpc)
elif mode == 'L':
    extract_left(path, file, desired_SNR, centre, pix_per_kpc)
    extract_1kpc_left(path, file, centre, pix_per_kpc)
elif mode == 'R':
    extract_right(path, file, desired_SNR, centre, pix_per_kpc)
    extract_1kpc_right(path, file, centre, pix_per_kpc)
