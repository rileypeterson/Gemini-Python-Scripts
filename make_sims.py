### User options
### Riley Peterson October 2017

### I need to credit this website for teaching me the inverse transform sampling
### http://www.nehalemlabs.net/prototype/blog/2013/12/16/how-to-do-inverse-transformation-sampling-in-scipy-and-numpy/
### Though ultimately I didn't use this, it was interesting and useful nevertheless

### This is here, so that make pairs works
import numpy as np

### Logistics
base_folder = "/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/" #directory where simulations will be made, created if it doens't exist
login_cl_path = "/Users/rpeterson/" #so that from pyraf import iraf is imported properly
psf_image = "/Users/rpeterson/Downloads/uds_2epoch_f160w_v0.3_psf.fits" #psf kernel that will be convolve with the galaxies
galfit_binary = "/net/frigg/Users/inger/bin/galfit" #path to GALFIT, if you can run GALFIT with "$galfit obj.in", just put galfit here



### Make input files for your galaxies
make_feedmes = 'yes' # if "yes", feedmes are made for the number of galaxies in number_of_sims. If running make_pairs I believe this needs to be "yes".
feedme_verbose = "no"

large_field_xdim,large_field_ydim = 1014,1014 # x,y ; the dimensions of your simulated image
mag_range = [16,28] # The next five parameters are ranges for the parameters, if make_interpolate_chain = "yes" and the appropriate variables are defined these are interpolated, else selected random from within the ranges
re_range = [2,40]
n_range = [0.5,5.5]
q_range = [0.05,1.0]
pa_range = [-90,90]

number_of_sims = 1037 # Should be even int, number of simulated galaxies in the image, if make_pairs="yes" and the maximum pixel distance is large it might not be possible to place this number, the script will inform the user if this happens
cutout_size = 400 # Should be even int, length size of the cutouts created (recommend large ~400 so that none of the light is accidentally clipped if it is a large galaxy)
sky_range = [0,0] #Range of sky values, I generally leave this as [0,0], works like the ranges above but is for the sky
mag_zp = 25.9463 #Magnitude zeropoint
#plt_scale_ratio=2.1375 #I don't think this is necessary anymore. For drizzling, since my flts had plate scale (0.13) and final drizzle was (0.06) I was trying to simulate galaxies in the flt using data from the drz, therefore their radii had to be shrunk. I figured out how to use the blot task to create an flt from the drz, so this happens then.
universal_seed = 6061944 ### can be None or int, super helpful for reproducibility. Same seed produces same results
    
make_interpolate_chain = "yes" # if "yes", values are interpolated from SExtractor output tables. Currently, needs to be "yes" if making pairs.
if make_interpolate_chain=="yes":
    #first parameter to be interpolated
    chain_type1 = "cubic"
    chain_plot1 = "yes"
    chain_bins1 = 40
    chain_data1 = "/Users/rpeterson/FP/GALFIT_Simulations1/SExtractor_missing_problem/actual_data_trial_field/UV105842F160W/UV105842F160W_phot/UV105842F160W_tab_all.fits" #Only works for fits files right now, might make for text later
    chain_column1 = "MAG_AUTO"
    chain_obj1 = "mag"
    
    do_chain2 = "yes"
    if do_chain2=="yes":
        #second parameter which is interpolated
        chain_type2 = "linear"
        chain_plot2 = "yes"
        chain_data2 = "/Users/rpeterson/FP/GALFIT_Simulations1/SExtractor_missing_problem/actual_data_trial_field/UV105842F160W/UV105842F160W_phot/UV105842F160W_tab_all.fits" #Only works for fits files right now, might make for text later
        chain_column2 = "FLUX_RADIUS"
        chain_obj2 = "re"
        #if "re" is chain_obj2, will assume that "ELLIPTICITY" exists within chain_data2, and will interpolate the circularized effective radius
        
    #chain_type3 = "linear"
    #chain_plot3 = "yes"
    #chain_data3 = "/Users/rpeterson/FP/UV-105842_id6p01/UV-105842_id6p01_phot/uv-105842_d6p-01-284.0-f160wtab_all.fits" #Only works for fits files right now, might make for text later
    #chain_column3 = "FLUX_RADIUS"
    #chain_obj3 = "re"
    
### Makes pairs of galaxies. First half are "parents" randomly placed, second half are "childs" determined by stepping through these intervals
### The last argument in each of these multiplied together should be half the number of sims
### For example pixel_distance_steps=np.linspace(5,55,10),luminosity_ratio_steps=np.linspace(0.5,2,4),re_ratio_steps=np.linspace(0.2,3,5),n_steps=np.linspace(1,4,2)
### My number of sims is 800, half of this is 400. 10*4*5*2=400
make_pairs="no"
if make_pairs=="yes":
    pixel_distance_steps=np.linspace(5,55,10) #These can also be custom non linear arrays
    luminosity_ratio_steps=np.linspace(0.5,2,4) #2 means child galaxy will be twice as bright
    interp_re="yes" # if "yes" the re value will be interpolated from the make_interpolate_chain
    if interp_re=="no":
        re_ratio_steps=np.linspace(0.2,3,5) #ratio is circularized effective radius, so it factors in axis ratio, 2 means child will be twice as large (circularized re)
    n_steps=np.linspace(1,4,2)
    
    
    
### Makes the galaxies by sending their feedme/input files to the terminal
make_galaxies = 'no' 
make_galaxies_verbose = "no"

### Convolves with the provided PSF using from astropy.convolution import convolve_fft
convolve = 'no'
convolution_verbose = "yes"

### Adds noise to individual cutouts
mknoise_cutouts = 'no' #Should probably be yes if adding to existing image
mknoise_cutouts_verbose = "yes"

if mknoise_cutouts=="yes":
    background = 0.136392
    gain = 6529.37744
    rdnoise = 4
    poisson = 'yes'

display_5 = 'no' ### if "yes", displays 5 galaxies at random and shows their feedme files


### Large field
create_large_field = 'yes' # if "yes", allows the large field to be created
display_large_field = "yes" 

if create_large_field=="yes":
    image_to_sim = "default" #should be full path to it, shouldn't be in base_folder. If not "default", the image provided will be used to determine possible coordinates (useful if image is rotated)
    image_to_sim_ext=0 #default is 0, if image_to_sim has multiple extensions this is the extensions number you want to simulate
    image_to_add_to="default" #should be full path to it. If not "default" galaxies will be added to this image
    sky_value = 0 # Should probably be zero, if nonzero it sets all the points in your simulated image that aren't zero in the original to sky_value
    
    add_mknoise_large1 = "no" # I doubt this will ever be used. Probably should be set to "no". Adds noise to a blank image.
    
    if add_mknoise_large1=="yes":
        background1 = 0.136392
        gain1=6529.37744
        rdnoise1=2
        poisson1="yes"
        
    add_models = "yes" #if "yes", models will be added to the image
    add_models_verbose = "yes" #my preference is to leave this as "yes"
    label = "yes" #if "yes" a region file will be created and if display_large_field="yes" the region file is displayed. Note their is a 3 pixel offset so that labels don't block the actual galaxy.
    label_verbose = "no" 
    save = "no" #if "yes", saves final output
    if add_models == "yes":
        which_to_add="_convolved" #options are "", "_convolved", "_noised" or "_convolved_noised" #the type of model that will be added
        how_many_to_add = 1037 # should be <= number_of_sims, number of models that will be added to the large field
        not_near_edge="yes" # if "yes", galaxies are prevented from going near the edge of the image according to the pixel distance clearance
        if not_near_edge=="yes":
            clearance=50

        ### Adds noise after galaxies have been placed in the image, should probably be no if adding to existing images
        add_mknoise_large2="yes" #If drizzling this should be "no", but the values below should be defined. if "yes", then drizzle should be set to "no", and noise is added using IRAF's mknoise.
        background2 = 0.623
        gain2=1632.34375/2.5 #mknoise gain
        rdnoise2=20 #mknoise read noise
        poisson2="yes" #mknoise add poisson noise
            
        drizzle="yes" #if "yes", drizzling is performed
        if drizzle=="yes":
            path_to_flts="/Users/rpeterson/Downloads/id6p01???_flt.fits" #Path to your flt images wildcard (*) or question mark (?) should be placed where necessary. E.g. path_to_flts="/Users/rpeterson/Downloads/id6p01???_flt.fits", such that glob.glob(path_to_flts) will get all your flts
            path_to_actual_driz="/Users/rpeterson/FP/do_everything/UV105842F160W/UV105842F160W_drz.fits"
            display_blots="yes" #Will display blotted (distorted flts)
            backgrounds=[0.6351,0.606,0.6373,0.606] #If "default" will use the background2 value above as the value for each flt, otherwise it is the list of specific values for each simulated flt e.g. backgrounds=[0.6351,0.606,0.6373,0.606]
            
            ###Originally I had the sky comparison be for user specified regions but this became way too complicated, instead the entire image is compared (histogram ranges are determined by range_sky=(im1_med-3*im1_std,im1_med+3*im1_std) see definition sky_compare for more details)
            
            compare_noised_blots_sky="yes" #Will display the blotted noised images
            #Astrodrizzle Parameters, default is to just make the finally drizzled image without sky subtraction and to clean intermediate images, the user can edit this within the actual code to suite their needs
            final_bits,final_pixfrac,final_fillval,final_scale,final_rot,resetbits=576,0.8,0,0.06,0,0

            compare_final_drz="yes"


    compare_sky="yes" #if "yes", compares sky of final_image for example drizzle="no" and add_mknoise_large2="yes" to image to sim. If image_to_sim is default makes a histogram of the simulate image
    #If you really want to compare specific regions use the commented part below, otherwise will use the whole image and clip the sky according to the median
    """
    if compare_sky=="yes":
        real_sky_lowerx,real_sky_lowery, real_sky_upperx,real_sky_uppery=489,690,605,765    ### Eight coordinates here, if any==None, then both images are display, and user is prompted for values
        sim_sky_lowerx,sim_sky_lowery, sim_sky_upperx,sim_sky_uppery=None,446,593,498 
    """


############################################################ Cross with caution        

### Imports
import random
import os
###Check for slash at end of base_folder
if base_folder[-1]!="/":
    print("Adding a slash to the end of the base_folder")
    base_folder=base_folder+"/"
if os.path.exists(base_folder)==False:
    os.mkdir(base_folder)
import subprocess
import time
from astropy.io import fits
from astropy.convolution import convolve_fft
import sys
#from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interpolate
from drizzlepac import astrodrizzle, ablot #If you don't have these and don't want to drizzle you can comment them out
from stwcs import wcsutil
import glob
import pandas as pd #in case I want to incorporate the master_stats.py definitions into here



os.chdir(login_cl_path)
from pyraf import iraf
iraf.artdata(_doprint=0)
os.chdir(base_folder)

plt_fig=1
plt.close("all")
np.random.seed(universal_seed)
random.seed(universal_seed)
chain_seed1 = universal_seed 
chain_seed2 = universal_seed
seed1=universal_seed
seed2=universal_seed
seed = universal_seed
if universal_seed is None:
    seed,seed1,seed2="INDEF"

### Definitions
def run(source, reference, output=None, interp='poly5', sinscl=1.0, stepsize=10, clobber=False):
#This script given by Warren Hack (hack@stsci.edu), STScI Science Software Branch, 14 Sept 2015
#Modified by Ryan Cole, 17 Sept. 2015
#Modified by Riley Peterson, 5 Dec. 2017
    """ Blot the source image (simple FITS format) to match the WCS of the reference image

    Parameters
    ==========
    source : string
        Filename of input/source image to be blotted
    reference : string
        Filename of image which the input image will be blotted to match
    output : string
        Filename of output blotted image
    clobber : boolean
        Specify whether or not to overwrite (delete then replace) the output file if a file
        with that same name already exists on disk
    
        
    ASSUMPTIONS:
      1. All files (source, reference and blotted image)are in simple FITS format.
      2. WCS information does not include NPOL or D2IM corrections 
            (simple FITS does not support such corrections)
      3. Source and reference files have the same units, and exptime can be set to 1 
            to match count rate units.
    """
    if output:
        # Insure output file does not already exist
        if os.path.exists(output) and not clobber:
            print 'Quitting... Output file already exists. Please rename or delete, then restart.'
            return
        if os.path.exists(output) and clobber:
            print 'Removing ',output,' to replace with new result...'
            os.remove(output)

    # Load data and WCS from "input" source image
    #My simulated images don't have WCS coordinates so this doesn't make sense and I edit it so that the source has the same wcs as the reference
    src = fits.getdata(source)
    src_wcs = wcsutil.HSTWCS(source)

    # Load WCS from reference image
    blot_wcs = wcsutil.HSTWCS(reference)
    out_wcs = blot_wcs.copy()
    print(out_wcs)
    
    # Define exptime to scale source drizzled image to reference image units
    exptime = 1

    # Blot src to reference
    outsci = ablot.do_blot(src, src_wcs, out_wcs, exptime, coeffs=True, 
                interp=interp, sinscl=sinscl, stepsize=stepsize) #,
    #            wcsmap=drizzlepac.wcs_functions.WCSMap)

    if output:           
        hdr = fits.getheader(source)
        #return hdr
        # Update WCS information of source image with 
        # WCS values from reference/output WCS
        wcs_hdr = blot_wcs.wcs2header(sip2hdr=True)
        for kw in wcs_hdr:
            hdr[kw]=wcs_hdr[kw]
            #hdr.update(kw,wcs_hdr[kw])
        hdr["VAFACTOR"]=blot_wcs.vafactor
        hdr["ORIENTAT"]=blot_wcs.orientat
        hdr["CTYPE1"]="RA---TAN-SIP"
        hdr["CTYPE2"]="DEC--TAN-SIP"
        #hdr.update('VAFACTOR',blot_wcs.vafactor)
        #hdr.update('ORIENTAT',blot_wcs.orientat)

        # Create output file here...
        phdu = fits.PrimaryHDU(data=outsci,header=hdr)
        # write out file now...
        phdu.writeto(output)

    return outsci    
    
###Analyzes sky of one or two images
def sky_compare(image1,x1,y1,x2,y2,ext1=0,image2=None,x1_2=None,y1_2=None,x2_2=None,y2_2=None,ext2=0,close="no",fig_num=100,clean2="no",clean1="no"):
    plt.figure(fig_num)
    if close=="yes":
        plt.close("all")
    im1_full_array=fits.open(image1,memmap=True)[ext1].data
    name1=image1.split("/")[-1]
    if clean2=="no" or image2 is None:
        im1_crop_array=im1_full_array[y1:y2,x1:x2].flat
        im1_mean,im1_std,im1_med=np.mean(im1_crop_array),np.std(im1_crop_array),np.median(im1_crop_array)
        if clean1=="yes":
            im1_crop_array=im1_crop_array[im1_med/3<im1_crop_array] #This is here to exclude bad pixels
            im1_crop_array=im1_crop_array[im1_crop_array<im1_med*1.66]
        im1_mean,im1_std,im1_med=np.mean(im1_crop_array),np.std(im1_crop_array),np.median(im1_crop_array)
        range_sky=(im1_med-3*im1_std,im1_med+3*im1_std)
        plt.hist(im1_crop_array,label=name1+" data "+str(len(im1_crop_array))+" mean:"+str(im1_mean)+" std:"+str(im1_std)+" median:"+str(im1_med),
                 alpha=0.5,bins=100,range=range_sky,normed=True) 
        print(name1+" mean,stddev,median: "+str(im1_mean)+" "+str(im1_std)+" "+str(im1_med))
    if image2 is not None:
        im2_full_array=fits.open(image2,memmap=True)[ext2].data
        name2=image2.split("/")[-1]
        im2_crop_array=im2_full_array[y1_2:y2_2,x1_2:x2_2].flat
        im2_mean,im2_std,im2_med=np.mean(im2_crop_array),np.std(im2_crop_array),np.median(im2_crop_array)
        if clean2=="yes":
            im2_crop_array=im2_crop_array[im2_med/3<im2_crop_array] #This is here to exclude bad pixels
            im2_crop_array=im2_crop_array[im2_crop_array<im2_med*1.66] #This is here to exclude bad pixels in the flts
            im2_mean,im2_std,im2_med=np.mean(im2_crop_array),np.std(im2_crop_array),np.median(im2_crop_array)
        range_sky=(im2_med-3*im2_std,im2_med+3*im2_std)
        plt.hist(im2_crop_array,label=image2.split("/")[-1]+" data "+str(len(im2_crop_array))+" mean:"+str(im2_mean)+" std:"+str(im2_std)+" median:"+str(im2_med),
                 alpha=0.5,bins=100,range=range_sky,normed=True)
        print(name2+" mean,stddev,median: "+str(im2_mean)+" "+str(im2_std)+" "+str(im2_med))
        if clean2=="yes":
            im1_crop_array=im1_full_array[y1:y2,x1:x2].flat
            im1_crop_array=im1_crop_array[im2_med/3<im1_crop_array] #This is here to exclude bad pixels
            im1_crop_array=im1_crop_array[im1_crop_array<im2_med*1.66] #This is here to exclude bad pixels in the flts
            im1_mean,im1_std,im1_med=np.mean(im1_crop_array),np.std(im1_crop_array),np.median(im1_crop_array)
            plt.hist(im1_crop_array,label=name1+" data "+str(len(im1_crop_array))+" mean:"+str(im1_mean)+" std:"+str(im1_std)+" median:"+str(im1_med),
                     alpha=0.5,bins=100,range=range_sky,normed=True)
            print(name1+" mean,stddev,median: "+str(im1_mean)+" "+str(im1_std)+" "+str(im1_med))
    plt.legend(loc="upper right")
    plt.xlabel("sky")
    plt.ylabel("number in bin (normalized)")  
    plt.show()
    
#Builds top of feedme file
def top(file1,output_name,cutout_size,mag_zp):
    file1.write("#================================================================================"+'\n')
    file1.write("# IMAGE PARAMETERS"+'\n')
    file1.write("B) "+output_name+"   # Output data image block"+'\n')
    file1.write("H) 1 "+str(cutout_size)+" 1 "+str(cutout_size)+" # Image region to fit (xmin xmax ymin ymax)"+'\n')
    file1.write("J) "+str(mag_zp)+"             # Magnitude photometric zeropoint"+'\n')
    file1.write("O) regular                # Display type (regular, curses, both)"+'\n')
    file1.write("P) 1                   # Create ouput only? (1=yes; 0=optimize)"+'\n')
    file1.write("S) 0                   # Modify/create objects interactively? )"+'\n')
    file1.write('\n')

#Builds component 1 randomly
def comp1(file1,cutout_size,mag_range,re_range,n_range,q_range,pa_range,plt_scale_ratio=1):
    mag=str(random.uniform(mag_range[0],mag_range[1]))
    re=str(random.uniform(re_range[0],re_range[1])/plt_scale_ratio)
    n=str(random.uniform(n_range[0],n_range[1]))
    q=str(random.uniform(q_range[0],q_range[1]))
    pa=str(random.uniform(pa_range[0],pa_range[1]))
    file1.write("#Simulated galaxy parameters"+'\n')
    file1.write(" 0) sersic                 #  Object type"+'\n')
    file1.write(" 1) "+str(cutout_size/2)+" "+str(cutout_size/2)+"  1 1  #  position x, y"+'\n')
    file1.write(" 3) "+mag+"     1          #  total magnitude"+'\n')
    file1.write(" 4) "+re+"      1          #      R_e [pixels]"+'\n')
    file1.write(" 5) "+n+"     1          #  exponent (de Vaucouleurs = 4)"+'\n')
    file1.write(" 9) "+q+"      1          #  axis ratio (b/a)"+'\n')
    file1.write("10) "+pa+"     1          #  position angle (PA) [deg: Up=0, Left=90]"+'\n')
    file1.write('\n')
    return mag,re,n,q,pa
    
#Builds sky according to sky range
def sky(file1,sky_range):
     file1.write("0) sky                    #    Object type"+'\n')
     file1.write("1) "+str(random.uniform(sky_range[0],sky_range[1]))+"     0       #  sky background"+'\n')
     file1.write('\n')
 
#Used for pairs approach to get the coordinates within a given distance and x,y    
def get_circle_coords(distance,x,y):
    xmin,xmax,ymin,ymax=x-distance,x+distance,y-distance,y+distance
    vals=[]
    for x1 in range(xmin,xmax+1):
        for y1 in range(ymin,ymax+1):
            distance_val=((x1-x)**2+(y1-y)**2)**0.5
            if distance_val<=distance:
                vals.append((x1,y1,distance_val))
    return vals



#Builds the feedme file for the interpolated values, this could probably made simpler
if make_interpolate_chain=="yes":
    def comp1_interp_chain(file1,cutout_size,mag_range,re_range,n_range,q_range,pa_range,chain_obj1,chain_val1,
                           chain_obj2,chain_val2,plt_scale_ratio=1,manual_overwrite="no",mag_val=None,re_val=None,n_val=None,q_val=None):
        file1.write("#Simulated galaxy parameters"+'\n')
        file1.write(" 0) sersic                 #  Object type"+'\n')
        file1.write(" 1) "+str(cutout_size/2+1)+" "+str(cutout_size/2+1)+"  1 1  #  position x, y"+'\n')
        mag=str(random.uniform(mag_range[0],mag_range[1]))
        re=str(random.uniform(re_range[0],re_range[1]))
        n=str(random.uniform(n_range[0],n_range[1]))
        ar=str(random.uniform(q_range[0],q_range[1]))
        pa=str(random.uniform(pa_range[0],pa_range[1]))
        ### do chain1
        if chain_obj1=="mag":
            mag=str(chain_val1)
        if chain_obj1=="re":
            re=str(chain_val1)
        if chain_obj1=="n":
            n=str(chain_val1)
        if chain_obj1=="ar":
            ar=str(chain_val1)
        if chain_obj1=="pa":
            pa=str(chain_val1)
        
        ### do chain2
        print("circularized effective radius",re)
        if chain_obj2 is not None and chain_val2 is not None:
            if chain_obj2=="mag":
                mag=str(chain_val2)
            if chain_obj2=="re":
                re=str(chain_val2/float(ar)**0.5) ###Because the chain_valued re is the circularized effective radius and we need the effective semi-major axis as input for GALFIT
            print("semi-major radius",re)
            if chain_obj2=="n":
                n=str(chain_val2)
            if chain_obj2=="ar":
                ar=str(chain_val2)
            if chain_obj2=="pa":
                pa=str(chain_val2)
        """
        ### do chain3
        if chain_obj3!=None and chain_array3!=None:
            if chain_obj3=="mag":
                mag=str(np.random.choice(chain_array3))
            if chain_obj3=="re":
                re=str(np.random.choice(chain_array3))
            if chain_obj3=="n":
                n=str(np.random.choice(chain_array3))
            if chain_obj3=="ar":
                ar=str(np.random.choice(chain_array3))
            if chain_obj3=="pa":
                pa=str(np.random.choice(chain_array3))
        """


        
        if manual_overwrite=="yes":
            mag,re,n,ar=str(mag_val),str(re_val),str(n_val),str(q_val)
        #mag,re,n,ar,pa="24.5993","4.9040","0.8374","1.0","-37.8909"
        #re=str(float(re)/plt_scale_ratio)
        file1.write(" 3) "+mag+"     1          #  total magnitude"+'\n')
        file1.write(" 4) "+re+"      1          #      R_e [pixels]"+'\n')
        file1.write(" 5) "+n+"     1          #  exponent (de Vaucouleurs = 4)"+'\n')
        file1.write(" 9) "+ar+"      1          #  axis ratio (b/a)"+'\n')
        file1.write("10) "+pa+"     1          #  position angle (PA) [deg: Up=0, Left=90]"+'\n')
        file1.write('\n')
        return mag,re,n,ar,pa

#Perform interpolation for chain1
def make_chain_vals1(chain_data1,chain_bins1,chain_seed1,chain_type1,chain_obj1,chain_column1,plt_fig=plt_fig):
    ###Open chain_data1 and get a histogram of the values
    chain_data1 = fits.open(chain_data1)[1].data[chain_column1]
    chain_hist1, chain_bins1 = np.histogram(np.sort(chain_data1), bins=chain_bins1, normed=True, density=True)
    
    ### There are two ways to proceed personally I think the second way is better
    ### Way #1 inverse transform sampling:
    
    chain_values1 = np.zeros(chain_bins1.shape)
    chain_values1[1:] = np.cumsum(chain_hist1*np.diff(chain_bins1))
    
    """
    ### There are two ways to proceed personally I think the second way is better
    ### Way #1 inverse transform sampling:
    epsilon=0
    for i in range(0,len(chain_values1)-1):
        if chain_values1[i]>=chain_values1[i+1]:
            chain_values1[i+1]=chain_values1[i+1]+1*10**-11+epsilon
            epsilon=epsilon+1*10**-12
    #print(chain_values1)
    #print(chain_bins1)
    chain_inv_cdf1 = interpolate.interp1d(chain_values1, chain_bins1, kind=chain_type1)
    np.random.seed(chain_seed1)
    chain_r1 = helper           #Should ensure that we make a random choice every time but from the same distribution
    mag_array = chain_inv_cdf1(chain_r1)
    plt.hist(chain_data1, bins=chain_bins1, alpha=0.4, normed=True, color='b',label="original data",zorder=2)
    plt.hist(chain_inv_cdf1(chain_r1), bins=chain_bins1, alpha=0.5, normed=True, color='g', label="interpolated data",zorder=1)
    plt.legend(loc="upper right")
    plt.xlabel(chain_obj1)
    plt.ylabel("Probability")
    plt.title("Chain1 Results")
    """
    
    ### Way #2 Interpolate PDF of data from histogram, finely sample from this to get
    ### a list of data values and their interpolated weight then randomly select from this distribution (with weighting)
    
    chain_linspace1 = np.linspace(np.min(chain_data1),np.max(chain_data1),chain_values1.size-1)
    chain_pre_interp1 = interpolate.interp1d(chain_linspace1, chain_hist1*np.diff(chain_bins1), kind=chain_type1, fill_value="extrapolate") #interpolates PDF
    x_chain=np.linspace(np.min(chain_data1),np.max(chain_data1),100000) #Creates discrete range of possible data values
    p = chain_pre_interp1(x_chain)/np.sum(chain_pre_interp1(x_chain))
    p[p<0]=0
    mag_array = np.random.choice(x_chain,size=number_of_sims,p=p/np.sum(p))
    chain_sim_data = np.random.choice(x_chain,size=100000,p=p/np.sum(p))
    plt.figure(plt_fig)
    plt.hist(mag_array,bins=chain_bins1,alpha=0.5,normed=True,label="simulated data sample size "+str(number_of_sims),color="r",zorder=3)
    plt.plot(chain_sim_data[chain_sim_data.argsort()],chain_pre_interp1(chain_sim_data)[chain_sim_data.argsort()]/np.diff(chain_bins1)[0], label="simulated data", color="g",zorder=1)
    plt.hist(chain_data1, bins=chain_bins1, alpha=0.5, normed=True, color='b',label="original data",zorder=2)
    #plt.plot(x,chain_pre_interp1(x))
    plt.legend(loc="upper right")
    plt.xlabel(chain_obj1)
    plt.ylabel("Probability")
    plt.title("Chain1 Results")

    return mag_array

def make_chain_func2(chain_data1,chain_column1,chain_array1,chain_obj1,chain_data2,chain_type2,chain_obj2,chain_column2,plt_fig=plt_fig):
    ###returns the function that interpolates between chain1 and chain2
    chain_data2_name=chain_data2
    chain_data2 = fits.open(chain_data2)[1].data[chain_column2]
    chain_data2 = chain_data2.astype("f8")
    if chain_obj2=="re":
        chain_data_q2 = fits.open(chain_data2_name)[1].data["ELLIPTICITY"]
        chain_data_q2 = (1-chain_data_q2.astype("f8"))**0.5 ###get sqrt of axis ratio
        chain_data2=np.multiply(chain_data_q2,chain_data2) ###multiply by FLUX_RADIUS to get a circularized effective ardius
    chain_data1 = fits.open(chain_data1)[1].data[chain_column1]
    chain_data1 = chain_data1.astype("f8") 
    chain_interp_func2 = interpolate.interp1d(chain_data1, chain_data2, kind=chain_type2, fill_value="extrapolate")
    chain_x2 = np.linspace(np.min(chain_data1),np.max(chain_data1),100000)
    plt.figure(plt_fig+1)
    plt.plot(chain_x2, chain_interp_func2(chain_x2),color='g',alpha=0.5,label="interpolated function",zorder=1)
    plt.scatter(chain_data1, chain_data2, label="original data",alpha=0.5,zorder=2, color='b')
    plt.scatter(chain_array1, chain_interp_func2(chain_array1), label="simulated data sample size "+str(number_of_sims), alpha=0.5, zorder=3, color='r')
    #df11=pd.DataFrame.from_dict(np.load("/Users/rpeterson/FP/GALFIT_Simulations1/SExtractor_missing_problem/inputs_dict.npy").item()).T.drop("seed")
    #foo=df11[~df11.id.isin(b.values())]   #Here b is SEx_to_sim as determined by master_stats.py
    #missed_ids,missed_mags,missed_res=foo["id"],foo["mag"].astype(float),np.multiply(foo["re"].astype(float),foo["q"].astype(float)**0.5)
    #print(len(missed_ids))
    #plt.scatter(missed_mags,missed_res,label="SExtractor misses", alpha=0.5, zorder=4, color='k')
    
    #[plt.annotate(obj[0],xy=(obj[1],obj[2])) for obj in zip(missed_ids,missed_mags,missed_res)]
    #plt.scatter(chain_val1, chain_val2, label="simulated data points")
    plt.legend(loc="upper right")
    plt.xlabel(chain_obj1)
    plt.ylabel(chain_obj2)
    plt.title("Chain2 Results")
    
    return chain_interp_func2
    
    
    
    

    
### Things that don't make sense
stop_button="no"
try:
    if how_many_to_add>number_of_sims:
        print("Warning:how_many_to_add>number_of_sims. Please set them equal or to values that make sense.")
        ask = raw_input("Continue? If 'y' then they will be set equal(y/n):")
        if ask=="y":
            how_many_to_add=number_of_sims    
        else:
            stop_button="yes"
except:
    pass
if stop_button=="yes":
    sys.exit()


### Prep for make_pairs
if make_pairs=="yes":
    pairs=[]
    if interp_re=="yes":
        re_ratio_steps=np.linspace(0,0,1)
    for a in pixel_distance_steps:
        for b in luminosity_ratio_steps:
            for c in re_ratio_steps:
                for d in n_steps:
                    pairs.append((a,b,c,d))
    if len(pairs)>number_of_sims/2:
        print("The length of pair values:"+str(len(pairs))+" is greater than half the number of sims. Either increase the number of sims or decrease the number of possible pair values.")
        sys.exit()
    
    
    
### Making input feedme files
initial_time = time.time()
chain_val1_list=[]
chain_val2_list=[]
inputs_dict_file=base_folder+"inputs_dict.npy"
k=1
inputs_dict=dict()
inputs_dict.update({"seed":universal_seed})
if make_feedmes=='yes':
    if make_interpolate_chain=="yes":
        ### Make chain values 1 and the interpolation function for chain values 2

        chain_array1=make_chain_vals1(chain_data1,chain_bins1,chain_seed1,chain_type1,chain_obj1,chain_column1,plt_fig=plt_fig)
        if do_chain2=="yes":
            chain_function2= make_chain_func2(chain_data1,chain_column1,chain_array1,chain_obj1,chain_data2,chain_type2,chain_obj2,chain_column2,plt_fig=plt_fig)
            
    for i in range(1,number_of_sims+1):
        if feedme_verbose=="yes":
            print("Creating feedme for galaxy number "+str(i)+" elapsed time is:"+str(round((time.time()-initial_time),4)))
        output_name = "galout_"+str(i)+"_model.fits"
        file1=open(base_folder+"obj_sim_"+str(i)+".in","w")
        top(file1,output_name,cutout_size,mag_zp)
        input_x,input_y=cutout_size/2,cutout_size/2
        if make_interpolate_chain=="yes":
            chain_val1=chain_array1[i-1]
            if do_chain2=="yes":
                chain_val2=chain_function2(chain_val1) ###Remember this is circularized effective radius
            else:
                chain_obj2=None
                chain_val2=None
            if make_pairs=="no" or i<=number_of_sims/2:
                input_mag,input_re,input_n,input_ar,input_pa=comp1_interp_chain(file1,cutout_size,mag_range,re_range,n_range,q_range,pa_range,chain_obj1,chain_val1,
                               chain_obj2=chain_obj2,chain_val2=chain_val2,plt_scale_ratio=1)
                inputs_dict[i]={"id":i,"x_cutout":input_x,"y_cutout":input_y,"mag":input_mag,"re":input_re,"n":input_n,"q":input_ar,"pa":input_pa}
            if make_pairs=="yes" and i>number_of_sims/2:
                mag_parent,re_parent,n_parent,q_parent=float(inputs_dict[k]["mag"]),float(inputs_dict[k]["re"]),float(inputs_dict[k]["n"]),float(inputs_dict[k]["q"])
                ###Want the radius to factor in q. Convert from a_eff to r_eff
                circ_re_parent=re_parent*q_parent**0.5
                ###Get random q (axis ratio) value for child
                q_child=random.uniform(q_range[0],q_range[1])
                ###Find r_eff for child based on the ratio in pairs
                circ_re_child=pairs[k-1][2]*circ_re_parent
                #circ_re_child=pairs[k-1][2]*circ_re_parent*(q_parent**0.5)/(q_child**0.5)
                ###Find child a_eff
                re_child=circ_re_child/(q_child**0.5)
                mag_child,n_child=-2.5*np.log10(pairs[k-1][1])+mag_parent,pairs[k-1][3]
                #print(n_child)
                if interp_re=="yes":
                    circ_re_child=chain_function2(mag_child)
                    re_child=circ_re_child/(q_child**0.5)
                input_mag,input_re,input_n,input_ar,input_pa=comp1_interp_chain(file1,cutout_size,mag_range,re_range,n_range,q_range,pa_range,chain_obj1,chain_val1,
                           chain_obj2,chain_val2,plt_scale_ratio=1,manual_overwrite="yes",mag_val=mag_child,re_val=re_child,n_val=n_child,q_val=q_child)
                inputs_dict[i]={"id":i,"x_cutout":input_x,"y_cutout":input_y,"mag":input_mag,"re":input_re,"n":input_n,"q":input_ar,"pa":input_pa}
                k=k+1
        else:
            mag_rand,re_rand,n_rand,q_rand,pa_rand=comp1(file1,cutout_size,mag_range,re_range,n_range,q_range,pa_range,plt_scale_ratio=1)
            inputs_dict[i]={"id":i,"x_cutout":input_x,"y_cutout":input_y,"mag":mag_rand,"re":re_rand,"n":n_rand,"q":q_rand,"pa":pa_rand}
        sky(file1,sky_range)
        
        file1.close()
    np.save(inputs_dict_file,inputs_dict)
    
    
    
### Runnning GALFIT to actually create models (in parallel)
command = galfit_binary 
#initial_time = time.time()

if make_galaxies=='yes':

    for i in range(1,number_of_sims+1):
        os.chdir(base_folder)
        name = base_folder+"obj_sim_"+str(i)+".in"
        subprocess.Popen([command,name])
        if make_galaxies_verbose=="yes":
            print("Making galaxy image for number "+str(i)+" elapsed time is:"+str(round((time.time()-initial_time),4)))
    """
    def run_galfit(name):
        subprocess.Popen([command,name])    
        print("Making galaxy image for number "+str(i)+" "+str(time.time()-initial_time))
    name_list=[]
    for i in range(1,number_of_sims+1):
        name_list.append(base_folder+"obj_sim_"+str(i)+".in")
    Pool().map(run_galfit,name_list)
Pool.close()
print("done")
"""
#initial_time = time.time()
if convolve=='yes':
    #ADDRESS NEED TO NORMALIZE PSF BEFORE CONVOLVING
    for i in range(1,number_of_sims+1):
        try:
            os.remove("galout_"+str(i)+"_model_convolved.fits")
        except:
            pass
        input_name = "galout_"+str(i)+"_model.fits"
        sim = fits.open(input_name)
        psf = fits.open(psf_image)
        new_image = convolve_fft(sim[0].data,psf[0].data,normalize_kernel=True) #You can also try boundary="wrap", this seems to be faster, but I'm not sure what the difference is  
        new_image = np.float32(new_image)
        hdu = fits.PrimaryHDU(new_image)
        hdu.writeto("galout_"+str(i)+"_model_convolved.fits")
        if convolution_verbose=="yes":
            print("Performing convolution for galaxy number "+str(i)+" elapsed time is:"+str(round((time.time()-initial_time),4)))


###Make noise for individual cutouts
if mknoise_cutouts=='yes':
    for i in range(1,number_of_sims+1):
        try:
            os.remove("galout_"+str(i)+"_model_convolved_noised.fits")
        except:
            pass
        try:
            os.remove("galout_"+str(i)+"_model_noised.fits")
        except:
            pass
        input_name = "galout_"+str(i)+"_model_convolved.fits"
        iraf.cd(base_folder)
        if convolve=="yes":
            iraf.mknoise(input_name,output="galout_"+str(i)+"_model_convolved_noised.fits",background=background,gain=gain,rdnoise=rdnoise,poisson=poisson,seed=seed)
        else:
            input_name = "galout_"+str(i)+"_model.fits"
            iraf.mknoise(input_name,output="galout_"+str(i)+"_model_noised.fits",background=background,gain=gain,rdnoise=rdnoise,poisson=poisson,seed=seed)
        if mknoise_cutouts_verbose=="yes":    
            print("Adding noise for galaxy number "+str(i)+" elapsed time is:"+str(round((time.time()-initial_time),4)))

###Display 5 galaxies
if display_5=='yes':
    listo = [random.randrange(1,number_of_sims) for i in range(1,6)]
    for i in listo:
        objects=[]
        for obj in ["galout_"+str(i)+"_model.fits",psf_image,"galout_"+str(i)+"_model_convolved.fits","galout_"+str(i)+"_model_noised.fits","galout_"+str(i)+"_model_convolved_noised.fits"]:
            if os.path.exists(obj):
                objects.append(obj)
        os.system("ds9 "+" -fits ".join(objects)+" &")
        os.system("open "+"obj_sim_"+str(i)+".in &")
        




if create_large_field=='yes':
    os.chdir(base_folder)
    iraf.cd(base_folder)
    iraf.protect("*saved.fits")
    #Might be better just to overwrite instead of deleting
    try:
        iraf.delete("full*")
    except:
        pass
    try:
        iraf.delete("*.reg")
    except:
        pass
    if image_to_sim=="default":
        if drizzle=="yes":
            print("If you are trying to create a simulated drizzle image_to_sim should be your real final drizzled image")
            sys.exit()
        n = np.ones([large_field_ydim,large_field_xdim])
        hdu = fits.PrimaryHDU(n)
        hdu.writeto('full_field_with_sky.fits')
    else:
        os.system("cp "+image_to_sim+" "+base_folder)
        os.rename(image_to_sim.split("/")[-1],"full_field_with_sky.fits")
    image_now = "full_field_with_sky.fits"
    hdulist=fits.open(image_now, mode="update")
    print("Making data points that were non zero in original image equal to sky_value in sim image")
    image_data=hdulist[0].data
    image_data[image_data!=0]=sky_value
    hdulist.flush()
    hdulist.close()
    if add_mknoise_large1=="no" and add_models=="no":
        final_large_image=image_now
    ###I dont expect mknoise1 to get any use
    if add_mknoise_large1=="yes" and add_models=="no":
        print("Running mknoise1 (before galaxies have been added)")
        iraf.mknoise(image_now,output=image_now.replace(".fits","_noise.fits"), background=background1, gain=gain1, rdnoise=rdnoise1, poisson=poisson1, seed=seed1)
        final_large_image=image_now.replace(".fits","_noise.fits")
    if add_mknoise_large1=="yes" and add_models=="yes":
        print("Running mknoise1 (before galaxies have been added)")
        iraf.mknoise(image_now,output=image_now.replace(".fits","_noise.fits"), background=background1, gain=gain1, rdnoise=rdnoise1, poisson=poisson1, seed=seed1)
    background_field = "full_field_with_sky.fits"
    noised_background_field = "full_field_with_sky_noise.fits"
    if os.path.exists(noised_background_field):
        os.system("cp "+noised_background_field+" full_field_with_sky_noise_galaxies.fits")
        image_to_add_gals = "full_field_with_sky_noise_galaxies.fits"
    else:
        os.system("cp "+background_field+" full_field_with_sky_galaxies.fits")
        image_to_add_gals = "full_field_with_sky_galaxies.fits"
        
    ###Add models
    if add_models == "yes":
        if image_to_sim=="default":
            image_data=np.ones([large_field_ydim,large_field_xdim])
        else:
            hdulist=fits.open(image_to_sim,mode="readonly")
            image_data=hdulist[0].data
            hdulist.close()
        ###Get list of possible coordinates
        i,j=np.nonzero(image_data)
        possible_coords = list(zip(i,j))
        if not_near_edge=="yes":
            ###Revise possible coordinates if we don't want galaxies near the edge
            non_edge_coords=[]
            print("Ensuring galaxies aren't near the edge, this might take a minute")
            for coord in possible_coords:
                try:
                    if image_data[coord[0]+clearance,coord[1]+clearance]!=0 and image_data[coord[0]-clearance,coord[1]+clearance]!=0 and image_data[coord[0]+clearance,coord[1]-clearance]!=0 and image_data[coord[0]-clearance,coord[1]-clearance]!=0 and coord[0]-clearance>0 and coord[1]-clearance>0:
                        non_edge_coords.append(coord)
                except:
                    pass
            possible_coords=non_edge_coords
        possible_coords_copy=possible_coords
        inputs_dict=np.load(inputs_dict_file).item()
        if image_to_add_to!="default":
            os.system("cp "+image_to_add_to+" "+base_folder)
            os.rename(image_to_add_to.split("/")[-1],"full_image_to_add_to.fits")
            image_to_add_gals="full_image_to_add_to.fits"
        random.seed(universal_seed)
        n=1
        i=1
        skipped_ids=[]
        master_break="no"
        skip="no"
        how_many_to_add_copy=how_many_to_add
        parents=how_many_to_add/2
        
        ###Start going through the galaxies
        while i<=how_many_to_add:
            added_image=fits.open(image_to_add_gals,mode="update")
            added_image_array=added_image[0].data
            cutout=fits.open("galout_"+str(i)+"_model"+which_to_add+".fits")    #Issue
            cutout_array=cutout[0].data
            ###I can't remember why I did this if statement here
            if len(possible_coords)>0:
                idx = random.randint(0,len(possible_coords)-1)
                y = possible_coords[idx][0]     ### This is opposite of what you expect
                x = possible_coords[idx][1]
            else:
                idx = random.randint(0,len(possible_coords_copy)-1)
                y = possible_coords_copy[idx][0]     ### This is opposite of what you expect
                x = possible_coords_copy[idx][1]
            inputs_dict[i].update({"x_global":x,"y_global":y})
            if make_pairs=="yes" and i<=how_many_to_add/2:
                ###We will need to make sure for the first half that nothing is placed within max interdistance of pairs
                ###This is kind of ugly, I will consider rewritting
                selected_coords_so_far=[(int(inputs_dict[id]["x_global"]),int(inputs_dict[id]["y_global"])) for id in sorted(inputs_dict.keys()) if id!="seed" and id<i]
                if len(selected_coords_so_far)>0:
                    start_time=time.time()
                    while True:
                        hit_here="no"
                        for coord in selected_coords_so_far:
                            distance_val=((coord[0]-x)**2+(coord[1]-y)**2)**0.5
                            try:
                                image_data[y,x]
                            except:
                                distance_val=0
                            if distance_val<=pairs[-1][0] or image_data[y,x]==0 or y<0 or x<0:
                                #print("within radius",distance_val,x,y,len(possible_coords))
                                idx = random.randint(0,len(possible_coords)-1)
                                y = possible_coords[idx][0]     ### This is opposite of what you expect
                                x = possible_coords[idx][1]
                                hit_here="yes"
                                break
                        #print(min([((coord1[0]-x)**2+(coord1[1]-y)**2)**0.5 for coord1 in selected_coords_so_far]),x,y)
                        if hit_here=="no":
                            #print(min([((coord1[0]-x)**2+(coord1[1]-y)**2)**0.5 for coord1 in selected_coords_so_far]),x,y)
                            break
                        if (time.time()-start_time)>1:
                            print("\nWarning: It is getting difficult to place parent galaxies. Please considering reducing the maximum spacing between pairs. Figuring out the possible coordinates...")
                            impossible_coords1=[]
                            for coord1 in selected_coords_so_far:
                                box=get_circle_coords(int(pairs[-1][0]),int(coord1[0]),int(coord1[1]))
                                impossible_coords=[(q[1],q[0]) for q in box]
                                [impossible_coords1.append(r) for r in impossible_coords]
                            #impossible_coords=[get_circle_coords(int(pairs[-1][0]),int(coord1[0]),int(coord1[1]),possible_coords) for coord1 in selected_coords_so_far]
                            #impossible_coords=sum(impossible_coords, [])
                            print("Figured out the impossible coordinates")
                            #impossible_coords1=sorted([(coor[1],coor[0]) for coor in impossible_coords])
                            possible_coords=sorted(list(set(possible_coords)-set(impossible_coords1)))
                            print("Figured out the new possible coordinates, of which there are: "+str(len(possible_coords)))
                            start_time=time.time()
                            if len(possible_coords)==0:
                                print("There aren't any places where its possible to place a parent galaxy. You should reduce the maximum distance if you want to add more parents. Moving on to child galaxies")
                                print("\nCould only add: "+str(i-1))
                                parents=i-1
                                how_many_to_add=2*(i-1)
                                time.sleep(10)
                                break
                            
                        #if hit_here=="yes" and (time.time()-start_time)>2:
                        #    coords_to_exclude=[get_circle_coords(int(pairs[-1][0]),int(inputs_dict[id]["x_global"]),int(inputs_dict[id]["y_global"]),possible_coords) for id in sorted(inputs_dict.keys()) if id!="seed" and id<=i]
                        #    coords_to_exclude=[(coords_to_exclude[l][1],coords_to_exclude[l][0]) for l in range(len(coords_to_exclude))]
                        #    possible_coords=list(set(possible_coords)-set(coords_to_exclude))
                        #    print(possible_coords)
            if make_pairs=="yes" and i>how_many_to_add/2:
                possible_child_coords=get_circle_coords(int(pairs[n-1][0]),int(inputs_dict[n]["x_global"]),int(inputs_dict[n]["y_global"]))
                possible_child_coords=[l for l in possible_child_coords if round(l[2])==int(pairs[n-1][0])]
                if not_near_edge=="yes":
                    ###This is disgusting, if this throws an error, put it in a try statement as I did above, it might fail because image_data might be out of bound
                    possible_child_coords=[l for l in possible_child_coords if round(l[2])==int(pairs[n-1][0]) and image_data[l[1]+clearance,l[0]+clearance]!=0 and image_data[l[1]-clearance,l[0]+clearance]!=0 and image_data[l[1]+clearance,l[0]-clearance]!=0 and image_data[l[1]-clearance,l[0]-clearance]!=0 and l[0]-clearance>0 and l[1]-clearance>0]
                idx = random.randint(0,len(possible_child_coords)-1)
                x = possible_child_coords[idx][0]
                y = possible_child_coords[idx][1]
                
                selected_coords_so_far=[(int(inputs_dict[id]["x_global"]),int(inputs_dict[id]["y_global"])) for id in sorted(inputs_dict.keys()) if id!="seed" and id<i and id!=i-how_many_to_add/2]
                start_time=time.time()
                while True:
                    skip="no"
                    hit_here="no"
                    for coord in selected_coords_so_far:
                        distance_val=((coord[0]-x)**2+(coord[1]-y)**2)**0.5
                        try:
                            image_data[y,x]
                        except:
                            distance_val=0
                            #print("hit child except")
                        if distance_val<=pairs[-1][0] or image_data[y,x]==0 or y<0 or x<0:
                            #print("within radius",distance_val,x,y)
                            idx = random.randint(0,len(possible_child_coords)-1)
                            x = possible_child_coords[idx][0]
                            y = possible_child_coords[idx][1]
                            hit_here="yes"
                            break
                    #print(min([((coord1[0]-x)**2+(coord1[1]-y)**2)**0.5 for coord1 in selected_coords_so_far]),x,y)
                    if (time.time()-start_time)>1:
                        start_time=time.time()
                        print("Warning: It is getting difficult to place child galaxies")
                        for indexo in range(len(possible_child_coords)):
                            hit_here="no"
                            x,y=possible_child_coords[indexo][0],possible_child_coords[indexo][1]
                            #print("Trying "+str(x)+" "+str(y))
                            for coord in selected_coords_so_far:
                                distance_val=((coord[0]-x)**2+(coord[1]-y)**2)**0.5
                                try:
                                    image_data[y,x]
                                except:
                                    distance_val=0
                                if distance_val<=pairs[-1][0] or image_data[y,x]==0 or y<0 or x<0:
                                    #print("within radius",distance_val,x,y)
                                    #print("Condition met for "+str(x)+" "+str(y))
                                    idx = random.randint(0,len(possible_child_coords)-1)
                                    x = possible_child_coords[idx][0]
                                    y = possible_child_coords[idx][1]
                                    hit_here="yes"
                                    break   ###ADDRESS
                            #print("Should be trying the next coordinates")
                            #time.sleep(0.5)
                        if hit_here=="yes":
                            skip="yes"
                            print("Skipping "+str(i)+" because there is no place for it to go")
                            #sys.exit()
                            break
                        
                    if hit_here=="no":
                        break
                    
                            
                
                n=n+1
                #coords_to_exclude=[get_circle_coords(int(pairs[-1][0]),int(inputs_dict[id]["x_global"]),int(inputs_dict[id]["y_global"]),possible_coords) for id in sorted(inputs_dict.keys()) if id!="seed" and id<=i]
                #coords_to_exclude=[(coords_to_exclude[l][1],coords_to_exclude[l][0]) for l in range(len(coords_to_exclude))]
                #possible_coords=list(set(possible_coords)-set(coords_to_exclude))
            if skip=="no":
                inputs_dict[i].update({"x_global":x,"y_global":y})
                if (cutout_size)%2!=0:
                    print("Cutout_size must be even")
                    break
                xmin,xmax = x-cutout_size/2-1,x+cutout_size/2-1
                ymin,ymax = y-cutout_size/2-1,y+cutout_size/2-1
                xmini,ymini,xmaxi,ymaxi=0,0,cutout_size,cutout_size
                xmin_org = xmin
                ymin_org = ymin
                #print(cutout_array.shape)
                #print(x,y)
                #print(xmin,xmax)
                #print(ymin,ymax)
                ### Could have done a try except statement here but it would omit galaxies cut off at edge of image
                if xmin<0:
                    print("Galaxy #"+str(i)+" was off edge"+" xmin was:"+str(xmin)+" new xmin:"+str(abs(xmin)))
                    xmini=abs(xmin)
                    #cutout_array=cutout_array[:,xmini:] 
                    #print(cutout_array.shape)
                    xmin=0
                if ymin<0:
                    print("Galaxy #"+str(i)+" was off edge"+" ymin was:"+str(ymin)+" new ymin:"+str(abs(ymin)))
                    ymini=abs(ymin)
                    #cutout_array=cutout_array[ymini:,:]
                    #print(cutout_array.shape)
                    ymin=0
                if xmax>large_field_xdim:
                    print("Galaxy #"+str(i)+" was off edge"+" xmax was:"+str(xmax)+" new xmax:"+str(large_field_xdim-xmin_org))
                    xmaxi=(xmax-x+1)+(large_field_xdim-x+1)
                    #cutout_array=cutout_array[:,:xmaxi]
                    #print(cutout_array.shape)
                    xmax=large_field_xdim
                if ymax>large_field_ydim:
                    print("Galaxy #"+str(i)+" was off edge"+" ymax was:"+str(ymax)+" new ymax:"+str(large_field_ydim-ymin_org))
                    ymaxi=(ymax-y+1)+(large_field_ydim-y+1)
                    #cutout_array=cutout_array[:ymaxi,:]   
                    #print(cutout_array.shape)
                    ymax=large_field_ydim
                    
                added_image_array[ymin:ymax,xmin:xmax]=added_image_array[ymin:ymax,xmin:xmax]+cutout_array[ymini:ymaxi,xmini:xmaxi]
                cutout.close()
                if add_models_verbose=="yes":
                    print("\nFinished adding galaxy #"+str(i)+" at "+str(x)+" "+str(y)+" to "+image_to_add_gals)
                if label=="yes":            ### Not best way to do this but putting it outside loop, would need to make lists, etc.
                    if label_verbose=="yes":
                        print("Making label for "+str(i))
                    file1=open("region_file.reg","a+")
                    file1.write("text "+str(x+3)+" "+str(y+3)+" # text ={"+str(i)+"}"+'\n')
                    file1.close()
            if skip=="yes":
                skipped_ids.append(i)
            i=i+1
        if make_pairs=="yes":
            for d in skipped_ids:
                del inputs_dict[d]
            #if it broke at 688, but you wanted to add 800, it deletes 689 through 800 cause they don't exist
            for d in range(how_many_to_add+1,how_many_to_add_copy+1):
                del inputs_dict[d]
        np.save(inputs_dict_file,inputs_dict)


        added_image_array[image_data==0]=0 #mask zero values so we don't get square cutouts on diagonal edges
        added_image.close()
        #hdulist.close()
    if add_mknoise_large2=="no" and add_models=="yes":
        final_large_image=image_to_add_gals
    if add_mknoise_large2=="yes" and add_models=="yes":
        print("Running mknoise2 (after galaxies have been added)")
        """
        Thought that convolving after adding instead of individual cutouts might have some effect, but it doesn't seem to change anything
        In comments are other things are tried to obtain the detection of faint galaxies
        os.rename(image_to_add_gals,image_to_add_gals.replace(".fits","_temp.fits"))
        sim = fits.open(image_to_add_gals.replace(".fits","_temp.fits"))
        psf = fits.open(psf_image)
        new_image = convolve_fft(sim[0].data,psf[0].data,boundary="wrap",normalize_kernel=True)
        new_image = np.float32(new_image)
        sim.close()
        hdu = fits.PrimaryHDU(new_image)
        try:
            os.remove(image_to_add_gals)
        except:
            pass
        hdu.writeto(image_to_add_gals,overwrite=True)
        """
        #os.rename(image_to_add_gals,image_to_add_gals.replace(".fits","_temp.fits"))
        #iraf.imexpr("(a/0.47)", image_to_add_gals, a=image_to_add_gals.replace(".fits","_temp.fits"))
        #for i in range(1,5):
        #    iraf.mknoise(image_to_add_gals,output=image_to_add_gals.replace(".fits","_noise2"+"_"+str(i)+".fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
        #    seed2=universal_seed*2
        #seed2=universal_seed
        #iraf.imcombine(image_to_add_gals.replace(".fits","_noise2"+"_*"+".fits"), image_to_add_gals.replace(".fits","_noise2.fits"), combine="average")
        iraf.mknoise(image_to_add_gals,output=image_to_add_gals.replace(".fits","_noise2.fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
        #os.rename(image_to_add_gals.replace(".fits","_noise2.fits"),image_to_add_gals.replace(".fits","_noise2_temp.fits"))
        #iraf.imexpr("(a)", image_to_add_gals.replace(".fits","_noise2.fits"), a=image_to_add_gals.replace(".fits","_noise2_temp.fits"))
        final_large_image=image_to_add_gals.replace(".fits","_noise2.fits")
        original_image_to_add_gals=fits.open(image_to_add_gals)
        mknoise2_file=fits.open(final_large_image,mode="update")
        mknoise2_file[0].data=np.nan_to_num(mknoise2_file[0].data)
        mknoise2_file[0].data[original_image_to_add_gals[0].data==0]=0
    #sys.exit()
    ###Thought this would free up memory, but I don't think its an issue
    added_image_array,chain_array1,cutout_array,image_data,j,non_edge_coords,possible_coords,possible_coords_copy=None,None,None,None,None,None,None,None
    num=100
    #sys.exit()
    if drizzle=="yes":
        print("Drizzling")
        
        #1. Blot to distorted frames of original flts
        list_of_original_flts=glob.glob(path_to_flts)
        
        #Create simple fits for run definition, assume SCI extension is 1
        [iraf.imcopy(obj+"[1]",base_folder+obj.split("/")[-1].replace(".fits","_simple.fits")) for obj in list_of_original_flts if os.path.exists(base_folder+obj.split("/")[-1].replace(".fits","_simple.fits"))==False]
        flts_list_simple=glob.glob("*_simple.fits") #Hopefully this is only thing will name simple
        
        #Create copies of the simulated image which will be blotted
        full_with_gals=image_to_add_gals
        [iraf.imcopy(full_with_gals+"["+str(image_to_sim_ext)+"]",full_with_gals.replace(".fits","_"+str(i)+".fits")) for i in range(1,len(list_of_original_flts)+1)]
        full_with_gals_list=glob.glob(full_with_gals.replace(".fits","_?.fits"))
        #noised_list=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_noise2_?.fits")
        
        #Blot the images
        objs=zip(flts_list_simple,full_with_gals_list)
        for obj in objs:
            sim,ref_flt=fits.open(obj[1]),fits.open(objs[0][0]) #Need to update the header of the simulated image so that it works, if we replace objs[0][0] with obj[0] the relative offsets in the frames won't be introduced.
            #Use the simulated data, but the header from the zero extension of the flts, took a bit of trial and error to figure this out
            drz_hdr=fits.open(path_to_actual_driz)
            updated_sim = fits.PrimaryHDU(np.float32(sim[0].data),drz_hdr[0].header) #By giving each sim the same reference header, offsets as they appear in the original flts will be introduced
            #The np.float32 is super necessary else the kernel dies???????????? Took a long time to figure that out
            updated_sim.writeto(obj[1],overwrite=True)
            run(obj[1], obj[0], output=obj[1].replace(".fits","_distorted.fits"), interp='poly5', sinscl=1.0, stepsize=10, clobber=True)
        blots_list=glob.glob("full_*_distorted.fits")
        if display_blots=="yes":
            os.system("ds9 full_*_distorted.fits &")
        
        
        #2. Obtain WCS Coordinates (I don't know if running tweak reg is necessary but here is how I did it)
        """
        #TweakREG
        #Run SExtractor 
        config_file="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/config.sex"
        pos_list_file=open("pos.list","w")
        astdriz_file=open("astdriz.list","w")
        for i in range(1,5):
            new_config_name="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/config"+str(i)+".sex"
            new_config_file=open(new_config_name,"w")
            [new_config_file.write(line.replace("astrodrizzle","sim"+str(i))) for line in open(config_file,"r").readlines()]
            new_config_file.close()
            os.system("sex sim_flt_"+str(i)+"_new.fits -c config"+str(i)+".sex")
            pos_list_file.write("sim_flt_"+str(i)+"_new.fits\n")
            astdriz_file.write("sim_flt_"+str(i)+"_new.fits "+"sim"+str(i)+"_tabext1.cat\n")
        pos_list_file.close()
        astdriz_file.close()
        tweakreg.TweakReg("@pos.list",catfile="astdriz.list",xcol=3,ycol=4,updatehdr="yes",
                          wcsname="HILO",searchrad=2,fluxunits="mag",rfluxunits="mag",threshold=3,clean=True)
        """


        #3. Add noise with different seed values
        
        #Use these if you need custom noise characteristics
        #backgrounds,ind=[0.71,0.681,0.71,0.68],0
        #backgrounds=[0.695 for b in backgrounds]
        #readnoises=[20.200,19.7999,19.900,20.100]
        #readnoises=[r/10 for r in readnoises]
        #gains=[652.9375*0.695,652.9375*0.6725,652.9375*0.67,652.9375*0.69]
        #gains=[2611.75*2.5 for i in gains]
        ind=0
        for blot in blots_list:
            if backgrounds!="default":
                background2=backgrounds[ind]
            #rdnoise2=readnoises[ind]
            #gain2=gains[ind]
            print(background2)
            iraf.mknoise(blot,output=blot.replace(".fits","_noise2.fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
            seed2=universal_seed*seed2
            ind=ind+1
        seed2=universal_seed
        
        if compare_noised_blots_sky=="yes":
            noised_blots_list=glob.glob("full_*_distorted_noise2.fits")
            os.system("ds9 full_*_distorted_noise2.fits &")
            os.system("ds9 *_simple.fits &")
            for obj in zip(flts_list_simple,noised_blots_list):
                shape=fits.open(obj[0])[0].data.shape
                x1_real,y1_real,x2_real,y2_real=1,1,shape[0],shape[1]
                x1_sim,y1_sim,x2_sim,y2_sim=1,1,shape[0],shape[1]
                sky_compare(obj[1],x1_real,y1_real,x2_real,y2_real,ext1=0,  #This is real and sim are switched, but it doesn't matter
                            image2=obj[0],x1_2=x1_sim,y1_2=y1_sim,x2_2=x2_sim,y2_2=y2_sim,ext2=0,
                            close="no",fig_num=num,clean1="yes",clean2="yes")
                num=num+1
        
        
        #4. Assembly the simulated image for astrodrizzle, build the other extensions
        iraf.delete("sim_flt_*.fits")
        noised_blots_list=glob.glob("full_*_distorted_noise2.fits")
        i=1
        for blot in noised_blots_list:
            #hdu_sim=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+".fits","update")
            hdu_1_new=fits.open(blot)
            #hdu_sim[1].data=hdu_1_new[0].data
            #hdu_sim.close()
            
            hdu_old_name=list_of_original_flts[i-1]
            hdu_old=fits.open(hdu_old_name)
            hdu_1_new[0].header["EXTVER"]=1
            hdu_old[3].data=np.zeros(hdu_old[3].data.shape,dtype="int")
            
            #Zero extension which is just header information is taken from the original flt
            #First extension is the SCI image which has the simulated data and its header (might use old header)
            #Third extension is the data quality flag image which since the simulated data doesn't have issues is all zeros as detailed above
            #If preferred the user could change hdu_old[1].header to hdu_1_new[0].header to use the new header, not sure it makes a difference
            hdulist = fits.HDUList([hdu_old[0],fits.ImageHDU(hdu_1_new[0].data,hdu_old[1].header,name="SCI"),hdu_old[3]])
            hdulist.writeto(base_folder+"sim_flt_"+str(i)+".fits")
            i=i+1
            
            
        #5. Run astrodrizzle, user can change if necessary
        astrodrizzle.AstroDrizzle(input=glob.glob(base_folder+"sim_flt_?.fits"),output="full",skysub=False,#wcskey="HILO", <--If using Tweakreg
                                  final_pixfrac=final_pixfrac,final_bits=final_bits,driz_cr=False,median=False,blot=False,driz_separate=False,
                                  context=True,resetbits=resetbits,clean=True,final_rot=final_rot,final_scale=final_scale,final_fillval=final_fillval)

        
        
        #5. Compare to actual drz if desired
        if compare_final_drz=="yes":
            os.system("ds9 full_drz_sci.fits -fits "+path_to_actual_driz+" &")
            shape=fits.open(path_to_actual_driz)[0].data.shape
            x1_real,y1_real,x2_real,y2_real=1,1,shape[0],shape[1]
            x1_sim,y1_sim,x2_sim,y2_sim=1,1,shape[0],shape[1]
            sky_compare("full_drz_sci.fits",x1_sim,y1_sim,x2_sim,y2_sim,ext1=0,
                        image2=path_to_actual_driz,x1_2=x1_real,y1_2=y1_real,x2_2=x2_real,y2_2=y2_real,ext2=0,
                        close="no",fig_num=num,clean1="yes",clean2="yes")
            num=num+1
        final_large_image="full_drz_sci.fits"
        
        
        #Copy to folder to be tested, this was just for my personal use.
        #The other things listed are artifacts of things I tried to get the drizzling to work
        """
        outputs=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_drz_*.fits")
        test_folder_outputs=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/testing_my_simulated_drizzled_image/do-everything-with-blot/do-everything-with-blot_*.fits")
        iraf.delete("/Users/rpeterson/FP/GALFIT_Simulations1/testing_my_simulated_drizzled_image/do-everything-with-blot/*.fits")
        for obj in zip(outputs,test_folder_outputs):
            os.system("cp "+obj[0]+" "+obj[1])

        sys.exit()
        
        sys.exit()
        
        
        x1_real:524

y1_real:689

x2_real:602

x2_real:738

x1_sim:598

y1_sim:570

x2_sim:653

y2_sim:605


        #4. 
        
        #Offsets, if necessary?
        
        x1,y1=418.119, 667.413
        x2,y2=421.533, 670.826
        x3,y3=446.567, 695.572
        x4,y4=450.014, 699.035
        offsets=dict()
        offsets[1],offsets[2],offsets[3],offsets[4]={"x":x1,"y":y1},{"x":x2,"y":y2},{"x":x3,"y":y3},{"x":x4,"y":y4}
        for i in range(1,5):
            file="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+".fits"
            offset_file_temp=file.replace(".fits","_offset_temp.fits")
            offset_file=file.replace(".fits","_offset.fits")
            iraf.delete(offset_file_temp)
            iraf.imcopy("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies.fits",offset_file_temp)
            iraf.imshift(offset_file_temp, offset_file, offsets[i]["x"]-offsets[1]["x"],offsets[i]["y"]-offsets[1]["y"])
            iraf.delete(offset_file_temp)
        
        #Do blotting
        os.system("cp "+image_to_sim+" "+base_folder)
        flts_list=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01???_flt.fits")
        [iraf.imcopy(obj+"[1]",obj.replace(".fits","_simple.fits")) for obj in flts_list if os.path.exists(obj.replace(".fits","_simple.fits"))==False]
        flts_list_simple=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01???_flt_simple.fits")
        full_with_gals="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies.fits"
        [iraf.imcopy(full_with_gals+"[1]",full_with_gals.replace(".fits","_"+str(i)+".fits")) for i in range(1,5)]
        full_with_gals_list=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_?.fits")
        noised_list=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_noise2_?.fits")
        for obj in zip(flts_list_simple,full_with_gals_list):
            run(obj[1], obj[0], output=obj[1].replace(".fits","_distorted.fits"), interp='poly5', sinscl=1.0, stepsize=10, clobber=True)
        #backgrounds,ind=[0.63739,0.61008,0.64024,0.60927],0
        #readnoises=[20.200,19.7999,19.900,20.100]
        #gains=[652.9375,652.9375,652.9375,652.9375]
        #gains=[g*2.5 for g in gains]
        #sys.exit()
        
        for image in glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_?_distorted.fits"):
            background2=backgrounds[ind]-0.01
            rdnoise2=readnoises[ind]
            gain2=gains[ind]
            if image_to_sim_ext!=0:
                iraf.mknoise(image,output=image.replace(".fits","_noise2.fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
            else:
                iraf.mknoise(image,output=image.replace(".fits","_noise2.fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
            seed2=universal_seed*seed2
            ind=ind+1
        
        seed2=universal_seed
        
        
        sky_compare("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_4_distorted_noise2.fits",708,102,823,174,
                    ext1=0,image2="/Users/rpeterson/Downloads/id6p01y4q_flt.fits",x1_2=124,y1_2=609,
                    x2_2=178,y2_2=644,ext2=1,close="no")
        sky_compare("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_3_distorted_noise2.fits",708,102,823,174,
                    ext1=0,image2="/Users/rpeterson/Downloads/id6p01y2q_flt.fits",x1_2=893,y1_2=283,
                    x2_2=957,y2_2=424,ext2=1,close="no",fig_num=101)
        sky_compare("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_2_distorted_noise2.fits",693,233,786,293,
                    ext1=0,image2="/Users/rpeterson/Downloads/id6p01y0q_flt.fits",x1_2=306,y1_2=232,
                    x2_2=386,y2_2=283,ext2=1,close="no",fig_num=103)
        sky_compare("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_1_distorted_noise2.fits",36,550,122,606,
                    ext1=0,image2="/Users/rpeterson/Downloads/id6p01xyq_flt.fits",x1_2=451,y1_2=686,
                    x2_2=548,y2_2=747,ext2=1,close="no",fig_num=104)
        #sys.exit()
        #Copy other extensions from flats, should glob this
        
        iraf.delete("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim*.fits")
        
        overwrite_ext_1="yes"
        
        for i in range(0,1):
            try:
                #I have no idea why this is necessary, but if this isn't here we get a segmentation fault from iraf??
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01xyq_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_1.fits["+str(i)+",overwrite]")
            except:
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01xyq_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_1.fits["+str(i)+",overwrite]")
            iraf.imcopy("/Users/rpeterson/Downloads/id6p01y0q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_2.fits["+str(i)+",overwrite]")
            iraf.imcopy("/Users/rpeterson/Downloads/id6p01y2q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_3.fits["+str(i)+",overwrite]")
            iraf.imcopy("/Users/rpeterson/Downloads/id6p01y4q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_4.fits["+str(i)+",overwrite]")
        #iraf.
        for i in range(1,4): #Make a DQ image of all zero, only need SCI,DQ and WCSCORR
            if i in [2]:
                continue
            if False:
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01y0q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_1.fits[append]")
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01y0q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_2.fits[append]")
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01y0q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_3.fits[append]")
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01y0q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_4.fits[append]")
            else:
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01xyq_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_1.fits[append]")
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01y0q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_2.fits[append]")
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01y2q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_3.fits[append]")
                iraf.imcopy("/Users/rpeterson/Downloads/id6p01y4q_flt.fits["+str(i)+"]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_4.fits[append]")
        
        iraf.copy("/Users/rpeterson/Downloads/id6p01xyq_flt.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_1.fits")
        iraf.copy("/Users/rpeterson/Downloads/id6p01y0q_flt.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_2.fits")
        iraf.copy("/Users/rpeterson/Downloads/id6p01y2q_flt.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_3.fits")
        iraf.copy("/Users/rpeterson/Downloads/id6p01y4q_flt.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_4.fits")
        
        if overwrite_ext_1=="yes":
            for i in range(1,5):
                #iraf.imarith("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted.fits","*",4,"/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted_x4.fits")

                
                
                image="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted.fits"
                iraf.mknoise(image,output=image.replace(".fits","_noise2.fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
                seed2=universal_seed*seed2
                #iraf.imcopy("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted_noise2.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+".fits[DQ,3,overwrite]")
                #hdu_sim=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+".fits","update")
                #hdu_new=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted_noise2.fits")
                #hdu_sim[1].data=hdu_new[0].data
                #hdu_sim.close()
                #sys.exit()
                iraf.imcopy("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted.fits[0]","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+".fits[SCI,1,overwrite]")
                
        
        #sys.exit()
        #Drizzle
        hdu_sim[1].data=hdu_new[0].data
hdu_sim=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_1.fits","update")

hdu_sim[1].data=hdu_new[0].data

hdu_sim
Out[522]: [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x1fd1c52d0>, <astropy.io.fits.hdu.image.ImageHDU object at 0x1f9259210>, <astropy.io.fits.hdu.image.ImageHDU object at 0x1fc155290>, <astropy.io.fits.hdu.image.ImageHDU object at 0x20ef07390>, <astropy.io.fits.hdu.image.ImageHDU object at 0x1835dd790>, <astropy.io.fits.hdu.image.ImageHDU object at 0x1f964b650>, <astropy.io.fits.hdu.table.BinTableHDU object at 0x1f4abd210>]

hdu_sim.close()
        #iraf.delete("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01???_flt.fits")
        #os.rename("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_1.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01xyq_flt.fits")
        #os.rename("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_2.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01y0q_flt.fits")
        #os.rename("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_3.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01y2q_flt.fits")
        #os.rename("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_4.fits","/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01y4q_flt.fits")
        
        actual_driz=fits.open("/Users/rpeterson/FP/do_everything/UV105842F160W/UV105842F160W_drz.fits")[0].data
        
        x_dimension,y_dimension=actual_driz.shape[1],actual_driz.shape[0]
        
        flts_list=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/id6p01???_flt.fits")
        #sys.exit()
        #/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_noise2_4_distorted.fits
        for i in range(1,5):
            #hdu_sim=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+".fits","update")
            hdu_1_new=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted.fits")
            #hdu_sim[1].data=hdu_1_new[0].data
            #hdu_sim.close()
            
            hdu_old_name=flts_list[i-1]
            hdu_old=fits.open(hdu_old_name)
            #hdu_new[0].header["EXTVER"]=1
            hdu_old[3].data=np.zeros(hdu_old[3].data.shape,dtype="int")
            hdulist = fits.HDUList([hdu_old[0],fits.ImageHDU(hdu_1_new[0].data*1.0355,hdu_old[1].header,name="SCI"),hdu_old[3]])
            hdulist.writeto("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+"_new.fits")
            print(fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+"_new.fits").info())
            print(hdu_old.info())
            
        #sys.exit()
        
        #TweakREG
        #Run SExtractor 
        config_file="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/config.sex"
        pos_list_file=open("pos.list","w")
        astdriz_file=open("astdriz.list","w")
        for i in range(1,5):
            new_config_name="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/config"+str(i)+".sex"
            new_config_file=open(new_config_name,"w")
            [new_config_file.write(line.replace("astrodrizzle","sim"+str(i))) for line in open(config_file,"r").readlines()]
            new_config_file.close()
            os.system("sex sim_flt_"+str(i)+"_new.fits -c config"+str(i)+".sex")
            pos_list_file.write("sim_flt_"+str(i)+"_new.fits\n")
            astdriz_file.write("sim_flt_"+str(i)+"_new.fits "+"sim"+str(i)+"_tabext1.cat\n")
        pos_list_file.close()
        astdriz_file.close()
        tweakreg.TweakReg("@pos.list",catfile="astdriz.list",xcol=3,ycol=4,updatehdr="yes",
                          wcsname="HILO",searchrad=2,fluxunits="mag",rfluxunits="mag",threshold=3,clean=True)
        #sys.exit()
        for i in range(1,5):
            hdu_sim=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+"_new.fits","update")
            image_for_mknoise="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted.fits"
            if image_to_sim_ext!=0:
                iraf.mknoise(image_for_mknoise,output=image_for_mknoise.replace(".fits","_noise2.fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
            else:
                iraf.mknoise(image_for_mknoise,output=image_for_mknoise.replace(".fits","_noise2.fits"), background=background2, gain=gain2, rdnoise=rdnoise2, poisson=poisson2, seed=seed2)
            seed2=universal_seed*seed2
            hdu_sim[1].data=fits.open(image_for_mknoise.replace(".fits","_noise2.fits"))[0].data
            hdu_sim.close()
            
            #ind=ind+1
            hdu_new=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_field_with_sky_galaxies_"+str(i)+"_distorted_noise2.fits")
            hdu_new_sim=fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+"_new.fits")
            hdu_old_name=flts_list[i-1]
            hdu_old=fits.open(hdu_old_name)
            hdu_new[0].header["EXTVER"]=1
            hdu_old[3].data=np.zeros(hdu_old[3].data.shape,dtype="int")
            hdu_new_sim[3].header["WCS_key"]="HILO"
            hdu_new_sim[3].header["WCS_KEY"]="HILO"
            for item in hdu_new_sim[3].header.keys():
                if len(item)>0:
                    print(item)
                    hdu_old[0].header[item]=hdu_new_sim[3].header[item]
                    hdu_old[3].header[item]=hdu_new_sim[3].header[item]
                    hdu_new[0].header[item]=hdu_new_sim[3].header[item]
            hdulist = fits.HDUList([hdu_old[0],fits.ImageHDU(hdu_new[0].data,hdu_new[0].header,name="SCI"),hdu_old[3],fits.BinTableHDU(hdu_new_sim[3].data,hdu_new_sim[3].header,"WCSCORR")])
            hdulist.writeto("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+"_new_new.fits")
            print(fits.open("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_"+str(i)+"_new_new.fits").info())
            print(hdu_old.info())
            
        seed2=universal_seed   
        #JUST ADD NOISE AFTER
        
        #sys.exit()
        FITS_rec([ ('OPUS', 1, 'O', '', 'id6p01y4q_w3m18525i', '', '',  150.25825241,  2.02019512,  557.,  557.,  -3.23250000e-05,  -1.74431000e-05,  -1.93362000e-05,   2.88303000e-05, 'RA---TAN', 'DEC--TAN',  0.,  0.,  0.,  0., 0, '', ''),
   ('IDC_w3m18525i', 1, '', '', 'id6p01y4q_w3m18525i', 'NOMODEL', 'NOMODEL',  150.26074236,  2.01972041,  507.,  507.,  -3.22642068e-05,  -1.73714625e-05,  -1.93409820e-05,   2.87600398e-05, 'RA---TAN', 'DEC--TAN',  0.,  0.,  0.,  0., 0, '', ''),
   ('', 0, '', '', '', '', '',    0.        ,  0.        ,    0.,    0.,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, '', '',  0.,  0.,  0.,  0., 0, '', ''),
   ('', 0, '', '', '', '', '',    0.        ,  0.        ,    0.,    0.,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, '', '',  0.,  0.,  0.,  0., 0, '', ''),
   ('', 0, '', '', '', '', '',    0.        ,  0.        ,    0.,    0.,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, '', '',  0.,  0.,  0.,  0., 0, '', ''),
   ('', 0, '', '', '', '', '',    0.        ,  0.        ,    0.,    0.,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, '', '',  0.,  0.,  0.,  0., 0, '', ''),
   ('', 0, '', '', '', '', '',    0.        ,  0.        ,    0.,    0.,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   0.00000000e+00, '', '',  0.,  0.,  0.,  0., 0, '', '')], 
  dtype=(numpy.record, [('WCS_ID', 'S40'), ('EXTVER', '>i2'), ('WCS_key', 'S1'), ('HDRNAME', 'S24'), ('SIPNAME', 'S24'), ('NPOLNAME', 'S24'), ('D2IMNAME', 'S24'), ('CRVAL1', '>f8'), ('CRVAL2', '>f8'), ('CRPIX1', '>f8'), ('CRPIX2', '>f8'), ('CD1_1', '>f8'), ('CD1_2', '>f8'), ('CD2_1', '>f8'), ('CD2_2', '>f8'), ('CTYPE1', 'S24'), ('CTYPE2', 'S24'), ('ORIENTAT', '>f8'), ('PA_V3', '>f8'), ('RMS_RA', '>f8'), ('RMS_Dec', '>f8'), ('NMatch', '>i4'), ('Catalog', 'S40'), ('DESCRIP', 'S128')]))
        astrodrizzle.AstroDrizzle(input=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/sim_flt_?_new.fits"),output="full",skysub=False,
                      wcskey="HILO",final_pixfrac=0.8,final_bits=576,driz_cr=False,median=False,blot=False,driz_separate=False,
                      context=True,resetbits=0,clean=True,final_rot=0,final_scale=0.06,final_fillval=0)
        #sys.exit()
        final_large_image="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_drz_sci.fits"
        #mknoise2_file=fits.open(final_large_image,mode="update")
        #mknoise2_file[0].data=np.nan_to_num(mknoise2_file[0].data)
        #mknoise2_file[0].data[actual_driz==0]=0
        #mknoise2_file.close()
        outputs=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_astrodrizzled/full_drz_*.fits")
        test_folder_outputs=glob.glob("/Users/rpeterson/FP/GALFIT_Simulations1/testing_my_simulated_drizzled_image/do-everything-with-blot/do-everything-with-blot_*.fits")
        iraf.delete("/Users/rpeterson/FP/GALFIT_Simulations1/testing_my_simulated_drizzled_image/do-everything-with-blot/*.fits")
        for obj in zip(outputs,test_folder_outputs):
            os.system("cp "+obj[0]+" "+obj[1])
        """
        
        
        
        
        
        
        
        

    ###Displays large field
    if display_large_field=="yes":
        if label=="yes" and os.path.exists("region_file.reg")==True:
            print("\nDisplaying labelled "+final_large_image)
            os.system("ds9 "+final_large_image+" -regions load region_file.reg &")
        if label=="no":
            print("\nDisplaying "+final_large_image)
            os.system("ds9 "+final_large_image+" &")
    ###Saves large field
    if save=="yes":
        save_image = final_large_image.replace(".fits","_saved.fits")
        if os.path.exists(save_image):
            choice = raw_input("Saved field already exists. Overwrite?(y/n):")
            if choice=="y":
                iraf.delete(save_image)
                iraf.imcopy(final_large_image,save_image)
            if choice=="n":
                pass
        else:
            iraf.imcopy(final_large_image,save_image)
    print("Based on your selected inputs your final field is: "+final_large_image)
    if make_pairs=="yes":
        print(str(parents)+" parent galaxies could be added")
        print("These childs were skipped "+str(skipped_ids))
    ###Now compares the entire image
    if compare_sky=="yes":
        #Before I thought it might be a good idea to plot specific regions, but I think it overcomplicates matters
        """
        print("Reminder: if you change seed or number of sims, you should pick a new sky region")
        if type(universal_seed)!=int:
            print("Should have universal_seed fixed")
        sky_coords=[real_sky_lowerx,real_sky_lowery,real_sky_upperx,real_sky_uppery,sim_sky_lowerx,sim_sky_lowery,sim_sky_upperx,sim_sky_uppery]
        sky_coords_dict=dict()
        sky_coords_names = "real_sky_lowerx,real_sky_lowery,real_sky_upperx,real_sky_uppery,sim_sky_lowerx,sim_sky_lowery,sim_sky_upperx,sim_sky_uppery".split(",")
        if any(isinstance(x, type(None)) for x in sky_coords):
            
            os.system("ds9 "+image_to_sim+" -fits "+final_large_image+" -lock frame image &")
            print("I recommend zooming to some region of sky which is shared by both, then you can just put same after picking the first")
            j=1
            for i in sky_coords_names:
                val=raw_input("Pick a coordinate for "+str(i)+":")
                sky_coords_dict[str(i)]=int(val)
                if j==4:
                    ask=raw_input("Make sim_sky_coords the same?(y/n):")
                    if ask=="y":
                        sim_sky_lowerx,sim_sky_lowery,sim_sky_upperx,sim_sky_uppery=sky_coords_dict["real_sky_lowerx"],sky_coords_dict["real_sky_lowery"],sky_coords_dict["real_sky_upperx"],sky_coords_dict["real_sky_uppery"]
                        break
                    else:
                        print("not making same")
                j=j+1
        else:
            sky_coords_dict["real_sky_lowerx"],sky_coords_dict["real_sky_lowery"],sky_coords_dict["real_sky_upperx"],sky_coords_dict["real_sky_uppery"]=real_sky_lowerx,real_sky_lowery,real_sky_upperx,real_sky_uppery
            sky_coords_dict["sim_sky_lowerx"],sky_coords_dict["sim_sky_lowery"],sky_coords_dict["sim_sky_upperx"],sky_coords_dict["sim_sky_uppery"]=sim_sky_lowerx,sim_sky_lowery,sim_sky_upperx,sim_sky_uppery
        """
        shape_sim=fits.open(final_large_image)[0].data.shape
        if image_to_sim!="default":
            sky_compare(final_large_image,1,1,shape_sim[0],shape_sim[1],ext1=0,
                        image2=image_to_sim,x1_2=1,y1_2=1,x2_2=shape_sim[0],y2_2=shape_sim[1],
                        ext2=image_to_sim_ext,close="no",fig_num=num,clean2="yes",clean1="yes")
        else:
            sky_compare(final_large_image,1,1,shape_sim[0],shape_sim[1],ext1=0,
                close="no",fig_num=num,clean2="no",clean1="yes")
            
            


                
try:                ### I don't need this but other versions of python ide's might
    plt.show()
except:
    pass
