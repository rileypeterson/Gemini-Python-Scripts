### File is created on September 21st by Riley Peterson
### Attempts combine the tasks in the GCP_Photometry Guide into one script

### TASKS
make_EXPTIME_1 = 'yes' #if "yes", changes EXPTIME keyword in the header of the drz, wht, and ctx image. Saves the original EXPTIME in important_parameters.txt
make_ncomb='no' #if "yes", makes the ncomb image. The pixel values in this image represent the number of stacked images which contributed to said pixel.
make_weight='no' #if "yes", creates the weight image. If a pixel in the ncomb image has a value >= nmin, it keeps the corresponding value in the original wht image. Otherwise, the value is changed to zero.
make_sigma='no' #if "yes", creates the sigma image. Formula for this is in the manual.
run_sex1='no'  #if "yes", runs SExtractor. Outputs the table of information whose columns are the parameters listed in the param file. Also outputs the aperture, segmentation, and background image.
auto_review_sex_output='no' #if "yes", opens SExtractor first run output in an excel file for you to review. I normally never set this to "yes". No subsequent tasks depend on this being "yes".
make_assoc='no' #if "yes", creates an association file which distinguishes target galaxies. Galaxies are selected based on MAG_AUTO, MAG_AUTO_ERR, and CLASS_STAR from SExtractor.
make_mask='no' #if "yes", creates the mask image according to the limiting magnitude.
make_mask_name="default" #On the first run this should probably be "default". If not "default" the string supplied here will be appended to the mask name. 
                         #e.g. make_mask_name="final" the mask will be base_name_mask_final.fits
final_mask_for_galfit="default" #full path to the mask that is to be used for GALFIT, if default uses base_name_mask.fits.

make_sex_sky="no" #if "yes", checks it itersky.cat exists, if it does (it won't exist until make_sky_plots="yes") then the changes in sky values are made according to the lists below 
if make_sex_sky=="yes":    
    make_sky_plots="no"   #Should probably be "yes" if first time running it, then can set to "no". This will make plots comparing SExtractor and iterstat sky values.
make_psf='no' #if "yes", makes the point spread function for the field according to the parameters below. Should have run SExtractor.
make_constraint="no" #if "yes", creates an absolute constraint file, this feature is not very robust. I recommend not using a constraint file (i.e. leave this as "no")

do_run_galfit="yes" #if "yes", the GALFIT tasks below can be executed
if do_run_galfit=="yes":
    changes_type="local" #1. As of December 18th, these three lines should not be editted.
    use_output_as_input="no" #2.
    method="stop_and_go" #3.
    
    make_cutouts="no" #if "yes", creates cutouts
    cutouts_list_to_make="default" #if "default" makes cutouts for target galaxies (non-zero ids). Else, this is a list of galaxies to be made. Can be useful for remaking the mask for specific galaxies. E.g. [230,233].
    
    ###Sersic1
    make_feedmes_sersic1="no" #if "yes", creates feedme files for the sersic1 fit, make sure PSF folder is in 2D folder
    n_start1=1.5 #start value for the sersic index
    
    run_sersic1="no" #if "yes", runs the sersic1 fits, once the sersic output are made this should probably be set to "no"
    run_sersic1_list="default" #if "default", the entire batch of sersic1 fits will be run. I think most of the time this should be "default". DO NOT use this list as a refit list.
    #the main use I see for providing a list here is that the user wants to check and make sure a few of the fits look good before running the entire batch, though they could also just ctrl+c to stop the script if it is "default".
    number_sent=10 ###The number of feedmes sent to the shell at a given time, do not make ridiculously high, probably should be number of cpu cores. I used 10.
    display_blocks_sersic1="no" ###if "yes", it displays the original,model, and residual as outputs are made, useful to look at the fit.log at the same time
    timer=20 ###Delay between showing output blocks, does not matter if display_blocks="no"

    ###Sersic1 Refits
    check_residuals_sersic1="no" #if "yes", cycles through the outputs that have been made so the user can check the residuals
    check_residuals_sersic1_list="default" #if "default", cycles through all the made outputs, otherwise just the list of ids provided here.
    plot_stats_sersic1="no" #if "yes", plots output statistics and suggests a refit list
    
    run_refits_sersic1="no" #if "yes", copies the obj.in and output files of cutouts in run_refits_sersic1_list to the refit folder so that the user can begin refitting these galaxies
    run_refits_sersic1_list=[2, 5, 30, 47, 62, 79, 88, 91, 113, 117, 145, 154, 155, 160, 162, 166, 175, 182, 186, 229, 235, 252, 263, 265, 271, 287, 289, 298, 367, 371, 372, 386, 399, 406, 412, 426, 440, 456, 474, 475, 481, 515, 524, 12, 130, 148, 183, 210, 3, 392, 410, 42, 487, 6, 60, 94] #Should NOT be "default", this is the list of galaxies that need to be refit
    
    enact_refit_edits_sersic1="yes" #if "yes", makes the changes to the longfeedme reflected in the summary_of_changes file (syntax is important for the summary file, see manual)
    combine_refits_and_original_to_final_folder_sersic1_and_make_table="yes" #if "yes", first copies the sersic1 original fits into a newly made final folder, then overwrites the refit galaxies by copying their output from the refit folder, creates txt and fits output results table


    ###Devauc
    make_feedmes_devauc1="no" #if "yes", makes the devauc obj.in files according to the sersic1 longfeedme file(edited for refits) 
    
    run_devauc1="no" #if "yes", runs all the devauc fits. number_sent and timer are defined above
    display_blocks_devauc1="no" #if "yes", displays outputs from the devauc fits as they are made
    

    ###Devauc1 Refits
    check_residuals_devauc1="no" #if "yes", checks the devauc outputs
    check_residuals_devauc1_list="default" #if "default", cycles through all the output, else the input should be a list (it will only cycle through these)
    plot_stats_devauc1="no" #if "yes", plots statistics from the devauc output
    
    run_refits_devauc1="no" #if "yes", the obj.in and output files (if available) from the list below are copied into the refit folder
    run_refits_devauc1_list=[] #Should NOT be "default", should be the list of galaxy ids that need to be refit
    
    enact_refit_edits_devauc1="no" #if "yes", enacts the changes in the summary file onto the longfeedme
    combine_refits_and_original_to_final_folder_devauc1_and_make_table="no" #if "yes", will compile the original and refits (overwrite the originals). Creates txt and fits outputs








    #Sersic2
    make_feedmes_sersic2="no" #if "yes", makes feedmes for galaxies which are eligible for sersic2 fits (see manual for eligibility criteria)
    
    run_sersic2="no" #if "yes", runs these galaxies using a starting n value of 4
    display_blocks_sersic2="no" #if "yes", will display outputs as they are made
    

    ###Sersic2 Refits
    check_residuals_sersic2="no" #same as above but for sersic2
    check_residuals_sersic2_list="default" #same as above but for sersic2
    plot_stats_sersic2="no" #same as above but for sersic2
    
    run_refits_sersic2="no" #same as above but for sersic2
    run_refits_sersic2_list=[144,209,397,406,423,436,537,721] #same as above but for sersic2
    
    enact_refit_edits_sersic2="no" #same as above but for sersic2
    use_sersic2_list=[144,209,397,406,423,721] #These two lists do not appear above, this is the list of galaxies for which the sersic2 fit performed better than the original sersic1 fit. If not in this list the galaxy will use the sersic1 output.
    remove_from_final_catalog_sersic2_list=[226,360,398] #These are typically problem galaxies for which neither the sersic1 or 2 fit were any good. These won't show up in the final catalog.
    combine_refits_and_original_to_final_folder_sersic2_and_make_table="no" #if "yes", copies from the sersic1 final folder then overwrites for the galaxies in use_sersic2_list. Creates txt and fits output tables
    



###INPUT PARAMETERS FOR TASKS

login_cl_dir = "/home/riley/" #Directory with your login.cl file. This is so when import pyraf from iraf is called, your login.cl will be recognized
base_folder = "/home/riley/Gemini/Gemini-Python-Scripts/Examples/example_1/do_ev_ex1/"  #Should not end with "/", CANNOT CONTAIN "_" AFTER THE LAST SLASH, I.E. BASE_NAME CANNOT NOT HAVE AN UNDERSCORE (IRONY). Now these syntactical nuances are corrected for by the script, so its not a huge deal.
drz_image = "do-everything_drz.fits"  #Should be like this UV108899id6p03F160W_drz.fits, I might implement this automatically in the future, but then I would need to guess which image is which (probably by looking for wht or drz or ctx in the name)
wht_image = "do-everything_wht.fits"  #Should be like this UV108899id6p03F160W_wht.fits
ctx_image = "do-everything_ctx.fits"  #Should be like this UV108899id6p03F160W_ctx.fits

### Mkncomb, mkweight, sigma image inputs
n_min = 3 #For mkweight minimum number of contributing stacked images for a valid pixel
read_noise = 20 #For sigma image creation, read noise
sqrt_FA = 0.4708 #For sigma image creation, variance reduction factor see Casertano et al. 2000
norm_factor = 4 #For sigma image creation, normalization factor. 
                #If max(wht_image)~original EXPTIME of wht image (see important_parameters.txt) then this should be 1
                #Else this should be ~(original EXPTIME of wht image)/n_images 
                #Basically the wht image can be used as a normalization for the sigma image if its values(in the center) are
                # roughly equal to the EXPTIME of the stack
                #For example my EXPTIME was ~2611 but my wht image max values were ~ 750. 750 * 4 ~ 2611
                #Basically creates an exposure time map (i.e. total EXPTIME for each individual pixels)
 

### SExtractor inputs
### IF PARAMETER IS NOT LISTED HERE IT DOESN'T GET INCLUDED IN YOUR .SEX FILE
### see your default.sex file in your SExtractor folder for reference, or the official documentation
### mine is here: /Users/rpeterson/anaconda/pkgs/sextractor-2.19.5-0/share/sextractor/default.sex
### if you have something you want to go off of, e.g. the default.sex file or another .sex file then set print_vars="yes"
### then make path_to_sex_config="/your/path/to/any/sex/file.sex"
print_vars="no" #When this is "yes", subsequent tasks do not run because the script exits (allowing you to paste in the SExtractor variables)
path_to_sex_config = "/Users/rpeterson/FP/do_everything/UV105842F160W/UV105842F160W_phot/UV105842F160W.sex"
### now you can paste this output below, and change things as needed, should all be strings, you are taking a gamble if you
### change a string to float or int by removing the quotes, so just leave them as strings
### note that CHECKIMAGE_TYPE,CHECKIMAGE_NAME,CATALOG_NAME,WEIGHT_IMAGE,WEIGHT_TYPE
### will be changed farther down, so that everything else runs smoothly, so it won't matter if you define these
### once you have it figured out set print_vars to "no"

### SExtractor variables
### It is best if _NAME variables are listed with absolute path
### e.g. FILTER_NAME = "/Users/rpeterson/anaconda/envs/astroconda/share/sextractor/default.conv"

ANALYSIS_THRESH = "2"
BACKPHOTO_THICK = "24"
BACKPHOTO_TYPE = "GLOBAL"
BACK_FILTERSIZE = "3"
BACK_SIZE = "128"
CATALOG_TYPE = "FITS_1.0"
CLEAN = "Y"
CLEAN_PARAM = "1.0"
DEBLEND_MINCONT = "0.01"
DEBLEND_NTHRESH = "16"
DETECT_MINAREA = "9"
DETECT_THRESH = "2"
DETECT_TYPE = "CCD"
FILTER = "Y"
FILTER_NAME = "/Users/rpeterson/anaconda/envs/astroconda/share/sextractor/default.conv"
FLAG_IMAGE = "flag.fits"
GAIN = "6529.37744"
GEOM_LEVELS = "24.0"
GEOM_TYPE = "SB"
MAG_GAMMA = "4.0"
MAG_ZEROPOINT = "25.9463"
MASK_TYPE = "CORRECT"
MEMORY_BUFSIZE = "8192"
MEMORY_OBJSTACK = "6000"
MEMORY_PIXSTACK = "600000"
MOFFAT_BETA = "1.0"
N_GEOM_LEVELS = "1.0"
N_MOFFAT_BETA = "1.0"
PARAMETERS_NAME = "/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_python_scripts/most_current/examples/example_1/MS.param"
PETROSIAN_ETA = "0.277217"
PHOT_APERTURES = "16.67"
PHOT_AUTOPARAMS = "2.5,3.5"
PHOT_FLUXFRAC = "0.5"
PIXEL_SCALE = "0.06"
SATUR_LEVEL = "1500.0"
SEEING_FWHM = "0.2"
STARNNW_NAME = "/Users/rpeterson/anaconda/envs/astroconda/share/sextractor/default.nnw"
VERBOSE_TYPE = "NORMAL"


### PSF inputs
psf_verify="no" ### daophot verify
psf_verbose="no" ### daophot verify
psf_apertures=12  ### photpars Aperture Size
psf_zmag=25.9463  ### photpars Should be same as above
psf_scale=1. ### datapars scale
psf_fwhmpsf=4 ### datapars fwhmpsf measure off stars [pixels]
psf_epadu = 2.5 #datapars epadu (electrons per adu?)
psf_readnoise = 2 #datapars readnoise
psf_itime=1  ### datapars exposure time typically 1
psf_sigma=0.004563  ### datapars sigma measure off of image background using imstat for a background region (use mean), will take fine tuning together with threshold    
psf_datamax="INDEF"  ### datapars datamax, if using saturated stars to build psf set to INDEF
psf_datamin="INDEF"### datapars datamin
psf_saturated="yes"  ### daopars saturated, if yes uses wings of saturated stars while making PSF 
psf_calgorithm="centroid" ### centerpars calgorithm
psf_cbox=20. ### centerpars cbox
psf_salgorithm="centroid" ### fitskypars salgorithm
psf_annulus=40.  ### fitskypars annulus
psf_dannulus=75 ### fitskypars dannulus 
psf_smooth="yes" ### fitskypars smooth
psf_matchrad=10. ### daopars matchrad
psf_psfrad=83 ### daopars psfrad, final psf will have this radius
psf_fitrad=10 ### daopars fitrad
psf_threshold=15 ### findpars threshold, will look for objects > sigma*threshold so this is important parameter together with sigma
max_psf = 200 ### maximum number of stars to include when making the psf for pstselect
star_thresh = 0.25 ### CLASS_STAR value threshold, if a star's CLASS_STAR is greater than star_thresh it contributes to PSF
limiting_mag = 24 ### If CLASS_STAR>star_thresh and MAG<limiting_mag then star will build psf, weeds out faint "stars"
datamax_sat = 40 ### datapars datamax, level for which stars are saturated, anything greater will still be used if psf_saturated='yes'
n_clean = 13 ### daopars nclean
psf_x_one = 1500
psf_y_one = 1500 ### output name of psf will be psf-psf_x_one-psf_y_one.fits, should be x,y center most location of region representing the PSF (if multiple PSFs)
psf_comparator="/Users/rpeterson/FP/UV-105842_id6p01_imaging/UV_105842_id6p01_psf_use/rotated_margalef_PSF.fits" #if not "default", will plot how your psf compares to this one (absolute path)

### Making association file inputs
assoc_mag_cut = 24.2 ### MAG_AUTO cut for association file
assoc_noise_cut = 0.04 ### MAG_AUTO_ERR cut for association file remember MAG_ERR=0.04 corresponds to S/N~25 (Ryan Cole)
assoc_class_star_cut = 0.8 ### CLASS_STAR cut for association file

### Mask inputs
mask_maglim = 25.2 ### maglim = min(magbright+6.2,magfaint+0.5), very important, objects fainter than this will be masked
add_to_mask_list=[] ### ids in this list will be added to the mask
remove_from_mask_list=[] ### ids in this list will be removed from the mask

### If you want add something to mask which isn't detected (therefore doesn't have an id #), see section 1.7.1 in the manual
### Perhaps I will implement this as a task in the future

### SEx sky inputs
change_to_iterstat_list=[] # list of ids whose sky values should be changed to the iterstat value
take_average_list=[127] #list of ids whose sky value should be the average of the iterstat and sextractor value
remove_entirely_list=[] #list of ids that are to be COMPLETELY removed from the list. I would do this if an object was just completely a misdetection, otherwise I would add to mask list.
# if all empty then no sky values are changed in the default SExtractor sky will be used.


### Constraint inputs
n_constraint=None ### Format is [lower_n_value, upper_n_value], not very robust, could use some more development and testing
re_constraint=None


### GALFIT inputs
min_cutout_radius = 125  ###Cutouts will be at least twice this distance on a side. E.g. min_cutout_radius = 125, my cutouts are 250 by 250
mag_zp = 25.9463 ###Magnitude zeropoint
plt_scale = 0.06 ###arcsec/pixel
fit_sky = 'no' ###Traditionally set to "no", if "yes" GALFIT will fit the sky, if "no" the sky value remains fixed
psf_size = 168 #size of convolution box for PSF
galfit_binary = "/net/frigg/Users/inger/bin/galfit" ### path to galfit executable, if you can run galfit from the terminal by typing "$galfit some_obj_in_file.in" then you can probably get away with this just being "galfit"











############################################################ Cross with caution        










if base_folder[-1]=="/":
    print("Removing last slash from base_folder")
    base_folder=base_folder[:-1]
if "_" in base_folder.split("/")[-1]:
    print("Removing underscores from base_name, replacing them with dashes")
    new_base_name=base_folder.split("/")[-1].replace("_","-")
    base_folder="/".join(base_folder.split("/")[:-1])+"/"+new_base_name
    





### Imports
import os
os.chdir(login_cl_dir)
from pyraf import iraf
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from astropy.io import fits
import time
import numpy as np
import subprocess
#import glob
import sys


final_mask_list=[]
### Definitions
def num_comp(dictionary,id):
    ###Determines the number of components for a given id within a longfeedme dictionary
    ###Implemented after a bug was discovered which did not allow more than 9 components for a cutout
    list_of_keys=[i  for i in dictionary[id].keys() if "id" in i]
    list_of_ints=[int(i.replace("id","")) for i in list_of_keys]
    return max(list_of_ints)

def copy_longfeedme_dictionary(dic_to_be_copied,folder_to_be_copied_to):
    ###Makes a copy of a longfeedme dictionary to a new folder if it doesn't already exist
    ###Perhaps this should overwrite
    folder_to_be_copied_to=folder_to_be_copied_to[:-1]
    final_name_of_dic=folder_to_be_copied_to+"/"+folder_to_be_copied_to.split("/")[-1]+"_longfeedme.npy"
    name_of_dic_to_be_copied=dic_to_be_copied.split("/")[-1]
    if os.path.exists(final_name_of_dic)==False:
        print(dic_to_be_copied,folder_to_be_copied_to)
        iraf.copy(dic_to_be_copied, folder_to_be_copied_to)
        print(folder_to_be_copied_to+"/"+name_of_dic_to_be_copied,final_name_of_dic)
        os.rename(folder_to_be_copied_to+"/"+name_of_dic_to_be_copied,final_name_of_dic)
    else:
        print(final_name_of_dic+" already exists.")
    return final_name_of_dic


def cat2longfeedme(catfile,long_out_path,dict_out_path,fit_type="sersic",n_start=1.5,min_cutout_radius=min_cutout_radius,final_mask_list=final_mask_list):
    ### Converts the .cat file from SExtractor to a long feedme file, creates a dictionary of this same information
    file = open(long_out_path,"w")
    mega_dict=dict()
    for line in open(catfile,"r").readlines():
        if line.split()[0]!="0" and int(line.split()[0]) not in final_mask_list:
            
            #Get values of main component
            id,mag,sky,q = int(line.split()[0]),float(line.split()[3]),float(line.split()[6]),(1-float(line.split()[14]))
            re,pa = float(line.split()[8])/q**0.5,float(line.split()[13])-90
            x_center,y_center = float(line.split()[9]),float(line.split()[10])
            x_center_int,y_center_int = int(x_center),int(y_center)
            if re > min_cutout_radius:
                radius = int(re)
            else:
                radius = min_cutout_radius
            #Determine cutout size
            start_x,end_x = x_center_int - radius,x_center_int + radius - 1
            start_y,end_y = y_center_int - radius,y_center_int + radius - 1
            x,y = x_center - start_x + 1,y_center - start_y + 1
            #Create main dictionary entry
            mega_dict[id]={"id1":id,"row1":id,"type1":fit_type,"x_global1":x_center,"y_global1":y_center,
                     "x_start1":start_x,"x_end1":end_x,"y_start1":start_y,"y_end1":end_y,"x_cutout1":x,"y_cutout1":y,
                     "mag1":mag,"re1":re,"n1":n_start,"q1":q,"pa1":pa}
            attribs = ["id1","row1","type1","x_global1","y_global1","x_start1","x_end1","y_start1","y_end1","x_cutout1","y_cutout1","mag1","re1","n1","q1","pa1"]
            orig_attribs = attribs
            k=1
            comp_num=2
            ind=1
            #Establish the other galaxies/objects in the cutout
            for line1 in open(catfile,"r").readlines():
                row_comp,mag_comp,x_comp,y_comp = k,float(line1.split()[3]),float(line1.split()[9]),float(line1.split()[10])
                if start_x<x_comp<end_x and start_y<y_comp<end_y and row_comp not in final_mask_list and row_comp!=id:
                    x_coord_comp,y_coord_comp = x_comp-start_x+1,y_comp-start_y+1
                    q_comp=(1-float(line1.split()[14]))
                    re_comp,pa_comp = float(line1.split()[8])/q_comp**0.5,float(line1.split()[13])-90
                    attribs=[attrib.replace(str(ind),str(comp_num)) for attrib in attribs]
                    ind=comp_num
                    nums=[id,row_comp,"sersic",x_comp,y_comp,start_x,end_x,start_y,end_y,x_coord_comp,y_coord_comp,mag_comp,re_comp,n_start,q_comp,pa_comp]
                    for i in range(len(attribs)):
                        mega_dict[id].update({attribs[i]:nums[i]})
                    comp_num=comp_num+1
                k=k+1
            mega_dict[id].update({"sky":sky})
            
            file.write(str(mega_dict[id]["id1"])+' \n')
            num_of_comps = num_comp(mega_dict,id)
            ind=1
            attribs = orig_attribs
            #Write the attributes to the file
            for i in range(0,num_of_comps):
                [file.write("    "+attribs[j]+" "+str(mega_dict[id][attribs[j]])+"\n") for j in range(len(attribs))]
                file.write("\n")
                attribs=[attrib.replace(str(ind),str(i+2)) for attrib in attribs]
                ind=i+2
            file.write("    sky "+str(mega_dict[id]["sky"])+"\n\n")
    file.close()
    #Save the dictionary
    np.save(dict_out_path,mega_dict)
    return mega_dict


def longfeedme2feedmes(input_longfeedme_dict,cutout_folder_path,out_folder_for_feedmes,list_to_make="default",psf_size=psf_size,fit_sky=fit_sky):
    ### Converts a longfeedme file to individual feedme files
    ### Assumes galfit will be run from out_folder_for_feedmes and that this folder is in two_D, eventually this might change
    try:
        os.mkdir(out_folder_for_feedmes)
    except:
        pass
    try:
        l_dict = np.load(input_longfeedme_dict).item()
    except:
        l_dict = input_longfeedme_dict
    iraf.cd(out_folder_for_feedmes)
    os.chdir(out_folder_for_feedmes)
    cutout_folder_path=cutout_folder_path.replace(two_D,"../")
    out_folder_for_feedmes=out_folder_for_feedmes.replace(two_D,"../")
    if list_to_make=="default":
        list_to_make=sorted(l_dict.keys())
    for id in list_to_make:
        file=open(out_folder_for_feedmes+"obj"+str(id)+".in","w")
        file.write("#================================================================================"+'\n')
        file.write("## IMAGE PARAMETERS for "+str(id)+'\n')
        file.write("A) "+cutout_folder_path+base_name+"_"+str(id)+"_drz.fits    # Input data image (FITS file)"+'\n')
        file.write("B) "+out_folder_for_feedmes+base_name+"_"+str(id)+"_out.fits    # Output data image block"+'\n')
        file.write("C) "+cutout_folder_path+base_name+"_"+str(id)+"_sigma.fits    # Sigma image name (made from data if blank or 'none')"+'\n')
        psf_distance=1000000000000000
        #Find closest PSF
        for psf in iraf.dir("../PSF/psf*",Stdout=1)[0].split():
            psf_x,psf_y = float(psf.split("-")[1]),float(psf.split("-")[2].replace(".fits",""))
            x_center_int,y_center_int=int(float(l_dict[id]["x_global1"])),int(float(l_dict[id]["y_global1"]))
            distance = ((x_center_int-psf_x)**2+(y_center_int-psf_y)**2)**0.5
            if distance < psf_distance:
                psf_distance=distance
                final_x=str(int(psf_x))
                final_y=str(int(psf_y))
        psf_name = "psf-"+final_x+"-"+final_y+".fits"    #ADDRESS
        file.write("D) ../PSF/"+psf_name+"    # Input PSF image and (optional) diffusion kernel"+'\n')
        file.write("E) 1                   # PSF fine sampling factor relative to data"+'\n')
        file.write("F) "+cutout_folder_path+base_name+"_"+str(id)+"_mask.fits"+"      # Bad pixel mask (FITS image or ASCII coord list)"+'\n')
        if make_constraint!="no":
            file.write("G) constraint_file.txt      # File with parameter constraints (ASCII file)"+'\n')
        xmax=int(l_dict[id]["x_end1"])-int(l_dict[id]["x_start1"])+1
        ymax=int(l_dict[id]["y_end1"])-int(l_dict[id]["y_start1"])+1
        file.write("H) 1 "+str(xmax)+" 1 "+str(ymax)+"    # Image region to fit (xmin xmax ymin ymax)"+'\n')
        file.write("I) "+str(psf_size)+' '+str(psf_size)+"           # Size of the convolution box (x y)"+'\n')
        file.write("J) "+str(mag_zp)+"            # Magnitude photometric zeropoint"+'\n')
        file.write("K) "+str(plt_scale)+" "+str(plt_scale)+"         # Plate scale (dx dy)   [arcsec per pixel]"+'\n')
        file.write("O) regular                # Display type (regular, curses, both)"+'\n')
        file.write("P) 0                   # Create ouput only? (1=yes; 0=optimize)"+'\n')
        file.write("""# INITIAL FITTING PARAMETERS 
        # 
        #   For object type, allowed functions are: sersic, nuker, 
        #                       expdisk, devauc, moffat, gaussian.""")
        
        file.write('\n'+'\n'+'\n'+'\n'+"### Beginning of Fitting for #"+str(id)+'\n'+'\n'+'\n')
        num_of_comps = num_comp(l_dict,id)
        for i in range(1,num_of_comps+1):
            i=str(i)
            file.write("#Entry for row number "+str(l_dict[id]["row"+i])+'\n')
            file.write("#Component #"+i+'\n')
            file.write(" 0) "+str(l_dict[id]["type"+i])+"                 #  Object type"+'\n')
            file.write(" 1) "+str(l_dict[id]["x_cutout"+i])+" "+str(l_dict[id]["y_cutout"+i])+"  1 1  #  position x, y"+'\n')
            file.write(" 3) "+str(l_dict[id]["mag"+i])+"     1          #  total magnitude"+'\n')
            file.write(" 4) "+str(l_dict[id]["re"+i])+"      1          #      R_e [pixels]"+'\n')
            if "sersic" in str(l_dict[id]["type"+i]):
                file.write(" 5) "+str(l_dict[id]["n"+i])+"      1          #  exponent (de Vaucouleurs = 4)"+'\n')
            if "devauc" in str(l_dict[id]["type"+i]):
                file.write(" 5) 4       0          #  exponent (de Vaucouleurs = 4)"+'\n')
            file.write(" 6) 0.0000      0          #     -----"+'\n')
            file.write(" 7) 0.0000      0          #     -----"+'\n')
            file.write(" 8) 0.0000      0          #     -----"+'\n')
            file.write(" 9) "+str(l_dict[id]["q"+i])+"      1          #  axis ratio (b/a)"+'\n')
            file.write("10) "+str(l_dict[id]["pa"+i])+"     1          #  position angle (PA) [deg: Up=0, Left=90]"+'\n')
            file.write(" Z) 0                      #  Output option (0 = residual, 1 = Don't subtract)"+'\n')
            file.write('\n')    
        file.write('\n'+"#Sky"+'\n')
        file.write(" 0) sky                    #    Object type"+'\n')
        if fit_sky=='yes':
            file.write(" 1) "+str(l_dict[id]["sky"])+"     1       #  sky background"+'\n')
        if fit_sky=='no':
            file.write(" 1) "+str(l_dict[id]["sky"])+"     0       #  sky background"+'\n')
        file.write(" Z) 0                      #  Output option (0 = residual, 1 = Don't subtract)"+'\n')
    file.close()




def do_make_cutouts(cutout_folder_path,longfeedme_dict,drz_image,sigma_image,mask_image,list_to_make="default"):
    ### Makes sci, sigma, and mask cutouts for the ids in list_to_make
    #print(list_to_make)
    try:
        mega_dict = np.load(longfeedme_dict).item()
    except:
        mega_dict = longfeedme_dict
    try:
        os.mkdir(cutout_folder_path)
    except:
        pass
    k=0
    if list_to_make=="default":
        list_to_make=sorted(mega_dict.keys())
    for i in list_to_make:
        print("Making cutouts for #"+str(i)+" Percent done:"+str(k*100/len(list_to_make)))
        start_x,end_x=mega_dict[i]["x_start1"],mega_dict[i]["x_end1"]
        start_y,end_y=mega_dict[i]["y_start1"],mega_dict[i]["y_end1"]
        id=str(mega_dict[i]["id1"])
        if start_x<1:   ###If the end of the cutout falls off the frame, than we make it 1. I thought there was going to be an issue here 
            start_x=1   ###because the cutout size parameter("H) 1 250 1 250 # Image region to fit (xmin xmax ymin ymax)") should obviously change as well, but
        if start_y<1:   ###I just tried running a larger region i.e. "H) 1 300 1 300 # Image region to fit (xmin xmax ymin ymax)" on 250 by 250 cutouts and the results seem the same.
            start_x=1   ###So it seems this isn't an issue and GALFIT adjusts for this. If it is, you could make a dictionary of cutout size in the def for line2feedme etc. email me if you need help rileypeterson21@gmail.com
                        ###The point is that if the end if off the edge of the image, the dictionary entries need to be updated else the global coordinates at the end might be off (since they are calculated using xstart and ystart)
            ###This might look like:
            ###for j in range(1,num_comp(mega_dict,i)):
                ###mega_dict[i]["y_start"+str(j)]=1
                ###etc... for x
                ###I would write this right now but I don't have time to test it.
        x_dim_max = int(iraf.imhead(drz_image,Stdout=1)[0].split("[")[1].split(",")[0])
        y_dim_max = int(iraf.imhead(drz_image,Stdout=1)[0].split("[")[1].split(",")[1].replace("]",""))  ###get dimensions of image :/
        if end_x>x_dim_max:
            end_x=x_dim_max
        if end_y>y_dim_max:
            end_y=y_dim_max
        if os.path.exists(cutout_folder_path+base_name+"_"+id+"_drz.fits"):
            os.remove(cutout_folder_path+base_name+"_"+id+"_drz.fits")
        if os.path.exists(cutout_folder_path+base_name+"_"+id+"_sigma.fits"):
            os.remove(cutout_folder_path+base_name+"_"+id+"_sigma.fits")
        if os.path.exists(cutout_folder_path+base_name+"_"+id+"_mask.fits"):
            os.remove(cutout_folder_path+base_name+"_"+id+"_mask.fits")
        print(start_x,end_x,start_y,end_y)
        iraf.imcopy(drz_image+"["+str(start_x)+":"+str(end_x)+","+str(start_y)+":"+str(end_y)+"]",cutout_folder_path+base_name+"_"+id+"_drz.fits",verbose="no")    ###Use imcopy to make cutouts
        iraf.imcopy(sigma_image+"["+str(start_x)+":"+str(end_x)+","+str(start_y)+":"+str(end_y)+"]",cutout_folder_path+base_name+"_"+id+"_sigma.fits",verbose="no")
        iraf.imcopy(mask_image+"["+str(start_x)+":"+str(end_x)+","+str(start_y)+":"+str(end_y)+"]",cutout_folder_path+base_name+"_"+id+"_mask.fits",verbose="no")
        k=k+1


def longfeedme2longfeedme_dict(longfeedme_txt,longfeedme_dict_out):
    ### Converts a longfeedme text file into a longfeedme dictionary
    longfeedme_dict=dict()
    txt=open(longfeedme_txt,"r").readlines()
    for line in txt:
        if line.startswith("    ")==False and line!="\n":
            id=int(line.split()[0])
            longfeedme_dict[id]={}
            ind=txt.index(line)+1
            while True:
                if ind>=len(txt)-1:
                    break
                if txt[ind]=="\n":
                    ind=ind+1
                if txt[ind].startswith("    "):
                    longfeedme_dict[id].update({txt[ind].split()[0]:txt[ind].split()[1]})
                    ind=ind+1
                else:
                    break
            
    np.save(longfeedme_dict_out,longfeedme_dict)
    return longfeedme_dict

def run_galfit(path_to_feedmes,path_to_outs,cutout_folder_path,fit="sersic1",input_longfeedme_dict=None,list_to_make=[],number_sent=10,display_blocks="yes",timer=20,verbose="yes",command=galfit_binary):
    ###Runs GALFIT in parallel
    ###A very arduous effort was made to correctly assess which galaxy ids crashed while this was running
    ###Unfortunately, a consistent method to accomplish this was not found, though I would welcome such a solution
    ###Crashed galaxies are listed at very end and written to file
    if list_to_make=="default":
        try:
            l_dict = np.load(input_longfeedme_dict).item()
        except:
            l_dict = input_longfeedme_dict
    list_to_make=sorted(l_dict.keys())
    iraf.cd(path_to_feedmes)
    os.chdir(path_to_feedmes)
    if len(list_to_make)<number_sent:    ###In the case you only have a few in your cat file
        number_sent=len(list_to_make)  
    try:
        iraf.delete(path_to_outs+"outfile.txt")
    except:
        pass
    outfile=open(path_to_outs+"outfile.txt","a+")    ###Create an output file for the GALFIT run
    i=0
    j=0
    k=0
    ind=0
    crash_count=0
    reset_time=timer
    out_files=[out for out in os.listdir(path_to_outs) if "_out.fits" in out and int(out.split("_")[1]) in list_to_make]
    out_times=[os.path.getmtime(path_to_outs+out) for out in os.listdir(path_to_outs) if "_out.fits" in out and int(out.split("_")[1]) in list_to_make]
    outs=zip(out_files,out_times)
    orig_outs=outs
    sent_files=[]
    created_fits=[]
    initial_time=time.time()
    while True:
        time.sleep(0.05)    #If this isn't here the script get a little wacky, I think this is because the os.listdir calls are sort of overloading the system
                            #Note this doesn't really make the script too much slower, because what you are really waiting on is the GALFIT output to finish, which takes the most time
                            #What could really be happening if this isn't here is some floating point error in the os.path.getmtime, perhaps rounding to the nearest minute is a good solution
        outfile_copy_file=open(path_to_outs+"outfile.txt","r")
        outfile_copy=outfile_copy_file.readlines()
        if j<number_sent and i<len(list_to_make):
            name="obj"+str(list_to_make[i])+".in"
            print("Sending "+name+" to the shell")
            process = subprocess.Popen([command,name],stdout=outfile)
            sent_files.append(name)
            i=i+1
            j=j+1
        if i==len(list_to_make):
            print("All of them have been sent, just waiting for them to finish")
            i=i+1
            
        try:
            new_out_files=[out for out in os.listdir(path_to_outs) if "_out.fits" in out and int(out.split("_")[1]) in list_to_make]
            new_out_times=[os.path.getmtime(path_to_outs+out) for out in os.listdir(path_to_outs) if "_out.fits" in out and int(out.split("_")[1]) in list_to_make]
            new_outs=zip(new_out_files,new_out_times)
        except:
            pass
        if len(list(set(new_outs)-set(orig_outs)))>ind:
            ind=ind+1
            new_outputs=[list(set(new_outs)-set(orig_outs))[f][0] for f in range(len(list(set(new_outs)-set(orig_outs)))) if list(set(new_outs)-set(orig_outs))[f][0] not in created_fits]
            for new_output in new_outputs:
                created_fits.append(new_output)
                print("\n"+new_output+" just finished\n")
                print("Elapsed time is: "+str(time.strftime("%H:%M:%S", time.gmtime(time.time()-initial_time)))+"\n")
                print(str(100*(len(created_fits)+int(str(outfile_copy).count("crash")))/len(list_to_make))+" percent done with "+fit+"\n")
                try:
                    if int(str(outfile_copy).count("crash"))>0:
                        print("Number that have crashed: "+str(str(outfile_copy).count("crash"))+"\n")
                    if int(str(outfile_copy).count("crash"))>crash_count:
                        crash_count=int(str(outfile_copy).count("crash"))
                        j=j-1
                except:
                    pass
                j=j-1
                if j<number_sent and i<len(list_to_make):
                    name="obj"+str(list_to_make[i])+".in"
                    print("Sending "+name+" to the shell")
                    process = subprocess.Popen([command,name],stdout=outfile)
                    sent_files.append(name)
                    i=i+1
                    j=j+1
            if i==len(list_to_make):
                print("All of them have been sent, just waiting for them to finish")
                i=i+1

    
        if display_blocks=="yes" and time.time()-reset_time>timer:
            safety_time=time.time()
            while time.time()-safety_time<timer:
                try:
                    model_out=fits.open(created_fits[k])[2].header
                    chisqr=model_out["CHI2NU"]
                    model_out.count("CHI2NU")
                    break
                except:
                    pass
                    break
            try:
                id=created_fits[k].split("_")[1]
                os.system("ds9 "+created_fits[k]+"[1] -fits "+created_fits[k]+"[2] -fits "+created_fits[k]+"[3] -fits "+cutout_folder_path+base_name+"_"+id+"_mask.fits &")
                k=k+1
                reset_time=time.time()
            except:
                pass
        outfile_copy_file.close()
        outfile_copy_file=open(path_to_outs+"outfile.txt","r")
        outfile_copy=outfile_copy_file.readlines()
        if len(list(set(new_outs)-set(orig_outs)))+int(str(outfile_copy).count("crash"))>=len(list_to_make):
            print("Done.")
            all_of_the_obj_ins=["obj"+str(integer)+".in" for integer in list_to_make]
            all_of_the_made_fits=["obj"+thing[0].split("_")[1]+".in" for thing in list(set(new_outs)-set(orig_outs))]
            real_crashed_list=sorted(list(set(all_of_the_obj_ins)-set(all_of_the_made_fits)))
            print("Crashed list should be "+str(real_crashed_list))
            crash_file=open(path_to_outs+"crashes.txt","w")
            crash_file.write(str(real_crashed_list))
            crash_file.close()
            break
        outfile_copy_file.close()
        
        
    if display_blocks=="yes":
        print("\nDisplaying rest of blocks")
        while True:
            if time.time()-reset_time>timer:
                try:
                    id=created_fits[k].split("_")[1]
                    os.system("ds9 "+created_fits[k]+"[1] -fits "+created_fits[k]+"[2] -fits "+created_fits[k]+"[3] -fits "+cutout_folder_path+base_name+"_"+id+"_mask.fits &")
                    k=k+1
                    reset_time=time.time()
                except:
                    print("Hopefully this is the end")
                    break



def output_dict2longfeedme(output_dict,SEx_in_dict,longfeedmename_path,list_to_make="default",longfeedme_outfit_type="sersic",n_start=1.5,min_cutout_radius=min_cutout_radius,final_mask_list=final_mask_list):
	###Converts an output dictionary to a longfeedme dictionary
    ###This is incomplete and not used, this would be useful for running using the output as input for the devauc fit etc.
    if list_to_make=="default":
        try:
            output_dict = np.load(output_dict).item()
        except:
            output_dict = output_dict
        list_to_make=sorted(output_dict.keys())
    new_dict=SEx_in_dict
    [new_dict[id].update(output_dict[id]) for id in list_to_make]
    return new_dict


def check_residuals(path_to_outputs,cutout_folder_path,list_to_check="default"):
    ### Checks residuals of original fits
    iraf.cd(path_to_outputs)
    os.chdir(path_to_outputs)
    out_files=[file for file in os.listdir(path_to_outputs) if "_out.fits" in file]
    if list_to_check!="default":
        out_files=[file for file in out_files if int(file.split("_")[1]) in list_to_check]
#    if start_id!="default" and list_to_check=="default":
#        out_files=[file for file in out_files if int(file.split("_")[1]) >= start_id]
    out_files_nums=sorted([int(file.split("_")[1]) for file in out_files])
    new_out_files=[]
    for id in out_files_nums:
        for file in out_files:
            if "_"+str(id)+"_" in file:
                new_out_files.append(file)
    print(new_out_files)
    ind=0
    while True:
        ask = raw_input("Press any key to continue(q to quit, int to go to specific galaxy):")
        if ask=="q":
            break
        try:
            if type(int(ask))==int:
                id=str(ask)
                for file1 in new_out_files:
                    if "_"+id+"_" in file1:
                        ind=new_out_files.index(file1)
                        file=file1
                        break
                if file!=file1:
                    print("\n\n\nCOULDN'T FIND OUTPUT FILE "+id+"\n\n\n")
                    
            else:
                id=new_out_files[ind].split("_")[1]
                file=new_out_files[ind]
        except:
            id=new_out_files[ind].split("_")[1]
            file=new_out_files[ind]
        try:
            os.system("ds9 "+file+"[1] -fits "+file+"[2] -fits "+file+"[3] -fits "+cutout_folder_path+base_name+"_"+id+"_mask.fits &")
            ind=ind+1
        except:
            print("This should never happen")
        try:
            print("The ID is: "+id)
            model_out=fits.open(file)[2].header
            print(model_out["C*MP_1"],)
            stats=model_out["1_*"]
            print(stats,)
            print(model_out["CHI2N*"],)
        except:
            print("This should never happen2")
    return model_out



def make_refits(obj_in_folder,input_cutout_folder,refit_folder_path,summary_file,list_to_refit=[]):
    ### Copies the obj.in and out.fits files (if they exist) to the refit folder so the refitting process can begin
    ###Slash or not at end is important for all inputs
    os.chdir(two_D)
    print("Performing sersic refits")
    try:
        os.mkdir(refit_folder_path)
    except:
        pass
    obj_files=[file for file in os.listdir(obj_in_folder) if "obj" and ".in" in file]
    obj_files=[file for file in obj_files if int(file.replace("obj","").replace(".in","")) in list_to_refit]
    obj_files_ints=sorted([int(file.replace("obj","").replace(".in","")) for file in obj_files])
    obj_files=["obj"+str(i)+".in" for i in obj_files_ints]
    for obj in obj_files:
        if os.path.exists(refit_folder_path+obj)==False:
            print(refit_folder_path+obj)
            os.system("cp "+obj_in_folder+obj+" "+refit_folder_path)
            try:
                if os.path.exists(refit_folder_path+base_name+"_"+obj.replace("obj","").replace(".in","")+"_out.fits")==False:
                    os.system("cp "+obj_in_folder+base_name+"_"+obj.replace("obj","").replace(".in","")+"_out.fits"+" "+refit_folder_path)
            except:
                print("here")
                print("Couldn't copy out file for "+obj)  ###Doesn't occur because error is disassociated when sent to shell
        
        ###I think I could put this in an else:, but I would need to test that
        txt=open(refit_folder_path+obj,"r").readlines()
        os.remove(refit_folder_path+obj)
        new_file=open(refit_folder_path+obj,"w")
        ###Here is where the line output destination is changed so that these "refits" don't replace the original, instead they land in the refit folder
        for line in txt:
            if line.startswith("B) ") and "/" in line:
                
                if "/" in line:
                    foo = line.split("/")[-1]
                    line="B) "+foo
                new_file.write(line)
                continue
            new_file.write(line)
        new_file.close()
    ind=0
    os.chdir(refit_folder_path)
    os.system("open "+obj_files[ind])
    os.system("open "+summary_file)
    while True:

        ask = raw_input("Press any key to continue(q to quit, int to go to specific galaxy, n for next):")
        
        if ask=="q":
            break
        if ask=="n":
            ind=ind+1
        try:
            if type(int(ask))==int:
                id=str(ask)
                for file1 in obj_files:
                    if "obj"+id+".in" in file1:
                        ind=obj_files.index(file1)
                        break
        except:
            pass
        os.system("open "+obj_files[ind])
        if os.path.exists(refit_folder_path+base_name+"_"+obj_files[ind].replace("obj","").replace(".in","")+"_out.fits"):
            file=base_name+"_"+obj_files[ind].replace("obj","").replace(".in","")+"_out.fits"
            id=obj_files[ind].replace("obj","").replace(".in","")
            
            try:
                os.system("ds9 "+file+"[1] -fits "+file+"[2] -fits "+file+"[3] -fits "+input_cutout_folder+base_name+"_"+id+"_mask.fits &")
            except:
                print("This should never happen")
            try:
                print("The ID is: "+id)
                model_out=fits.open(file)[2].header
                print(model_out["C*MP_1"],)
                stats=model_out["1_*"]
                print(stats,)
                print(model_out["CHI2N*"],)
            except:
                print("This should never happen2")
                
                
def find_crashed(obj_in_folder):
    ###Assumes outs are in obj_in_folder
    outs_ints=[int(file.split("_")[1]) for file in os.listdir(obj_in_folder) if "out.fits" in file]
    in_ints=[int(obj.replace("obj","").replace(".in","")) for obj in os.listdir(obj_in_folder) if "obj" and ".in" in obj]
    crashed=sorted(set(in_ints)-set(outs_ints))
    print(outs_ints)
    print(in_ints)
    print(crashed)
    return crashed

def make_not_deblended_additions_to_mask(initial_mask,further_deblended_segm,output_mask,ids_in_deblended_mask,start_value_for_ids_in_mask=1000):
    ###Makes a new mask which includes originally undetected objects (see manual)
    ###To keep track of input for documentation purposes
    name=output_mask.split("/")[-1]
    in_file=open(output_mask.replace(name,"in_file.txt"),"w")
    in_file.write("intial_mask:"+initial_mask+"\n")
    in_file.write("further_deblended_segm:"+further_deblended_segm+"\n")
    in_file.write("output_mask:"+output_mask+"\n")
    in_file.write("ids_in_deblended_mask:"+ids_in_deblended_mask+"\n")
    in_file.write("start_value_for_ids_in_mask:"+start_value_for_ids_in_mask+"\n")
    in_file.close()
    hdulist = fits.open(further_deblended_segm)
    deblended_mask_data = hdulist[0].data
    start_id=start_value_for_ids_in_mask
    for id in ids_in_deblended_mask:
        deblended_mask_data[deblended_mask_data==id]=start_id+id
    deblended_mask_data[deblended_mask_data<start_value_for_ids_in_mask]=0
    print(deblended_mask_data)
    initial_mask = fits.open(initial_mask)
    initial_mask_data = initial_mask[0].data
    deblended_mask_data=deblended_mask_data+initial_mask_data
    fits.writeto(output_mask,deblended_mask_data,overwrite=True)
    os.system("ds9 "+output_mask+" &")




def longfeedme_dict2longfeedme(longfeedme_dict,path_name_of_output_longfeedme):
    ###Converts a longfeedme dictionary to a longfeedme.txt file
    try:
        dic=np.load(longfeedme_dict).item()
    except:
        dic=longfeedme_dict
    list_to_make=sorted(dic.keys())
    file=open(path_name_of_output_longfeedme,"w")
    for id in list_to_make:
        file.write(str(id)+" \n")
        num_of_comps = num_comp(dic,id) 
        attribs = ['id', 'row', 'type', 'x_global', 'y_global', 'x_start', 'x_end', 'y_start', 'y_end', 'x_cutout', 'y_cutout', 'mag', 're', 'n', 'q', 'pa']
        for comp in range(1,num_of_comps+1):
            for attrib in attribs:
                file.write("    "+attrib+str(comp)+" "+str(dic[id][attrib+str(comp)])+"\n")
            file.write("\n")
        file.write("    "+"sky "+str(dic[id]["sky"])+"\n")
    file.close()
    return path_name_of_output_longfeedme

def enact_summary_of_changes_to_longfeedme(summary_of_changes_txt,long_feedme_dict,path_name_of_output_longfeedme,changes_type="local"):
    ###Enacts the changes in the summary file on the longfeedme file
    ###This would be much more effectively written and MUCH CLEANER if I had used a dictionary
    ###Perhaps I will rewrite this in the future
    ###Not super confident this will always do what you expect, please check the output to make sure
    ###Like I said this would be way simpler with a dictionary and then writing to file
    try:
        dic=np.load(long_feedme_dict).item()
    except:
        dic=long_feedme_dict
    attribs = ['id', 'row', 'type', 'x_global', 'y_global', 'x_start', 'x_end', 'y_start', 'y_end', 'x_cutout', 'y_cutout', 'mag', 're', 'n', 'q', 'pa']
    txt = open(summary_of_changes_txt,"r").readlines()
    ###This is added because when Laura ran it, she accidentally had a space and new line which caused problems, thus we just strip each line
    txt = [line1.strip() for line1 in txt]
    if changes_type=="local":
        for line in txt:
            if line.startswith("#")==False and len(line)>1:
                if line.split()[0]=="delete":
                    id,to_remove=int(line.split()[1]),int(line.split()[2])
                    num_comps=num_comp(dic,id)
                    if to_remove==1:
                        print("\nRemoving object #"+str(id))
                        del dic[id]
                        continue
                    else:
                        print("\nRemoving component #"+str(to_remove)+" from postage stamp #"+str(id))
                        for attrib in attribs:
                            break_button="no"
                            del dic[id][attrib+str(to_remove)] 
                            if to_remove<num_comps:
                                print("Downshifting the other components so they fill the gap left by component #"+str(to_remove))
                                for j in range(to_remove+1,num_comps+1):
                                    for attrib in attribs:
                                        dic[id][attrib+str(j-1)]=dic[id][attrib+str(j)]
                                for attrib in attribs:
                                    del dic[id][attrib+str(num_comps)]
                                    break_button="yes"
                            if break_button=="yes":
                                break
                if line.split()[0].startswith("id"):
                    id=int(line.split()[1])
                    idx=txt.index(line)
                    print("\nAdding the following entry:")
                    for attrib in attribs:
                        name,value=txt[idx].split()[0],txt[idx].split()[1]
                        print(name+" "+value)
                        dic[id].update({name:value})
                        idx=idx+1
                try:
                    id,attrib,old,new=int(line.split()[0]),line.split()[1],line.split()[2],line.split()[3]
                    print("\nChanging "+str(id)+" component #"+attrib[-1]+" to have "+str(attrib)+" = "+str(new))
                    dic[id].update({attrib:new})
                except:
                    pass
    if changes_type=="global":
        print(dic.keys())
        for line in txt:
            if line.startswith("#")==False and len(line)>1:
                if line.split()[0].startswith("id"):
                    id=int(line.split()[1])
                    idx=txt.index(line)
                    print("\nAdding the following entry:")
                    mini_dic=dict()
                    mini_dic[id]={"dummy":"2"}
                    for attrib in attribs:
                        value=txt[idx].split()[1]
                        print(attrib+" "+value)
                        mini_dic[id].update({attrib:value})
                        idx=idx+1
                    x_global,y_global=float(mini_dic[id]["x_cutout"])+int(mini_dic[id]["x_start"])-1,float(mini_dic[id]["y_cutout"])+int(mini_dic[id]["y_start"])-1
                    for id1 in sorted(dic.keys()):
                        x_start,x_end,y_start,y_end=float(dic[id1]["x_start1"]),float(dic[id1]["x_end1"]),float(dic[id1]["y_start1"]),float(dic[id1]["y_end1"])
                        if x_start<x_global<x_end and y_start<y_global<y_end:
                            num_comps=num_comp(dic,id1)
                            mini_dic[id]["x_cutout"],mini_dic[id]["y_cutout"]=x_global-x_start+1,y_global-y_start+1
                            mini_dic[id]["id"]=dic[id1]["id1"]
                            mini_dic[id]["x_start"],mini_dic[id]["x_end"],mini_dic[id]["y_start"],mini_dic[id]["y_end"]=int(x_start),int(x_end),int(y_start),int(y_end)
                            mini_dic[id]["x_global"],mini_dic[id]["y_global"]=x_global,y_global
                            distance_thresh=3
                            for j in range(1,num_comps+1):
                                distance=((float(mini_dic[id]["x_cutout"])-float(dic[id1]["x_cutout"+str(j)]))**2+(float(mini_dic[id]["y_cutout"])-float(dic[id1]["y_cutout"+str(j)]))**2)**0.5
                                if distance<distance_thresh:
                                    print("Instead of adding a new component, the old values of component #"+str(j)+" are being changed")
                                    print("Also maintaining the original row number "+str(dic[id1]["row"+str(j)])+" instead of changing it to "+str(mini_dic[id]["row"]))
                                    mini_dic[id]["row"]=dic[id1]["row"+str(j)]
                                    for attrib in attribs:
                                        dic[id1].update({attrib+str(j):mini_dic[id][attrib]})
                                    print("Check id number "+str(id1))
                                    break
                                if distance>distance_thresh and j==num_comps:
                                    print("\nAdding the entry above to "+str(id1))
                                    for attrib in attribs:
                                        dic[id1].update({attrib+str(num_comps+1):mini_dic[id][attrib]})
                                    break

                                    
                            
                try:
                    id1,attrib1,old,new1=int(line.split()[0]),line.split()[1],line.split()[2],line.split()[3]
                    if "x_cutout" in attrib1:
                        x_global=float(new1)+int(dic[id1]["x_start1"])-1
                    if "y_cutout" in attrib1:
                        y_global=float(new1)+int(dic[id1]["y_start1"])-1
                    for attrib in attribs:
                        if attrib in line.split()[1]:
                            number=attrib1.replace(attrib,"")
                            row=dic[id1]["row"+number]
                            id_list=dic.keys()
                            for id in id_list:
                                num_comps=num_comp(dic,id)
                                for j in range(1,num_comps+1):
                                    if dic[id]["row"+str(j)]==row:
                                        if "x_cutout" in attrib:
                                            new1=x_global-int(dic[id]["x_start1"])+1
                                        if "y_cutout" in attrib:
                                            new1=y_global-int(dic[id]["y_start1"])+1                                            
                                        dic[id].update({attrib+str(j):new1})
                                        print("\nChanging "+str(id)+" component #"+str(j)+" to have "+str(attrib)+" = "+str(new1))
                except:
                    pass

        for line in txt:
            if line.startswith("#")==False and len(line)>1:
                if line.split()[0]=="delete":
                    ###Need to get the row
                    id1,to_remove=int(line.split()[1]),int(line.split()[2])
                    try:
                        row=dic[id1]["row"+str(to_remove)]
                    except:
                        print("Object already removed")
                        continue
                    if to_remove==1:
                        print("\nRemoving object #"+str(id1))
                        del dic[id1]
                    id_list=sorted(dic.keys())
                    for id in id_list:
                        num_comps=num_comp(dic,id)
                        for j in range(1,num_comps+1):
                            break_button="no"
                            if dic[id]["row"+str(j)]==row:
                                print("\nRemoving object #"+str(row)+" from postage stamp #"+str(id))
                                for attrib in attribs:
                                    del dic[id][attrib+str(j)]
                                if j<num_comps:
                                    print("Downshifting the other components so they fill the gap left by component #"+str(j))
                                    for n in range(j+1,num_comps+1):
                                        for attrib in attribs:
                                            dic[id][attrib+str(n-1)]=dic[id][attrib+str(n)]
                                    for attrib in attribs:
                                        del dic[id][attrib+str(num_comps)]
                                        break_button="yes"
                            if break_button=="yes":
                                break
    longfeedme_dict2longfeedme(dic,path_name_of_output_longfeedme)
    

preamble_for_changes=['#Summarizes changes made to refit your_fit_type_here galaxies from GALFIT\n',
 '#These changes will be made to the your_long_feedme_here\n',
 '#Format is as follows:\n',
 '#ID of galaxy, with explanation for why it needs refitting\n',
 '#ID of central galaxy   parameter to change   old value   new value\n',
 '#If adding a galaxy to the feedme file, create a new entry as you would in the longfeedme file\n',
 '#If removing a galaxy write \xe2\x80\x9cdelete\xe2\x80\x9d  ID of central galaxy  Component number to be deleted\n',
 '#Additional comments and a new line, here are some examples:\n',
 '\n',
 '#90, several large errors on x,y,mag,q, and pa. Residual doesn\xe2\x80\x99t look good\n',
 '\n',
 '#21 mag1 21 23 (would be unhashed)\n',
 '#21 re2 5.14043432066 10 (would be unhashed)\n',
 '\n',
 '#Completely fixes the problem, residual and errors look great now\n',
 '\n',
 '#398, there is an extremely faint galaxy which GALFIT can\xe2\x80\x99t recognize\n',
 '\n',
 '#delete 398 5 (would be unhashed)\n',
 '\n',
 '#Fit works well now that GALFIT is not trying to fit this object, could also probably mask\n',
 '\n',
 '#799, need to add a star which is just outside the cutout (assume original has 2 objects)\n',
 '#Best way to do this is copy and paste a previous entry for the same id from the longfeedme and then edit\n',
 '\n',
 '#id3 799 (would be unhashed)\n',
 '#row3 2000 (would be unhashed) #make large enough that no other object will have this row\n',
 '#type3 psf (would be unhashed)\n',
 '#x_global3 1903.66259766\t (would be unhashed) #Only necessary to be accurate if changes_type="global"\n',
 '#y_global3 277.49807739 (would be unhashed) #same\n',
 '#x_start3 1753 (would be unhashed) \n',
 '#x_end3 2002  (would be unhashed)\n',
 '#y_start3 61  (would be unhashed)\n',
 '#y_end3 310  (would be unhashed)\n',
 '#x_cutout3 -5 (would be unhashed)\n',
 '#y_cutout3 190 (would be unhashed)\n',
 '#mag3 19 (would be unhashed)\n',
 '#re3 2 (would be unhashed)\n',
 '#n3 1.5 (would be unhashed)\n',
 '#q3 0.99 (would be unhashed)\n',
 '#pa3 -45 (would be unhashed)\n',
 '\n',
 '#Fit looks great now\n',
 '#If changes_type is global all these changes are implemented throughout the longfeedme, not just in the cutout which you refit\n\n',
 '##########################################################################']



def make_output_dict(folder_to_output_cutouts,list_to_make="default"):
    ###Scraps the output from the 2nd extensions of the out file and makes a dictionary
    output_dict=dict()
    if list_to_make=="default":
        list_of_outs=[file for file in os.listdir(folder_to_output_cutouts) if "_out.fits" in file]
    else:
        list_of_outs=[file for file in os.listdir(folder_to_output_cutouts) if "_out.fits" in file and int(file.split("_")[1]) in list_to_make]
    for out in list_of_outs:
        model_out=fits.open(folder_to_output_cutouts+out)[2].header
        id = int(model_out["DATAIN"].split("_")[1])
        output_dict[id]={"id1":id}
        num_comps=len(model_out["COMP*"])-1
        m=model_out
        for i in range(1,num_comps+1):
            i=str(i)
            comp_type,x,x_e,y,y_e=m["COMP_"+i],m[i+"_XC"].split()[0],m[i+"_XC"].split()[2],m[i+"_YC"].split()[0],m[i+"_YC"].split()[2]
            mag,mag_e=m[i+"_MAG"].split()[0],m[i+"_MAG"].split()[2]
            if comp_type=="sersic":
                n,n_e=m[i+"_N"].split()[0],m[i+"_N"].split()[2]
            if comp_type=="devauc":
                n,n_e="4","0"
            if comp_type!="psf":
                r,r_e=m[i+"_RE"].split()[0],m[i+"_RE"].split()[2]
                q,q_e,pa,pa_e=m[i+"_AR"].split()[0],m[i+"_AR"].split()[2],m[i+"_PA"].split()[0],m[i+"_PA"].split()[2]
                all_values=[id,comp_type,x,x_e,y,y_e,mag,mag_e,r,r_e,n,n_e,q,q_e,pa,pa_e]
                values_names=["id","type","x_cutout","x_cutout_err","y_cutout","y_cutout_err",
                              "mag","mag_err","re","re_err","n","n_err","q","q_err","pa","pa_err"]
                zipped=zip(values_names,all_values)
                for j in range(len(zipped)):
                    name,value=zipped[j][0]+i,zipped[j][1]
                    output_dict[id].update({name:value})
            if comp_type=="psf":
                all_values=[id,comp_type,x,x_e,y,y_e,mag,mag_e]
                values_names=["id","type","x_cutout","x_cutout_err","y_cutout","y_cutout_err",
                              "mag","mag_err"]
                zipped=zip(values_names,all_values)
                for j in range(len(zipped)):
                    name,value=zipped[j][0]+i,zipped[j][1]
                    output_dict[id].update({name:value})
        output_dict[id].update({"chi2nu":m["CHI2NU"],"sky":m[str(num_comps+1)+"_SKY"].replace("]","").replace("[",""),"num_comps":num_comps-1})
    np.save(folder_to_output_cutouts+"output_dict.npy",output_dict)
    return output_dict

def output_dict2tab(output_dict,SEx_in_dict,cutout_dict,input_dict,output_tab_path,list_to_make="default",plt_scale=plt_scale):
    ###Converts the output dictionary made above into the final table
    ###I copied the math straight in from the other scripts but I should double check this is correct
    try:
        dic=np.load(output_dict).item()
    except:
        dic=output_dict
    try:
        SEx_in_dict=np.load(SEx_in_dict).item()
    except:
        pass
    try:
        input_dict=np.load(input_dict).item()
    except:
        pass
    try:
        cutout_dict=np.load(cutout_dict).item()
    except:
        pass
    #Create the table text file 
    file=open(output_tab_path,"w")
    if list_to_make=="default":
        list_to_make=sorted(dic.keys())
    
    #utilize output_dict to get results of GALFIT
    for id in list_to_make:
        print(id)
        mtot,mtot_err,aeff,aeff_err=dic[id]["mag1"],dic[id]["mag_err1"],dic[id]["re1"],dic[id]["re_err1"]
        ar,ar_err,n,n_err,x,x_err,y,y_err=dic[id]["q1"],dic[id]["q_err1"],dic[id]["n1"],dic[id]["n_err1"],dic[id]["x_cutout1"],dic[id]["x_cutout_err1"],dic[id]["y_cutout1"],dic[id]["y_cutout_err1"]
        eps,eps_err,pa,pa_err=1-float(ar),ar_err,dic[id]["pa1"],dic[id]["pa_err1"]
        if mtot_err!="nan":
            reff,reff_err=float(aeff)*np.sqrt(float(ar))*plt_scale,np.sqrt(float(ar)*(float(aeff_err)*plt_scale)**2+float(ar_err)**2*(float(aeff)*plt_scale)**2/(4*float(ar)))
            logre,logre_err=np.log10(float(reff)),(float(reff_err)/float(reff))*(1./np.log(10.))
            mueff,mueff_err=float(mtot)+2.5*np.log10(2.*3.14159*float(reff)**2),np.sqrt(float(mtot_err)**2+4.71529*(float(reff_err)/float(reff))**2)
            ###ADDRESS Need to check those
        if mtot_err=="nan":
            print("nan for "+str(id))
            reff,reff_err=float(aeff)*np.sqrt(float(ar))*plt_scale,"nan"
            logre,logre_err=np.log10(float(reff)),"nan"
            mueff,mueff_err=float(mtot)+2.5*np.log10(2.*3.14159*float(reff)**2),"nan"
        sky,chisqr,comps=dic[id]["sky"],dic[id]["chi2nu"],dic[id]["num_comps"]
    
        #utilize the SEx_in_dict to figure out the input SExtractor values
        SEx_mag,SEx_aeff,SEx_q,SEx_pa,SEx_x,SEx_y=SEx_in_dict[id]["mag1"],SEx_in_dict[id]["re1"],SEx_in_dict[id]["q1"],SEx_in_dict[id]["pa1"],SEx_in_dict[id]["x_global1"],SEx_in_dict[id]["y_global1"]
        SEx_aeff=str(float(SEx_aeff)*plt_scale)
    
    
        #utilize input_dict because we need it for flag #1 to figure out input re
        flag=0
        for i in range(2,comps+2):
            x_comp,y_comp=dic[id]["x_cutout"+str(i)],dic[id]["y_cutout"+str(i)]
            if np.sqrt((float(x_comp)-float(x))**2+(float(y_comp)-float(y))**2)<float(input_dict[id]["re1"]):
                flag=flag+1
                break
        if float(chisqr)>3:
            flag=flag+2
        for i in range(2,comps+2):
            mag_err_comp=dic[id]["mag_err"+str(i)]
            mag_comp=dic[id]["mag"+str(i)]
            if float(mag_err_comp)/float(mag_comp)>0.5:
                flag=flag+4
                break
        for i in range(2,comps+2):
            try:
                re_err_comp=dic[id]["re_err"+str(i)]
            except:
                print("Is psf so no re")
                continue
            re_comp=dic[id]["re"+str(i)]
            if float(re_err_comp)/float(re_comp)>0.5:
                flag=flag+8
                break
        if abs(float(mtot)-float(SEx_mag))>1:
            flag=flag+16
        #omitting flag 32
        if float(n)>4.5:
            flag=flag+64
            
        #utilize cutout_dict to get global position
        x_global,y_global=float(cutout_dict[id]["x_start1"])+float(x)-1,float(cutout_dict[id]["y_start1"])+float(y)-1
        master_list=[id,mtot,mtot_err,aeff,aeff_err,reff,reff_err,logre,logre_err,mueff,mueff_err,n,n_err,
                    ar,ar_err,eps,eps_err,pa,pa_err,x_global,x_err,y_global,y_err,sky,chisqr,comps,flag,
                    SEx_mag,SEx_aeff,SEx_q,SEx_pa,SEx_x,SEx_y]
        master_list=[str(attrib) for attrib in master_list]
        rounding_list=[None,3,3,2,2,3,3,3,3,3,3,2,2,3,3,3,3,1,1,2,2,2,2,4,2,None,None,3,3,3,1,2,2]
        for i in range(len(rounding_list)):
            if rounding_list[i] is not None:
                f="{:."+str(rounding_list[i])+"f}"
                master_list[i]=f.format(float(master_list[i]))
        master_list=[format(str(i),'>10') for i in master_list]
        file.write("".join(master_list)+'\n')
    file.close()


def make_fits(outstats_file):
    ###Converts the text file made above into a fits table
    os.chdir(two_D)
    iraf.cd(two_D)
    try:
        name = outstats_file.replace(".txt",".fits")
    except:
        pass
    cd_file = open("outtab.cd","w")
    columns=['NUMBERI', 'MTOT_SER', 'E_MTOT_SER', 'AE_SER', 'E_AE_SER',
             'RE_SER', 'E_RE_SER', 'LRE_SER', 'E_LRE_SER', 'MUE_SER',
             'E_MUE_SER', 'N_SER', 'E_N_SER', 'ARATIO_SER', 'E_ARATIO_SER',
             'EPS_SER', 'E_EPS_SER', 'PA_SER', 'E_PA_SER', 'X_SER', 'E_X_SER',
             'Y_SER', 'E_Y_SER', 'SKY_SER', 'CHISQ_SER', 'NCOMP_SER', 'FLAG_SER',
             'MAGin_SER', 'AEFFin_SER', 'ARATIOin_SER', 'PAin_SER', 'Xin_SER', 'Yin_SER']
    types=['i', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r',
             'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r',
             'r', 'r', 'r', 'r', 'r', 'i', 'i', 'r', 'r', 'r',
             'r', 'r', 'r']

    units=["", "", "", 'pixels', 'pixels', 'arcseconds', 'arcseconds',
           "", "", "", "", "", "", "", "", "", "",
           'degrees', 'degrees', "", "", "", "", "", "", "",
           "", "", "", "", "", "", ""]
    disps = ['I4', 'F6.3', 'F6.3', 'F6.2', 'F6.2', 'F6.3', 'F6.3', 'F6.3',
             'F6.3', 'F6.3', 'F6.3', 'F4.2', 'F4.2', 'F5.3', 'F5.3', 'F5.3',
             'F5.3', 'F5.1', 'F5.1', 'F7.2', 'F7.2', 'F7.2', 'F7.2', 'F8.4',
             'F7.2', 'I4', 'I3', 'F7.3', 'F7.3', 'F7.3', 'F7.1', 'F7.2', 'F7.2']
    if "devauc" in name:
        columns=[i.replace("SER","DEV") for i in columns]
    for i in range(len(columns)):
        cd_file.write(columns[i]+" "+types[i]+" "+disps[i]+" "+units[i]+" \n")
    cd_file.close()
    temp_file_for_nans,temp_file_for_nans_name,temp_file_readlines=open(outstats_file.replace(".txt","_temp.txt"),"w"),outstats_file.replace(".txt","_temp.txt"),open(outstats_file,"r").readlines()
    [temp_file_for_nans.write(line.replace("nan","INDEF")) for line in temp_file_readlines]
    temp_file_for_nans.close()
    iraf.delete(name)
    iraf.tcreate(name,"outtab.cd",temp_file_for_nans_name)
    iraf.delete(temp_file_for_nans_name)
    iraf.delete("outtab.cd")


def combine_refits(refit_folder,original_folder,output_folder,final_dict,overwrite="no",make_changes_list=[]):
    ###Combines the original fits with the refits into the final folder
    try:
        os.mkdir(output_folder)
    except:
        pass
    try:
        dic=np.load(final_dict).item()
    except:
        dic=final_dict
    final_list=sorted(dic.keys())
    
    for id in final_list:
        if overwrite=="no":
            if os.path.exists(output_folder+"obj"+str(id)+".in")==False and os.path.exists(output_folder+base_name+"_"+str(id)+"_out.fits")==False:
                os.chdir(original_folder)
                os.system("cp -rf obj"+str(id)+".in "+output_folder)
                os.system("cp -rf "+base_name+"_"+str(id)+"_out.fits "+output_folder)
                os.chdir(refit_folder)
                os.system("cp -rf obj"+str(id)+".in "+output_folder)
                os.system("cp -rf "+base_name+"_"+str(id)+"_out.fits "+output_folder)
            else:
                print("Already in the final folder")
                break
        if overwrite=="yes":
            os.chdir(original_folder)
            os.system("cp -rf obj"+str(id)+".in "+output_folder)
            os.system("cp -rf "+base_name+"_"+str(id)+"_out.fits "+output_folder)
            if id in make_changes_list:
                os.chdir(refit_folder)
                os.system("cp -rf obj"+str(id)+".in "+output_folder)
                os.system("cp -rf "+base_name+"_"+str(id)+"_out.fits "+output_folder)


def plot_stats(folder_with_outs):
    ###Plots statistics of an output folder
    iraf.cd(folder_with_outs)
    os.chdir(folder_with_outs)
    list_of_outs=[file for file in os.listdir(folder_with_outs) if "_out.fits" in file]    
    j=1
    id_list,x_list,x_e_list=[],[],[]
    y_list,y_e_list,n_list=[],[],[]
    n_e_list,r_list,r_e_list=[],[],[]
    q_list,q_e_list,pa_list=[],[],[]
    mag_list,mag_e_list,pa_e_list,refit_list=[],[],[],[]
    for out in list_of_outs:
        m=fits.open(out)[2].header
        i="1"
        id=out.split("_")[1]
        comp_type,x,x_e,y,y_e=m["COMP_"+i],m[i+"_XC"].split()[0],m[i+"_XC"].split()[2],m[i+"_YC"].split()[0],m[i+"_YC"].split()[2]
        mag,mag_e=m[i+"_MAG"].split()[0],m[i+"_MAG"].split()[2]
        if comp_type=="sersic":
            n,n_e=m[i+"_N"].split()[0],m[i+"_N"].split()[2]
            n_list.append(n)
            n_e_list.append(n_e)
        r,r_e=m[i+"_RE"].split()[0],m[i+"_RE"].split()[2]
        q,q_e,pa,pa_e=m[i+"_AR"].split()[0],m[i+"_AR"].split()[2],m[i+"_PA"].split()[0],m[i+"_PA"].split()[2]
        id_list.append(id),x_list.append(x),x_e_list.append(x_e)
        y_list.append(y),y_e_list.append(y_e),r_list.append(r)
        r_e_list.append(r_e),q_list.append(q),q_e_list.append(q_e)
        pa_list.append(pa),pa_e_list.append(pa_e),mag_list.append(mag),mag_e_list.append(mag_e)

    plt.figure(j)
    plt.scatter(x_list,x_e_list)
    [plt.annotate(id_list[l],xy=(x_list[l],x_e_list[l])) for l in range(len(id_list))]
    plt.xlabel("x"),plt.ylabel("x_err")
    plt.figure(j+1)
    plt.scatter(y_list,y_e_list)
    [plt.annotate(id_list[l],xy=(y_list[l],y_e_list[l])) for l in range(len(id_list))]
    plt.xlabel("y"),plt.ylabel("y_err")
    plt.figure(j+2)
    plt.scatter(mag_list,mag_e_list)
    [plt.annotate(id_list[l],xy=(mag_list[l],mag_e_list[l])) for l in range(len(id_list))]
    plt.xlabel("mag"),plt.ylabel("mag_err")
    plt.figure(j+3)
    if comp_type=="sersic":
        plt.scatter(n_list,n_e_list)
        [plt.annotate(id_list[l],xy=(n_list[l],n_e_list[l])) for l in range(len(id_list))]
        plt.xlabel("n"),plt.ylabel("n_err")
        plt.figure(j+4)
    plt.scatter(r_list,r_e_list)
    [plt.annotate(id_list[l],xy=(r_list[l],r_e_list[l])) for l in range(len(id_list))]
    plt.xlabel("r"),plt.ylabel("r_err")
    plt.figure(j+5)
    plt.scatter(q_list,q_e_list)
    [plt.annotate(id_list[l],xy=(q_list[l],q_e_list[l])) for l in range(len(id_list))]
    plt.xlabel("q"),plt.ylabel("q_err")
    plt.figure(j+6)
    plt.scatter(pa_list,pa_e_list)
    [plt.annotate(id_list[l],xy=(pa_list[l],pa_e_list[l])) for l in range(len(id_list))]
    plt.xlabel("pa"),plt.ylabel("pa_err")
    for n in range(len(x_e_list)):
        if float(x_e_list[n])>=1:
            print("#"+id_list[n]+", has x_err >= 1")
            refit_list.append(id_list[n])
        if float(y_e_list[n])>=1:
            print("#"+id_list[n]+", has y_err >= 1")
            refit_list.append(id_list[n])
        if float(mag_e_list[n])>=1:
            print("#"+id_list[n]+", has mag_err >= 1")
            refit_list.append(id_list[n])
        if comp_type=="sersic":
            if float(n_list[n])>=6:
                refit_list.append(id_list[n])
                print("#"+id_list[n]+", has n >= 6")
            if float(n_e_list[n])>=1:
                refit_list.append(id_list[n])
                print("#"+id_list[n]+", has n_err >= 1")
        if float(q_e_list[n])>=1:
            refit_list.append(id_list[n])
            print("#"+id_list[n]+", has q_err >= 1")
        if float(r_e_list[n])>=1:
            refit_list.append(id_list[n])
            print("#"+id_list[n]+", has r_err >= 1")
        if x_e_list[n]=="nan":
            print("#"+id_list[n]+", has nan error")
            refit_list.append(id_list[n])
        if x_e_list[n]=="nan":
            print("#"+id_list[n]+", has nan error")
            refit_list.append(id_list[n])
        if float(mag_e_list[n])<=0.0001:
            print("#"+id_list[n]+", has low error (mag_err<=0.0001)")
            refit_list.append(id_list[n])
        if id_list[n] in refit_list:
            print("\n")
    refit_list=[int(n) for n in refit_list]
    print(sorted(list(set(refit_list))))
    print("You need to add ones that crashed")
    iraf.cd(two_D)
    os.chdir(two_D)








base_name = base_folder.split('/')[-1]
plt.close("all")
### Begins script
### Makes sure we are in base_folder
if os.path.exists(base_folder)==False:
    os.mkdir(base_folder)

### In case you supply the full path to images, they get renamed accordingly
os.chdir(base_folder)
iraf.cd(base_folder)
zip(["drz","wht","ctx"],[drz_image,wht_image,ctx_image])
for im in zip(["drz","wht","ctx"],[drz_image,wht_image,ctx_image]):
    if os.path.exists(base_folder+base_name+"_"+im[0]+".fits")==False:
        prior_name=im[1].split("/")[-1]
        os.system("cp "+im[1]+" .")
        os.rename(prior_name,base_name+"_"+im[0]+".fits")
drz_image = base_name+"_drz.fits"  
wht_image = base_name+"_wht.fits"  
ctx_image = base_name+"_ctx.fits"  

### Primitively checks that there are only 3 fits files in base_folder, I don't know if this is really necessary
if sum([i.count('.fits') for i in os.listdir(base_folder)])!=3:
    print("WARNING: Are your ctx, wht, and drz images in your base_folder??? Should only be those in there.")
    time.sleep(10)


### If make_EXPTIME_1 is 'yes' this hedits them to 1 if it hasn't been done already, tested done
if make_EXPTIME_1=='yes':
    if os.path.exists(base_folder+"/important_parameters.txt")==False:
        important_parameters = open("important_parameters.txt","w")
        for line in iraf.imhead(base_folder+"/*.fits",long="yes",Stdout=1):
            if line.startswith("FILENAME"):
                important_parameters.write("Original EXPTIME:"+"\n")
                important_parameters.write(line+'\n')
            if line.startswith("EXPTIME"):
                important_parameters.write(line+'\n'+'\n')
        important_parameters.close()
    EXP_TIME_list=[]
    for line in iraf.imhead(base_folder+"/*.fits",long="yes",Stdout=1):
        if line.startswith("EXPTIME"):
            EXP_TIME_list.append(float(line.split()[2]))
    if make_EXPTIME_1=='yes' and sum(EXP_TIME_list)!=3: #dirty check, if ctx,drz,wht are all EXPTIME=1 then this should be 3
        print("Changing EXPTIME to 1"+'\n'+'\n')
        iraf.hedit("*.fits","EXPTIME",1,verify="no",show="yes",update="yes")

    

### If directory structure already exists, skip this, tested done

if os.path.exists(base_folder+"/"+base_name+"_imaging")==False or os.path.exists(base_folder+"/"+base_name+"_phot")==False or os.path.exists(base_folder+"/"+base_name+"_2D")==False:

    ### Create directory structure
    ### Should have three images: drz w/ background, wht, and ctx in a folder base_folder
    ### Creates base_folder_2D, base_folder_phot, and base_folder_imaging within the base_folder
    if os.path.exists(base_folder+"/"+base_name+"_2D")==False:
        os.mkdir(base_folder+"/"+base_name+"_2D")
    if os.path.exists(base_folder+"/"+base_name+"_phot")==False:
        os.mkdir(base_folder+"/"+base_name+"_phot")
    if os.path.exists(base_folder+"/"+base_name+"_imaging")==False:
        os.mkdir(base_folder+"/"+base_name+"_imaging")
    ### Copies files to these directories
    os.system("cp "+base_folder+"/*.fits"+" "+base_folder+"/"+base_name+"_2D")
    os.system("cp "+base_folder+"/*.fits"+" "+base_folder+"/"+base_name+"_phot")
    os.system("cp "+base_folder+"/*.fits"+" "+base_folder+"/"+base_name+"_imaging")

### Nicknames these folders for easier access
phot=base_folder+"/"+base_name+"_phot/"
two_D=base_folder+"/"+base_name+"_2D/"
imaging=base_folder+"/"+base_name+"_imaging/"



### Makes ncomb, weight, and sigma image 
if make_ncomb=='yes':
    iraf.cd(imaging)
    os.chdir(imaging)
    if os.path.exists(base_folder+"/"+base_name+"_imaging/"+base_name+"_ncomb.fits"):
        print("Ncomb image has already been made. Now deleting and remaking."+'\n'+'\n')
    print('\n'+'\n'+'\n'+"Making ncomb image..."+'\n'+'\n')
    def convert(binary):
        binary = bin(binary).replace('0b','')
        binary = list(str(binary))
        binary = [int(x) for x in binary]
        return sum(binary)
    try:
        iraf.delete(base_name+"_ncomb.fits")
    except:
        pass
    iraf.imcopy(ctx_image,base_name+"_ncomb.fits")
    hdulist = fits.open(base_name+"_ncomb.fits")
    ncomb_data = hdulist[0].data
    #for i in range(0,len(ncomb_data)):
    #    for j in range(0,len(ncomb_data[i])):
    #        ncomb_data[i][j]=convert(ncomb_data[i][j])
    #print(np.max(ncomb_data))
    for i in range(0,int(np.max(ncomb_data))+1):
        ncomb_data[ncomb_data==i]=convert(i)
    hdulist.writeto(base_name+"_ncomb.fits",overwrite=True)
    print("Done making ncomb image. Displaying...")
    os.system("ds9 "+base_name+"_ncomb.fits &")
if make_weight=='yes':
    if os.path.exists(base_folder+"/"+base_name+"_imaging/"+base_name+"_weight.fits"):
        print("Weight image has already been made. Now deleting and remaking."+'\n'+'\n')
    print('\n'+'\n'+'\n'+"Making weight image..."+'\n'+'\n')
    try:
        iraf.delete(imaging+"*_weight.fits")
        #iraf.delete(imaging+"n_*")
    except:
        pass
    iraf.cd(imaging) ### IMPORTANT
    iraf.imexpr("a>=b?c:0",base_name+"_weight.fits",base_name+"_ncomb.fits",n_min,wht_image)
    #iraf.mkweight(imaging+base_name+"_weight.fits",imaging+base_name+"_ncomb.fits",imaging+wht_image,nmin=n_min)
    #iraf.imreplace(imaging+base_name+"_weight.fits", 0, lower=0.9999999, upper=1.00000001) #Changes 1 values to zero, this helped me for the SExtractor run
    #iraf.delete(imaging+"n_*")
    print("Done making weight image. Displaying...")
    os.system("ds9 "+imaging+base_name+"_weight.fits &")
if make_sigma=='yes':
    if os.path.exists(base_folder+"/"+base_name+"_imaging/"+base_name+"_sigma.fits"):
        print("Sigma image has already been made. Now deleting and remaking."+'\n'+'\n')
    print('\n'+'\n'+'\n'+"Making sigma image..."+'\n'+'\n')
    os.chdir(imaging)
    iraf.cd(imaging)
    iraf.delete(imaging+"*_norm.fits")
    iraf.delete(imaging+"*_sigma.fits") 
    
    ###Perhaps instead of the wht_image * norm_factor this should be weight_image * norm_factor
    
    iraf.imarith(imaging+wht_image,"*",norm_factor,imaging+base_name+"_norm.fits") #This is to create "exposure map" which will be the normalization for sigma image  
    iraf.imexpr("sqrt((a)/c+b*(d/c)**2 )*e",imaging+base_name+"_sigma.fits",imaging+drz_image,imaging+base_name+"_ncomb.fits",imaging+base_name+"_norm.fits",read_noise,sqrt_FA)
    print("Done making sigma image. Displaying...")
    os.system("ds9 "+imaging+base_name+"_sigma.fits &")




### Runs SExtractor (first time)
try:
    del CHECKIMAGE_TYPE,CHECKIMAGE_NAME,CATALOG_NAME,WEIGHT_IMAGE,WEIGHT_TYPE
except:
    pass
### These are automatically chosen so that they have the right names for the rest of the script
CHECKIMAGE_TYPE = "SEGMENTATION,APERTURES,BACKGROUND"
CHECKIMAGE_NAME = base_name+"_segm.fits,"+base_name+"_aper.fits,"+base_name+"_background.fits"
CATALOG_NAME = base_name+"_tab_all.fits"
WEIGHT_IMAGE = imaging+base_name+"_weight.fits"
WEIGHT_TYPE = "MAP_WEIGHT"

if run_sex1=='yes':
    os.chdir(phot)
    iraf.cd(phot)
    if print_vars=="yes":
        ### This is kind of shaky, but is supposed to print all the variables within a .sex config file
        txt = open(path_to_sex_config,"r").readlines()
        for line in txt:
            if line.startswith("#")==False and line[0].isupper() and len(line)>0:
                if "#" in line:
                    name,entry=line.strip("\n").split()[0],line.strip("\n").split()[1:]
                    idx=entry.index([i for i in entry if "#" in i][0])
                    comment = " ".join(entry[idx:])
                    entry = "".join(entry[:idx])
                    print(name+' = "'+entry+'"     '+comment)
                else:
                    print(line.split()[0]+' = "'+line.split()[1]+'"') #in the case you have one of my ugly formatted ones
                
        sys.exit()
    txt=""
    file1=open(base_name+".sex","w")
    file1.write("#SExtractor Configuration File\n")
    file1.write("#Not formatted beautifully, but will get the job done\n")
    for item in dir():
        if item.isupper() and "=" not in str(item):
            file1.write(str(item)+" "+str(eval(item))+'\n')
    file1.close()
        
    
    ### Creates runsex.csh and runs first SExtractor run
    file2=open("runsex.csh","w")
    file2.write("#!/bin/csh"+"\n"+"\n")
    file2.write("sex "+drz_image+" -c "+base_name+".sex")
    file2.close()
    ### Might be better to use Popen and direct the stdout and stderr to a file (like a do with GALFIT)
    os.system("source runsex.csh")
    SEx_tab_all = phot+base_name+"_tab_all.fits"
    #os.system("ds9 "+phot+base_name+"_aper.fits &")
    ### Reviews output of SExtractor first run if auto_review_sex_output is yes
    if auto_review_sex_output=='yes' and os.path.exists(phot+base_name+"_tab_all.fits"):
        #os.system("ds9 "+phot+base_name+"_aper.fits &")
        image = fits.open(base_name+"_tab_all.fits", memmap=True)
        image_data = Table(image[1].data)
        image_data.write("SExtractor_output.csv", overwrite=True)
        os.system("open SExtractor_output.csv")




if make_assoc=="yes":
    os.chdir(phot)
    iraf.cd(phot)
    image = fits.open(base_name+"_tab_all.fits")
    image_data = Table(image[1].data)
    ### Future note to self, when writing a longer script with multiple plots elect an integer like z, declare it at the beginning and then
    ### Update figures as they are made e.g. plt.figure(z+1), z=z+1
    
    ### Makes the MAG_AUTO vs MAG_AUTO_ERR and MAG_AUTO vs. CLASS_STAR plots
    plt.figure(1)
    plt.scatter(image_data["MAG_AUTO"], image_data["CLASS_STAR"], c="k", marker="x")
    plt.axhline(y=assoc_class_star_cut, color='r')
    plt.axvline(x=assoc_mag_cut, color='r')
    plt.xlabel("MAG_AUTO")
    plt.ylabel("CLASS_STAR")
    plt.figure(2)
    plt.scatter(image_data["MAG_AUTO"], image_data["MAGERR_AUTO"], c="k", marker="x")
    plt.axhline(y=assoc_noise_cut, color='r')
    plt.axvline(x=assoc_mag_cut, color='r')
    plt.xlabel("MAG_AUTO")
    plt.ylabel("MAGERR_AUTO")
    file = open(base_name+"_ALLsex.cat","w")
    iraf.delete("region_file.reg")
    file1=open("region_file.reg","a+")
    ### All these columns need to be in the param file
    columns=["NUMBER","FLUX_AUTO","FLUXERR_AUTO","MAG_AUTO","MAGERR_AUTO","KRON_RADIUS","BACKGROUND","ISOAREA_IMAGE","FLUX_RADIUS","X_IMAGE","Y_IMAGE","ALPHA_J2000","DELTA_J2000","THETA_IMAGE","ELLIPTICITY","FWHM_IMAGE","FLAGS","CLASS_STAR"]
    try:
        reduced_data = image_data[columns]
    except:
        print("Do you have "+" ".join(columns)+" in your .param file?")
        sys.exit()
    select_data = image_data["MAG_AUTO","MAGERR_AUTO","CLASS_STAR"]
    k=0
    for i in range(0,len(image_data)):
        if select_data[i]["MAGERR_AUTO"]<assoc_noise_cut and select_data[i]["MAG_AUTO"]<assoc_mag_cut and select_data[i]["CLASS_STAR"]<assoc_class_star_cut:
            file.write(str(np.array(reduced_data[i])).replace(")","").replace("(","").replace(",","")+'\n')
            file1.write("text "+str(reduced_data[i]["X_IMAGE"]+3)+" "+str(reduced_data[i]["Y_IMAGE"]+3)+" # text ={"+str(i+1)+"}"+'\n')
            k=k+1
        else:
            ### If it doesn't meet the criteria it gets a zero for its id
            reduced_data[i][0]=0
            file.write(str(np.array(reduced_data[i])).replace(")","").replace("(","").replace(",","")+'\n')
            #print(i)
    file.close()
    file1.close()
    os.system("ds9 "+base_name+"_aper.fits"+" -regions load region_file.reg &")
    print("There are "+str(k)+" that meet the criteria")
    magbright,magfaint=min(select_data["MAG_AUTO"]),max(select_data["MAG_AUTO"])
    maglimit_guess=min(magbright+6.2,magfaint+0.5)
    print("min(magbright+6.2,magfaint+0.5)="+str(maglimit_guess))

    
    
    
if make_mask=='yes':
    os.chdir(phot)
    iraf.cd(phot)
    file1 = open(base_name+"_ALLsex.cat","r")
    txt1 = file1.readlines()
    unmask_list=[]
    total_list=[]
    k=1
    for line in txt1:
        if float(line.split()[3])<=mask_maglim:
            unmask_list.append(int(k))
        total_list.append(k)
        k=k+1
    for i in add_to_mask_list:
        unmask_list.remove(i)
    for i in remove_from_mask_list:
        unmask_list.append(i)
    unmask_list=sorted(unmask_list)
    hdulist = fits.open(phot+base_name+"_segm.fits")
    segm_data = hdulist[0].data
    drz_image_data = fits.open(phot+drz_image)[0].data
    
    ### The background value of 9 here is hard coded, could cause problems. I'm going to change this right now
    ### Gets a value which isn't in the unmask list and 
    ### sets this value to be the background
    backgrnd_val=sorted(list(set(range(1,max(unmask_list)))-set(unmask_list)))[-1]
    
    
    segm_data[drz_image_data==0]=backgrnd_val  
    for i in unmask_list:
        segm_data[segm_data==i]=0    
    if make_mask_name=="default":
        hdulist.writeto(phot+base_name+"_mask.fits",overwrite=True)
        os.system("ds9 "+phot+base_name+"_mask.fits &")
    else:
        hdulist.writeto(phot+base_name+"_mask_"+make_mask_name+".fits",overwrite=True)
        os.system("ds9 "+phot+base_name+"_mask_"+make_mask_name+".fits &")
    final_mask_list=sorted(set(total_list)-set(unmask_list))
    if make_mask_name=="default":
        mask_file=open(phot+base_name+"_mask_list.txt","w")
    else:
        mask_file=open(phot+base_name+"_mask_"+make_mask_name+"_list.txt","w")
    image = fits.open(base_name+"_tab_all.fits")
    image_data = Table(image[1].data)
    [mask_file.write(str(id)+" "+str(image_data[id-1]["X_IMAGE"])+" "+str(image_data[id-1]["Y_IMAGE"])+"\n") for id in final_mask_list]
    mask_file.close()
    
    
    
        
    
    
    
    
    
if make_sex_sky=='yes':
    #This could use some improvement, I think I will rewrite it in the future
    os.chdir(phot)
    iraf.cd(phot)
    ### Not every iraf installation has this packages, couldn't get these on my vm at home
    iraf.stsdas(_doprint=0)
    iraf.hst_calib(_doprint=0)
    iraf.nicmos(_doprint=0)
    if os.path.exists(base_name+"_itersky.cat"):
        print("\nYou can make changes to the sky values easily by giving the appropriate id's to change_to_iterstat_list and take_average_list")
        cat = open(base_name+"_ALLsex.cat","r").readlines()
        iter_cat = open(base_name+"_itersky.cat","r").readlines()
        final_cat = open(base_name+"_ALLsex_final.cat","w")
        for line in cat:
            if int(line.split()[0]) not in change_to_iterstat_list and int(line.split()[0]) not in take_average_list and int(line.split()[0]) not in remove_entirely_list:
                final_cat.write(line)
            if int(line.split()[0]) in change_to_iterstat_list:
                print(line.split()[0]+" is in change_to_iterstat_list")
                for line1 in iter_cat:
                    if line1.split()[0]==line.split()[0]:
                        print("replacing "+line.split()[6]+" with "+line1.split()[1])
                        final_cat.write(line.replace(line.split()[6],line1.split()[1]))
            if int(line.split()[0]) in take_average_list:
                print(line.split()[0]+" is in take_average_list")
                for line1 in iter_cat:
                    if line1.split()[0]==line.split()[0]:
                        avg=np.mean([float(line.split()[6]),float(line1.split()[1])])
                        print("replacing "+line.split()[6]+" with "+str(avg))
                        final_cat.write(line.replace(line.split()[6],str(avg)))
            if int(line.split()[0]) in remove_entirely_list:
                print(line.split()[0]+" is in remove_entirely_list")
                pass
        final_cat.close()
    if make_sky_plots=="yes":
        image_data = fits.open(drz_image,memmap=True)[0].data
        cat = open(base_name+"_ALLsex.cat","r").readlines()
        cat_iter = open(base_name+"_itersky.cat","w")
        id = [int(line.split()[0]) for line in cat]
        sky_values = [float(line.split()[6]) for line in cat]
        x_values = [int(round(float(line.split()[9]))) for line in cat]
        y_values = [int(round(float(line.split()[10]))) for line in cat]
        k=0
        j=0
        pp = PdfPages("sky_results.pdf")
        i=0
        l=0
        fig, axs = plt.subplots(3, 3, sharex=True, tight_layout=False)
        f = plt.figure(1)
        plt.plot([1,2],[3,4])
        plt.close(f)
        #len(sky_values)
        while k<len(sky_values):
            if id[k]!=0:
                try:
                    x,y=x_values[k],y_values[k]
                    ### Don't like that this is hardcoded could cause errors at the edge of the image
                    drz_row_sum = np.sum(image_data[y-125:y+125,x-25:x+25],1)/(50) #sums row values e.g. a = np.array([[1,2],[3,4]]), np.sum(a,1) --> [3,7] not [4,6]
                    drz_col_sum = np.sum(image_data[y-25:y+25,x-125:x+125],0)/(50) #pancakes the array
                    #print(drz_col_sum)
                    #print(drz_row_sum)
                    axs[i,l].plot(np.arange(1,251),drz_row_sum,color="b",label="rows",linewidth=0.5)
                    axs[i,l].plot(np.arange(1,251),drz_col_sum,color="k",label="columns",linewidth=0.5)
                    axs[i,l].axhline(sky_values[k],color="g",label="SExtractor",linewidth=0.5)
                    
                    xmin,xmax,ymin,ymax=x-125,x+125,y-125,y+125
                    print(sky_values[k])
                    print(xmin,xmax,ymin,ymax)
                    print(id[k])
                    cutout_iter = drz_image+str("[")+str(xmin)+":"+str(xmax)+","+str(ymin)+":"+str(ymax)+"]"
                    initial_time=time.time()
                    iterstatus = float(iraf.iterstat(cutout_iter,maxiter=20,nsigrej=2.4,Stdout=1)[-1].split()[5].replace("mode=",""))
                    print("Time to run iterstat:"+str(time.time()-initial_time))
                    cat_iter.write(str(id[k])+" "+str(iterstatus)+"\n")
                    axs[i,l].axhline(iterstatus,color='r',label="iterstat",linewidth=0.5)
                    axs[i,l].set_ylim([iterstatus*0.8,iterstatus*1.2])
                    axs[i,l].set_yticks(np.round((np.linspace(sky_values[k]*0.8,sky_values[k]*1.2,10)),decimals=2))
                    axs[i,l].tick_params(axis="y",labelsize=6)
                    axs[i,l].minorticks_on()
                    axs[i,l].legend(loc="lower left",prop={'size': 4})
                    axs[i,l].annotate("id="+str(id[k]),xycoords="axes fraction", xy=(0.05,0.85))
                    i=i+1
                    
                    if i>2:
                        i=0
                        l=l+1
                        if l>2:
                            pp.savefig(fig)
                            #plt.close(fig)
                            l=0
                            i=0
                            j=j+1
                            fig, axs = plt.subplots(3, 3, sharex=True, tight_layout=False)
                    y=y-85
                    drz_row_sum = np.sum(image_data[y-25:y+25,x-125:x+125],0)/(50) #sums row values e.g. a = np.array([[1,2],[3,4]]), np.sum(a,1) --> [3,7] not [4,6]
                    y=y+2*85
                    drz_row1_sum = np.sum(image_data[y-25:y+25,x-125:x+125],0)/(50) #pancakes the array
                    #print(drz_row1_sum)
                    #print(drz_row_sum)
                    axs[i,l].plot(np.arange(1,251),drz_row_sum,color="b",label="rows below",linewidth=0.5)
                    axs[i,l].plot(np.arange(1,251),drz_row1_sum,color="k",label="rows above",linewidth=0.5)
                    axs[i,l].axhline(sky_values[k],color="g",label="SExtractor",linewidth=0.5)
                    #y=y-85
                    #xmin,xmax,ymin,ymax=x-125,x+125,y-125,y+125
                    #cutout_iter = drz_image+str("[")+str(xmin)+":"+str(xmax)+","+str(ymin)+":"+str(ymax)+"]"
                    axs[i,l].axhline(iterstatus,color='r',label="iterstat",linewidth=0.5)
                    axs[i,l].legend(loc="lower left",prop={'size': 4})
                    #print(iterstatus)
                    axs[i,l].set_ylim([sky_values[k]*0.8,sky_values[k]*1.2])
                    axs[i,l].set_yticks(np.round((np.linspace(sky_values[k]*0.8,sky_values[k]*1.2,10)),decimals=2))
                    axs[i,l].tick_params(axis="y",labelsize=6)
                    axs[i,l].minorticks_on()
                    axs[i,l].annotate("id="+str(id[k]),xycoords="axes fraction", xy=(0.05,0.85))
                    print(id[k],str(100*k/len(sky_values))+" percent done.")
                    i=i+1
                    if i>2:
                        i=0
                        l=l+1
                        if l>2:
                            pp.savefig(fig)
                            #plt.close(fig)
                            l=0
                            i=0
                            j=j+1
                            fig, axs = plt.subplots(3, 3, sharex=True, tight_layout=False)
                except:
                    print("hit except")
            k=k+1
        pp.savefig(fig)
        cat_iter.close()
        #plt.close(fig)
    

            
    
        pp.close()






if make_psf=='yes':
    ### Makes crude PSF, see manual
    print("Making PSF...")
    iraf.cd(two_D)
    try:
        os.mkdir(two_D+"PSF")
    except:
        pass
    try:
        iraf.copy(drz_image,"PSF")
    except:
        pass
    iraf.cd("PSF")
    ### Load iraf packages
    iraf.digiphot(_doprint=0)
    iraf.daophot(_doprint=0)
    
    ### Set variables for PSF
    iraf.daophot.verify=psf_verify
    iraf.daophot.verbose=psf_verbose
    iraf.photpars.apertures=psf_apertures
    iraf.photpars.zmag=psf_zmag 
    iraf.datapars.scale=psf_scale
    iraf.datapars.fwhmpsf=psf_fwhmpsf
    #iraf.datapars.ccdread=psf_ccdread
    #iraf.datapars.gain=psf_gain
    iraf.datapars.epadu=psf_epadu
    iraf.datapars.readnoise=psf_readnoise
    iraf.datapars.itime=psf_itime
    iraf.datapars.sigma=psf_sigma
    iraf.datapars.datamax=psf_datamax
    iraf.datapars.datamin=psf_datamin
    iraf.daopars.saturated=psf_saturated
    iraf.centerpars.calgorithm=psf_calgorithm
    iraf.centerpars.cbox=psf_cbox
    iraf.fitskypars.salgorithm=psf_salgorithm
    iraf.fitskypars.annulus=psf_annulus
    iraf.fitskypars.dannulus=psf_dannulus
    iraf.fitskypars.smooth=psf_smooth
    iraf.daopars.matchrad=psf_matchrad
    iraf.daopars.psfrad=psf_psfrad
    iraf.daopars.fitrad=psf_fitrad
    iraf.findpars.threshold=psf_threshold
    

    iraf.delete("*.coo.*")
    iraf.delete("*.mag.*")
    iraf.delete("*.pst.*")
    iraf.delete("*.psg.*")
    iraf.delete("*.psf.*")
    iraf.delete("*.sub.*")
    iraf.delete("*.nrj.*")
    iraf.delete("*.nst.*")
    iraf.delete("psf-"+str(psf_x_one)+"-"+str(psf_y_one)+".fits")
    iraf.delete("*.reg")
    drz_name=drz_image.split(".")[0]
    print("""This is loosely based off of the instructions found here: 
        http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?psf.hl 
        This could definitely be improved""")
    print("\nRunning daofind...")
    iraf.daofind(drz_name,"default")
    iraf.phot(drz_name,"default","default")
    print("\nRunning phot...")
    iraf.pstselect(drz_name,"default","default",maxnpsf=max_psf)
    print("\nRunning pstselect...")
    if os.path.exists(phot+base_name+"_tab_all.fits"):
        SEx_output_fits = phot+base_name+"_tab_all.fits"
        image = fits.open(SEx_output_fits, memmap=True)
        image_data = Table(image[1].data)
        class_star_list = list(image[1].data["CLASS_STAR"])
        mag_list = list(image[1].data["MAG_AUTO"])
        x_list,y_list = list(image[1].data["X_IMAGE"]),list(image[1].data["Y_IMAGE"])
        file0=open(drz_name+".pst.1","r").readlines()
        file1=open(drz_name+".pst.2","w")
        reg = open("reg.reg","w")
        ### This verifies with SExtractor CLASS_STAR classification
        for line in file0:
            if line.startswith("#"):
                file1.write(line)
            try:
                int(line.split()[0])
                x_pst,y_pst=float(line.split()[1]),float(line.split()[2])
                final_distance=100000
                #print(len(class_star_list),len(x_list))
                for i in range(0,len(x_list)):
                    distance=((x_pst-float(x_list[i]))**2+(y_pst-float(y_list[i]))**2)**0.5
                    if distance<final_distance:
                        final_distance=distance
                        final_class_star=float(class_star_list[i])
                        final_mag=float(mag_list[i])
                if final_class_star>star_thresh and final_mag<limiting_mag:
                    file1.write(line)
                    reg.write("circle "+str(line.split()[1])+" "+str(line.split()[2])+" "+str(15)+'\n')
            except:
                pass
        file1.close()
        reg.close()
        print("\nSelected PSF stars with CLASS_STAR > "+str(star_thresh)+" and MAG_AUTO < "+str(limiting_mag)+" and wrote them to .pst.2")
    print("\nRunning psf...")
    iraf.datapars.datamax=datamax_sat
    iraf.daopars.nclean=n_clean
    iraf.psf(drz_name,"default","default","default","default","default",interactive="no")
    iraf.nstar(drz_name,drz_name+".psg.1","default","default","default")
    iraf.substar(drz_name,"default","","default","default")
    iraf.seepsf(drz_name+".psf.1.fits","psf-"+str(psf_x_one)+"-"+str(psf_y_one),magnitude=mag_zp)
    print("\nDisplaying...")
    
    os.system("ds9 "+drz_image+" -regions load 'reg.reg' -fits "+drz_name+".sub.1.fits -regions load 'reg.reg' &")
    #iraf.tvmark(2,"STDIN",mark="circle",radii=15,color=207,Stdin=daofind)
    #pstselect = iraf.pdump(drz_name+".pst.2","XCENTER,YCENTER","yes",Stdout=1)
    #iraf.tvmark(1,"STDIN",mark="circle",radii=25,color=208,Stdin=pstselect)
    #x=(iraf.imstat(drz_name+"sub.1.fits",fields="min,max",format="no", Stdout=1))
    psf_name="psf-"+str(psf_x_one)+"-"+str(psf_y_one)+".fits"
    os.system("ds9 "+psf_name+" &")
    psf_data = fits.open(psf_name)[0].data
    psf_dim=psf_data.shape[0]
    plt.figure(100)
    plt.plot(np.arange(1,psf_dim+1),np.sum(psf_data,0)/psf_dim)
    axes = plt.gca()
    ### Hard coded limits but should be ok
    axes.set_ylim([-1*10**-5,1*10**-5]) #Might not work for everyone
    if psf_comparator!="default":
        psf_data = fits.open(psf_comparator)[0].data
        psf_dim=psf_data.shape[0]
        plt.plot(np.arange(1,psf_dim+1),np.sum(psf_data,0)/psf_dim)
    axes = plt.gca()
    axes.set_ylim([-1*10**-5,1*10**-5])
    plt.title("average zoomed in at base")
    plt.figure(101)
    psf_data = fits.open(psf_name)[0].data
    psf_dim=psf_data.shape[0]
    plt.plot(np.arange(1,psf_dim+1),np.sum(psf_data,0)/psf_dim)
    if psf_comparator!="default":
        psf_data = fits.open(psf_comparator)[0].data
        psf_dim=psf_data.shape[0]
        plt.plot(np.arange(1,psf_dim+1),np.sum(psf_data,0)/psf_dim)
    plt.title("average of rows")
    plt.figure(102)
    psf_data = fits.open(psf_name)[0].data
    psf_dim=psf_data.shape[0]
    
    ### I assume this is the central row, haven't had time to think if this is true always
    psf_mid_row=int(np.floor(psf_dim/2))
    plt.plot(np.arange(1,psf_dim+1),psf_data[:,psf_mid_row])
    plt.title("central row")
    if psf_comparator!="default":
        psf_data = fits.open(psf_comparator)[0].data
        psf_dim=psf_data.shape[0]
        plt.plot(np.arange(1,psf_dim+1),psf_data[:,psf_mid_row])
    #Wrote some of this kind of hastily, but it seems to work for me
    
# login.cl is here /Users/rileypeterson/.iraf/login.cl

#ascii.read(drz_name+".coo.1")["XCENTER","YCENTER"]





#GALFITTING
    
    
    
    
    


    

        
if make_constraint=="no":
    os.chdir(two_D)
    iraf.cd(two_D)
    try:
        iraf.delete("constraint_file.txt")
        print("deleted constraint_file.txt")
    except:
        pass


if make_constraint=="yes":
    
    ### This isn't so robust... very simple to add to though
    os.chdir(two_D)
    iraf.cd(two_D)
    try:
        iraf.delete("constraint_file.txt")
    except:
        pass
    file=open("constraint_file.txt","a")
    for i in range(2,16):
        if n_constraint!=None:
            file.write(str(i)+"    n    "+str(n_constraint[0])+" to "+str(n_constraint[1])+'\n')
        if re_constraint!=None:
            file.write(str(i)+"    re    "+str(re_constraint[0])+" to "+str(re_constraint[1])+'\n')
    print("made constraint_file.txt")
    file.close()














if do_run_galfit=='yes':
    os.chdir(two_D)
    iraf.cd(two_D)
    
    
    ###Get mask list from final_mask_for_galfit
    if final_mask_for_galfit=="default":
        final_mask_for_galfit=phot+base_name+"_mask.fits"
    mask_data=fits.open(final_mask_for_galfit)[0].data
    final_mask_list=list(set(mask_data[np.nonzero(mask_data)]))
    
    ###Copy final.cat file to the two_D folder
    if os.path.exists(base_name+"_ALLsex_final.cat")==False:
        iraf.copy(phot+base_name+"_ALLsex_final.cat",two_D)
        
    ###Establish names
    final_cat=base_name+"_ALLsex_final.cat"
    sersic1_longfeedme_txt=final_cat.replace(".cat","_longfeedme.txt")
    sersic1_longfeedme_dict=final_cat.replace(".cat","_longfeedme.npy") #Not full path
    sersic1_base_folder=two_D+base_name+"_sersic1/"
    sersic1_refit_folder=two_D+base_name+"_sersic1_refits/"
    cutout_folder=two_D+"galfit-input-cutouts/"
    sersic1_refit_summary=sersic1_refit_folder+base_name+"_sersic1_refits_summary.txt"
    sersic1_refit_longfeedme_with_changes=sersic1_refit_folder+base_name+"_sersic1_longfeedme_w_refits.txt"
    sersic1_refit_longfeedme_with_changes_dict=sersic1_refit_longfeedme_with_changes.replace(".txt",".npy")
    sersic1_final_folder=two_D+base_name+"_sersic1_final/"
    sersic1_final_outtab=two_D+base_name+"_sersic1_final/"+base_name+"_sersic1_output_table_final.txt"
    sersic1_longfeedme_txt_new=sersic1_base_folder+sersic1_longfeedme_txt.replace(".txt","_sersic1.txt").replace("_ALLsex_final","")
    sersic1_longfeedme_dict_new=sersic1_base_folder+sersic1_longfeedme_dict.replace(".npy","_sersic1.npy").replace("_ALLsex_final","")
    

    
    ###If the dictionary doesn't exist then we need to make it, the longfeedme dictionary is the most essential part of the operation
    if os.path.exists(sersic1_longfeedme_dict)==False:
        cat_dict=cat2longfeedme(final_cat,sersic1_longfeedme_txt,sersic1_longfeedme_dict,n_start=n_start1,min_cutout_radius=min_cutout_radius,final_mask_list=final_mask_list)
    
    ###Make cutouts using the dictionary created in the previous step
    if make_cutouts=="yes":
        do_make_cutouts(cutout_folder,sersic1_longfeedme_dict,drz_image,imaging+base_name+"_sigma.fits",final_mask_for_galfit,list_to_make=cutouts_list_to_make)   
    
    ###Make the individual feedme files for sersic1
    if make_feedmes_sersic1=="yes":
        longfeedme2feedmes(sersic1_longfeedme_dict,cutout_folder,sersic1_base_folder,list_to_make="default",psf_size=psf_size,fit_sky=fit_sky)
        ###If either the longfeedme txt or longfeedme dictionary haven't been added to the _sersic1 folder then make a copy, rename, and add them.
        ### This was eventually written as a definition copy_longfeedme_dictionary
        os.chdir(two_D)
        iraf.cd(two_D)
        if os.path.exists(sersic1_base_folder+sersic1_longfeedme_txt)==False:
            iraf.copy(sersic1_longfeedme_txt, sersic1_base_folder)
            os.rename(sersic1_base_folder+sersic1_longfeedme_txt, sersic1_longfeedme_txt_new)
        if os.path.exists(sersic1_base_folder+sersic1_longfeedme_dict)==False:
            iraf.copy(sersic1_longfeedme_dict, sersic1_base_folder)
            os.rename(sersic1_base_folder+sersic1_longfeedme_dict, sersic1_longfeedme_dict_new)
    
    ###Execute GALFIT on the feedmes in sersic1
    if run_sersic1=="yes":
        run_galfit(sersic1_base_folder,sersic1_base_folder,cutout_folder,fit="sersic1",input_longfeedme_dict=sersic1_longfeedme_dict,
                   list_to_make=run_sersic1_list,number_sent=number_sent,display_blocks=display_blocks_sersic1,timer=timer,verbose="yes",command=galfit_binary)    
    
    ###Evaluate which need to be refit, this is not really necessary if display_blocks="yes", as it pretty much does the same
    if check_residuals_sersic1=="yes":
        model_out=check_residuals(sersic1_base_folder,cutout_folder,list_to_check=check_residuals_sersic1_list)
        
    ###Plot statistics of the out.fits files in the sersic1 folder, this is a useful definition as it just takes the folder as the arg
    if plot_stats_sersic1=="yes":
        plot_stats(sersic1_base_folder)
        
    ###Refitting Sersic1   
    if run_refits_sersic1=="yes":
        if os.path.exists(sersic1_refit_summary)==False:
            try:
                os.mkdir(sersic1_refit_folder)
            except:
                pass
            file2=open(sersic1_refit_summary,"w")
            [file2.write(line) for line in preamble_for_changes]
            file2.close()
        make_refits(sersic1_base_folder,cutout_folder,sersic1_refit_folder,sersic1_refit_summary,list_to_refit=run_refits_sersic1_list)
    
    ###Enacts the changes on the longfeedme according to the summary of changes file
    if enact_refit_edits_sersic1=="yes":
        enact_summary_of_changes_to_longfeedme(sersic1_refit_summary,sersic1_longfeedme_dict_new,sersic1_refit_longfeedme_with_changes,changes_type=changes_type)
        longfeedme2longfeedme_dict(sersic1_refit_longfeedme_with_changes,sersic1_refit_longfeedme_with_changes_dict)
        
    ###Combines refits and puts the versions in the final folder
    if combine_refits_and_original_to_final_folder_sersic1_and_make_table=="yes":
        combine_refits(sersic1_refit_folder,sersic1_base_folder,sersic1_final_folder,sersic1_refit_longfeedme_with_changes_dict)
        #First make a dictionary then create the text for the dictionary
        output_dict_sersic1=make_output_dict(sersic1_final_folder,list_to_make="default")
        ###Assumes sersic1_longfeedme.npy is dictionary values from SExtractor else the last columns of the rows aren't the values from SExtractor obviously
        output_dict2tab(output_dict_sersic1,sersic1_longfeedme_dict_new,sersic1_longfeedme_dict_new,sersic1_longfeedme_dict_new,sersic1_final_outtab,list_to_make="default",plt_scale=plt_scale)
        make_fits(sersic1_final_outtab)
        if os.path.exists(sersic1_final_folder+base_name+"_sersic1_final_longfeedme.npy")==False:
            iraf.copy(sersic1_refit_longfeedme_with_changes_dict,sersic1_final_folder)
            os.rename(sersic1_final_folder+base_name+"_sersic1_longfeedme_w_refits.npy",sersic1_final_folder+base_name+"_sersic1_final_longfeedme.npy")
        if os.path.exists(sersic1_final_folder+base_name+"_sersic1_final_longfeedme.txt")==False:
            iraf.copy(sersic1_refit_longfeedme_with_changes,sersic1_final_folder)
            os.rename(sersic1_final_folder+base_name+"_sersic1_longfeedme_w_refits.txt",sersic1_final_folder+base_name+"_sersic1_final_longfeedme.txt")
        
        
    sersic1_final_longfeedme_dict=sersic1_final_folder+base_name+"_sersic1_final_longfeedme.npy"
    sersic1_final_longfeedme_txt=sersic1_final_folder+base_name+"_sersic1_final_longfeedme.txt"    
        
        
        
        
        
    

    #Devauc
    #Names
    devauc1_base_folder=two_D+base_name+"_devauc1/"
    devauc1_refit_folder=two_D+base_name+"_devauc1_refits/"
    devauc1_refit_summary=devauc1_refit_folder+base_name+"_devauc1_refits_summary.txt"
    devauc1_refit_longfeedme_with_changes=devauc1_refit_folder+base_name+"_devauc1_longfeedme_w_refits.txt"
    devauc1_refit_longfeedme_with_changes_dict=devauc1_refit_longfeedme_with_changes.replace(".txt",".npy")
    devauc1_final_folder=two_D+base_name+"_devauc1_final/"
    devauc1_final_outtab=two_D+base_name+"_devauc1_final/"+base_name+"_devauc1_output_table_final.txt"


    ###Currently these are the only options
    if use_output_as_input=="no" and method=="stop_and_go" and make_feedmes_devauc1=="yes":
        try:
            os.mkdir(devauc1_base_folder)
        except:
            pass
        devauc1_longfeedme_dict=copy_longfeedme_dictionary(sersic1_final_longfeedme_dict,devauc1_base_folder)
        devauc_dic=np.load(devauc1_longfeedme_dict).item()
        [devauc_dic[id].update({"type1":"devauc"}) for id in sorted(devauc_dic.keys())]
        [devauc_dic[id].update({"n1":"4.0"}) for id in sorted(devauc_dic.keys())]
        np.save(devauc1_longfeedme_dict,devauc_dic)
        devauc1_longfeedme_txt=longfeedme_dict2longfeedme(devauc1_longfeedme_dict,devauc1_longfeedme_dict.replace(".npy",".txt"))
    
    ###Makes devauc feedmes
    if make_feedmes_devauc1=="yes":
        longfeedme2feedmes(devauc1_longfeedme_dict,cutout_folder,devauc1_base_folder,list_to_make="default",psf_size=psf_size,fit_sky=fit_sky)
    
    ###Runs devauc feedmes
    if run_devauc1=="yes":
            run_galfit(devauc1_base_folder,devauc1_base_folder,cutout_folder,fit="devauc1",input_longfeedme_dict=devauc1_longfeedme_dict,
                       list_to_make="default",number_sent=number_sent,display_blocks=display_blocks_devauc1,timer=timer,verbose="yes",command=galfit_binary)    
    
    ###Checks residual of devauc
    if check_residuals_devauc1=="yes":
        model_out=check_residuals(devauc1_base_folder,cutout_folder,list_to_check=check_residuals_devauc1_list)
    
    ###Plots devauc stats
    if plot_stats_devauc1=="yes":
        plot_stats(devauc1_base_folder)

    ###Creates refit folder for devauc
    if run_refits_devauc1=="yes":
        if os.path.exists(devauc1_refit_summary)==False:
            try:
                os.mkdir(devauc1_refit_folder)
            except:
                pass
            file2=open(devauc1_refit_summary,"w")
            [file2.write(line) for line in preamble_for_changes]
            file2.close()
        make_refits(devauc1_base_folder,cutout_folder,devauc1_refit_folder,devauc1_refit_summary,list_to_refit=run_refits_devauc1_list)
    
    ###Enacts changes on longfeedme according to summary file
    if enact_refit_edits_devauc1=="yes":
        enact_summary_of_changes_to_longfeedme(devauc1_refit_summary,devauc1_longfeedme_dict,devauc1_refit_longfeedme_with_changes)
        longfeedme2longfeedme_dict(devauc1_refit_longfeedme_with_changes,devauc1_refit_longfeedme_with_changes_dict)
    
    ###Combines refits and originals in final folder, also makes table (txt and fits)
    if combine_refits_and_original_to_final_folder_devauc1_and_make_table=="yes":
        combine_refits(devauc1_refit_folder,devauc1_base_folder,devauc1_final_folder,devauc1_refit_longfeedme_with_changes_dict)
        #First make a dictionary then create the text for the dictionary
        output_dict_devauc1=make_output_dict(devauc1_final_folder,list_to_make="default")
        ###Assumes sersic1_longfeedme.npy is dictionary values from SExtractor else the last columns of the rows aren't the values from SExtractor obviously
        output_dict2tab(output_dict_devauc1,sersic1_longfeedme_dict_new,sersic1_longfeedme_dict_new,devauc1_longfeedme_dict,devauc1_final_outtab,list_to_make="default",plt_scale=plt_scale)
        make_fits(devauc1_final_outtab)
        devauc_final_dict=copy_longfeedme_dictionary(devauc1_refit_longfeedme_with_changes_dict,devauc1_final_folder)
        longfeedme_dict2longfeedme(devauc_final_dict,devauc_final_dict.replace(".npy",".txt"))
    
    devauc1_final_longfeedme_dict=devauc1_final_folder+base_name+"_devauc1_final_longfeedme.npy"
    devauc1_final_longfeedme_txt=devauc1_final_folder+base_name+"_devauc1_final_longfeedme.txt"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #Sersic2 
    sersic2_base_folder=two_D+base_name+"_sersic2/"
    sersic2_refit_folder=two_D+base_name+"_sersic2_refits/"
    sersic2_refit_summary=sersic2_refit_folder+base_name+"_sersic2_refits_summary.txt"
    sersic2_refit_longfeedme_with_changes=sersic2_refit_folder+base_name+"_sersic2_longfeedme_w_refits.txt"
    sersic2_refit_longfeedme_with_changes_dict=sersic2_refit_longfeedme_with_changes.replace(".txt",".npy")
    sersic2_final_folder=two_D+base_name+"_sersic2_final/"
    sersic2_final_outtab=two_D+base_name+"_sersic2_final/"+base_name+"_sersic2_output_table_final.txt"
    sersic2_longfeedme_dict=sersic2_base_folder+base_name+"_sersic2_longfeedme.npy"
    sersic2_longfeedme_txt=sersic2_longfeedme_dict.replace(".npy",".txt")
    sersic2_final_longfeedme_dict=sersic2_final_folder+base_name+"_sersic2_final_longfeedme.npy"
    sersic2_final_longfeedme_txt=sersic2_final_folder+base_name+"_sersic2_final_longfeedme.txt"
     
    
    
    
    ###First thing we need to do is figure out which are eligible for a sersic2 fit
    ###Get the necessary information from the output dictionaries
    
    
    
    
    ### Again this is currently the only mode
    if use_output_as_input=="no" and method=="stop_and_go" and make_feedmes_sersic2=="yes":
        try:
            os.mkdir(sersic2_base_folder)
        except:
            pass
        ###Load the dictionaries 
        sersic1_final_longfeedme_dict2=make_output_dict(sersic1_final_folder,list_to_make="default")
        devauc1_final_longfeedme_dict2=make_output_dict(devauc1_final_folder,list_to_make="default")
        ###Get entries which they both have, otherwise how can we compare the chi square and sersic values?
        sersic2_list_of_keys=sorted([i for i in sersic1_final_longfeedme_dict2.keys() if i in devauc1_final_longfeedme_dict2.keys()])
        sersic2_longfeedme_dic=np.load(sersic1_final_longfeedme_dict).item()
        for id in sersic2_longfeedme_dic.keys():
            if id not in sersic2_list_of_keys:
                del sersic2_longfeedme_dic[id]
        ###Select better chisqr in Devauc fit and final sersic1 n > 3
        for id in sersic2_list_of_keys:
            if float(sersic1_final_longfeedme_dict2[id]["chi2nu"])>float(devauc1_final_longfeedme_dict2[id]["chi2nu"]) and float(sersic1_final_longfeedme_dict2[id]["n1"])>3.0:
                print("id is:"+str(id), "sersic chi2nu is:"+str(float(sersic1_final_longfeedme_dict2[id]["chi2nu"])), "devauc chi2nu is:"+str(float(devauc1_final_longfeedme_dict2[id]["chi2nu"])),
                      "n from sersic fit is:"+str(float(sersic1_final_longfeedme_dict2[id]["n1"])))
                sersic2_longfeedme_dic[id].update({"n1":"4.0"})
            else:
                del sersic2_longfeedme_dic[id]
        np.save(sersic2_longfeedme_dict,sersic2_longfeedme_dic)
        longfeedme_dict2longfeedme(sersic2_longfeedme_dic,sersic2_longfeedme_txt)
    
    
    ###Make feedmes for eligible sersic2 values
    if make_feedmes_sersic2=="yes":
        longfeedme2feedmes(sersic2_longfeedme_dict,cutout_folder,sersic2_base_folder,list_to_make="default",psf_size=psf_size,fit_sky=fit_sky)
        
    ###Run sersic2
    if run_sersic2=="yes":
            run_galfit(sersic2_base_folder,sersic2_base_folder,cutout_folder,fit="devauc1",input_longfeedme_dict=sersic2_longfeedme_dict,
                       list_to_make="default",number_sent=number_sent,display_blocks=display_blocks_sersic2,timer=timer,verbose="yes",command=galfit_binary)    
    
    ###Check residuals, here you will want to compare to sersic1 probably (can use the output txt or fits table)
    if check_residuals_sersic2=="yes":
        model_out=check_residuals(sersic2_base_folder,cutout_folder,list_to_check=check_residuals_sersic2_list)
    
    ###Plot the stats of sersic2
    if plot_stats_sersic2=="yes":
        plot_stats(sersic2_base_folder)

    ###Select galaxies that have the possibility of being improved by the sersic2 fit
    if run_refits_sersic2=="yes":
        if os.path.exists(sersic2_refit_summary)==False:
            try:
                os.mkdir(sersic2_refit_folder)
            except:
                pass
            file2=open(sersic2_refit_summary,"w")
            [file2.write(line) for line in preamble_for_changes]
            file2.close()
        make_refits(sersic2_base_folder,cutout_folder,sersic2_refit_folder,sersic2_refit_summary,list_to_refit=run_refits_sersic2_list)
    
    ###Enact the changes to the longfeedme
    if enact_refit_edits_sersic2=="yes":
        enact_summary_of_changes_to_longfeedme(sersic2_refit_summary,sersic2_longfeedme_dict,sersic2_refit_longfeedme_with_changes)
        longfeedme2longfeedme_dict(sersic2_refit_longfeedme_with_changes,sersic2_refit_longfeedme_with_changes_dict)
    
    
    ###Combine with sersic1 based on use_sersic2_list, create the final sersic catalog
    if combine_refits_and_original_to_final_folder_sersic2_and_make_table=="yes":
        try:
            os.mkdir(sersic2_final_folder)
        except:
            pass
        ###First I need to merge the refit dictionary with the sersic1 final dictionary for the ids in use_sersic2_list
        sersic1_final_longfeedme_dic=np.load(sersic1_final_longfeedme_dict).item()
        sersic2_final_longfeedme_dic=sersic1_final_longfeedme_dic
        sersic2_refit_longfeedme_with_changes_dic=np.load(sersic2_refit_longfeedme_with_changes_dict).item()
        for id in use_sersic2_list:
            sersic2_final_longfeedme_dic[id]=sersic2_refit_longfeedme_with_changes_dic[id]
        for id in remove_from_final_catalog_sersic2_list:
            del sersic2_final_longfeedme_dic[id]
        np.save(sersic2_final_longfeedme_dict,sersic2_final_longfeedme_dic)
        longfeedme_dict2longfeedme(sersic2_final_longfeedme_dic,sersic2_final_longfeedme_txt)
        combine_refits(sersic2_refit_folder,sersic1_final_folder,sersic2_final_folder,sersic2_final_longfeedme_dict,overwrite="yes",make_changes_list=use_sersic2_list)
        output_dict_sersic2=make_output_dict(sersic2_final_folder,list_to_make="default")
        ###Assumes sersic1_longfeedme.npy is dictionary values from SExtractor else the last columns of the rows aren't the values from SExtractor obviously  <--- if you run the script as is this isn't a big deal cause they are the SExtractor values
        output_dict2tab(output_dict_sersic2,sersic1_longfeedme_dict_new,sersic1_longfeedme_dict_new,sersic2_final_longfeedme_dict,sersic2_final_outtab,list_to_make="default",plt_scale=plt_scale)
        make_fits(sersic2_final_outtab)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    

    

    
    
    
    
    
    
    
    
    
    
    
    
    