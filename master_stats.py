#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:16:05 2017

@author: rpeterson
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.table import Table
from astropy.io import fits, ascii
import pandas as pd

def make_stats_dict(path_to_input_dict,SEx_tab_all,path_to_output_fits_table=None,plt_scale=0.06):
    #Every galaxy in the simulated will have an entry
    #Every key will be prefixed with either sim_, sext_, or gal_
    #Inputs will be the input simulated dictionary, SExtractor tab_all.fits file, and optional GALFIT output

    sext = Table(fits.open(SEx_tab_all)[1].data)
    sim = np.load(path_to_input_dict).item()
    sim.pop("seed")
    orig_attribs=['y_cutout', 'n', 'q', 're', 'x_global', 'pa', 'mag', 'y_global', 'x_cutout']
    new_sim=dict()
    for id in sim.keys():
        new_sim[id]={"sim_id":id}
        for attrib in orig_attribs:
            new_sim[id].update({"sim_"+attrib:sim[id][attrib]})
    sim=new_sim
    master_stats=sim
    sim = pd.DataFrame.from_dict(sim)
    #First thing is to give every simulated galaxy its SExtractor id
    #Go through and find the closest SExtractor object
    sim_xs,sim_ys,sim_ids=list(sim.loc["sim_x_global"].astype(int)),list(sim.loc["sim_y_global"].astype(int)),list(sim.loc["sim_id"].astype(int))
    sext_xs,sext_ys,sext_ids=list(sext["X_IMAGE"]),list(sext["Y_IMAGE"]),list(sext["NUMBER"])
    sim_coords,sext_coords=zip(sim_xs,sim_ys,sim_ids),zip(sext_xs,sext_ys,sext_ids)
    sim_to_sext=dict()
    for sim_coord in sim_coords:
        min_distance=1000
        x,y,id=sim_coord[0],sim_coord[1],sim_coord[2]
        for sext_coord in sext_coords:
            x1,y1,id1=sext_coord[0],sext_coord[1],sext_coord[2]
            distance=((x1-x)**2+(y1-y)**2)**0.5
            if distance<min_distance:
                min_distance=distance
                sim_to_sext[id]={"sext_id":id1,"distance":min_distance} 
    SEx_to_sim=dict()
    for id in sim_to_sext.keys():
        sext_id,distance=sim_to_sext[id]["sext_id"],sim_to_sext[id]["distance"]
        for id1 in sim_to_sext.keys():
            sext_id1,distance1=sim_to_sext[id1]["sext_id"],sim_to_sext[id1]["distance"]
            if sext_id==sext_id1 and distance1<=distance:
                distance=distance1
                SEx_to_sim[sext_id1]=id1
                
            
    #return SEx_to_sim
    [master_stats[SEx_to_sim[i]].update({"sext_id":i}) for i in SEx_to_sim.keys()]            
    [master_stats[SEx_to_sim[obj[0]]].update({"sext_mag":obj[1]}) for obj in zip(list(sext["NUMBER"]),list(sext["MAG_AUTO"])) if obj[0] in SEx_to_sim.keys()] #This is digusting, but it works, it will work for any SExtractor data value
    [master_stats[SEx_to_sim[obj[0]]].update({"sext_mag_err":obj[1]}) for obj in zip(list(sext["NUMBER"]),list(sext["MAGERR_AUTO"])) if obj[0] in SEx_to_sim.keys()]
    [master_stats[SEx_to_sim[obj[0]]].update({"sext_re":obj[1]}) for obj in zip(list(sext["NUMBER"]),list(sext["FLUX_RADIUS"])) if obj[0] in SEx_to_sim.keys()]
    [master_stats[SEx_to_sim[obj[0]]].update({"sext_x_global":obj[1]}) for obj in zip(list(sext["NUMBER"]),list(sext["X_IMAGE"])) if obj[0] in SEx_to_sim.keys()] 
    [master_stats[SEx_to_sim[obj[0]]].update({"sext_y_global":obj[1]}) for obj in zip(list(sext["NUMBER"]),list(sext["Y_IMAGE"])) if obj[0] in SEx_to_sim.keys()]
    [master_stats[SEx_to_sim[obj[0]]].update({"sext_pa":obj[1]-90}) for obj in zip(list(sext["NUMBER"]),list(sext["THETA_IMAGE"])) if obj[0] in SEx_to_sim.keys()]
    [master_stats[SEx_to_sim[obj[0]]].update({"sext_q":1-obj[1]}) for obj in zip(list(sext["NUMBER"]),list(sext["ELLIPTICITY"])) if obj[0] in SEx_to_sim.keys()] 
    #That's pretty much all the useful parameters I can think of
    #Move on the GALFIT portion
    if path_to_output_fits_table is not None:
        gal = Table(fits.open(path_to_output_fits_table)[1].data)
        params=['MTOT_', 'E_MTOT_', 'AE_', 'E_AE_', 'RE_', 'E_RE_', 'N_', 'E_N_', 'ARATIO_', 'E_ARATIO_',
                'PA_', 'E_PA_', 'X_', 'E_X_', 'Y_', 'E_Y_', 'SKY_', 'CHISQ_']
        try:
            gal["MTOT_SER"]
            fit_type="SER"
        except:
            fit_type="DEV"
        #wish I could just use wildcard *, but I don't think its possible
        params=[i+fit_type for i in params]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_id":obj[0]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[0]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_mag":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[0]]))] 
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_mag_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[1]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_ae":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[2]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_ae_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[3]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_recirc":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[4]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_recirc_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[5]]))]        
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_n":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[6]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_n_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[7]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_q":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[8]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_q_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[9]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_pa":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[10]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_pa_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[11]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_x_global":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[12]]))] 
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_x_global_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[13]]))] 
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_y_global":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[14]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_y_global_err":obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[15]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_sky":obj[1]-90}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[16]]))]
        [master_stats[SEx_to_sim[obj[0]]].update({"gal_chisqr":1-obj[1]}) for obj in zip(list(gal["NUMBERI"]),list(gal[params[17]]))]
        
    #get true intergalactic distance
    for id in master_stats.keys():
        x,y=master_stats[id]["sim_x_global"],master_stats[id]["sim_y_global"]
        min_distance=1000
        for id1 in master_stats.keys():
            x1,y1=master_stats[id1]["sim_x_global"],master_stats[id1]["sim_y_global"]
            distance=((x1-x)**2+(y1-y)**2)**0.5
            if 0<distance<min_distance:
                min_distance=distance
                master_stats[id].update({"sim_distance":min_distance})
                
    #get SExtractor perceived distance
    for id in master_stats.keys():
        try:
            x,y=master_stats[id]["sext_x_global"],master_stats[id]["sext_y_global"]
        except:
            pass
        min_distance=1000
        for id1 in master_stats.keys():
            try:
                x1,y1=master_stats[id1]["sext_x_global"],master_stats[id1]["sext_y_global"]
            except:
                pass
            distance=((x1-x)**2+(y1-y)**2)**0.5
            if 0<distance<min_distance:
                min_distance=distance
                master_stats[id].update({"sext_distance":min_distance})
                
    #logre for sims
    [master_stats[id].update({"sim_logre":np.log10(float(master_stats[id]["sim_re"])*np.sqrt(float(master_stats[id]["sim_q"]))*plt_scale)}) for id in master_stats.keys()]
    #logre for GALFIT
    if path_to_output_fits_table is not None:
        GALFIT_ids=[SEx_to_sim[i] for i in list(gal["NUMBERI"])]
        
        [master_stats[id].update({"gal_logre":np.log10(float(master_stats[id]["gal_ae"])*np.sqrt(float(master_stats[id]["gal_q"]))*plt_scale)}) for id in GALFIT_ids]
       
        [master_stats[id].update({"gal_logre_err":np.sqrt(float(master_stats[id]["gal_q"])*(float(master_stats[id]["gal_q_err"])*plt_scale)**2+float(master_stats[id]["gal_q_err"])**2*(float(master_stats[id]["gal_ae"])*plt_scale)**2/(4*float(master_stats[id]["gal_q"])))}) for id in GALFIT_ids]
    
    """
    aeff,aeff_err=out[id]["re1"],out[id]["re_err1"]
    ar,ar_err=out[id]["q1"],out[id]["q_err1"]
    plt_scale=0.06
    reff,reff_err=float(aeff)*np.sqrt(float(ar))*plt_scale,np.sqrt(float(ar)*(float(aeff_err)*plt_scale)**2+float(ar_err)**2*(float(aeff)*plt_scale)**2/(4*float(ar)))
    logre,logre_err=np.log10(float(reff)),(float(reff_err)/float(reff))*(1./np.log(10.))
    aeff_sim=float(sim[sim_id]["re"])
    ar_sim=float(sim[sim_id]["q"])
    reff_sim=float(aeff_sim)*np.sqrt(float(ar_sim))*plt_scale
    logre_sim=np.log10(float(reff_sim))
    x,y,yerr=logre_sim,(logre-logre_sim),logre_err
    print(x,y,yerr)
    plt.errorbar(x,y,yerr=yerr,c="b",marker="x",ms=4)
    """
    
    return pd.DataFrame.from_dict(master_stats),SEx_to_sim,master_stats
def get_misses(sext_df,full_df):
    parents_with_misses=[]
    df_length=len(full_df.T)
    print(df_length)
    for id in range(1,df_length/2+1):
        try:
            gal_or_sext_df[id],sext_df[id+df_length/2]
        except:
            try:
                sext_df[id]
                print(id,"missing")
                parents_with_misses.append(id)
            except:
                try:
                    sext_df[id+df_length/2]
                    print(id+df_length/2,"missing")
                    parents_with_misses.append(id+df_length/2)
                except:
                    pass
    df_of_misses=sext_df.filter(np.unique(parents_with_misses),axis=1)
    return df_of_misses


def plots_sim_stats(dataframe_from_make_stats_dict,plot_sext="yes",plot_gal="yes",plot_sext_misses="yes",plot_gal_misses="yes"):
    df=dataframe_from_make_stats_dict
    plt.close("all")
    fig, ax = plt.subplots(3, 3,squeeze=False)
    ul,u,ur,l,c,r,bl,b,br=[ax[i,j] for i in range(0,3) for j in range(0,3)]
    fig.subplots_adjust(hspace=.25)
    possible_patches=[]
    
    #Get the dataframes from objects which have GALFIT and SExtractor values exclusively
    if plot_sext=="yes": df_sext=df.T.dropna(subset=["sext_id"],how="any").T
    if plot_gal=="yes": df_galfit=df.T.dropna().T
    #if plot_sext_misses=="yes": df_sext_misses=get_misses(df_sext,df)
    #if plot_gal_misses=="yes": df_galfit_misses=df_sext_misses.filter(df_galfit.columns,axis=1)

    

    #Descriptive stats, upper left
    ul.set_title("Some stats about the simulation");ul.axis([0, 10, 0, 10])
    where=ul
    sim_num=len(df.columns)
    where.text(1,9,"# of simulated objects: "+str(sim_num))
    sim_mag_max,sim_mag_min=round(max(df.loc["sim_mag"].dropna().astype(float)),3),round(min(df.loc["sim_mag"].dropna().astype(float)),3)
    where.text(1,6,"simulated faintest,brightest: "+str(sim_mag_max)+", "+str(sim_mag_min))
    blue_patch=mpatches.Patch(color='blue', label='Simulated data is in blue')
    possible_patches.append(blue_patch)
    if plot_sext=="yes":
        sext_num=len(df.loc["sext_id"].dropna())
        ul.text(1,8,"# of reliable SExtracted detections: "+str(sext_num))
        sext_mag_max,sext_mag_min=round(max(df.loc["sext_mag"].dropna()),3),round(min(df.loc["sext_mag"].dropna()),3)
        where.text(1,5,"SExtractor faintest,brightest: "+str(sext_mag_max)+", "+str(sext_mag_min))
        yellow_patch=mpatches.Patch(color='yellow', label='SExtractor data is in yellow')
        possible_patches.append(yellow_patch)
    if plot_gal=="yes":
        gal_num=len(df.loc["gal_id"].dropna())
        where.text(1,7,"# of GALFITted objects: "+str(gal_num))
        gal_mag_max,gal_mag_min=round(max(df.loc["gal_mag"].dropna()),3),round(min(df.loc["gal_mag"].dropna()),3)
        where.text(1,4,"GALFIT faintest,brightest: "+str(gal_mag_max)+", "+str(gal_mag_min))
        red_patch=mpatches.Patch(color='red', label='GALFIT data is in red')
        possible_patches.append(red_patch)
    where.legend(handles=possible_patches,loc=8,fontsize=14)
    
    #Magnitude histograms
    where=u
    where.set_title("Magnitude Distribution");where.set_xlabel("Magnitude");where.set_ylabel("Number in bin")
    range_hist=(sim_mag_min,sim_mag_max)
    where.hist(df.loc["sim_mag"].dropna().astype(float),bins=40,alpha=0.8,color="b",normed=False,range=range_hist)
    if plot_sext=="yes": where.hist(df.loc["sext_mag"].dropna(),bins=40,alpha=0.8,color="y",normed=False,range=range_hist)
    if plot_gal=="yes": where.hist(df.loc["gal_mag"].dropna(),bins=40,alpha=0.8,color="r",normed=False,range=range_hist)
    
    #Magnitude errors against true magnitude
    where=ur
    where.set_title("Output Magnitude against True Magnitude");where.set_xlabel("Input Magnitude");where.set_ylabel("Output Magnitude - Input Magnitude")
    if plot_sext=="yes": where.errorbar(df_sext.T.sim_mag.astype(float),df_sext.T.sext_mag-df_sext.T.sim_mag.astype(float),yerr=df_sext.T.sext_mag_err,fmt="s",c="y",marker="x",ms=4)
    #if plot_sext_misses=="yes": where.errorbar(df_sext_misses.T.sim_mag.astype(float),df_sext_misses.T.sext_mag-df_sext_misses.T.sim_mag.astype(float),yerr=df_sext_misses.T.sext_mag_err,fmt="s",c="k",marker="x",ms=4)
    if plot_gal=="yes": where.errorbar(df_galfit.T.sim_mag.astype(float),df_galfit.T.gal_mag-df_galfit.T.sim_mag.astype(float),yerr=df_galfit.T.gal_mag_err,fmt="s",c="r",marker="x",ms=4)
    #if plot_gal=="yes": where.errorbar(df_galfit_misses.T.sim_mag.astype(float),df_galfit_misses.T.gal_mag-df_galfit_misses.T.sim_mag.astype(float),yerr=df_galfit_misses.T.gal_mag_err,fmt="s",c="k",marker="x",ms=4)

    #Magnitude errors against angular distance (True)
    where=l
    where.set_title("Output Magnitude against True Angular Distance");where.set_xlabel("Angular Distance (True) [pixels]");where.set_ylabel("Output Magnitude - Input Magnitude")
    if plot_sext=="yes": where.errorbar(df_sext.T.sim_distance.astype(float),df_sext.T.sext_mag-df_sext.T.sim_mag.astype(float),yerr=df_sext.T.sext_mag_err,fmt="s",c="y",marker="x",ms=4)
    #if plot_sext_misses=="yes": where.errorbar(df_sext_misses.T.sim_distance.astype(float),df_sext_misses.T.sext_mag-df_sext_misses.T.sim_mag.astype(float),yerr=df_sext_misses.T.sext_mag_err,fmt="s",c="k",marker="x",ms=4)
    if plot_gal=="yes": where.errorbar(df_galfit.T.sim_distance.astype(float),df_galfit.T.gal_mag-df_galfit.T.sim_mag.astype(float),yerr=df_galfit.T.gal_mag_err,fmt="s",c="r",marker="x",ms=4)
    #if plot_gal=="yes": where.errorbar(df_galfit_misses.T.sim_distance.astype(float),df_galfit_misses.T.gal_mag-df_galfit_misses.T.sim_mag.astype(float),yerr=df_galfit_misses.T.gal_mag_err,fmt="s",c="k",marker="x",ms=4)

    #Same but as a function of SExtractor perceived intergalactic distance
    where=c
    where.set_title("Output Magnitude against SExtractor Angular Distance");where.set_xlabel("Angular Distance (SExtractor) [pixels]");where.set_ylabel("Output Magnitude - Input Magnitude")
    if plot_sext=="yes": where.errorbar(df_sext.T.sext_distance.astype(float),df_sext.T.sext_mag-df_sext.T.sim_mag.astype(float),yerr=df_sext.T.sext_mag_err,fmt="s",c="y",marker="x",ms=4)
    #if plot_sext_misses=="yes": where.errorbar(df_sext_misses.T.sext_distance.astype(float),df_sext_misses.T.sext_mag-df_sext_misses.T.sim_mag.astype(float),yerr=df_sext_misses.T.sext_mag_err,fmt="s",c="k",marker="x",ms=4)
    if plot_gal=="yes": where.errorbar(df_galfit.T.sext_distance.astype(float),df_galfit.T.gal_mag-df_galfit.T.sim_mag.astype(float),yerr=df_galfit.T.gal_mag_err,fmt="s",c="r",marker="x",ms=4)
    #if plot_gal=="yes": where.errorbar(df_galfit_misses.T.sext_distance.astype(float),df_galfit_misses.T.gal_mag-df_galfit_misses.T.sim_mag.astype(float),yerr=df_galfit_misses.T.gal_mag_err,fmt="s",c="k",marker="x",ms=4)

    
    #logre errors against true logre
    where=r
    if plot_gal=="yes": where.set_title("log(re) Error for GALFIT against True log(re)");where.set_xlabel("Input log(re) [arcseconds]");where.set_ylabel("GALFIT log(re) - Input log(re)")
    if plot_gal=="yes": where.errorbar(df_galfit.T.sim_logre.astype(float),df_galfit.T.gal_logre-df_galfit.T.sim_logre.astype(float),yerr=df_galfit.T.gal_logre_err,fmt="s",c="r",marker="x",ms=4)
    
    #logre errors against angular distance (True)
    where=bl
    if plot_gal=="yes": where.set_title("log(re) Error for GALFIT against Angular Distance");where.set_xlabel("Angular Distance (True) [pixels]");where.set_ylabel("GALFIT log(re) - Input log(re)")
    if plot_gal=="yes": where.errorbar(df_galfit.T.sim_distance.astype(float),df_galfit.T.gal_logre-df_galfit.T.sim_logre.astype(float),yerr=df_galfit.T.gal_logre_err,fmt="s",c="r",marker="x",ms=4)
    #if plot_gal=="yes": where.errorbar(df_galfit_misses.T.sim_distance.astype(float),df_galfit_misses.T.gal_logre-df_galfit_misses.T.sim_logre.astype(float),yerr=df_galfit_misses.T.gal_logre_err,fmt="s",c="k",marker="x",ms=4)

    #logre errors against angular distance (SExtractor)
    where=b
    if plot_gal=="yes": where.set_title("log(re) Error for GALFIT against SExtractor Intergalactic Distance");where.set_xlabel("Intergalactic Distance (SExtractor) [pixels]");where.set_ylabel("GALFIT log(re) - Input log(re)")
    if plot_gal=="yes": where.errorbar(df_galfit.T.sext_distance.astype(float),df_galfit.T.gal_logre-df_galfit.T.sim_logre.astype(float),yerr=df_galfit.T.gal_logre_err,fmt="s",c="r",marker="x",ms=4)
    
    
    #sersic index errors against true sersic index
    where=br
    if plot_gal=="yes": where.set_title("Sersic Index Error for GALFIT against True Sersic Index");where.set_xlabel("Input n");where.set_ylabel("GALFIT n - Input n")
    if plot_gal=="yes": where.errorbar(df_galfit.T.sim_n.astype(float),df_galfit.T.gal_n-df_galfit.T.sim_n.astype(float),yerr=df_galfit.T.gal_n_err,fmt="s",c="r",marker="x",ms=4)
    
    
    
    """
    #r.errorbar(df_galfit.T.sim_distance.astype(float),df_galfit.T.gal_logre-df_galfit.T.sim_logre.astype(float),yerr=df_galfit.T.gal_logre_err,fmt="s",c="r",marker="x",ms=4)
    #r.set_title("Magnitude Error for GALFIT");r.set_xlabel("Intergalactic Distance (SExtractor) [pixels]");r.set_ylabel("GALFIT Magnitude - Input Magnitude")
    fig, ax = plt.subplots(3, 3,squeeze=False)
    ul,u,ur,l,c,r,bl,b,br=[ax[i,j] for i in range(0,3) for j in range(0,3)]
    fig.subplots_adjust(hspace=.25)
    #Descriptive stats
    sim_num,sext_num,gal_num=len(df.columns),len(df.loc["sext_id"].dropna()),len(df.loc["gal_id"].dropna())
    ul.set_title("Some stats about the simulation");ul.text(1,9,"# of simulated objects: "+str(sim_num));ul.axis([0, 10, 0, 10])
    ul.text(1,8,"# of reliable SExtracted detections: "+str(sext_num));ul.text(1,7,"# of GALFITted objects: "+str(gal_num))
    sim_mag_max,sim_mag_min=round(max(df.loc["sim_mag"].dropna().astype(float)),3),round(min(df.loc["sim_mag"].dropna().astype(float)),3)
    ul.text(1,6,"simulated faintest,brightest: "+str(sim_mag_max)+", "+str(sim_mag_min))
    sext_mag_max,sext_mag_min=round(max(df.loc["sext_mag"].dropna()),3),round(min(df.loc["sext_mag"].dropna()),3)
    ul.text(1,5,"SExtractor faintest,brightest: "+str(sext_mag_max)+", "+str(sext_mag_min))
    gal_mag_max,gal_mag_min=round(max(df.loc["gal_mag"].dropna()),3),round(min(df.loc["gal_mag"].dropna()),3)
    ul.text(1,4,"GALFIT faintest,brightest: "+str(gal_mag_max)+", "+str(gal_mag_min))
    blue_patch,red_patch,yellow_patch = mpatches.Patch(color='blue', label='Simulated data is in blue'),mpatches.Patch(color='red', label='GALFIT data is in red'),mpatches.Patch(color='yellow', label='SExtractor data is in yellow')
    ul.legend(handles=[blue_patch,yellow_patch,red_patch],loc=8,fontsize=14)
    
    where=u
    where.scatter(df_galfit.T.sim_mag.astype(float),np.multiply(df_galfit.T.sim_re.astype(float),df_galfit.T.sim_q.astype(float)**0.5),c="b",marker="x",s=16)
    where.errorbar(df_galfit.T.gal_mag.astype(float),np.multiply(df_galfit.T.gal_ae.astype(float),df_galfit.T.gal_q.astype(float)**0.5),yerr=df_galfit.T.gal_ae_err,fmt="s",c="r",marker="x",ms=4)
    where.set_title("Circularized Effective Radius vs. Magnitude");where.set_xlabel("Magnitude");where.set_ylabel("Circularized Effective Radius")

    #u.hist(bins=20, weights=np.ones_like(df[df.columns[0]]) * 100. / len(df))
    #input_mags=[float(input[id]["mag"]) for id in input.keys() if id!="seed"]
    """
    return df

    #plt.scatter(input_mags,output["MAGin_SER"],color="b")
    #plt.scatter(input_mags,output["MTOT_SER"],color="r")

    
a,b,c=make_stats_dict("/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_python_scripts/most_current/examples/example_1/simulations/inputs_dict.npy",
                "/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_python_scripts/most_current/examples/example_1/do-everything/do-everything_phot/do-everything_tab_all.fits",
                path_to_output_fits_table="/Users/rpeterson/FP/GALFIT_Simulations1/galfitsim1_python_scripts/most_current/examples/example_1/do-everything/do-everything_2D/do-everything_sersic1_final/do-everything_sersic1_output_table_final.fits")
df=plots_sim_stats(a,plot_sext="yes")
#a.loc["mag4"].dropna().astype(float).hist(bins=40)

def plot_gal_stats(final_output_table):
    data=Table(fits.open(final_output_table)[1].data)
    df=data.to_pandas()
    fig,axs=plt.subplots(2,2,squeeze=False)
    axs[0,0].set_ylabel("Number",fontsize=14),axs[1,0].set_ylabel("Number",fontsize=14)
    attribs=["MTOT_SER","LRE_SER","MUE_SER","N_SER"]
    names=["Absolute Magnitude Distribution","Log($r_e$) Distribution","$<\mu>_e$ Distribution","S"+r"$\acute{e}$"+"rsic Index Distribution"]
    xlabels=["$mag_{AB}$","log($r_e$) [arcsec]","$<\mu>_e$","n"]
    for obj in zip(attribs,names,axs.flat,xlabels):
        df[obj[0]].hist(normed=False,grid=False,bins=40,ax=obj[2])
        obj[2].set_title(obj[1],fontsize=14)
        obj[2].set_xlabel(obj[3],fontsize=14)
        obj[2].set_ylim(0,30)
    plt.subplots_adjust(hspace = 0.25)
        
    #ax=df.MTOT_SER.hist(normed=True,grid=False,bins=40,xlabel="hi")
    #ax.set_title("MAGNITUDE DISTRIBUTION")
    
    return data