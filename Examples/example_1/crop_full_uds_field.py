"""
01/26/2018
Made this script to crop the UDS WFC3 f160w field which I downloaded.
Arbitrarily selected at region x:19100:21100, y:4600:6600
"""
from astropy.io import fits
full_uds_drz_image=fits.open(r"C:\Users\Riley Peterson\Downloads\hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits",memmap=True)
full_uds_drz_data=full_uds_drz_image[0].data
crop_uds_drz_data=full_uds_drz_data[4600:6600,19100:21100]
fits.writeto("uds_example_drz.fits",crop_uds_drz_data,full_uds_drz_image[0].header,overwrite=True)
full_uds_wht_image=fits.open(r"C:\Users\Riley Peterson\Downloads\hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_wht.fits",memmap=True)
full_uds_wht_data=full_uds_wht_image[0].data
crop_uds_wht_data=full_uds_wht_data[4600:6600,19100:21100]
fits.writeto("uds_example_wht.fits",crop_uds_wht_data,full_uds_wht_image[0].header,overwrite=True)