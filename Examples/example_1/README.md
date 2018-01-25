01/24/2018  
Beginning to put tools in place to make an example. This example will focus on do_everything.py, but perhaps I will demonstrate the simulations as well.  

Steps I took to create this example:  
Download f160w drz and wht extensions from https://archive.stsci.edu/pub/hlsp/candels/uds/uds-tot/v1.0/  
These files are ~1.5GB. Github has a file size limit of 100MB so I'm going to crop a portion of these images so that I'll be able to push them to this folder.  

"""
The CANDELS overview/data reference papers:  
• Grogin et al., 2011 ApJS, CANDELS Survey and Observational Design, http://arxiv.org/abs/1105.3753  
• Koekemoer et al., 2011, ApJS, CANDELS HST Imaging Data Products, http://arxiv.org/abs/1105.3754  
"""  

Once I've cropped these images, I'm going to try and create a model PSF from stars in the image...hopefully there are some...  
Creating a good PSF is somewhat complicated. I may resort to using: http://www.stsci.edu/hst/wfc3/analysis/PSF
