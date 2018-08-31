# Gemini-Python-Scripts
Python scripts developed while working as an intern at the Gemini Observatory in Hilo, Hawaii.  
Manual.pdf - Users manual for Software  
GALFIT Wrapper - do_everything_v0.py  
Simulation Software -  make_sims.py  
Plotting and Analysis of Simulations - master_stats.py  

Works in progress:  
-Provide a full example using publicly available data  
-Resolve discrepancies with the inverse variance map and the current method for creating the sigma image  
-Consider eliminating the dependency of pyraf in favor of astropy modules  
-Consider simplifying/revising/automating the refitting process  

#### Status Updates/History:

Sunday March 25, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/648af5ecbbb6f1e69ad9b20906a22330b7f2d2fb#diff-66d1da86b9755f7c032b1c07934308a6)  
- Work on [Issue #1](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/1) on branch [update_sims](https://github.com/rileypeterson/Gemini-Python-Scripts/tree/update_sims)  

Sunday March 26, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/5e67ccca59d264e82f6cfa54f6f8817d60be6965#diff-66d1da86b9755f7c032b1c07934308a6)
- Continue working on [Issue #1](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/1) on branch [update_sims](https://github.com/rileypeterson/Gemini-Python-Scripts/tree/update_sims)
- Next step is to fill out comments, rename to full_config, and adjust make_sims.py to accomodate the config parser
*NTS* use fallback option:  
```
print(cfg.get('Section1', 'monster', fallback='No such things as monsters.'))
    # -> "No such things as monsters."
```
Sunday April 1st, 2018 : 
Still working on issue #1. Unfortunately fallback is only available in 3.2+ :(

Monday April 2nd, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/92e17cbee351d4e3da148e49d2f1e080a8bffa10)
- Continue working on [Issue #1](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/1) on branch [update_sims](https://github.com/rileypeterson/Gemini-Python-Scripts/tree/update_sims)
- Since they're closely related I have also started working on [Issue #2](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/2)

Friday/Saturday April 6th/7th, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/911b8bd065346d8a95d52edc73e6031d955c914b)
- Continue working on [Issue #1](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/1) on branch [update_sims](https://github.com/rileypeterson/Gemini-Python-Scripts/tree/update_sims) and [Issue #2](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/2)  
- I'm about to begin refactoring make_feedmes
- The configuration file is pretty much set up, it just needs to be connected to the variables within the script
- Started some basic error catching/logging, but also want to eventually create a test file

Wednesday April 11th, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/58aa0d2225ba23fe610aeb59170cfe778a47e5ce)
- Continue working on [Issue #1](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/1) on branch [update_sims](https://github.com/rileypeterson/Gemini-Python-Scripts/tree/update_sims) and [Issue #2](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/2)  
- Created some additional logging features and realize I need to give some thought about how to organize this

Monday April 16th, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/e22f9ac2c946831e6814098534f6b7733e826e4c)
- Continue working on the same issues. Next will be to refactor make_feedmes.

Sunday/Monday April 22nd/23rd, 2018 :  
- Refactored the first interpolation. Slowly but surely...

Saturday/Sunday June 2nd/3rd, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/50e3d7d5b4cb15c78c17ab39b0d4297ff72ca08c)
- Travelled to Norway recently and have been busy with work
- Finished the second interpolation, refactored code so it *could* be imported (i.e. shove things into if __name__ == "__main__":), trying to keep the logger and config out of functions, pretty sure I'll have the interpolating/creating feedme files done by my next commit (or at least that's a goal to strive for), decided it probably makes a lot of sense to establish a class for galaxies (galaxies are objects why not use OOP): 
```
class Galaxy(object):
    def __init__(self, x, y, mag, re, n, q, pa):
        # ID number associated with
        # Maybe some like file path (to feedme, fits, etc.)
        # Cutout dimensions
        # Galaxy position in cutout (x,y). Position in larger image.
        # Could write methods to generate the sersic profile or add the cutout to a larger image
        self.mag = mag
        self.x = x
        self.y = y
        self.re = re
        self.n = n
        self.q = q
        self.pa = pa
        ...
can use vars(Galaxy(...)) to capture these attributes and save them off into dictionaries or fits tables (probably)
```
Monday June 4th, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/b907c65e64699a1316522a8c8b88a8b22970f39d)
- I think the galaxy class was a great idea, I wrote a function which generates galaxies objects (uses interpolated values if supplied, otherwise random between the range).  
- I think it'd be nice to not require the config file, so I'm trying to write functions that make it optional. It'd be nice to just run `python make_sims.py` and have it generate a small field of ~10 galaxies. I'll go back and do that at the end.
- Microsoft acquired github today, and it makes me nervous.

Tuesday June 12th, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/1a37b63d95fffefddf21a77dfc4a8db4e445372d)
- Got the feedme/input files being created. Logging working as expected. Configuration file working as expected. Putting issue [#10](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/10) on hold for now. Plots being generated as expected. Next is to run the feedmes (easy) and the convolutions.

Sunday July 22nd, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/15a43facf98bd4fd73baa38958229512f23ac679)
- Finally getting back to this project... work has kept me busy and there's a lot going on this summer. Created the actual galaxy fits models and also have the convolution working. Added a TODO for this [line](https://github.com/rileypeterson/Gemini-Python-Scripts/blob/update_sims/make_sims.py#L1204) since I'm still not exactly clear on the difference between "fill" and "wrap" for the [boundary](http://docs.astropy.org/en/stable/api/astropy.convolution.convolve_fft.html) arg... I recall "wrap" is faster. 
- Next will be creating a mknoise function like iraf's and applying it to individual cutouts. Will require some testing to make sure I get it right. I'm trying to [eliminate the need for iraf/pyraf](https://github.com/rileypeterson/Gemini-Python-Scripts/issues/9) since it is somewhat of a pain to install/setup and it is also archaic. If I can achieve the same result (of [mknoise](http://iraf.net/irafhelp.php?val=artdata.mknoise&help=Help+Page)) then I could also port the code to python 3 (I believe that's the only barrier).
        
Thursday August 30th, 2018 : [Commit](https://github.com/rileypeterson/Gemini-Python-Scripts/commit/87a66e912ca936cbc5849df5a0baa147ea6d2c96)
- Working on replacing `mknoise` found how iraf does this [here](https://github.com/iraf/iraf-v216/blob/9590f45760a4791f3305407fb51c87f1282b32be/noao/artdata/t_mknoise.x#L565).
