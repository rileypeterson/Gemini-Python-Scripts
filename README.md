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

Sunday/Monday April 22nd and 23rd, 2018 :  
- Refactored the first interpolation. Slowly but surely...
