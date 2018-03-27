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

## Status Updates/History:

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
