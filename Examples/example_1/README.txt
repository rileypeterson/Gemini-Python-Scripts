This is example_1
I’m going to make 800 simulated galaxies, they will be pairs, they will be astrodrizzled.

If you run make_sims.py as is it will make the exact same thing I did.

Not for do-everything with regard to the mask limit I got:

min(magbright+6.2,magfaint+0.5)=23.4791652679

This doesn’t make sense to implement because I have assoc_mag_cut = 24.2

So I want to at least be fitting 24.2 and less

The GCP guide says ideally we want to fit ~1mag below the faintest sample member, this isn’t reasonable either. I’m going to set my mask limit to 25.2.

For some reason two stars were placed into my image, I think this is due to the convolution. 

do-everything_525_out.fits just finished

Elapsed time is: 01:04:44

99 percent done with sersic1

Number that have crashed: 12

Done.
Crashed list should be ['obj12.in', 'obj130.in', 'obj148.in', 'obj183.in', 'obj210.in', 'obj3.in', 'obj392.in', 'obj410.in', 'obj42.in', 'obj487.in', 'obj6.in', 'obj60.in', 'obj94.in']

I have had a sudden change in thinking, I think it would be better if I increased the deb lending threshold and such for SExtractor, therefore more near neighbors will be recovered. Because looking at the residuals it is clear that most of the problems are because SExtractor didn’t detect the neighbor.

Actually, I don’t have time to go back and do this. Instead I will just add entries to the longfeedme

I returned the values to their original 16 and 0.01 and ran it so its back to the original

Refit list
[2, 5, 30, 47, 62, 79, 88, 91, 113, 117, 145, 154, 155, 160, 162, 166, 175, 182, 186, 229, 235, 252, 263, 265, 271, 287, 289, 298, 367, 371, 372, 386, 399, 406, 412, 426, 440, 456, 474, 475, 481, 515, 524, 12, 130, 148, 183, 210, 3, 392, 410, 42, 487, 6, 60, 94]

In the real data it isn’t always clear if there are two galaxies right next to each other or if the galaxy is two component. Since for the real data, I did one component fits, I did the same for the simulated. So if pairs looked like they were two components then I fit them with a single component fit (as if I don’t actually know they are separate galaxies). However, if it was clear that they were separate galaxies then I added another entry to fit the undetected neighbor. Hope that makes sense. Very rushed right now…


