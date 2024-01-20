## What is this?
    - I made the code for my bachelor thesis and won't update it unless I have very good reasons to do so
    - this code is able to analyse galaxy-data from the AMUSED website using the mpdaf package
    - the main goal was to get the radial velocities
    - i used the galaxy-data from here https://amused.univ-lyon1.fr/project/UDF/HUDF/browse
    - probably most of the well working galaxies are listed in "Galaxiensammlung"

## Feel free to write me if you need help or have questions
    - you should be able to contact me under etienne.beier@uni-potsdam.de
    - if I left the university by then you can probably still contact me by that mail or by etienne.beier@alumni.uni-potsdam.de

## How to use

### enter the file name in the definition of "glxy1"
    - radialmap1.3.py will work by default with 000855
    - other galaxies will need the following adjustments
      
### per default the programm will fit to the H-Alpha line
    - in the definition of "obs", "ha" can be switched to "hb" or "o3"
    - those are short for H-Beta and OIII emission lines
    - if a galaxy fails fit ha, try one of the other and it might work better
 
### the program has a simple filter to reduce noise: "flux_filter"
    - this will eliminate all light below given value
    - depending on the galaxy the value needs to be adjusted
 
### the programm has a filter to clean up velocity curves from to high values with the "arz"-factor
    - this is the biggest possible velocity shift between two neigbour pixels in km/s
    - for the most galaxies it should be smaller than te given 100km/s

### in ver rare cases the fit for the velocity curves has to be adjusted
    - the parameters can be found at pcov and p0 definition
    - they are normed to the values of 000855
    - the best way to adjust them is try and error
