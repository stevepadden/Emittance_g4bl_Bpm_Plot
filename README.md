# Emittance_g4bl_Bpm_Plot
Plots emittance from ASCII formatted G4Beamline (g4bl) beam position monitors within a folder.

Step 1: Generate all BPM's in ASCII format. (BLTrackfile)
Step 2: Open up python script and direct to that folder in the "dir" element
Step 3: Remove any existing BPM images from the folder (windows specific)

This will generate 2 images for each BPM, one for x x' and one for y y' - note this relies on the pariaxal approximation.

Some key variables are defined at the top of this code, ensure they make sense for your current setup.


