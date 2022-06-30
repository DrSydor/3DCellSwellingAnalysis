# 3D Cell Swelling Analysis
Analysis of 3D cell volumes during cell swelling assays. There are two main scripts here: one for the actual analysis of the timelapse images (cell_swelling_analysis.py) and a second to plot the data and create xz and yz cross sections of each of the cells included in the analysis (cell_swelling_plots.py).


## cell_swelling_analysis.py

This script automatically determines the 3D volume of cells in a timecourse
series. It uses the Napari viewer as an interface so the user can manually
select the cells for segmentation and analysis. The script then generates
.csv files containing the volumes of all the selected cells and an image of
the cell masks for each time point. These outputs can then be used directly in
a separate script that will generate a plot of the normalized data and images
of each of the selected cells.

When running the script, run it in two parts: The first part will bring up
Napari and allows for the selection of the cells. The second part will conduct
the analysis using the points selected in the first part.

This script was largely based on the 3D cell analysis script from
https://jni.github.io/i2k-skimage-napari/lectures/2_segmentation_and_regionprops.html

Images should be pre-processed in ImageJ as follows:
1. Import data directly into ImageJ using the Bio-Formats Importer
    (Plug-ins -> Bio-Formats -> Bio-Formats Importer). The color channels were
    separated upon import and only the channel that gives the full cytosolic
    volume was used for analysis.
2. Images converted to 8-bit
3. The images were photobleach corrected using histogram matching
    (Image -> Adjust -> Bleach Correction)
4. Convert to 32 bit float (Optional)
5. Run ROF Denoise with a theta = 50 (Optional)
6. Covert to 8 bit image
7. Save image as .tif


## cell_swelling_plot.py

This script takes the .csv files output from the cell_swelling_analysis.py
script and collects and normalizes the data to the pre-treatment (-1 sec)
timepoint. It then plots the average of the data, along with the individual
datapoints. It also creates images of the selected cells showing the maximum
instensity projection (xy) and single yz and xz slices.
