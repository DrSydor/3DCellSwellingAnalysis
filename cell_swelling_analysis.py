'''
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

'''


# PART I- Manusl definition of cells via Napari

import napari
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from skimage import (io, morphology, measure, filters, segmentation)

cells = io.imread("202206221-01_WTHenlewDynasore.tif") # Change file name as needed
times = [-1, 0, 68, 136, 204, 271, 1239, 2139, 3039, 3939] # in seconds
img_scale = [0.5, 0.244, 0.244] # Image pixel dimensions

viewer = napari.Viewer()

timepoint = 6 # using a mid-to-late timepoint to select the cells
cells_test = cells[timepoint]
edges = filters.scharr(cells_test)

# Using a guassian filter instead of denoising works better. Lower sigma is best.
denoised = filters.gaussian(cells_test, sigma=2)

# Otsu method was previously determined to work the best.
mean_thresholded = denoised > filters.threshold_otsu(denoised)

width = 30

remove_holes = morphology.remove_small_holes(mean_thresholded, width ** 3)
remove_objects = morphology.remove_small_objects(remove_holes, width ** 3)

viewer.add_image(
    remove_objects,
    name='cleaned',
    scale=img_scale,
    opacity=0.3)

labels = measure.label(remove_objects)

viewer.add_labels(
    labels,
    scale=img_scale,
    opacity=0.5)

viewer.dims.ndisplay = 2
viewer.dims.set_point(0, 30 * img_scale[0])

points = viewer.add_points(
    name='interactive points',
    scale=img_scale,
    ndim=3,
    size=4,
    n_dimensional=True)

points.mode = 'add' # now, we annotate the centers of the nuclei in your image

# !!!!!!!!!!!!!!!!!!!!Run scrip up to this point first!!!!!!!!!!!!!!!!!!!!!!



# PART II: Cell 3D volume analysis

marker_locations = points.data # The above points will be used for all time points

for time in range(len(times)):

    timepoint = time # change timepoint as needed
    cells_test = cells[timepoint]

    edges = filters.scharr(cells_test)

    # Using a guassian filter instead of denoising works better. Lower sigma is best.
    denoised = filters.gaussian(cells_test, sigma=2)

    # Otsu method was previously determined to work the best.
    mean_thresholded = denoised > filters.threshold_otsu(denoised)

    width = 30

    remove_holes = morphology.remove_small_holes(mean_thresholded, width ** 3)

    remove_objects = morphology.remove_small_objects(remove_holes, width ** 3)

    markers = np.zeros(cells_test.shape, dtype=np.uint32)
    marker_indices = tuple(np.round(marker_locations).astype(int).T)
    markers[marker_indices] = np.arange(len(marker_locations)) + 1
    markers_big = morphology.dilation(markers, morphology.ball(5))

    segmented = segmentation.watershed(
        edges,
        markers_big,
        mask=remove_objects)

    segmented_padded = np.pad(
        segmented,
        ((1, 1), (0, 0), (0, 0)),
        mode='constant',
        constant_values=0)

    # The [1:-1] below is to remove the padding added above
    interior_labels = segmentation.clear_border(segmented_padded)[1:-1]

    regionprops = measure.regionprops(interior_labels,
                                      intensity_image=cells_test)

    info_table = pd.DataFrame(
        measure.regionprops_table(
            interior_labels,
            intensity_image=cells_test,
            properties=['label', 'slice', 'area', 'centroid'])
    )

    time = str(times[timepoint])
    info_table.to_csv(f'{time}.csv')

    segmented_maxproj = np.max(segmented, axis = 0)
    plt.imshow(segmented_maxproj)

    for row in range(len(info_table.index)):
        plt.text(info_table['centroid-2'].iloc[row],
                 info_table['centroid-1'].iloc[row],
                 str(info_table['label'].iloc[row]))

    plt.savefig(f'Cell masks- {time} sec.png')
    plt.close()
