'''
This script takes the .csv files output from the cell_swelling_analysis.py
script and collects and normalizes the data to the pre-treatment (-1 sec)
timepoint. It then plots the average of the data, along with the individual
datapoints. It also creates images of the selected cells showing the maximum
instensity projection (xy) and single yz and xz slices.
'''

from matplotlib import pyplot as plt
import pandas as pd
import glob

from skimage import io
import numpy as np
from skimage import exposure
from scipy import ndimage
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

times = [-1, 0, 68, 136, 204, 271, 1239, 2139, 3039, 3939] # in seconds

path = r'D:\20220621- Cell swelling\WTwDynasore' # use your path
all_files = glob.glob("*.csv")

#Generates an initial dataframe using the rows from the initial (-1.csv) file
df = pd.read_csv('-1.csv', index_col=0, usecols= [1], header=0)

for time in times:
    filename = str(time) + ".csv"
    df2 = pd.read_csv(filename, index_col=0, usecols= [1,3], header=0)
    df[str(time)] = df2['area']

df_norm_og = pd.DataFrame()

for time in times:
    #col_name = str(time)
    df_norm_og[time] = df[str(time)]/df['-1']

# Bad cells; cells that were poorly segmented:
bad_cells = [2,5,6,16,18,20,21,22]

df_norm = df_norm_og.drop(bad_cells)

av_col = df_norm.mean(axis = 0).to_frame()
stdev_col = df_norm.std().tolist()

fig, ax = plt.subplots(figsize = (6,4))

ax.plot(av_col, 'bo-', label = 'Average')
ax.errorbar(av_col.index, av_col[0], yerr = stdev_col)

for n in range(len(df_norm.index)):
    n += 1
    if n in bad_cells:
        pass
    else:
        ax.scatter(times, df_norm.T[n], s = 5, c = 'black')

ax.set_ylabel('Normalized cell volume', fontsize = 14)
ax.set_xlabel('Time (sec)', fontsize = 14)

plt.title('WT Henle + Dynasore Cells', fontsize = 16)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.savefig('WT Henle w Dynasore cells.pdf', dpi = 300)
plt.show()

# Create a datraframe with all of the coordinates of the target cells
# centroid-1 = y; centroid-2 = x
df_coord = pd.read_csv('-1.csv', index_col=0, usecols= [1,5,6], header=0)

cells = io.imread("202206221-01_WTHenlewDynasore.tif")
img_scale = [0.5, 0.244, 0.244] # Image pixel dimensions

# Maximum intensity projection
cells_maxproj = np.max(cells, axis = 1) # axis = 1 gives a z projection

# Contrast cells
vmin, vmax = np.quantile(cells_maxproj, q=(0.001, 1))
cells_proj_contrasted = exposure.rescale_intensity(
    cells_maxproj,
    in_range=(vmin, vmax),
    out_range=np.float32)


for n in range(len(df_norm_og.index)):
    axis_yz = int(df_coord.iloc[n]['centroid-1'])
    axis_xz = int(df_coord.iloc[n]['centroid-2'])

    n += 1
    if n in bad_cells:
        pass
    else:

        # yz slice
        #axis_yz = 200 # specify axis here

        cells_yz = cells[:, :, axis_yz, :]

        vmin, vmax = np.quantile(cells_yz, q=(0.05, 1))
        cells_yz_contrasted = exposure.rescale_intensity(
            cells_yz,
            in_range=(vmin, vmax),
            out_range=np.float32)

        # yz slice
       #axis_xz = 325 # specify axis here

        cells_xz = cells[:, :, :, axis_xz]

        vmin, vmax = np.quantile(cells_xz, q=(0.05, 1))
        cells_xz_contrasted = exposure.rescale_intensity(
            cells_xz,
            in_range=(vmin, vmax),
            out_range=np.float32)

        timepoints = [-1, 0, 68, 136, 271, 1239, 3039, 3939]

        def cross_section_plots(image, imageyz, imagexz):
            fig = Figure(figsize=(6,4), dpi=600)
            canvas = FigureCanvasAgg(fig)

            widths = [1,5.1]
            heights = [1,5.1]

            axs = fig.add_gridspec(ncols=2, nrows=2, width_ratios = widths, height_ratios =heights)
            ax1 = fig.add_subplot(axs[1,1])
            ax1.imshow(image, cmap = 'gray')
            ax1.axis('off')
            ax1.plot([0,cells_maxproj.shape[2]], [axis_yz,axis_yz], 'c-', lw=2)
            ax1.plot([axis_xz,axis_xz], [0,cells_maxproj.shape[1]], 'y-', lw=2)
            #axs['BottomLeft'].set_title('yz plane')
            rotated_xz = ndimage.rotate(imagexz, -90)
            ax2 = fig.add_subplot(axs[1,0])
            ax2.imshow(rotated_xz, cmap = 'gray')
            ax2.axis('off')
            #axs['TopRight'].set_title('xz plane')
            ax3 = fig.add_subplot(axs[0,1])
            ax3.imshow(imageyz, cmap = 'gray')
            ax3.axis('off')

            fig.subplots_adjust(wspace= -0.8,
                                hspace = 0.05)
            fig.tight_layout()
            canvas.draw()
            buf = canvas.buffer_rgba()
            X = np.asarray(buf)
            return X

        fig = plt.figure(constrained_layout=True, figsize = (5,7), dpi = 600)
        fig.suptitle('WT Henle Cells + Dynasore- Cell Swelling', fontsize=14)
        axs = fig.subplot_mosaic([['TopTopLeft', 'TopTopRight'],['TopLeft', 'TopRight'],['MidLeft', 'MidRight'], ['BottomLeft', 'BottomRight']],
                                  gridspec_kw={'width_ratios':[4, 4]})

        axs['TopTopLeft'].set_title('Pre-treatment')
        axs['TopTopLeft'].imshow(cross_section_plots(cells_proj_contrasted[0],
                                                   cells_yz_contrasted[0],
                                                   cells_xz_contrasted[0]))
        axs['TopTopLeft'].axis('off')
        axs['TopTopRight'].set_title(f'{timepoints[1]} sec')
        axs['TopTopRight'].imshow(cross_section_plots(cells_proj_contrasted[1],
                                                   cells_yz_contrasted[1],
                                                   cells_xz_contrasted[1]))
        axs['TopTopRight'].axis('off')
        axs['TopLeft'].set_title(f'{timepoints[2]} sec')
        axs['TopLeft'].imshow(cross_section_plots(cells_proj_contrasted[2],
                                                  cells_yz_contrasted[2],
                                                  cells_xz_contrasted[2]))
        axs['TopLeft'].axis('off')
        axs['TopRight'].set_title(f'{timepoints[3]} sec')
        axs['TopRight'].imshow(cross_section_plots(cells_proj_contrasted[3],
                                                   cells_yz_contrasted[3],
                                                   cells_xz_contrasted[3]))
        axs['TopRight'].axis('off')
        axs['MidLeft'].set_title(f'{timepoints[4]} sec')
        axs['MidLeft'].imshow(cross_section_plots(cells_proj_contrasted[5],
                                                  cells_yz_contrasted[5],
                                                  cells_xz_contrasted[5]))
        axs['MidLeft'].axis('off')
        axs['MidRight'].set_title(f'{timepoints[5]} sec')
        axs['MidRight'].imshow(cross_section_plots(cells_proj_contrasted[6],
                                                   cells_yz_contrasted[6],
                                                   cells_xz_contrasted[6]))
        axs['MidRight'].axis('off')
        axs['BottomLeft'].set_title(f'{timepoints[6]} sec')
        axs['BottomLeft'].imshow(cross_section_plots(cells_proj_contrasted[8],
                                                     cells_yz_contrasted[8],
                                                     cells_xz_contrasted[8]))
        axs['BottomLeft'].axis('off')
        axs['BottomRight'].set_title(f'{timepoints[7]} sec')
        axs['BottomRight'].imshow(cross_section_plots(cells_proj_contrasted[9],
                                                      cells_yz_contrasted[9],
                                                      cells_xz_contrasted[9]))
        axs['BottomRight'].axis('off')

        filename = 'WTwDynasore Images- Cell ' + str(n)

        plt.savefig(filename + '.pdf')
        plt.savefig(filename + '.png')
