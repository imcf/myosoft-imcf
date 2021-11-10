# this is a python rewrite of the original ijm published at 
# https://github.com/Hyojung-Choo/Myosoft/blob/Myosoft-hub/Scripts/central%20nuclei%20counter.ijm

# IJ imports
# TODO: are the imports RoiManager and ResultsTable needed when using the services?
from ij import IJ, WindowManager as wm
from ij.plugin import Duplicator, RoiEnlarger, RoiScaler
from trainableSegmentation import WekaSegmentation
from de.biovoxxel.toolbox import Extended_Particle_Analyzer

# Bio-formats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions

# python imports
import time
import os

#@ String (visibility=MESSAGE, value="<html><b> Welcome to Myosoft! </b></html>") msg1
#@ File (label="Select fiber-ROIs zip-file", style="file") roi_zip
#@ File (label="Select image file", description="select your image") path_to_image
#@ File (label="Select directory for output", style="directory") output_dir
#@ Boolean (label="close image after processing", description="tick this box when using batch mode", value=False) close_raw
#@ String (visibility=MESSAGE, value="<html><b> channel positions in the hyperstack </b></html>") msg5
#@ Integer (label="Fiber staining 1 channel number (0=n.a.)", style="slider", min=0, max=5, value=1) fiber_channel_1
#@ Integer (label="Fiber staining 2 channel number (0=n.a.)", style="slider", min=0, max=5, value=2) fiber_channel_2
#@ Integer (label="Fiber staining 3 channel number (0=n.a.)", style="slider", min=0, max=5, value=3) fiber_channel_3
#@ Integer (label="minimum fiber intensity (0=auto)", description="0 = automatic threshold detection", value=0) min_fiber_intensity_1
#@ Integer (label="minimum fiber intensity (0=auto)", description="0 = automatic threshold detection", value=0) min_fiber_intensity_2
#@ Integer (label="minimum fiber intensity (0=auto)", description="0 = automatic threshold detection", value=0) min_fiber_intensity_3
#@ ResultsTable rt
#@ RoiManager rm


def fix_ij_options():
    """put IJ into a defined state
    """
    # disable inverting LUT
    IJ.run("Appearance...", " menu=0 16-bit=Automatic")
    # set foreground color to be white, background black
    IJ.run("Colors...", "foreground=white background=black selection=red")
    # black BG for binary images and pad edges when eroding
    IJ.run("Options...", "black pad")
    # set saving format to .txt files
    IJ.run("Input/Output...", "file=.txt save_column save_row")
    # ============= DON’T MOVE UPWARDS =============
    # set "Black Background" in "Binary Options"
    IJ.run("Options...", "black")
    # scale when converting = checked
    IJ.run("Conversions...", "scale")


def fix_ij_dirs(path):
    """use forward slashes in directory paths

    Parameters
    ----------
    path : string
        a directory path obtained from dialogue or script parameter

    Returns
    -------
    string
        a more robust path with forward slashes as separators
    """

    fixed_path = str(path).replace("\\", "/")
    # fixed_path = fixed_path + "/"

    return fixed_path


def open_image_with_BF(path_to_file):
    """ use Bio-Formats to opens the first image from an image file path

    Parameters
    ----------
    path_to_file : string
        path to the image file

    Returns
    -------
    ImagePlus
        the first imp stored in a give file
    """
    options = ImporterOptions()
    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE)
    options.setAutoscale(True)
    options.setId(path_to_file)
    imps = BF.openImagePlus(options) # is an array of ImagePlus

    return imps[0]


def fix_BF_czi_imagetitle(imp):
    image_title = os.path.basename( imp.getShortTitle() )
    image_title = image_title.replace(".czi", "")
    image_title = image_title.replace(" ", "_")
    image_title = image_title.replace("_-_", "")
    image_title = image_title.replace("__", "_")
    image_title = image_title.replace("#", "Series")

    return image_title


def clear_ij_roi_manager(rm):
    """delete all ROIs from the RoiManager

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    """
    rm.runCommand('reset')


def get_threshold_from_method(imp, channel, method):
    """returns the threshold value of chosen IJ AutoThreshold method in desired channel

    Parameters
    ----------
    imp : ImagePlus
        the imp from which to get the threshold value
    channel : integer
        the channel in which to get the treshold
    method : string
        the AutoThreshold method to use

    Returns
    -------
    list
        the upper and the lower threshold (integer values)
    """
    imp.setC(channel) # starts at 1
    ip = imp.getProcessor()
    ip.setAutoThreshold(method + " dark")
    lower_thr = ip.getMinThreshold()
    upper_thr = ip.getMaxThreshold()
    ip.resetThreshold()

    return lower_thr, upper_thr


def measure_in_all_rois( imp, channel, rm ):
    """measures in all ROIS on a given channel of imp all parameters that are set in IJ "Set Measurements"

    Parameters
    ----------
    imp : ImagePlus
        the imp to measure on
    channel : integer
        the channel to measure in. starts at 1.
    rm : RoiManager
        a reference of the IJ-RoiManager
    """
    imp.setC(channel)
    rm.runCommand(imp,"Deselect")
    rm.runCommand(imp,"Measure")


def change_all_roi_color( rm, color ):
    """change the color of all ROIs in the RoiManager

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    color : string
        the desired color. e.g. "green", "red", "yellow", "magenta" ...
    """
    number_of_rois = rm.getCount()
    for roi in range( number_of_rois ):
        rm.select(roi)
        rm.runCommand("Set Color", color)


def change_subset_roi_color( rm, selected_rois, color ):
    """change the color of selected ROIs in the RoiManager

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    selected_rois : array
        ROIs in the RoiManager to change
    color : string
        the desired color. e.g. "green", "red", "yellow", "magenta" ...
    """
    rm.runCommand("Deselect")
    rm.setSelectedIndexes(selected_rois)
    rm.runCommand("Set Color", color)
    rm.runCommand("Deselect")


def show_all_rois_on_image(rm, imp):
    """shows all ROIs in the ROiManager on imp

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    imp : ImagePlus
        the imp on which to show the ROIs
    """
    imp.show()
    rm.runCommand(imp,"Show All")


def save_all_rois(rm, target):
    """save all ROIs in the RoiManager as zip to target path

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    target : string
        the path in to store the ROIs. e.g. /my-images/resulting_rois.zip
    """
    rm.runCommand("Save", target)


def save_selected_rois( rm, selected_rois, target ):
    """save selected ROIs in the RoiManager as zip to target path

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    selected_rois : array
        ROIs in the RoiManager to save
    target : string
        the path in to store the ROIs. e.g. /my-images/resulting_rois_subset.zip
    """
    rm.runCommand("Deselect")
    rm.setSelectedIndexes(selected_rois)
    rm.runCommand("save selected", target)
    rm.runCommand("Deselect")


def select_positive_fibers( imp, channel, rm, min_intensity ):
    """For all ROIs in the RoiManager, select ROIs based on intensity measurement in given channel of imp.
    See https://imagej.nih.gov/ij/developer/api/ij/process/ImageStatistics.html

    Parameters
    ----------
    imp : ImagePlus
        the imp on which to measure
    channel : integer
        the channel on which to measure. starts at 1
    rm : RoiManager
        a reference of the IJ-RoiManager
    min_intensity : integer
        the selection criterion (here: intensity threshold)
    
    Returns
    -------
    array
        a selection of ROIs which passed the selection criterion (are above the threshold)
    """
    imp.setC(channel)
    all_rois = rm.getRoisAsArray()
    selected_rois = []
    for i, roi in enumerate(all_rois):
        imp.setRoi(roi)
        stats = imp.getStatistics()
        if stats.mean > min_intensity:
            selected_rois.append(i)

    return selected_rois 


def open_rois_from_zip( rm, path ):
    """open RoiManager ROIs from zip and adds them to the RoiManager

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    path : string
        path to the ROI zip file
    """
    rm.runCommand("Open", path)


def preset_results_column( rt, column, value):
    """pre-set all rows in given column of the IJ-ResultsTable with desired value

    Parameters
    ----------
    rt : ResultsTable
        a reference of the IJ-ResultsTable
    column : string
        the desired column. will be created if it does not yet exist
    value : string or float or integer
        the value to be set
    """
    for i in range( rt.size() ):
        rt.setValue(column, i, value)

    rt.show("Results")


def add_results( rt, column, row, value ):
    """adds a value in desired rows of a given column

    Parameters
    ----------
    rt : ResultsTable
        a reference of the IJ-ResultsTable
    column : string
        the column in which to add the values
    row : array
        the row numbers in which too add the values.
    value : string or float or integer
        the value to be set
    """
    for i in range( len( row ) ):
        rt.setValue(column, row[i], value)

    rt.show("Results")


def enhance_contrast( imp ):
    """use "Auto" Contrast & Brightness settings in each channel of imp

    Parameters
    ----------
    imp : ImagePlus
        the imp on which to change C&B
    """
    for channel in range( imp.getDimensions()[2] ):
        imp.setC(channel + 1) # IJ channels start at 1
        IJ.run(imp, "Enhance Contrast", "saturated=0.35")


def renumber_rois(rm):
    """rename all ROIs in the RoiManager according to their number

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    """
    number_of_rois = rm.getCount()
    for roi in range( number_of_rois ):
        rm.rename( roi, str(roi + 1) )


def setup_defined_ij(rm, rt):
    """set up a clean and defined Fiji user environment

    Parameters
    ----------
    rm : RoiManager
        a reference of the IJ-RoiManager
    rt : ResultsTable
        a reference of the IJ-ResultsTable
    """
    fix_ij_options()
    rm.runCommand('reset')
    rt.reset()
    IJ.log("\\Clear")


execution_start_time = time.time()
setup_defined_ij(rm, rt)

# open image using Bio-Formats
path_to_image = fix_ij_dirs(path_to_image)
raw = open_image_with_BF(path_to_image)

# get image info
raw_image_calibration = raw.getCalibration()
raw_image_title = fix_BF_czi_imagetitle(raw)

# take care of paths and directories
input_rois_path = fix_ij_dirs( roi_zip )
output_dir = fix_ij_dirs(output_dir) + "/2c_fibertyping"

if not os.path.exists( str(output_dir) ):
    os.makedirs( str(output_dir) )

# open ROIS and show on image
open_rois_from_zip( rm, str(input_rois_path) )
change_all_roi_color(rm, "blue")
show_all_rois_on_image( rm, raw )

# update the log for the user
IJ.log( "Now working on " + str(raw_image_title) )
if raw_image_calibration.scaled() == False:
    IJ.log("Your image is not spatially calibrated! Size measurements are only possible in [px].")
IJ.log( " -- settings used -- ")
IJ.log( "Selected fiber-ROIs zip-file = " + str(input_rois_path) )
IJ.log( "Fiber staining 1 channel number = " + str(fiber_channel_1) )
IJ.log( "Fiber staining 2 channel number = " + str(fiber_channel_2) )
IJ.log( "Fiber staining 3 channel number = " + str(fiber_channel_3) )
IJ.log( " -- settings used -- ")

# measure size & shape,
IJ.run("Set Measurements...", "area perimeter shape feret's redirect=None decimal=4")
IJ.run("Clear Results", "")
measure_in_all_rois( raw, fiber_channel_1, rm )

# loop through the fiber channels, check if positive, add info to results table
all_fiber_channels = [fiber_channel_1, fiber_channel_2, fiber_channel_3]
all_min_fiber_intensities = [min_fiber_intensity_1, min_fiber_intensity_2, min_fiber_intensity_3]
roi_colors = ["green", "orange", "red"]
all_fiber_subsets =[ [], [], [] ]

for index, fiber_channel in enumerate(all_fiber_channels):
    if fiber_channel > 0:
        preset_results_column( rt, "channel " + str(fiber_channel) + " positive (" + roi_colors[fiber_channel-1] + ")", "NO" )
        if all_min_fiber_intensities[index] == 0:
            all_min_fiber_intensities[index] = get_threshold_from_method(raw, fiber_channel, "Mean")[0]
        IJ.log( "fiber channel " + str(fiber_channel) + " intensity threshold: " + str(all_min_fiber_intensities[index]) ) 
        positive_fibers = select_positive_fibers( raw, fiber_channel, rm, all_min_fiber_intensities[index] )
        all_fiber_subsets[index] = positive_fibers
        if len(positive_fibers) > 0:
            change_subset_roi_color(rm, positive_fibers, roi_colors[index])
            save_selected_rois( rm, positive_fibers, output_dir + "/" + raw_image_title + "_positive_fiber_rois_c" + str( fiber_channel ) + ".zip")
            add_results( rt, "channel " + str(fiber_channel) + " positive (" + roi_colors[fiber_channel-1] + ")", positive_fibers, "YES")

# single positive
positive_c1 = all_fiber_subsets[0]
positive_c2 = all_fiber_subsets[1]
positive_c3 = all_fiber_subsets[2]
# double positive
positive_c1_c2 = list( set(all_fiber_subsets[0]).intersection(all_fiber_subsets[1]) )
positive_c1_c3 = list( set(all_fiber_subsets[0]).intersection(all_fiber_subsets[2]) )
positive_c2_c3 = list( set(all_fiber_subsets[1]).intersection(all_fiber_subsets[2]) )
# triple positive
positive_c1_c2_c3 = list( set(positive_c1_c2).intersection(all_fiber_subsets[2]) )

# update ROI color & results table for double and triple positives
if len(positive_c1_c2) > 0:
    preset_results_column( rt, "channel 1,2 positive (magenta)", "NO" )
    change_subset_roi_color(rm, positive_c1_c2, "magenta")
    save_selected_rois( rm, positive_c1_c2, output_dir + "/" + raw_image_title + "_positive_fiber_rois_c1_c2.zip")
    add_results( rt, "channel 1,2 positive (magenta)", positive_c1_c2, "YES")

if len(positive_c1_c3) > 0:
    preset_results_column( rt, "channel 1,3 positive (yellow)", "NO" )
    change_subset_roi_color(rm, positive_c1_c3, "yellow")
    save_selected_rois( rm, positive_c1_c3, output_dir + "/" + raw_image_title + "_positive_fiber_rois_c1_c3.zip")
    add_results( rt, "channel 1,3 positive (yellow)", positive_c1_c3, "YES")

if len(positive_c2_c3) > 0:
    preset_results_column( rt, "channel 2,3 positive (cyan)", "NO" )
    change_subset_roi_color(rm, positive_c2_c3, "cyan")
    save_selected_rois( rm, positive_c2_c3, output_dir + "/" + raw_image_title + "_positive_fiber_rois_c2_c3.zip")
    add_results( rt, "channel 2,3 positive (cyan)", positive_c2_c3, "YES")

if len(positive_c1_c2_c3) > 0:
    preset_results_column( rt, "channel 1,2,3 positive (white)", "NO" )
    change_subset_roi_color(rm, positive_c1_c2_c3, "white")
    save_selected_rois( rm, positive_c1_c2_c3, output_dir + "/" + raw_image_title + "_positive_fiber_rois_c1_c2_c3.zip")
    add_results( rt, "channel 1,2,3 positive (white)", positive_c1_c2_c3, "YES")

# save all results together
save_all_rois( rm, output_dir + "/" + raw_image_title + "_all_fiber_type_rois_color-coded.zip" )
rt.save(output_dir + "/" + raw_image_title + "_fibertyping_results.csv")

# dress up the original image, save a overlay-png, present original to the user
raw.show()
show_all_rois_on_image( rm, raw )
raw.setDisplayMode(IJ.COMPOSITE)
enhance_contrast( raw )
IJ.run("From ROI Manager", "") # ROIs -> overlays so they show up in the saved png
qc_duplicate = raw.duplicate()
IJ.saveAs(qc_duplicate, "PNG", output_dir + "/" + raw_image_title + "_fibertyping")
qc_duplicate.close()
wm.toFront( raw.getWindow() )
IJ.run("Remove Overlay", "")
raw.setDisplayMode(IJ.GRAYSCALE)
show_all_rois_on_image( rm, raw )
total_execution_time_min = (time.time() - execution_start_time) / 60.0
IJ.log("total time in minutes: " + str(total_execution_time_min))
IJ.log( "~~ all done ~~" )
IJ.selectWindow("Log")
IJ.saveAs("Text", str(output_dir + "/" + raw_image_title + "_fibertyping_Log"))
if close_raw == True:
    raw.close()