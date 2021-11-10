# this is a python rewrite of the original ijm published at
# https://github.com/Hyojung-Choo/Myosoft/blob/Myosoft-hub/Scripts/central%20nuclei%20counter.ijm

# IJ imports
# TODO: are the imports RoiManager and ResultsTable needed when using the services?
from ij import IJ, WindowManager as wm
from ij.plugin import Duplicator, RoiEnlarger, RoiScaler
from trainableSegmentation import WekaSegmentation
from de.biovoxxel.toolbox import Extended_Particle_Analyzer
from ij.measure import ResultsTable

# Bio-formats imports
from loci.plugins import BF
from loci.plugins.in import ImporterOptions

# python imports
import time
import os

#@ String (visibility=MESSAGE, value="<html><b> Welcome to Myosoft - identify fibers! </b></html>") msg1
#@ File (label="Select directory with classifiers", style="directory") classifiers_dir
#@ File (label="Select directory for output", style="directory") output_dir
#@ File (label="Select image file", description="select your image")  path_to_image
#@ Boolean (label="close image after processing", description="tick this box when using batch mode", value=False) close_raw
#@ String (visibility=MESSAGE, value="<html><b> Morphometric Gates </b></html>") msg2
#@ Integer (label="Min Area [um²]", value=10) minAr
#@ Integer (label="Max Area [um²]", value=6000) maxAr
#@ Float (label="Min Circularity", value=0.5) minCir
#@ Float (label="Max Circularity", value=1) maxCir
#@ Float (label="Min solidity", value=0.0) minSol
#@ Float (label="Max solidity", value=1) maxSol
#@ Integer (label="Min perimeter [um]", value=5) minPer
#@ Integer (label="Max perimeter [um]", value=300) maxPer
#@ Integer (label="Min min ferret [um]", value=0.1) minMinFer
#@ Integer (label="Max min ferret [um]", value=100) maxMinFer
#@ Integer (label="Min ferret AR", value=0) minFAR
#@ Integer (label="Max ferret AR", value=8) maxFAR
#@ Float (label="Min roundess", value=0.2) minRnd
#@ Float (label="Max roundess", value=1) maxRnd
#@ String (visibility=MESSAGE, value="<html><b> Expand ROIS to match fibers </b></html>") msg3
#@ Float (label="ROI expansion [microns]", value=1) enlarge
#@ String (visibility=MESSAGE, value="<html><b> channel positions in the hyperstack </b></html>") msg5
#@ Integer (label="Membrane staining channel number", style="slider", min=1, max=5, value=1) membrane_channel
#@ Integer (label="Fiber staining (MHC) channel number (0=skip)", style="slider", min=0, max=5, value=3) fiber_channel
#@ Integer (label="minimum fiber intensity (0=auto)", description="0 = automatic threshold detection", value=0) min_fiber_intensity
#@ Integer (label="sub-tiling to economize RAM", style="slider", min=1, max=8, value=4) tiling_factor

#@ RoiManager rm
#@ ResultsTable rt


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


def preprocess_membrane_channel(imp):
    """apply myosoft pre-processing steps for the membrane channel

    Parameters
    ----------
    imp : ImagePlus
        a single channel image of the membrane staining
    """
    IJ.run(imp, "Enhance Contrast", "saturated=0.35")
    IJ.run(imp, "Apply LUT", "")
    IJ.run(imp, "Enhance Contrast", "saturated=1")
    IJ.run(imp, "8-bit", "")
    IJ.run(imp, "Invert", "")
    IJ.run(imp, "Convolve...", "text1=[-1.0 -1.0 -1.0 -1.0 -1.0\n-1.0 -1.0 -1.0 -1.0 0\n-1.0 -1.0 24.0 -1.0 -1.0\n-1.0 -1.0 -1.0 -1.0 -1.0\n-1.0 -1.0 -1.0 -1.0 0] normalize")


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


def apply_weka_model(model_path, imp, tiles_per_dim):
    """apply a pretrained WEKA model to an ImagePlus

    Parameters
    ----------
    model_path : string
        path to the model file
    imp : ImagePlus
        ImagePlus to apply the model to
    tiles_per_dim : integer
        tiles the imp to save RAM

    Returns
    -------
    ImagePlus
        the result of the WEKA segmentation. One channel per class.
    """
    segmentator = WekaSegmentation()
    segmentator.loadClassifier( model_path )
    result = segmentator.applyClassifier( imp, [tiles_per_dim, tiles_per_dim], 0, True ) #ImagePlus imp, int[x,y,z] tilesPerDim, int numThreads (0=all), boolean probabilityMaps

    return result


def process_weka_result(imp):
    """apply myosoft pre-processing steps for the imp after WEKA classification to prepare it
    for ROI detection with the extended particle analyzer

    Parameters
    ----------
    imp : ImagePlus
        a single channel (= desired class) of the WEKA classification result imp
    """
    IJ.run(imp, "8-bit", "")
    IJ.run(imp, "Median...", "radius=3")
    IJ.run(imp, "Gaussian Blur...", "sigma=2")
    IJ.run(imp, "Auto Threshold", "method=MaxEntropy")
    IJ.run(imp, "Invert", "")


def delete_channel(imp, channel_number):
    """delete a channel from target imp

    Parameters
    ----------
    imp : ImagePlus
        the imp from which to delete target channel
    channel_number : integer
        the channel number to be deleted. starts at 0.
    """
    imp.setC(channel_number)
    IJ.run(imp, "Delete Slice", "delete=channel")


def run_extended_particle_analyzer( imp, eda_parameters ):
    """identifies ROIs in target imp using the extended particle analyzer of the BioVoxxel toolbox
    with given parameters

    Parameters
    ----------
    imp : ImagePlus
        the image on which to run the EPA on. Should be 8-bit thresholded
    eda_parameters : array
        all user defined parameters to restrict ROI identification
    """
    epa = Extended_Particle_Analyzer()
    epa.readInputImageParameters(imp)
    epa.setDefaultParameterFields()

    # expose all parameters explicitly
    epa.usePixel = False
    epa.usePixelForOutput = False
    epa.Area = str(eda_parameters[0]) + "-" + str(eda_parameters[1])
    epa.Extent = "0.00-1.00"
    epa.Perimeter = str(eda_parameters[2]) + "-" + str(eda_parameters[3])
    epa.Circularity = str(eda_parameters[4]) + "-" + str(eda_parameters[5])
    epa.Roundness = str(eda_parameters[6]) + "-" + str(eda_parameters[7])
    epa.Solidity = str(eda_parameters[8]) + "-" + str(eda_parameters[9])
    epa.Compactness = "0.00-1.00"
    epa.AR = "0-Infinity"
    epa.FeretAR = str(eda_parameters[10]) + "-" + str(eda_parameters[11])
    epa.EllipsoidAngle = "0-180"
    epa.MaxFeret = "0-Infinity"
    epa.MinFeret = str(eda_parameters[12]) + "-" + str(eda_parameters[13])
    epa.FeretAngle = "0-180"
    epa.COV = "0.00-1.00"
    epa.Output = "Nothing"
    epa.Redirect = "None"
    epa.Correction = "None"
    epa.Reset = False
    epa.DisplayResults = False
    epa.ClearResults = False
    epa.Summarize = False
    epa.AddToManager = True
    epa.ExcludeEdges = False
    epa.IncludeHoles = False

    epa.defineParticleAnalyzers()
    epa.particleAnalysis( imp.getProcessor(), imp, imp.getTitle() )


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


def enlarge_all_rois( amount_in_um, rm, pixel_size_in_um ):
    """enlarges all ROIs in the RoiManager by x scaled units

    Parameters
    ----------
    amount_in_um : float
        the value by which to enlarge in scaled units, e.g 3.5
    rm : RoiManager
        a reference of the IJ-RoiManager
    pixel_size_in_um : float
        the pixel size, e.g. 0.65 px/um
    """
    amount_px = amount_in_um / pixel_size_in_um
    all_rois = rm.getRoisAsArray()
    rm.reset()
    for roi in all_rois:
        enlarged_roi = RoiEnlarger.enlarge(roi, amount_px)
        rm.addRoi(enlarged_roi)


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

print rt.size()

# open image using Bio-Formats
path_to_image = fix_ij_dirs(path_to_image)
raw = open_image_with_BF(path_to_image)

# get image info
raw_image_calibration = raw.getCalibration()
raw_image_title = fix_BF_czi_imagetitle(raw)
print("raw image title: ", str(raw_image_title))

# take care of paths and directories
output_dir = fix_ij_dirs(output_dir) + "/" + str(raw_image_title) + "/1_identify_fibers"
print("output_dir: ", str(output_dir))

if not os.path.exists( str(output_dir) ):
    os.makedirs( str(output_dir) )

classifiers_dir = fix_ij_dirs(classifiers_dir)
primary_model = classifiers_dir + "/" + "primary.model"
secondary_model = classifiers_dir + "/" + "secondary_central_nuclei.model"

# update the log for the user
IJ.log( "Now working on " + str(raw_image_title) )
if raw_image_calibration.scaled() == False:
    IJ.log("Your image is not spatially calibrated! Size measurements are only possible in [px].")
IJ.log( " -- settings used -- ")
IJ.log( "area = " + str(minAr) + "-" + str(maxAr) )
IJ.log( "perimeter = " + str(minPer) + "-" + str(maxPer) )
IJ.log( "circularity = " + str(minCir) + "-" + str(maxCir) )
IJ.log( "roundness = " + str(minRnd) + "-" + str(maxRnd) )
IJ.log( "solidity = " + str(minSol) + "-" + str(maxSol) )
IJ.log( "feret_ar = " + str(minFAR) + "-" + str(maxFAR) )
IJ.log( "min_feret = " + str(minMinFer) + "-" + str(maxMinFer) )
IJ.log( "ROI expansion [microns] = " + str(enlarge) )
IJ.log( "Membrane channel = " + str(membrane_channel) )
IJ.log( "MHC positive fiber channel = " + str(fiber_channel) )
IJ.log( "sub-tiling = " + str(tiling_factor) )
IJ.log( " -- settings used -- ")

# image (pre)processing and segmentation (-> ROIs)
membrane = Duplicator().run(raw, membrane_channel, membrane_channel, 1, 1, 1, 1) # imp, firstC, lastC, firstZ, lastZ, firstT, lastT
preprocess_membrane_channel(membrane)
weka_result1 = apply_weka_model(primary_model, membrane, tiling_factor )
delete_channel(weka_result1, 1)
weka_result2 = apply_weka_model(secondary_model, weka_result1, tiling_factor )
delete_channel(weka_result2, 1)
weka_result2.setCalibration(raw_image_calibration)
process_weka_result(weka_result2)
IJ.saveAs(weka_result2, "Tiff", output_dir + "/" + raw_image_title + "_all_fibers_binary")
eda_parameters = [minAr, maxAr, minPer, maxPer, minCir, maxCir, minRnd, maxRnd, minSol, maxSol, minFAR, maxFAR, minMinFer, maxMinFer]
raw.show() # EPA will not work if no image is shown
run_extended_particle_analyzer(weka_result2, eda_parameters)

# modify rois
rm.hide()
raw.hide()
enlarge_all_rois( enlarge, rm, raw_image_calibration.pixelWidth )
renumber_rois(rm)
save_all_rois( rm, output_dir + "/" + raw_image_title + "_all_fiber_rois.zip" )

# check for positive fibers
if fiber_channel > 0:
    if min_fiber_intensity == 0:
        min_fiber_intensity = get_threshold_from_method(raw, fiber_channel, "Mean")[0]
        IJ.log( "automatic intensity threshold detection: True" )

    IJ.log( "fiber intensity threshold: " + str(min_fiber_intensity) )
    change_all_roi_color(rm, "blue")
    positive_fibers = select_positive_fibers( raw, fiber_channel, rm, min_fiber_intensity  )
    change_subset_roi_color(rm, positive_fibers, "magenta")
    save_selected_rois( rm, positive_fibers, output_dir + "/" + raw_image_title + "_mhc_positive_fiber_rois.zip")

# measure size & shape, save
IJ.run("Set Measurements...", "area perimeter shape feret's redirect=None decimal=4")
IJ.run("Clear Results", "")
measure_in_all_rois( raw, membrane_channel, rm )

rt = ResultsTable.getResultsTable("Results")

print rt.size()

if fiber_channel > 0:
    print rt.size()
    preset_results_column( rt, "MHC Positive Fibers (magenta)", "NO" )
    print rt.size()
    add_results( rt, "MHC Positive Fibers (magenta)", positive_fibers, "YES")
    print rt.size()

rt.save(output_dir + "/" + raw_image_title + "_all_fibers_results.csv")
print "saved the all_fibers_results.csv"
# dress up the original image, save a overlay-png, present original to the user
rm.show()
raw.show()
show_all_rois_on_image( rm, raw )
raw.setDisplayMode(IJ.COMPOSITE)
enhance_contrast( raw )
IJ.run("From ROI Manager", "") # ROIs -> overlays so they show up in the saved png
qc_duplicate = raw.duplicate()
IJ.saveAs(qc_duplicate, "PNG", output_dir + "/" + raw_image_title + "_all_fibers")
qc_duplicate.close()
wm.toFront( raw.getWindow() )
IJ.run("Remove Overlay", "")
raw.setDisplayMode(IJ.GRAYSCALE)
show_all_rois_on_image( rm, raw )
total_execution_time_min = (time.time() - execution_start_time) / 60.0
IJ.log("total time in minutes: " + str(total_execution_time_min))
IJ.log( "~~ all done ~~" )
IJ.selectWindow("Log")
IJ.saveAs("Text", str(output_dir + "/" + raw_image_title + "_all_fibers_Log"))
if close_raw == True:
    raw.close()