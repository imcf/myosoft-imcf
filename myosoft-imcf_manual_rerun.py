from ij import IJ, WindowManager as wm
from ij.gui import WaitForUserDialog
import os

#@ ImagePlus raw
#@ RoiManager rm
#@ ResultsTable rt
#@ File (label="Select directory for output", style="directory") output_dir
#@ Integer (label="Measure in this channel", style="slider", min=1, max=5, value=1) measurement_channel


def fix_ij_dirs(path):
    """use forward slashes in directory paths

    Parameters
    ----------
    path : string
        a directory path obtained from dialogue or script parameter

    Returns
    -------
    string
        a more robust path with forward slashes as 
    """

    fixed_path = str(path).replace("\\", "/")
    fixed_path = fixed_path + "/"

    return fixed_path


def fix_BF_czi_imagetitle(imp):
    image_title = os.path.basename( imp.getTitle() )
    image_title = image_title.replace(".czi", "")
    image_title = image_title.replace(" ", "_")
    image_title = image_title.replace("_-_", "")
    image_title = image_title.replace("__", "_")
    image_title = image_title.replace("#", "Series")

    return image_title


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


def extract_color_of_all_rois(rm):
    """get the RGB color of ROIs in the RoiManager and match it to the colors name string

    Parameters
    ----------
    rm : RoiManager
        the IJ-RoiManager

    Returns
    -------
    array
        an array containing the corresponding color name string for each roi in the ROiManager
    """
    rgb_color_lookup = {
    -65536: "red", 
    -65281: "magenta", 
    -16711936: "green", 
    -256: "yellow", 
    -1: "white", 
    -16776961: "blue", 
    -16777216: "black", 
    -14336: "orange", 
    -16711681: "cyan" 
    }

    all_rois = rm.getRoisAsArray()
    roi_colors = []
    for roi in all_rois:
        if roi.getStrokeColor() == None:
            roi_colors.append(rgb_color_lookup[roi.getColor().getRGB()])
        else:
            roi_colors.append(rgb_color_lookup[roi.getStrokeColor().getRGB()])

    return roi_colors


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


def add_results_to_resultstable( rt, column, values ):
    """add values to the ResultsTable starting from row 0 of a given column

    Parameters
    ----------
    rt : ResultsTable
        a reference of the IJ-ResultsTable
    column : string
        the column in which to add the values
    values : array
        tarray with values to be added
    """
    for i in range( len( values ) ):
        rt.setValue(column, i, values[i])

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


output_dir = fix_ij_dirs(output_dir)
rt.reset()
raw_image_title = fix_BF_czi_imagetitle(raw)
renumber_rois(rm)
save_all_rois( rm, output_dir + "manual_rerun_all_fiber_rois_color-coded.zip" )
roi_colors = extract_color_of_all_rois(rm)
WaitForUserDialog("Choose measurements", "Set measurements in Analyze > Set Measurements, then click OK").show()
measure_in_all_rois(raw, measurement_channel, rm)
add_results_to_resultstable(rt, "ROI color", roi_colors )
rt.save(output_dir + "manual_rerun_results.csv")

# dress up the original image, save a overlay-png, present original to the user
raw.show()
show_all_rois_on_image( rm, raw )
raw.setDisplayMode(IJ.COMPOSITE)
enhance_contrast( raw )
IJ.run("From ROI Manager", "") # ROIs -> overlays so they show up in the saved png
qc_duplicate = raw.duplicate()
IJ.saveAs(qc_duplicate, "PNG", output_dir + raw_image_title + "_manual_rerun")
qc_duplicate.close()
wm.toFront( raw.getWindow() )
IJ.run("Remove Overlay", "")
raw.setDisplayMode(IJ.GRAYSCALE)
show_all_rois_on_image( rm, raw )