# Myosoft-IMCF

IMCF-adaptation of Myosoft, a Fiji script that identifies muscle fibers in
images of sections.

Original publication: <https://doi.org/10.1371/journal.pone.0229041>

Original code: <https://github.com/Hyojung-Choo/Myosoft/tree/Myosoft-hub>

## `1_identify_fibers.py`

- Will identify all fibers based on the membrane staining using WEKA pixel
  classification, filter them according to the morphometric gates and save the
  corresponding ROIs.
- Will now also save the WEKA segmentation as a binary so it can be edited
  manually. If you do so, you need to run the "extended particle analyzer"
  manually as well to choose & apply the morphometric gates.
- Can be run in batch.

## `2a_identify_MHC_positive_fibers.py`

- Allows to manual re-run the MHC positive fiber detection. Useful in case you
  would like to re-run detection with a manual threshold for an image.

## `2b_central_nuclei_counter.py`

- Will identify centralized nuclei given a ROI-zip together with its
  corresponding image.
- Identification is based on the same logic as before incorporating the
  information of a MHC staining channel.
- The ROI color code is annotated in the results table.

## `2c_fibertyping.py`

- Identifies positive fibers in up to 3 channels given a ROI-zip together with
  its corresponding image.
- Includes identification of double and triple positive combinations.
- The ROI color code is annotated in the results table.

## `3_manual_rerun.py`

- Requires an already open image with an already populated ROI manager.
- Allows to manually select measurement parameters and the measurement channel.
- Extracts the ROI color code and stores it in the result table.

All scripts store resulting ROI-zips, logs, result tables and overview PNGs.

A potential workflow could look like this:

1. Run script 1) over night in batch mode on as many images as desired.
2. You can potentially manually curate the resulting ROIs now, or directly move
   on to the next step.
3. Run either script 2b) or 2c), depending on the assay.
4. With the results open, manually edit the ROIs and run script 3) for the final
   result.
