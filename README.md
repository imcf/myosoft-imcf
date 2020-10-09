# myosoft-imcf

imcf-adaptation of Myosoft, a Fiji script that identifies muscle fibers in images of sections

original publication: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0229041
original code: https://github.com/Hyojung-Choo/Myosoft/tree/Myosoft-hub


1) myosoft-imcf_identify_fibers.py
- Will identify all fibers based on the membrane staining using WEKA pixel classification, filter them according to the morphometric gates and save the corresponding ROIs
- will now also save the WEKA segmentation as a binary so it can be edited manually. If you do so, you need to run the "extended particle analyzer" manually as well to choose & apply the morphometric gates.
- can be run in batch

2a) myosoft-imcf_central_nuclei_counter.py
- will identify centralized nuclei given a ROI-zip together with its corresponding image
- identification is based on the same logic as before incorporating the information of a MHC staining channel
- the ROI color code is annotated in the results table

2b) myosoft-imcf_fibertyping.py
- identifies positive fibers in up to 3 channels given a ROI-zip together with its corresponding image
- includes identification of double and triple positive combinations
- the ROI color code is annotated in the results table

3) myosoft-imcf_manual_rerun.py
- requires an already open image with an already populated ROI manager
- allows to manually select measurement parameters and the measurement channel
- extracts the ROI color code and stores it in the result table

All scripts  store resulting ROI-zips, logs, Result tables and overview pngs.

A potential workflow would look like this:

1. Run script 1) over night in batch mode on as many images as desired
2. you can potentially manually curate the resulting ROIs now, or directly move on to the next step
3. run either script 2a) or 2b), depending on the assay
4. with the results open, manually edit the ROIs and run script 3) for the final result
