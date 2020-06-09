# nuclear foci
A macro to automatically quantify the foci (e.g. DNA damage) in immunofluorescence images.

Input images are read with BioFormats and should have at least two channels: a foci channel and a nuclear marker channel (e.g. DAPI)
Quantification is done in 2D. 3D images are first projected.
Currently only one channel at the time can be quantified.

Nuclear segmentation is primarily done automatically, with an option to manually edit the segmentation.

The macro supports processing multiple files in a folder. Results files can be appended using the macro 'Append_results_files.ijm'.
A simple colocalization analysis can be done using the macro 'foci_colocalization.ijm'.
