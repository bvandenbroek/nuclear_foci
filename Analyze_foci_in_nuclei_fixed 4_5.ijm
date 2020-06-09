/*
 * Macro to count damage foci in one channel taking nucleus information from a second channel.
 * Required input: an image (stack) with at least these two channels.
 * 
 * Multi-file support: the user should first select a file in the to-be-analyzed directory.
 * and set the correct settings in a dialog. Then the user is asked whether to analyze all the
 * files (with the same extension) in the directory and whether the macro should pause after each image.
 * 
 * For time series another macro is available.
 *  
 * Bram van den Broek, Netherlands Cancer Institute, 2012-2020
 * b.vd.broek@nki.nl
 * 
 * version 3, February 2014:
 * - possibility for threshold bias for nuclei segmentation
 * - added total nucleus intensity measurement
 * - possibility to do the analysis on 'focussed' images (best slice selected for every nucleus)
 * - huge speed improvement, especially for large images
 * - fixed small mistake in automatic foci threshold
 * - fixed mistake in stddev_foci_intensity
 * - slightly improved foci background subtraction
 * - 
 * - not background-correcting foci now actually works
 * - foci statistics plot is now also sorted on nr. of foci and intensity
 * 
 * version 3.2, March 2014:
 * - background and stddev estimation in foci channel is now correctly measured. SNR should be set a little higher
 * - Log window displays more info and is also saved
 * - small bugfixes
 * 
 * version 3.3, February 2015:
 * - files are saved in a "results" directory
 * 
 * version 3.4, April 2015:
 * - Added the possibility to shrink the nuclear ROIs, useful in case of cytoplasmic or other staining very close to the nucleus
 * 
 * version 3.5, May 2015:
 * - The distance of the foci to the nuclear perimeter is also calculated
 * - Possibility to manually merge and delete ROIs
 * 
 * version 3.6, June 2015:
 * - Fixed an error in the measurement of the channel that is being excluded.
 * 
 * version 3.7, July 2015:
 * - Fixed Java bug about sometimes getting NaN when applying getResult
 * 
 * version 3.8, March 2016:
 * - CHANGED: threshold is now calculated for every nucleus seperately. Use with care, for highly variable cell brightnesses.
 * 
 * version 3.9 September 2017, For Apostolos:
 * - Removed distance measuring of foci to the nucleus edge
 * - Changed the order of background subtraction and get_mean_and_stddev
 * - Several speed improvements (multi-measure in get_mean_and_stddev, and not displaying while counting foci).
 * 
 * version 3.91 April 2018:
 * - Fixed a mistake regarding sigma_large (and sigma_small) to be constant in stead of changing with the foci_diameter
 * - Some other small things
 * 
 * version 3.93 April 2018:
 * - Added a variable for the auto local thresholding of nuclei
 * 
 * version 3.98 October 2018:
 * - Fixed issue with newer versions of ImageJ (>1.52e or so)
 * - Added local thresholding option in the dialog. It's becoming ugly now, but will all be much better in the next incarnation...
 * 
 * version 4.00 November 2018:
 * - Added measurement of nuclear shape factors (circularity, roundness, solidity; all in 2D)
 * - decreased the effect of the estimated foci size by a factor 1.5, to better match the actual foci size.
 * - files are saved in a folder named "results_ch"+ch_foci.
 * 
 *  version 4.02 May 2019:
 *  - Show DAPI channel in analyzed image by default
 *  
 *  version 4.1 December 2019: Improved nuclear segmentation (for large images)
 *  - Changed some things regarding nuclear segmentation (median filtering is now properly done for glabal thresholding)
 *  - Changed nucleus median filter size to 5 (was 3)
 *  - Disabled background subtraction in nuclei channel (new default)
 *  
 *  version 4.2 Januari 2020: Foci intensity counting in addional channel
 *  - New in results table: In the foci mask, mean and total intensity of additional channel (if selected)
 *  - Changed back nucleus median filter size to 3 (was not succesful for some images)
 *  - No differential threshold if threshold_scaling_coeff == -1 (was -999 or so).
 *  
 *  version 4.3 Januari 2020: Included percentile measurements (lower and upper) for the additional channel
 *  2020-01-27 - Doesn't work properly yet!
 *  
 *  version 4.4 Januari 2020:
 *  - Saves .tsv files now
 *  - Trying to play around with background subtraction/outlier removal
 *  
 *  version 4.5 May 2020:
 *  - Finally ditched the BioFormats windowless importer
 *  - Improved nuclei editing. N.B. Editing now ends with space bar
 */


requires("1.52n");

saveSettings();
setOption("BlackBackground", true);
//setOption("DisablePopupMenu", true);
run("Conversions...", " ");
roiManager("Show None");

var exclude_edges = false;
var SNR = 10;
var pause=false;
var clean=true;
var large_dots=true;
var current_slice = 1;
var median_radius_nuclei = 5;	//Was set at 2!
var lower_limit = 0;
var analyze_all = false;
var pause_after_file = false;

var auto_threshold_nuclei = true;
var auto_local_th_nuc = false;

	//LOCAL THRESHOLD PARAMETERS - CHANGE AS SEEN FIT
	var alth_radius = 100;
	var alth_parameter_1 = -5;

var bgsubtr_nuclei_segmentation = false;
var rolling_ball_radius = 100;

var watershed = true;		//separate touching nuclei
var speed_mode = true;
var remove_pixel_calibration = false;
var logarithm = false;				//this may help if the nuclei signal is very rough
var threshold_scaling_coeff = -1;	//-1 for no differential thresholding, 0.5 for the square root of the difference with the average threshold of all nuclei. 
var lowerPercentile = 0.0;
var upperPercentile = 1.0;


var nuclei_found = true;
var current_image_nr = 0;
var start_time = 0;
var run_time = 0;
var th_nuc_min = 0;
var th_nuc_max = 255;
var threshold_nuc_bias = 0;		//manual absolute bias of the autothreshold method
var outliers_threshold = 20;	//threshold for foci outlier removal before retreiving median and stddev of nucleus - TO DO: measure it!
var outliers_diameter = 5;		//Remove outliers in the nuclei channel before segmentation
var avg_median;					//average of outlier-removed median of all nuclei (ok, a bit dubble but you can never be too sure) 
var avg_stddev;					//average of outlier-removed stddev of all nuclei
var background;					//background value (everything but nuclei)
var threshold;
var merged_image;
var slices=0;


//settings from dialog
var ch_nuclei = 1;
var ch_foci = 2;
var nuclei_method_array = newArray("automatically select slice with sharpest nuclei", "maximum intensity projection", "sum slices", "manually select slice");
var nuclei_method = "automatically select slice with sharpest nuclei";	//default method
var Min_Nucleus_Size = 4;	//apparently this is always in pixels - check it!
var Max_Nucleus_Size = 40;
foci_projection_method_array = newArray("maximum intensity projection", "sum slices", "smart z-stack projection (recommended for widefield)", "best slice for every nucleus");
var foci_projection_method = "maximum intensity projection";	//default
var foci_diameter = 3.0;	//estimated size of foci (used in outlier removal and background subtraction)
var min_foci_size = 0.0;	//minimum foci size in pixels
var focus_nuclei = false;	//select slices of the nucleus signal with slice numbers found in the foci channel
var choice_additional_nucleus_channel = false;
var ch_measure_additional = 3;
var choice_exclude_nucleus_channel = false;
var ch_exclude = 4;
var foci_background_correct = true;
var manual_foci_threshold = 0;
choice_exclude_threshold_array = newArray("below", "above");
var choice_exclude_threshold = "below";	//default
var exclude_threshold = 0;
var manually_edit_nuclei = true;
var shrink_nuclei = false;		//erode the segmented nuclei selections. Useful if the cytoplasm also has staining
var shrink_nuclei_size = 0;		//size of erosion of the nucleus, in units (e.g. microns)

var verbose = false;

var sigma_large;	//size of Gaussian blur filter (large) for foci background subtraction; adjusted after the dialog
var sigma_small;	//size of Gaussian blur filter (small) for foci background subtraction; adjusted after the dialog


if(nImages>0) run("Close All");
path = File.openDialog("Select a File");
run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
run("Bio-Formats Macro Extensions");
dir = getDirectory("image");
file_name = getInfo("image.filename");
rename(file_name);

Stack.getDimensions(width, height, channels, slices, frames);
if (bitDepth()!=24 && channels==1) exit("Multi-channel image required.");
if (bitDepth()==24) run("Make Composite");	//convert RGB images
if (bitDepth()==16) bits=16;
else bits=8;
run("32-bit");
Stack.getDimensions(width, height, channels, slices, frames);
if(remove_pixel_calibration==true) run("Properties...", "unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1");	//For Martijn's wrong pixel calibration
getPixelSize(unit, pw, ph, pd);
setLocation(0, 0);


extension_length=(lengthOf(file_name)- lengthOf(File.nameWithoutExtension)-1);
extension = substring(file_name, (lengthOf(file_name)-extension_length));
file_list = getFileList(dir); //get filenames of directory

//make a list of images in that folder with 'extension' as extension.
j=0;
image_list=newArray(file_list.length);	//Dynamic array size doesn't work on some Macs, so first make image_list the maximal size and then trim.
for(i=0; i<file_list.length; i++){
	if (endsWith(file_list[i],extension)) {
		image_list[j] = file_list[i];
		j++;
	}
}
image_list = Array.trim(image_list, j);	//Trimming the array of images

print("\\Clear");
print("Directory contains "+file_list.length+" files, of which "+image_list.length+" "+extension+" files.");
print("\n");

//---------CONFIG FILE INITIATIONS
tempdir = getDirectory("temp");
config_file = tempdir+"\\foci_macro_config.txt";
if (File.exists(config_file)) {
	config_string = File.openAsString(config_file);
	config_array = split(config_string,"\n");
	if (config_array.length==25 || config_array.length==23) {
		ch_nuclei = parseInt(config_array[0]);
		ch_foci = parseInt(config_array[1]);

		Min_Nucleus_Size = parseFloat(config_array[2]);
		Max_Nucleus_Size = parseFloat(config_array[3]);
		auto_threshold_nuclei = parseInt(config_array[4]);
		threshold_nuc_bias = parseFloat(config_array[5]);
		auto_local_th_nuc = parseInt(config_array[6]);
		exclude_edges = parseInt(config_array[7]);
		shrink_nuclei = parseInt(config_array[8]);
		shrink_nuclei_size = parseFloat(config_array[9]);		
		manually_edit_nuclei = parseInt(config_array[10]);
		
		foci_diameter = parseFloat(config_array[11]);
		min_foci_size = parseFloat(config_array[12]);
		SNR = parseFloat(config_array[13]);
		foci_background_correct=parseInt(config_array[14]);
		manual_foci_threshold = config_array[15];

		choice_additional_nucleus_channel = parseInt(config_array[16]);
		ch_measure_additional = parseInt(config_array[17]);
		choice_exclude_nucleus_channel = parseInt(config_array[18]);
		ch_exclude = parseInt(config_array[19]);
		choice_exclude_threshold = config_array[20];
		exclude_threshold = parseInt(config_array[21]);

		verbose = parseInt(config_array[22]);

		if(config_array.length==25) nuclei_method = config_array[23];
		if(config_array.length==25) foci_projection_method = config_array[24];
	}
}



//---------OPEN DIALOG
Dialog.createNonBlocking("Options");
	Dialog.addSlider("nuclei channel nr", 1, channels, ch_nuclei);
	Dialog.addSlider("foci channel nr", 1, channels, ch_foci);
	//Dialog.setInsets(0, 20, 0);
	if(slices>1) Dialog.addChoice("Select method for nuclei z-projection", nuclei_method_array, nuclei_method);
	/*
	if(unit=="microns") {
		Dialog.addSlider("Mininum nucleus diameter ("+unit+")", 1, 100, Min_Nucleus_Size);
		Dialog.addSlider("Maximum nucleus diameter ("+unit+")", 1, 100, Max_Nucleus_Size);	
	}
	else {
		Dialog.addSlider("Mininum nucleus diameter ("+unit+")", 1, 1000, Min_Nucleus_Size);
		Dialog.addSlider("Maximum nucleus diameter ("+unit+")", 1, 1000, Max_Nucleus_Size);	
	}
	*/
	Dialog.addNumber("Nucleus diameter between", Min_Nucleus_Size, 0, 2, "and");
	Dialog.setInsets(-28, 70, 0);
	Dialog.addNumber("", Max_Nucleus_Size, 0, 2, unit);
	Dialog.addCheckbox("Automatic nuclei segmentation, with threshold bias", auto_threshold_nuclei);
	Dialog.setInsets(-22, 150, 0);
	Dialog.addNumber("", threshold_nuc_bias, 1, 4, "(-"+pow(2,bits)-1+" to "+pow(2,bits)-1+")");
	Dialog.addCheckbox("If automatic segmentation, use local thresholding (adjust parameters in the source code)", auto_local_th_nuc);
	Dialog.addCheckbox("Exclude nuclei on edge of image", exclude_edges);
	Dialog.addCheckbox("Shrink segmented nuclei with", shrink_nuclei);
	Dialog.setInsets(-22, 20, 22);
	Dialog.addNumber("", shrink_nuclei_size, 1, 4, unit);
	Dialog.setInsets(-22, 20, 0);
	Dialog.addCheckbox("Manually edit segmented nuclei?", manually_edit_nuclei);
	Dialog.setInsets(15, 0, 5);
	if(slices>1) Dialog.addChoice("Projection method for foci z-stack", foci_projection_method_array, foci_projection_method);
	Dialog.addNumber("Estimated foci diameter", foci_diameter, 1, 2, "pixels");
	Dialog.addNumber("Minimum foci diameter", min_foci_size, 1, 2, "pixels");
	Dialog.addNumber("Signal-to-Noise Ratio of foci", SNR, 1, 2, "");
	Dialog.addCheckbox("Correct foci background noise", foci_background_correct);
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Minimum threshold for foci", manual_foci_threshold, 0, 4, "(0  to "+pow(2,bits)-1+")");
	//Dialog.setInsets(0,20,20);
	Dialog.addMessage("Additional Options");
	//Dialog.addCheckbox("Measure nucleus intensity", choice_nucleus_intensity);
	Dialog.addCheckbox("Measure mean intensity (of projection) in channel", choice_additional_nucleus_channel);
	Dialog.setInsets(-23, 125, 0);
	Dialog.addNumber("", ch_measure_additional, 0, 1, "");
	Dialog.addCheckbox("Exclude certain nuclei in channel", choice_exclude_nucleus_channel);
	Dialog.setInsets(-23, 125, 0);
	Dialog.addNumber("", ch_exclude, 0, 1, "");
	Dialog.setInsets(0, 0, 0);
	Dialog.addChoice("with signal", choice_exclude_threshold_array, choice_exclude_threshold);
	Dialog.setInsets(-22, 125, 0);
	Dialog.addNumber("threshold", exclude_threshold, 0, 1, "gray levels"); //, pow(2,bits)-1, exclude_threshold);
	//Dialog.setInsets(0,20,0);
	Dialog.addCheckbox("Verbose mode (show intermediate screens during processing)", verbose);
	//Dialog.addHelp(url) For later use
Dialog.show;
	ch_nuclei=Dialog.getNumber();					//0
	ch_foci=Dialog.getNumber();						//1

	if(slices>1) nuclei_method = Dialog.getChoice();//21	//only when it is a z-stack
	Min_Nucleus_Size = Dialog.getNumber();			//2
	Max_Nucleus_Size = Dialog.getNumber();			//3
	auto_threshold_nuclei = Dialog.getCheckbox();	//4
	threshold_nuc_bias = Dialog.getNumber();		//5
	auto_local_th_nuc = Dialog.getCheckbox();		//6
	exclude_edges = Dialog.getCheckbox();			//7
	shrink_nuclei = Dialog.getCheckbox();			//8
	shrink_nuclei_size = Dialog.getNumber();		//9
	manually_edit_nuclei = Dialog.getCheckbox();	//10

	if(slices>1) foci_projection_method = Dialog.getChoice();//24	//only when it is a z-stack
	if (slices>1 && foci_projection_method=="best slice for every nucleus") focus_nuclei=true;
	foci_diameter = Dialog.getNumber();				//11
	min_foci_size = Dialog.getNumber();				//12
	SNR = Dialog.getNumber();						//13
	foci_background_correct=Dialog.getCheckbox();	//14
	manual_foci_threshold = Dialog.getNumber();		//15
	
	choice_additional_nucleus_channel = Dialog.getCheckbox();//16
	ch_measure_additional = Dialog.getNumber();		//17
	if(choice_additional_nucleus_channel==true) if(ch_measure_additional<1 || ch_measure_additional>channels) exit("channel "+ch_measure_additional+" for measuring intensity does not exist. Exiting macro...");
	choice_exclude_nucleus_channel = Dialog.getCheckbox();	//18
	ch_exclude = Dialog.getNumber();				//19
	if (choice_exclude_nucleus_channel==true) if(ch_exclude<1 || ch_exclude>channels) exit("channel "+ch_exclude+" for exclusion of nuclei does not exist. Exiting macro...");
	choice_exclude_threshold = Dialog.getChoice();	//20
	exclude_threshold = Dialog.getNumber();			//21

	verbose=Dialog.getCheckbox();					//22

if(manual_foci_threshold==0) manual_foci_threshold=1;	//needs to be at least 1

sigma_large = 3*foci_diameter;	//size of Gaussian blur filter (large) for foci background subtraction
sigma_small = 0.2*foci_diameter;				//size of Gaussian blur filter (small) for foci background subtraction

//Define save paths
newdir= dir+"\\results_ch"+ch_foci+"\\";
if(!File.exists(newdir)) File.makeDirectory(newdir);
extension_index = indexOf(file_name, extension)-1;				//index of file extension
results_file = newdir+"\\"+substring(file_name,0,extension_index)+"_results_foci.tsv";	//name of results file
merged_image_file = newdir+"\\"+File.nameWithoutExtension+"_analyzed";
foci_statistics_file = newdir+"\\"+File.nameWithoutExtension+"_foci_statistics";
foci_statistics_ordered_file = newdir+"\\"+File.nameWithoutExtension+"_foci_statistics_ordered";
ROIs_file = newdir+"\\"+File.nameWithoutExtension+"_ROIs.zip";
log_file = newdir+"\\log.txt";


//---------SAVE SETTINGS IN CONFIG FILE
save_config_file();

//enquire if all files in the directory should be analyzed if the directory contains more than one files with the same extension.
if(image_list.length>1){
	analyze_all=getBoolean("Analyze all "+image_list.length+" "+extension+" files in this directory with these settings?");
	if(analyze_all==true) pause_after_file=getBoolean("Pause after each file?");
}

start_time=getTime();

//START OF DO...WHILE LOOP FOR ANALYZING ALL IMAGES IN A DIRECTORY
do{

roiManager("Reset");
run("Clear Results");

if(analyze_all==true) {
	run("Close All");
	file_name = image_list[current_image_nr];	//retrieve file name from image list
	Ext.openImagePlus(dir+file_name);		//open file using LOCI Bioformats plugin
	//convert to uint16 in case of int16 calibration (y=-32768+x)
	calibration = calibrate(0);	//retreive calibration (weird method, but ok)
	if(calibration!=0)
	run("Calibrate...", "function=None");	
	run("Subtract...", "value="+-calibration+" stack");
	resetMinAndMax();
	if(remove_pixel_calibration==true) run("Properties...", "unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1");
	rename(file_name);
	print("current image: "+file_name);
	extension_index = indexOf(file_name, extension)-1;				//index of file extension
	results_file = newdir+"\\"+substring(file_name,0,extension_index)+"_results_foci.tsv";	//name of results file
	merged_image_file = newdir+"\\"+substring(file_name,0,extension_index)+"_analyzed";//name of analyzed image
	ROIs_file = newdir+"\\"+substring(file_name,0,extension_index)+"_ROIs.zip";			//name of file with ROIs
	Stack.getDimensions(width, height, channels, slices, frames);
	current_image_nr++;
}



//////////////////////////////////////////
// TEMPORARY STATEMENT ABOUT PIXEL SIZE //
//////////////////////////////////////////
//run("Properties...", "channels=2 slices=1 frames=1 unit=micron pixel_width=0.2145240 pixel_height=0.2145240 voxel_depth=0.2500000 frame=[1 sec]");



if(verbose==false) setBatchMode(true);

run("Split Channels");
selectWindow("C"+ch_nuclei+"-"+file_name);

//---------METHODS FOR NUCLEI DETECTION
if (slices==1) {};
else if (nuclei_method=="automatically select slice with sharpest nuclei") {
	run("Duplicate...", "title=edges_nuclei duplicate");
	run("Set Measurements...", "area mean standard median stack limit redirect=None decimal=3");
	run("Find Edges", "stack");
	var stdev_array = newArray(slices);
	for(i=0;i<slices;i++) {
		setSlice(i+1);
		run("Measure");
		stdev_array[i] = getResult("StdDev",i);	//measure stddev of edges image; sharpest image has largest stddev.
	}
	Array.getStatistics(stdev_array, min, max);
	index_max = newArray();
	index_max = indexOfArray(stdev_array, max);
	close();
	selectWindow("C"+ch_nuclei+"-"+file_name);
	setSlice(index_max[0]+1);	//select slice with largest standard deviation
	run("Duplicate...", "title=nuclei");
	message("selected slice: "+index_max[0]+1);
}
else if (nuclei_method=="manually select slice") {
	selectWindow("C"+ch_nuclei+"-"+file_name);
	setBatchMode("show");
	waitForUser("Select slice for analysis of nuclei");
	if (verbose==false) setBatchMode("hide");
	current_slice = getSliceNumber();
	setSlice(current_slice);
	run("Duplicate...", "title=nuclei");
}
else if (nuclei_method=="maximum intensity projection") {
	run("Z Project...", "projection=[Max Intensity]");
}
else if (nuclei_method=="sum slices") {
	run("Z Project...", "projection=[Sum Slices]");
}
rename("nuclei");



//---------MAIN PART

segment_nuclei();
if (manually_edit_nuclei==true) edit_ROIs("nuclei");
selectWindow("segmented_nuclei");
if (shrink_nuclei==true) {
	for(i=0;i<roiManager("count");i++) {
		roiManager("Select",i);
		run("Enlarge...", "enlarge=-"+shrink_nuclei_size);
		roiManager("Rename", "nucleus_"+i+1);
		roiManager("Update");
	}
}

var mean_intensity_exclude_channel = newArray(roiManager("count"));
if(choice_exclude_nucleus_channel == true) exclude_nuclei();

//create list of ROI coordinates
var ROI_x = newArray(roiManager("count"));	//containers for selection locations
var ROI_y = newArray(roiManager("count"));
for(i=0;i<roiManager("count");i++) {
	roiManager("Select",i);
	getSelectionBounds(x, y, ROI_width, ROI_height);
	ROI_x[i]=x;
	ROI_y[i]=y;
}

if(nuclei_found==true) {
	selectWindow("C"+ch_foci+"-"+file_name);
	// Quickfix for Luuk: duplicate so that you can measure intensity in the foci channel as well
	run("Duplicate...", "title=C"+ch_foci+"-"+file_name+" duplicate");
	message("Foci before background subtraction");
	if (slices>1) {
		z_project("C"+ch_foci+"-"+file_name, foci_projection_method);
		rename("foci");
	}
	if (focus_nuclei==true) {	//then also find best slice for the nuclei, but taking the foci channel as input. It's calculating it again, but is fast ayway.
		z_project("C"+ch_nuclei+"-"+file_name, foci_projection_method);
		rename("nuclei_slice_corresponding_with_foci");
	}
	else rename("foci");

	//data containers
	var nucleus_area = newArray(roiManager("count"));
	var nucleus_sum_intden = newArray(roiManager("count"));
	var nucleus_sum_stddev = newArray(roiManager("count"));
	if(focus_nuclei==true) var nucleus_slice_intensity = newArray(roiManager("count"));
	if(focus_nuclei==true) var nucleus_slice_stddev = newArray(roiManager("count"));
	var	nucleus_circ = newArray(roiManager("count"));
	var	nucleus_round = newArray(roiManager("count"));
	var	nucleus_solidity = newArray(roiManager("count"));

	var mean_intensity_additional_channel = newArray(roiManager("count"));		//mean intensity in the additional channel
	var total_intensity_additional_channel = newArray(roiManager("count"));		//total intensity in the additional channel
	var mean_foci_intensity_additional_channel = newArray(roiManager("count"));	//mean intensity in the foci mask in the additional channel
	var total_foci_intensity_additional_channel = newArray(roiManager("count"));	//total intensity in the foci mask in the additional channel

	var foci_count = newArray(roiManager("count"));			//number of foci per nucleus
	var foci_median = newArray(roiManager("count"));		//background signal of each nucleus
	var foci_stddev = newArray(roiManager("count"));		//standard deviation of the background of each nucleus

	//data containers per nucleus, up to a maximum of 10000 foci per nucleus.
	var foci_area = newArray(9999);					//area per focus
	var foci_mean_intensity = newArray(9999);		//mean intensity per focus
	var foci_integrated_density = newArray(9999);	//total intensity per focus
	var foci_distance = newArray(9999);		 		//distance of foci to nuclear perimeter

	var avg_foci_area = newArray(9999);				//mean of foci area in nucleus
	var stddev_foci_area = newArray(9999);			//stddev of foci area in nucleus
	var stderr_foci_area = newArray(9999);			//std error of foci area in nucleus
	var avg_foci_intensity = newArray(9999);		//mean of foci intensity in nucleus
	var stddev_foci_intensity = newArray(9999);		//stddev of foci intensity in nucleus
	var stderr_foci_intensity = newArray(9999);		//std error of foci intensity in nucleus
	var avg_foci_intden = newArray(9999);			//mean of total intensity per focus in nucleus (=area*intensity)
	var stddev_foci_intden = newArray(9999);		//stddev of total intensity per focus in nucleus
	var stderr_foci_intden = newArray(9999);		//std error of total intensity per focus in nucleus
	var total_foci_intensity = newArray(9999);		//integrated foci intensity in nucleus
	var stddev_total_foci_intensity = newArray(9999);//stddev of integrated foci intensity in nucleus
	var stderr_total_foci_intensity = newArray(9999);//std error of integrated foci intensity in nucleus
	var avg_foci_distance = newArray(9999);			//mean of distance of foci to nuclear perimeter
	var stddev_foci_distance = newArray(9999);		//syddev of distance of foci to nuclear perimeter
	var stderr_foci_distance = newArray(9999);		//stderr of distance of foci to nuclear perimeter

	get_nuclei_area_and_intensity();
	if(choice_additional_nucleus_channel == true) measure_additional_channel();
	selectWindow("foci");
	if(foci_background_correct == true) background_subtract("foci");
	else run("Duplicate...", "title=foci_original");	//Creating a duplicate is necessary because the macro needs this image later
	selectWindow("foci");
	get_mean_and_stddev();
	
	roiManager("Show None");

	run("Hyperstack...", "title=foci_statistics type=32-bit display=Color width="+roiManager("count")+" height=1 channels=4 slices=1 frames=1");

	if(speed_mode==true) {
		newImage("foci Maxima", "32-bit black", width, height, 1);
	}
	if (verbose==true) setBatchMode(true); //always use batch mode here
	analyze_foci();
	if (verbose==true) setBatchMode(false);

//	print("Applied threshold: "+threshold);
	selectWindow("foci_statistics");
	if(verbose==true) {
		setBatchMode("show");
		run("In [+]");
		run("In [+]");
		run("In [+]");
	}

	display_merged_image();
	if(verbose==false) cleanup();
//setBatchMode(true);
	order_foci_statistics();
	handle_results();

	if(verbose==false) {
		close("foci_statistics");
		close("foci_statistics_ordered");
	}				
	if(verbose==false) setBatchMode(false);
	//show_distributions();
	if(verbose==false) setBatchMode(true);
}

if(analyze_all==true && current_image_nr<image_list.length) {
	if(pause_after_file==true) {
		if (verbose==false) setBatchMode(false);
		selectWindow(merged_image);
		waitForUser("Inspect results and click OK to continue with the next file");
		if (verbose==false) setBatchMode(true);
	}
	if(pause_after_file==true) run("Close All");
	nuclei_found = true; //reset for next loop
}
if(analyze_all==true && current_image_nr==image_list.length) {
	selectWindow("Log");
	saveAs("text",log_file);
	run_time=round((getTime()-start_time)/1000);
	showMessage("Finished analyzing "+image_list.length+" files in "+run_time+" seconds.");
}
if(analyze_all==false) {
	run_time=round((getTime()-start_time)/1000);
	print("Finished in "+run_time+" seconds.");
}

//END OF DO...WHILE LOOP
} while(analyze_all==true && current_image_nr<image_list.length);

restoreSettings();


/*
 * Other approach - to be implemented
 * ----------------------------------
 * difference of Gaussians (in 3D)
 * add slices below and above (necessary for objects 'touching' the upper or lower plane
 * 3D object counter
 * * run("3D Objects Counter", "threshold=17 slice=5 min.=16 max.=12582912 exclude_objects_on_edges centres_of_masses statistics summary");
 * For inspection: merge max intensity projection of binary of centroids with extended depth image of original/background-subtracted foci stack
 * Table gives lots of info!
 * 
 */


//---------FUNCTIONS

function segment_nuclei() {
	run("Duplicate...", "title=segmented_nuclei duplicate");
	//run("Properties...", "channels=[] slices=[] frames=[] unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000 frame=0 origin=0,0");
	selectWindow("segmented_nuclei");
	run("32-bit");
	run("Remove Outliers...", "radius="+outliers_diameter+" threshold="+outliers_threshold+" which=Bright");
	run("Enhance Contrast", "saturated=0.35");
	if(logarithm==true) run("Log", "stack");
	selectWindow("segmented_nuclei");
//waitForUser("after median");
	run("Duplicate...", "title=nuclei_before_segmentation duplicate");
	resetMinAndMax();
	selectWindow("segmented_nuclei");
	resetMinAndMax();
	if(auto_threshold_nuclei==true)	{
//waitForUser("before background subtraction");
		if(bgsubtr_nuclei_segmentation == true) run("Subtract Background...", "rolling="+rolling_ball_radius);
//waitForUser("before 8-bit");
		if(auto_local_th_nuc==true) {
			run("Conversions...", "scale");
//setMinAndMax(0,4095);
			run("8-bit");
//waitForUser("before median filter");
			run("Median...", "radius="+median_radius_nuclei);
//waitForUser("after median filter");
			run("Auto Local Threshold", "method=Mean radius="+alth_radius+" parameter_1="+alth_parameter_1+" parameter_2=0 white");
			//setBatchMode("show");
		}
		else {
//waitForUser("before median filter");
			run("Median...", "radius="+median_radius_nuclei);
//waitForUser("after median filter");
			setAutoThreshold("Otsu dark");
//waitForUser("thresholded");
			getThreshold(min,max);
			resetThreshold();
			setThreshold(min+threshold_nuc_bias,max);
			run("Convert to Mask", "background=Dark black");
		}
	}
	else {
		if(bgsubtr_nuclei_segmentation == true) run("Subtract Background...", "rolling="+rolling_ball_radius);
		run("Median...", "radius="+median_radius_nuclei);
		setAutoThreshold("Otsu dark");
		getThreshold(min,max);
		resetThreshold();
		setThreshold(min+threshold_nuc_bias,max);
		min_old=min;
		selectWindow("segmented_nuclei");
		setBatchMode("show");
		run("Threshold...");
		setAutoThreshold("Otsu dark");
		selectWindow("Threshold");
		waitForUser("Set threshold for segmentation of nuclei and press OK");
		selectWindow("segmented_nuclei");
		if(verbose==false) setBatchMode("hide");
		getThreshold(min,max);
		threshold_nuc_bias=d2s(min-min_old,1);
		store_bias = getBoolean("Store threshold bias ("+threshold_nuc_bias+") in config file for later use?");
		if (store_bias==true) save_config_file();
		run("Convert to Mask", "  black");
	}
	run("Fill Holes");
	if(watershed==true) run("Watershed");
	setThreshold(127, 255);
	//run("Set Measurements...", "area mean standard centroid stack redirect=None decimal=3");
	if(exclude_edges==true) run("Analyze Particles...", "size="+Min_Nucleus_Size*Min_Nucleus_Size*PI/4+"-"+Max_Nucleus_Size*Max_Nucleus_Size*PI/4+" circularity=0.33-1.00 show=Nothing display exclude add");
	else run("Analyze Particles...", "size="+Min_Nucleus_Size*Min_Nucleus_Size*PI/4+"-"+Max_Nucleus_Size*Max_Nucleus_Size*PI/4+" circularity=0.33-1.00 show=Nothing display add");
	if(roiManager("count")==0) {
		nuclei_found=false;
		print("No nuclei found!");
	}
	else print(roiManager("count")+" nuclei found");
	resetThreshold();
//	run("Duplicate...", "title=distance_map");
//	run("Distance Map");
}


function exclude_nuclei() {
	selectWindow("C"+ch_exclude+"-"+file_name);
	if(slices>1) run("Z Project...", " projection=[Max Intensity]");
	rename("exclude_channel");
	run("Set Measurements...", "mean redirect=None decimal=3");
	run("Clear Results");
	ROIs_to_delete = newArray(roiManager("count"));

	j=0;
	for(i=0;i<roiManager("count");i++) {
		roiManager("select",i);
		//List.setMeasurements();
		//mean_intensity_exclude_channel[i] = List.getValue("Mean");
		run("Measure");
		intensity = getResult("Mean",i);
		//mean_intensity_exclude_channel[i] = getResult("Mean",i);
		//print("length: "+mean_intensity_exclude_channel.length);
		if(verbose==true) print("Mean intensity in exclude channel "+ch_exclude+": "+intensity);
		if(choice_exclude_threshold=="below") {
			if(intensity < exclude_threshold) {
				ROIs_to_delete[j]=i;
				j++;
			}
			else mean_intensity_exclude_channel[i-j] = intensity;	//measure the signal in the non-deleted nuclei
		}
		else if(choice_exclude_threshold=="above") {
			if(intensity > exclude_threshold) {
				ROIs_to_delete[j]=i;
				j++;
			}
			else mean_intensity_exclude_channel[i-j] = intensity;	//measure the signal in the non-deleted nuclei
		}
	}
	ROIs_to_delete = Array.trim(ROIs_to_delete, j);
	if(ROIs_to_delete.length>0) {
		roiManager("Select",ROIs_to_delete);
		roiManager("Delete");
	}
	print(ROIs_to_delete.length+" nuclei are being excluded. "+roiManager("Count")+" nuclei will be analyzed.");
	if(roiManager("Count")==0) nuclei_found=false;
	run("16-bit");	//for merging later
}


function get_nuclei_area_and_intensity() {
	if(slices>1) selectWindow("C"+ch_nuclei+"-"+file_name);
	else selectWindow("nuclei");
	run("Set Measurements...", "area mean shape standard integrated redirect=None decimal=3");
	if(slices>1) {
		run("Z Project...", " projection=[Sum Slices]");
		rename("nuclei_summed");
	}
	for(i=0;i<roiManager("count");i++) {
		roiManager("Select", i);
		List.setMeasurements();
		nucleus_area[i] = List.getValue("Area");
		nucleus_sum_intden[i]=List.getValue("IntDen")/10000;
		nucleus_sum_stddev[i]=List.getValue("StdDev");
		nucleus_circ[i]=List.getValue("Circ.");
		nucleus_round[i]=List.getValue("Round");
		nucleus_solidity[i]=List.getValue("Solidity");
	}
	run("Select None");
	if(focus_nuclei==true) {
		selectWindow("nuclei_slice_corresponding_with_foci");
		for(i=0;i<roiManager("count");i++) {
			roiManager("Select", i);
			List.setMeasurements();
			nucleus_slice_intensity[i] = List.getValue("Mean");
			nucleus_slice_stddev[i] = List.getValue("StdDev");
		}
	}
	run("Select None");
}


function measure_additional_channel() {
	selectWindow("C"+ch_measure_additional+"-"+file_name);
	if(slices>1) run("Z Project...", " projection=[Sum Slices]");	//Before version 4.2 this was 'Max Projection'
	rename("additional_channel");
//	run("16-bit");
	run("Set Measurements...", "mean redirect=None decimal=3");
	for(i=0;i<roiManager("count");i++) {
		roiManager("select",i);
		setPercentileThreshold(lowerPercentile, upperPercentile);
		List.setMeasurements("limit");
		mean_intensity_additional_channel[i] = List.getValue("Mean");
		total_intensity_additional_channel[i] = List.getValue("IntDen");
		//print("Mean intensity in channel "+ch_measure_additional+": "+mean_intensity_additional_channel[i]);
	}
}

//Because of local SNR check in 'find maxima' dim large 'foci' at the perimeter of the cell are picked more frequently
//Therefore an average nucleus-wide 'global' SNR threshold is used to measure the foci.
function get_mean_and_stddev() {
		//get median and stddev values of the foci channel (excluding foci (outliers)) of current nucleus
		selectWindow("foci");

		if (verbose==true) setBatchMode(true); //always use batch mode here
		run("Enhance Contrast", "saturated=0.2");

		showStatus("calculating background and standard deviation in all nuclei...");
if(foci_background_correct == false) {
		run("Duplicate...", "title=foci_outliers_removed duplicate");
		//remove outliers (foci) and then threshold to retrieve median and stddev of background in nucleus
		setBatchMode("show");	//must be visible, or remove outliers doesn't work for unknown reason.
		run("Remove Outliers...", "radius="+1.5*foci_diameter+" threshold="+outliers_threshold+" which=Bright slice");	//to get a better estimate of the stddev

		if(verbose == false) setBatchMode("hide");
}
		run("Set Measurements...", "mean standard median limit redirect=None decimal=3");
		run("Clear Results");
/*
		for(i=0;i<roiManager("count");i++) {
			selectWindow("foci_outliers_removed");
			roiManager("Select", i);
			setAutoThreshold("Otsu");			//or else try ("mean")
			//List.setMeasurements();	//List.setMeasurements cannot be limited to the thresholded pixels
			run("Measure");
			run("Select None");
			resetThreshold();
			//foci_median[i]=List.getValue("Median");	//background signal of each nucleus in the foci channel
			//foci_stddev[i]=List.getValue("StdDev"); 	//stddev of the background in each nucleus
			//foci_median[i]=getResult("Median");			//background signal of each nucleus in the foci channel
			foci_median[i]=getResult("Mean");			//background signal of each nucleus in the foci channel
			foci_stddev[i]=getResult("StdDev"); 		//stddev of the background in each nucleus
			print("nucleus "+i+1+": "+foci_median[i]+" +- "+foci_stddev[i]);
		}
*/
//Replaced by Multi measure (much faster)
		roiManager("Select All");
		roiManager("multi-measure");
		for(i=0;i<roiManager("count");i++) {
			foci_median[i]=getResult("Mean",i);	//Mean because Median after background subtraction is not integer any more
			foci_stddev[i]=getResult("StdDev",i);
		}
		run("Select None");
		resetThreshold();
		Array.getStatistics(foci_median, median_min, median_max, median_mean, median_stdDev);
		Array.getStatistics(foci_stddev, stddev_min, stddev_max, stddev_mean, stddev_stdDev);
		//Array.print(foci_median);
		//Array.print(foci_stddev);
		avg_median = median_mean;	//put in global variable
		avg_stddev = stddev_mean;	//put in global variable		
//		Array.sort(foci_stddev);
//		avg_stddev = foci_stddev[floor(foci_stddev.length/2)];	//take the median standard deviation
		print("Median nuclear background intensity in foci channel: "+d2s(avg_median,1)+" +- "+d2s(median_stdDev,1));
		print("Mean standard deviation of the background: "+d2s(avg_stddev,1)+" +- "+d2s(stddev_stdDev,1));

		//retreive background (offset) value (mean of all frames)
		showStatus("retrieving background signal from image...");
		run("Duplicate...", "title=background duplicate");
		run("Set Measurements...", "mean standard median limit redirect=None decimal=3");
		setAutoThreshold("Li bright");
		run("Measure");	//List.setMeasurements cannot be limited to the thresholded pixels
		background=getResult("Median");
		showStatus("");
		print("Median background signal (outside nuclei): "+d2s(background,1));
		run("Select None");
		selectWindow("foci_outliers_removed");
		print("Calculated 'average threshold' for foci: ("+SNR+" * "+d2s(avg_stddev,1)+") + "+d2s(background,1)+" = "+d2s(SNR*avg_stddev+background,1));
		if (verbose==true) setBatchMode(false);
		//resetThreshold();
}


function background_subtract(image) { //method: Difference of Gaussians (good for small spots)
	showStatus("subtracting background from foci...");
	selectWindow(image);
	
	if (verbose==true) setBatchMode(true); //always use batch mode here
	setBatchMode("show");
	run("Enhance Contrast", "saturated=0.2");
	run("Duplicate...", "title=foci_outliers_removed duplicate");
	run("Remove Outliers...", "radius="+1.5*foci_diameter+" threshold="+outliers_threshold+" which=Bright slice");	//to get a better estimate of the stddev
	if(verbose == false) setBatchMode("hide");

	run("Duplicate...", "title="+image+"_large_blur duplicate range=[]");
	run("Gaussian Blur...", "sigma="+sigma_large+" stack");
	selectWindow(image);
	rename(image+"_original");
	run("Duplicate...", "title="+image+"_small_blur duplicate range=[]");
	run("Gaussian Blur...", "sigma="+sigma_small+" stack");
	imageCalculator("Subtract stack", image+"_small_blur",image+"_large_blur");
	selectWindow(image+"_large_blur");
	run("Close");
	selectWindow(image+"_small_blur");
	rename(image);
	run("Enhance Contrast", "saturated=0.1");
	getMinAndMax(min, max);
	setMinAndMax(0, max);
	showStatus("");
	message("background-subtracted image");
	setBatchMode("show");
}


function z_project(stack, method) {
	selectWindow(stack);
	if (method=="maximum intensity projection") {
		run("Z Project...", " projection=[Max Intensity]");
	}
	else if (method=="sum slices") {
		run("Z Project...", " projection=[Sum Slices]");
	}
	else if (method=="smart z-stack projection (recommended for widefield)") {
		Extended_Depth_of_Field();
		selectWindow("C"+ch_foci+"-"+file_name+"_focused");
	}
	else if (method=="best slice for every nucleus") {
		run("Duplicate...", "title=best_slices");
		run("Select All");
		setBackgroundColor(0,0,0);
		run("Clear", "slice");	//create black image but retain pixel size, units, bitdepth, etc.
		focus_image = getTitle();
		if (verbose==true) setBatchMode(true); //always use batch mode here
		setBatchMode("show");
		for(i=0;i<roiManager("count");i++) {
			get_best_z_slice(i, "C"+ch_foci+"-"+file_name, stack, focus_image, "StdDev", 0);	//stack does not need to be twice the same, e.g. "C"+ch_foci+"-"+file_name
		}
		run("Select None");
		if (verbose==true) setBatchMode(false);
	}
}


function get_best_z_slice(nucleus_nr, image_measure, image_copy_from, image_paste_to, measure_type, start_index) {
	selectWindow(image_measure);
	roiManager("Select", start_index+nucleus_nr);
	var array = newArray(slices);

	//measure mean intensity and standard deviation of every slice and return the largest value
	for(j=0;j<slices;j++) {
		setSlice(j+1);
		List.setMeasurements();
		array[j] = List.getValue(measure_type);	//best slice has largest mean/stdDev/...
	}
	Array.getStatistics(array, min, max);
	index_max = newArray();
	index_max = indexOfArray(array, max);

	//print("cell "+i+": slice "+index_max[0]+1);
	showStatus("cell "+i+": slice "+index_max[0]+1);

	//copy the (whole) nucleus to the new image
	selectWindow(image_copy_from);
	setSlice(index_max[0]+1);		//select correct slice
	resetMinAndMax();
	roiManager("Select", nucleus_nr);
	run("Copy");
	selectWindow(image_paste_to);
	roiManager("Select", nucleus_nr);
	run("Paste");
}


function analyze_foci() {
	selectWindow("foci");
//	setBatchMode("show");
	setBatchMode("hide");	//Hide the image, because it is faster for large images
	
	run("Set Measurements...", "  area mean standard min integrated redirect=None decimal=3");
	message("starting analysis of "+roiManager("count")+" nuclei...");
	
	current_image_height = 1;	//height of foci_statistics (=max nr of foci in movie)
	
	//find mean and integrated intensity of foci in current nucleus
	for(i=0;i<roiManager("count");i++) {
		showStatus("detecting foci in nucleus "+i+"/"+roiManager("count"));
		showProgress(i/roiManager("count"));

		//get rid of surrounding nuclei
		selectWindow("foci");
		roiManager("Select", i);
		if(speed_mode==false) {
			run("Create Mask");
			rename("Mask_nucleus_"+i+1);
			run("Divide...", "value=255");
			run("16-bit");			//make mask 16-bit

			imageCalculator("Multiply", "Mask_nucleus_"+i+1, "foci" );
		}
		else {
			//create distance map for current nucleus
			selectWindow("segmented_nuclei");
			roiManager("Select", i);
			run("Duplicate...", "title=distance_map_nucleus_"+i+1);
//			run("Distance Map");
			//create foci image for current nucleus
			selectWindow("foci");
			roiManager("Select", i);
			run("Duplicate...", "title=foci_nucleus_"+i+1);
			run("Duplicate...", "title=Mask_nucleus_"+i+1);	//not used, only to prevent crashes later on

			selectWindow("Mask_nucleus_"+i+1);	//this is actually the foci image
		}
		//set Threshold for foci
//		if(foci_background_correct == true) threshold=maxOf(SNR*avg_stddev,manual_foci_threshold);
//		else threshold=maxOf(SNR*avg_stddev+background,manual_foci_threshold);
		
		avg_threshold = SNR*avg_stddev+background;
		if(threshold_scaling_coeff != -1) {
			if(foci_stddev[i]-avg_stddev>0) sqrt_diff = pow((SNR*foci_stddev[i]-avg_stddev),threshold_scaling_coeff);
			else sqrt_diff = -pow(abs(SNR*foci_stddev[i]-avg_stddev),threshold_scaling_coeff);
			if(foci_background_correct == true) threshold=maxOf(avg_threshold+sqrt_diff,manual_foci_threshold);
			else threshold=maxOf(avg_threshold+sqrt_diff+background,manual_foci_threshold);
		}
		else {
			if(foci_background_correct == true) threshold=maxOf(avg_threshold,manual_foci_threshold);
			else threshold=maxOf(avg_threshold+background,manual_foci_threshold);
		}
		setThreshold(threshold,65535);	//set foci threshold as a number times the average STD of all nuclei in this image, or the manually defined minimum
		if(verbose==true && threshold_scaling_coeff != -1) print("Diff: "+sqrt_diff+", Threshold "+i+1+": "+threshold);
//setBatchMode("show");
//waitForUser("after thresholding "+i+", with threshold at "+SNR*avg_stddev);
		run("Find Maxima...", "noise=0 output=[Segmented Particles] above");
		//First time 'Analyze Particles' to remove foci that are smaller than the size threshold
		setThreshold(1,65535);
		run("Analyze Particles...", "size="+PI/4*min_foci_size*min_foci_size+"-Infinity pixel circularity=0.00-1.00 show=Masks");
		rename("Maxima_nucleus_"+i+1);
		roiManager("Select",i);
		setSelectionLocation(0,0);
		setBackgroundColor(255, 255, 255);
		run("Clear Outside");
		setBackgroundColor(0, 0, 0);
		run("Divide...", "value=255");	//create mask
		run("16-bit");			//make mask 16-bit
//setBatchMode("show");
//waitForUser("Maxima_nucleus_"+i+1);
		if(speed_mode==false) imageCalculator("Multiply", "Maxima_nucleus_"+i+1,"foci");	//multiply mask with foci image
		else {
//			roiManager("Select",i);
//			setSelectionLocation(0,0);
			run("Copy");
			run("Select None");
			//create foci distance to nuclear perimeter image
			imageCalculator("Multiply create", "Maxima_nucleus_"+i+1,"distance_map_nucleus_"+i+1);
			rename("foci_distance_map_nucleus_"+i+1);
//setBatchMode("show");
//waitForUser("distance map "+i+1);
			//create foci intensity image
			imageCalculator("Multiply", "Maxima_nucleus_"+i+1,"foci_nucleus_"+i+1);
		}
		rename("detected_foci_nucleus_"+i+1);
//setBatchMode("show");
//waitForUser("detected foci "+i);
		close("Mask_nucleus_"+i+1);

		//Second time 'Analyze Particles' to get info on foci in nucleus (Note BvdB: should be doable in a single pass, right?)
		selectWindow("detected_foci_nucleus_"+i+1);
		run("Properties...", "channels=1 slices=1 frames=1 unit="+unit+" pixel_width="+pw+" pixel_height="+ph+" voxel_depth="+pd);	//set to original units/values in stead of pixels
		setThreshold(threshold,65535);	//include all pixels above the threshold

//		run("Analyze Particles...", "size="+PI/4*min_foci_size*min_foci_size+"-Infinity pixel circularity=0.00-1.00 show=Nothing display clear");
		run("Analyze Particles...", "size="+PI/4*min_foci_size*min_foci_size+"-Infinity pixel circularity=0.00-1.00 show=Masks display clear");
		resetThreshold();

//***** Stuff voor Inge de Krijger: measure the additional channel only in the detected foci
		if (choice_additional_nucleus_channel==true) {
			run("Invert LUT");
	//setBatchMode("show");
			run("Divide...", "value=255");
			selectWindow("additional_channel");
			roiManager("select",i);
			run("Duplicate...", "title=additional_channel_nucleus_"+i+1);
	//setBatchMode("show");
			imageCalculator("Multiply 32-bit", "additional_channel_nucleus_"+i+1,"Mask of detected_foci_nucleus_"+i+1);
			setThreshold(1, 65535);
			run("NaN Background");
			List.setMeasurements;
			close();
	//setBatchMode("show");
			mean_foci_intensity_additional_channel[i] = List.getValue("Mean");
			total_foci_intensity_additional_channel[i] = List.getValue("IntDen");
		}
//*****


		//write foci info to foci statistics image
		foci_count[i]=nResults;		//number of detected foci in current nucleus

		//Array.sort(foci_area);
		//Array.sort(foci_mean_intensity);
		//Array.sort(foci_integrated_density);

		selectWindow("foci_statistics");
		if(nResults>current_image_height) {
			run("Canvas Size...", "width="+roiManager("count")+" height="+nResults+" position=Top-Center zero");
			current_image_height=nResults;
		}
		Stack.setChannel(1);		//channel 1: foci area
		for(j=0;j<nResults;j++) {
			//setResult("Mean",j, getResult("Mean",j)-background);	//currently disabled because foci are measured on background-subtracted foci image 
			foci_area[j] = getResult("Area",j);
			setPixel(i,j,foci_area[j]);
		}
		Stack.setChannel(2);		//channel 2: mean intensity
		for(j=0;j<nResults;j++) {
			foci_mean_intensity[j] = getResult("Mean",j);
			setPixel(i,j,foci_mean_intensity[j]);
		}
		Stack.setChannel(3);		//channel 3: integrated intensity (=area*(raw)mean)
		for(j=0;j<nResults;j++) {
			foci_integrated_density[j] = getResult("RawIntDen",j);
			setPixel(i,j,foci_integrated_density[j]);
		}
/*
		//Third time 'Analyze Particles' on foci distance map to retreive distances to nuclear perimeter
		selectWindow("foci_distance_map_nucleus_"+i+1);
		setThreshold(1,65535);
		run("Analyze Particles...", "size="+PI/4*min_foci_size*min_foci_size+"-Infinity pixel circularity=0.00-1.00 show=Nothing display clear");
		resetThreshold();

		//write distance info to foci statistics image
		selectWindow("foci_statistics");
		Stack.setChannel(4);		//channel 4: distance to nuclear perimeter
		for(j=0;j<nResults;j++) {
			foci_distance[j] = getResult("Mean",j);
			setPixel(i,j,foci_distance[j]*pw);
		}
		Stack.setChannel(3);
*/		
		//Calculate averages and stddevs of foci in current nucleus
		temp_array = Array.slice(foci_area,0,nResults);	//perform statistics only over the non-zero entries
		Array.getStatistics(temp_array, area_min, area_max, area_mean, area_stdDev);
		avg_foci_area[i] = area_mean;
		stddev_foci_area[i] = area_stdDev;
		stderr_foci_area[i] = area_stdDev/sqrt(nResults);

		temp_array = Array.slice(foci_mean_intensity,0,nResults);
		Array.getStatistics(temp_array, intensity_min, intensity_max, intensity_mean, intensity_stdDev);
		avg_foci_intensity[i] = intensity_mean;
		stddev_foci_intensity[i] = intensity_stdDev;
		//stderr_foci_intensity[i] = intensity_stdDev/sqrt(nResults);
/*
		temp_array = Array.slice(foci_distance,0,nResults);
		Array.getStatistics(temp_array, distance_min, distance_max, distance_mean, distance_stdDev);
		avg_foci_distance[i] = distance_mean*pw;
		stddev_foci_distance[i] = distance_stdDev*pw;
*/
		// The following outcommented part somehow gives slightly different results for integrated density than what it should be.
		
		//temp_array = Array.slice(foci_integrated_density,0,nResults);
		//Array.getStatistics(temp_array, intden_min, intden_max, intden_mean, intden_stdDev);
		//avg_foci_intden[i] = intden_mean;
		//stddev_foci_intden[i] = intden_stdDev;
		//stderr_foci_intden[i] = intden_stdDev/sqrt(nResults);

		// This is easier and faster
		avg_foci_intden[i] = intensity_mean*area_mean;

		total_foci_intensity[i] = avg_foci_intden[i]*foci_count[i];		//total foci intensity in nucleus 
		//stddev_total_foci_intensity[i] = stddev_foci_intden[i]*foci_count[i];	//stddev of total foci intensity in nucleus 
		//stderr_total_foci_intensity[i] = stderr_foci_intden[i]*foci_count[i];	//std error of total foci intensity in nucleus 

		if(i>0) {
			if(speed_mode==false) imageCalculator("Add", "detected_foci_nucleus_"+i+1, "detected_foci_nucleus_"+i);
			close("detected_foci_nucleus_"+i);	//close previous nucleus
		}

		if(speed_mode==true) {
			selectWindow("foci Maxima");
			roiManager("Select", i);
			run("Paste");
		}

		//intermediate cleanup
		close("Mask_nucleus_"+i+1+" Segmented");
		close("Mask of detected_foci_nucleus_"+i+1);
		close("foci_nucleus_"+i+1);
		close("distance_map_nucleus_"+i+1);
		close("foci_distance_map_nucleus_"+i+1);
		if (choice_additional_nucleus_channel==true) close("additional_channel_nucleus_"+i+1);
	}
//	if(speed_mode==false) {
		selectWindow("detected_foci_nucleus_"+i);
		rename("detected_foci");
//	}
//	else {
//	}
	
	resetThreshold();
	
	selectWindow("foci_statistics");
	Stack.setChannel(1);
	resetMinAndMax();
	run("Magenta Hot");
	Stack.setChannel(2);
	run("Enhance Contrast", "saturated=0.35");
	run("Orange Hot");
	Stack.setChannel(3);
	resetMinAndMax();
	run("Cyan Hot");
	Stack.setChannel(4);
	resetMinAndMax();
	run("Red Hot");
	Stack.setChannel(3);

	if(speed_mode==false) {
		selectWindow("detected_foci");
		setThreshold(1,65535);
		run("Create Mask");
		resetThreshold();
		rename("foci Maxima");
		selectWindow("foci");
		setOption("Show All", true);
	}
}


function display_merged_image() {
	selectWindow("foci Maxima");
	run("Select None");
	run("Multiply...", "value=127");
	setMinAndMax(0, 127);
	run("16-bit");
	selectWindow("foci_original");
	run("16-bit");
	selectWindow("nuclei");
	run("16-bit");
	if(focus_nuclei==true) {
		selectWindow("nuclei");
		rename("nuclei_projection_used_in_segmentation");	//as used in segmentation
		setOption("Show All", true);
		setBatchMode("show");
		selectWindow("nuclei_slice_corresponding_with_foci");
		run("16-bit");
		rename("nuclei");
	}
	if (choice_additional_nucleus_channel==true) {
		selectWindow("additional_channel");
		run("16-bit");
	}
	if (choice_exclude_nucleus_channel==true && choice_additional_nucleus_channel==false) run("Merge Channels...", "c1=[foci Maxima] c2=foci_original c3=nuclei c4=exclude_channel gray=*None* create keep");
	else if (choice_exclude_nucleus_channel==false && choice_additional_nucleus_channel==true) run("Merge Channels...", "c1=[foci Maxima] c2=foci_original c3=nuclei c4=additional_channel gray=*None* create keep");
	else if (choice_exclude_nucleus_channel==true && choice_additional_nucleus_channel==true) run("Merge Channels...", "c1=[foci Maxima] c2=foci_original c3=nuclei c4=additional_channel c5=exclude_channel gray=*None* create keep");
	else run("Merge Channels...", "c1=[foci Maxima] c2=foci_original c3=nuclei gray=*None* create keep");
	getDimensions(width, height, channels2, slices2, frames2);
	run("Properties...", "channels="+channels2+" slices="+slices2+" frames="+frames2+" unit="+unit+" pixel_width="+pw+" pixel_height="+ph+" voxel_depth="+pd+" origin=0,0");	//set to original units/values in stead of pixels
	if (choice_exclude_nucleus_channel==true && choice_additional_nucleus_channel==false) {
		Stack.setChannel(4);
		run("Grays");
	}
	if (choice_exclude_nucleus_channel==false && choice_additional_nucleus_channel==true) {
		Stack.setChannel(4);
		run("Magenta");
	}
	if (choice_exclude_nucleus_channel==true && choice_additional_nucleus_channel==true) {
		Stack.setChannel(4);
		run("Magenta");
		Stack.setChannel(5);
		run("Grays");
	}
	Stack.setChannel(1);
	run("Red");
	setMinAndMax(0, 127);
	Stack.setChannel(2);
	run("Green");
	resetMinAndMax();
	run("Enhance Contrast", "saturated=0.05");
	Stack.setChannel(3);
	run("Blue");
	Stack.setChannel(2);
	//resetMinAndMax();

	Stack.setActiveChannels("11100");	//do not show the other channels
	
	rename("Merged final image");
	
	setBatchMode("show");
	roiManager("Show All with labels");
	run("Channels Tool...");
}

function cleanup() {
	if(clean==true) {
		for(i=1;i<=channels;i++) {
			if(isOpen("C"+i+"-"+file_name)) {
				close("C"+i+"-"+file_name);
			}
		}
		close("nuclei");
		close("foci");
		close("foci_original");
	}
}


function order_foci_statistics() {
	selectWindow("foci_statistics");
	run("Duplicate...", "title=foci_statistics_ordered duplicate");
	setBackgroundColor(0, 0, 0);
	run("Select All");
	run("Clear", "stack");
	if(verbose==true) setBatchMode("show");
	selectWindow("foci_statistics");
	Stack.setChannel(3);	//total foci intensity

	rankPosArr = Array.rankPositions(foci_count);	//determine rank array of nuclei in terms of nr of foci
	ranks = Array.rankPositions(rankPosArr);

	for(i=1;i<=4;i++) {
		foci_statistics_pixel_operations(i, ranks);
	}
	
	selectWindow("foci_statistics");
	run("Flip Vertically", "stack");
	
	selectWindow("foci_statistics_ordered");
	run("Select None");
	run("Rotate... ", "angle=180 grid=1 interpolation=None stack");
	//run("Log", "stack");
	//run("Divide...", "value="+log(10)+" stack");
	run("In [+]");
	run("In [+]");
	run("In [+]");
	//resetMinAndMax();
	Stack.setChannel(1);				//area
	resetMinAndMax();
	run("Enhance Contrast", "saturated=1");
	Stack.setChannel(2);				//intensity
	setMinAndMax(threshold, threshold*3);
	//run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(3);				//integrated density
	setMinAndMax(0, threshold*100);
}


function foci_statistics_pixel_operations(channel, ranks) {
	selectWindow("foci_statistics_ordered");
	Stack.setChannel(channel);
	selectWindow("foci_statistics");
	Stack.setChannel(channel);
	for(x=0;x<foci_count.length;x++) {		//for every nucleus (x)
		selectWindow("foci_statistics");
		foci_int = newArray(foci_count[x]);		//create array to hold foci intensity
		for(y=0;y<foci_count[x];y++) {			//for every row (y)
			foci_int[y] = getPixel(x,y);			//get intensity of nucleus x, focus y
		}
		foci_int_sorted = Array.sort(foci_int);		//sort array
		Array.reverse(foci_int);
		selectWindow("foci_statistics_ordered");
		for(y=0;y<foci_count[x];y++) {
			setPixel(ranks[x],y,foci_int_sorted[y]);//write intensity x,y in the ordered image
		}
	}
}

function handle_results() {
	run("Clear Results");
	run("Input/Output...", "file=.tsv use_file copy_row copy_column save_column");	//do not copy row numbers, create tsv file
	for(i=0;i<roiManager("Count");i++) {
		setResult("nucleus_nr", i, i+1);
		setResult("nuc_area ("+unit+"^2)", i, nucleus_area[i]);
		setResult("nuc_total_int (x10^4)", i, nucleus_sum_intden[i]);
		setResult("nuc_total_stddev", i, nucleus_sum_stddev[i]);
		if(slices>1 && focus_nuclei==true) setResult("nuc_slice_intensity", i, nucleus_slice_intensity[i]);
		if(slices>1 && focus_nuclei==true) setResult("nuc_slice_stddev", i, nucleus_slice_stddev[i]);
		setResult("nuc_circ",i, nucleus_circ[i]);
		setResult("nuc_round",i, nucleus_round[i]);
		setResult("nuc_solidity",i, nucleus_solidity[i]);

		setResult("nr_of_foci", i, foci_count[i]);
		setResult("mean_foci_intensity", i, avg_foci_intensity[i]);
		setResult("stddev_foci_intensity", i, stddev_foci_intensity[i]);
		setResult("total_intensity ", i, total_foci_intensity[i]);
		setResult("mean_foci_area ("+unit+"^2)", i, avg_foci_area[i]);
		setResult("stddev_foci_area", i, stddev_foci_area[i]);
		setResult("area_coverage (%)", i, avg_foci_area[i]*foci_count[i]/nucleus_area[i]*100);
//		setResult("mean_foci_distance ("+unit+")", i, avg_foci_distance[i]);
//		setResult("stddev_foci_distance", i, stddev_foci_distance[i]);
		if (choice_additional_nucleus_channel==true) setResult("mean_int_ch"+ch_measure_additional, i, mean_intensity_additional_channel[i]);
		if (choice_additional_nucleus_channel==true) setResult("mean_foci_int_ch"+ch_measure_additional, i, mean_foci_intensity_additional_channel[i]);
		if (choice_additional_nucleus_channel==true) setResult("total_int_ch"+ch_measure_additional, i, total_intensity_additional_channel[i]);
		if (choice_additional_nucleus_channel==true) setResult("total_foci_int_ch"+ch_measure_additional, i, total_foci_intensity_additional_channel[i]);

				resetThreshold();

		if (choice_exclude_nucleus_channel==true) setResult("mean_int_ch"+ch_exclude, i, mean_intensity_exclude_channel[i]);
	}
	updateResults();

	Array.getStatistics(foci_count, min, max, mean, stdDev);
	print(""+mean*foci_count.length+" foci detected in "+roiManager("Count")+" nuclei (average foci per nucleus: "+mean+")");
	print("\n");
	
	selectWindow("Results");
	saveAs("results", results_file);
	//File.close(results_file); doesn't work(?)
	
	selectWindow("foci_statistics");
	saveAs("tiff", foci_statistics_file);
	selectWindow("foci_statistics_ordered");
	saveAs("tiff", foci_statistics_ordered_file);
	selectWindow("Merged final image");
	saveAs("tiff", merged_image_file);
	merged_image=getTitle();
	roiManager("Select All");
	roiManager("Save", ROIs_file);
}


function show_distributions() {
	if(analyze_all==false) {
		//run("Distribution...", "parameter=[nuc_area ("+unit+"^2)] or=15 and=Min_Nucleus_Size-Max_Nucleus_Size");
		run("Distribution...", "parameter=[nuc_total_int (x10^4)] or=40 and=0-20");
		Plot.getValues(xpoints, ypoints);
		minima = Array.findMinima(ypoints, 10);
		print("Minima found at positions:");
		for(n=0;n<minima.length;n++) {
			print(minima[n]/2);
		}
	//	run("Distribution...", "parameter=foci_count or=10 and=0-0");
	//	run("Distribution...", "parameter=avg_max_int or=10 and=0-0");
	//	run("Distribution...", "parameter=total_int or=10 and=0-0");
	}
}


function indexOfArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (d2s(array[a],3)==d2s(value,3)) {
			count++;
		}
	}
	if (count>0) {
		indices=newArray(count);
		count=0;
		for (a=0; a<lengthOf(array); a++) {
			if (array[a]==value) {
				indices[count]=a;
				count++;
			}
		}
		return indices;
	}
}


function message(txt) {
if (pause==true){
	waitForUser(txt);
	}
}


function save_config_file() {
	config_file = File.open(tempdir+"\\foci_macro_config.txt");
	print(config_file, ch_nuclei);
	print(config_file, ch_foci);
	print(config_file, Min_Nucleus_Size);
	print(config_file, Max_Nucleus_Size);
	print(config_file, auto_threshold_nuclei);
	print(config_file, threshold_nuc_bias);
	print(config_file, auto_local_th_nuc);
	print(config_file, exclude_edges);
	print(config_file, shrink_nuclei);
	print(config_file, shrink_nuclei_size);
	print(config_file, manually_edit_nuclei);
	print(config_file, foci_diameter);
	print(config_file, min_foci_size);
	print(config_file, SNR);
	print(config_file, foci_background_correct);
	print(config_file, manual_foci_threshold);
	print(config_file, choice_additional_nucleus_channel);
	print(config_file, ch_measure_additional);
	print(config_file, choice_exclude_nucleus_channel);
	print(config_file, ch_exclude);
	print(config_file, choice_exclude_threshold);
	print(config_file, exclude_threshold);
	print(config_file, verbose);
	if(slices>1) print(config_file, nuclei_method);
	if(slices>1) print(config_file, foci_projection_method);
	File.close(config_file);
}


//Manual editing of detected nuclei
function edit_ROIs(image1) {
	shift=1;
	ctrl=2; 
	rightButton=4;
	alt=8;
	leftButton=16;
	insideROI = 32;
	
	flags=-1;
	//x2=-1; y2=-1; z2=-1; flags2=-1;
	
	selectWindow(image1);
	roiManager("Show All without labels");
	setOption("DisablePopupMenu", true);
	setBatchMode(true);
	resetMinAndMax();
	setBatchMode("show");
	color_ROIs();
	print("\\Clear");
	print("Delete, combine and draw new ROIs. \n- CTRL + leftclick deletes a ROI.\n- SHIFT + leftclick selects ROIs, rightclick merges selected ROIs. \n- Draw new ROIs using the magic wand or the freehand tool and press 't' to add. \n- Press space bar when finished editing.\n");
	showMessage("Delete, combine and draw new ROIs.","\n- CTRL + leftclick deletes a ROI.\n- SHIFT + leftclick selects ROIs, rightclick merges selected ROIs. \n- Draw new ROIs using the magic wand or the freehand tool and press 't' to add. \n- Press space bar when finished editing.\n\n\nThis information is also printed in the log window.");
	print("Starting editing "+roiManager("count")+" ROIs...");
	
	setTool("freehand");
	roiManager("Show All Without Labels");
	setOption("DisablePopupMenu", true);
	
	nROIs = roiManager("Count");
	
	while(!isKeyDown("space")) {		//exit by pressing space bar
		getCursorLoc(x, y, z, flags);
		//print("\\Update:"+flags);
		if(flags==17 || flags==18)	{	//(de)select multiple ROIs with shift-leftclick; delete ROI with rightclick
			for(i=0;i<roiManager("Count");i++) {
				roiManager("Select",i);
				if(Roi.contains(x, y)==true) {
				selected = Roi.getProperty("selected");
					//click to select a single ROI
					if(flags==17 && selected==false) {		//select ROI
						//print("selecting ROI "+i);
						Roi.setStrokeColor("red");
						Roi.setProperty("selected",true);
					}
					else if(flags==17 && selected==true) {	//shift-leftclick: deselect ROI
						//print("deselecting ROI "+i);
						Roi.setStrokeColor("cyan");
						//Roi.setFillColor("1900ffff");
						Roi.setProperty("selected",false);
					}
					else if(flags==18) {	//ctrl-leftclick: delete ROI
						roiManager("Delete");
						for(j=0;j<roiManager("Count");j++) {	//deselect all ROIs and rename
							roiManager("Select",j);
							roiManager("Rename", "ROI "+j);
						}
					}
				}
			}
			roiManager("Deselect");
			run("Select None");
			updateDisplay();
		}
	
		if(flags==4) {	//right button: combine selected ROIs
			selected_ROI_array = newArray(roiManager("Count"));	//create array with indices of selected ROIs
			j=0;
			for(i=0;i<roiManager("Count");i++) {
				roiManager("select",i);
				selected = Roi.getProperty("selected");
				if(selected==true) {
					selected_ROI_array[j] = i;
					j++;
					//print(j);
				}
			}
			//check if more than one ROI is selected. If yes, combine the selected ROIs and update the list
			selected_ROI_array = Array.trim(selected_ROI_array,j);
			//print(selected_ROI_array.length + " ROIs selected");
			if(selected_ROI_array.length > 1) {
				roiManager("Select",selected_ROI_array);
				roiManager("Combine");
				roiManager("Update");
				to_delete_array = Array.copy(selected_ROI_array);								//selecting and deleting redundant ROIs
				to_delete_array = Array.slice(selected_ROI_array,1,selected_ROI_array.length);	//create array without the first element
				roiManager("Deselect");
				roiManager("select", to_delete_array);
				roiManager("Delete");
				roiManager("Select",selected_ROI_array[0]);
				run("Enlarge...", "enlarge=1 pixel");			//remove wall between ROIs by enlarging and shrinking with 1 pixel
				run("Enlarge...", "enlarge=-1 pixel");
				roiManager("Update");
				
				setKeyDown("none");
				
				color_ROIs();
			}
		}
	
	
		if(nROIs!=roiManager("Count")) {	//change in the number of ROIs 
			run("Select None");
			color_ROIs();
			nROIs = roiManager("Count");
		}
	
		else wait(50);
	}	//end of while loop
	
	//Deselect and rename all ROIs once more
	color_ROIs();
}


function color_ROIs() {
	run("Remove Overlay");

	for(j=0;j<roiManager("Count");j++) {	//fill all ROIs
		roiManager("Select",j);
		roiManager("Rename", "ROI "+j+1);
		Roi.setProperty("selected",false);
		//Roi.setFillColor("1900ffff");	//10% cyan fill
	}
	roiManager("Deselect");
	if(roiManager("count")>0) run("From ROI Manager");	//Add overlay containing the ROI fill
	roiManager("Select All");
	roiManager("Set Color", "cyan");
	roiManager("Deselect");
	roiManager("Show All");
	updateDisplay();
}


function setPercentileThreshold(lowerPercentile, upperPercentile) {
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	
	lowerTotal = 0;
	upperTotal = nPixels;
	lower=0;
	while (lowerTotal < nPixels*lowerPercentile) {
		lowerTotal += histogram[lower];
		lower++;
	}
	upper=histogram.length-1;
	while (upperTotal > nPixels*upperPercentile) {
		upperTotal -= histogram[upper];
		upper--;
	}
	setThreshold(lower,upper);
}



function Extended_Depth_of_Field() {

	radius=3;
	
	//Get start image properties
	w=getWidth();
	h=getHeight();
	d=nSlices();
	source=getImageID();
	origtitle=getTitle();
	rename("tempnameforprocessing");
	sourcetitle=getTitle();

	//Generate edge-detected image for detecting focus
	run("Duplicate...", "title=["+sourcetitle+"_Heightmap] duplicate range=1-"+d);
	heightmap=getImageID();
	heightmaptitle=getTitle();
	run("Find Edges", "stack");
	run("Maximum...", "radius="+radius+" stack");

	//Alter edge detected image to desired structure
	run("32-bit");
	for (x=0; x<w; x++) {
		showStatus("Creating focused image from stack...");
		showProgress(x/w);
		for (y=0; y<h; y++) {
			slice=0;
			max=0;
			for (z=0; z<d; z++) {
				setZCoordinate(z);
				v=getPixel(x,y);
				if (v>=max) {
					max=v;
					slice=z;
				}
			}
			for (z=0; z<d; z++) {
				setZCoordinate(z);
				if (z==slice) {
					setPixel(x,y,1);
				} else {
					setPixel(x,y,0);
				}
			}
		}
	}
	run("Gaussian Blur...", "sigma="+radius+" stack");
	
	//Generation of the final image
	
	//Multiply modified edge detect (the depth map) with the source image
	run("Image Calculator...", "image1="+sourcetitle+" operation=Multiply image2="+heightmaptitle+" create 32-bit stack");multiplication=getImageID();
	//Z project the multiplication result
	run("Z Project...", "start=1 stop="+d+" projection=[Sum Slices]");
	//Some tidying
	rename(origtitle+"_focused");
	selectImage(heightmap);
	close();
	selectImage(multiplication);
	close();
	selectImage(source);
	rename(origtitle);
	
}