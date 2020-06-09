if(nImages>0) run("Close All");
roiManager("Reset");
run("Clear Results");
run("Set Measurements...", "area mean min integrated limit redirect=None decimal=3");

analyzedSuffix = "_analyzed.tif";

var nrNuclei = 0;
var total_foci = 0;
var total_coloc = 0;
var current_nucleus = 0;

dir_1 = getDirectory("Choose the result folder containing the first analyzed foci files.");
dir_2 = getDirectory("Choose the result folder containing the second analyzed foci files.");
ch1 = File.getName(dir_1);
ch2 = File.getName(dir_2);

dirParent = File.getParent(dir_1);
colocDir = dirParent+File.separator+"coloc_results";
print(dirParent);
print(colocDir);
if(!File.exists(colocDir)) File.makeDirectory(colocDir);

file_list = getFileList(dir_1); //get filenames of directory 1

//make a list of images with 'extension' as extension.
j=0;
image_list=newArray(file_list.length);	//Dynamic array size doesn't work on some systems, so first make image_list the maximal size and then trim.
for(i=0; i<file_list.length; i++){
	if (matches(file_list[i],".*analyzed.*")) {
		image_list[j] = file_list[i];
		j++;
	}
}
image_list = Array.trim(image_list, j);	//Trimming the array of images


print("\\Clear");
print("Directory contains "+file_list.length+" files, of which "+image_list.length+" analyzed files.");
print("\n");

//colocTable = "Colocalization results";
//Table.create(colocTable);
// (in results table)

setBatchMode(true);
for(f=0;f<image_list.length;f++) {

	run("Close All");

	nrNuclei = 0;
	total_foci = 0;
	total_coloc = 0;
	current_nucleus = 0;
	roiManager("Reset");
	run("Clear Results");
	print("-----------------------------------");
	print("file: "+image_list[f]);
	ROIs_file = replace(image_list[f],"_analyzed.tif","_ROIs.zip");
	roiManager("Open", dir_1+ROIs_file);
	roiManager("Show None");
	nrNuclei = roiManager("count");
	
	//Compute overlap in the images
	open(dir_1+image_list[f]);
	rename(image_list[f]+"_1");
	run("Remove Overlay");
	getPixelSize(unit, pw, ph, pd);
	if(unit == "micron" || unit == "microns") unit = "µm";
	
	run("Duplicate...", "duplicate channels=1");
	rename("foci_mask_1");
	run("Grays");
	setThreshold(64,255);
	run("Analyze Particles...", "add");
	run("32-bit");
	setThreshold(64,255);
	run("NaN Background");
	resetThreshold();
	List.setMeasurements();
	area_1 = List.getValue("Area");
	setMinAndMax(0, 127);
	
	open(dir_2+image_list[f]);
	rename(image_list[f]+"_2");
	run("Remove Overlay");
	run("Duplicate...", "duplicate channels=1");
	rename("foci_mask_2");
	run("Grays");
	run("32-bit");
	setThreshold(64,255);
	run("NaN Background");
	resetThreshold();
	List.setMeasurements();
	area_2 = List.getValue("Area");
	setMinAndMax(0, 127);
	total_foci = roiManager("count") - nrNuclei;
	
	//merge chanels for visualization
	run("Merge Channels...", "c1=foci_mask_1 c2=foci_mask_2 create keep ignore");
	rename(image_list[f]+"_colocalized_foci");
	Stack.setChannel(2);
	run("Green");
	Stack.setChannel(1);
	run("Red");

	Stack.setChannel(2);
	run("Add Slice", "add=channel");
	Stack.setChannel(3);
	run("Blue improved");
	setBatchMode("show");

	close(image_list[f]+"_1");
	close(image_list[f]+"_2");
	
	nrColoc = newArray(nrNuclei);
	nrNotColoc = newArray(nrNuclei);
	fociColocArea = newArray(nrNuclei);
	Array.fill(nrColoc,0);
	Array.fill(nrNotColoc,0);
	Array.fill(fociColocArea,0);

	//Loop over all foci of ch1 to check whether it colocalizes with any foci pixels in ch2
	for(i=nrNuclei;i<roiManager("count");i++) {
		selectWindow("foci_mask_2");
		roiManager("Select", i);
		Roi.getBounds(x, y, width, height);
		//run("Measure");
		List.setMeasurements();
		max = List.getValue("Max");
		area = List.getValue("Area");
		for(j=0;j<nrNuclei;j++) {
			roiManager("Select", j);
			if(Roi.contains(x+width/2, y+height/2)) current_nucleus = j;	//determine to which nucleus the focus belongs
		}
		if (max == 127) {
			nrColoc[current_nucleus]+=1;	//colocalized ch1 foci in current nucleus
			total_coloc+=1;					//all colocalized ch1 foci in the image
			fociColocArea[current_nucleus]+=area;
		}	
		else {
			nrNotColoc[current_nucleus]+=1;
		}
	}
	
	//Colorize nuclei (white) by number of colocalized foci and save data to results window
	run("Clear Results");
	for(i=0;i<nrNuclei;i++) {
		selectWindow(image_list[f]+"_colocalized_foci");
		roiManager("Select", i);
		color = nrColoc[i]/(nrColoc[i]+nrNotColoc[i]);
		changeValues(-1,1,color);

		//Measure total foci area for ch1
		selectWindow("foci_mask_1");
		roiManager("Select", i);
		List.setMeasurements();
		totalFociArea_ch1 = List.getValue("Area");

		//Measure total foci area for ch2
		selectWindow("foci_mask_2");
		roiManager("Select", i);
		List.setMeasurements();
		totalFociArea_ch2 = List.getValue("Area");


		setResult("cell nr", i, i+1);
		setResult("total nr foci "+ch1+"", i, nrColoc[i]+nrNotColoc[i]);
		setResult("nr foci coloc "+ch1+"->"+ch2+"", i, nrColoc[i]);
		setResult("coloc % "+ch1+"->"+ch2+"", i, 100*nrColoc[i]/(nrColoc[i]+nrNotColoc[i]) );				//% of coloc "+ch1+" foci with ch2
		setResult("total foci area "+ch1+" ("+unit+"^2)", i, totalFociArea_ch1);						//Area of all "+ch1+" foci
		setResult("total foci area "+ch2+" ("+unit+"^2)", i, totalFociArea_ch2);						//Area of all "+ch1+" foci
		setResult("coloc foci area "+ch1+" / total_"+ch1+" ("+unit+"^2)", i, fociColocArea[i]);				//Area of "+ch1+" foci that colocalizes with foci in ch2
		setResult("coloc foci area % "+ch1+" / total_"+ch1+"", i, 100*fociColocArea[i]/totalFociArea_ch1);		//Area% of colocalized foci in ch1 (over total ch1 foci area)
		setResult("coloc foci area % "+ch1+" / total_"+ch2+"", i, 100*fociColocArea[i]/totalFociArea_ch2);		//Area% of colocalized foci in ch1 (over total ch2 foci area)
	}
	selectWindow(image_list[f]+"_colocalized_foci");
	setMinAndMax(0,1);
	run("Select None");
	
	//calculate colocalized area of all foci
	imageCalculator("AND create", "foci_mask_1","foci_mask_2");
	rename("Overlap_mask");
	run("Grays");
	run("Select None");
	setThreshold(64,255);
	List.setMeasurements("limit");
	area_1_AND_2 = List.getValue("Area");
	close("Overlap_mask");
	close("foci_mask_1");
	close("foci_mask_2");

	//Save results
	saveAs("Results", colocDir + File.separator + substring(image_list[f], 0, lengthOf(image_list[f])-lengthOf(analyzedSuffix))+"_coloc_results.tsv");
	run("RGB Color");
	saveAs("Tiff",  colocDir + File.separator + substring(image_list[f], 0, lengthOf(image_list[f])-lengthOf(analyzedSuffix))+"_coloc_results_(RGB).tif");
	close();
	
	//print results
	print(""+total_coloc+" / "+total_foci+" "+ch1+" foci (red) show colocalization with "+ch2+" foci (green) ("+d2s(total_coloc/total_foci*100,1)+" %)");
	print("Total area foci "+ch1+": " +area_1+ " "+unit+"^2");
	print("Total area foci "+ch2+": " +area_2+ " "+unit+"^2");
	print("Area of overlap: " +area_1_AND_2+ " "+unit+"^2");
	print("Area% of colocalization (2->1): " + d2s(area_1_AND_2/area_2*100,1) + "%");
	print("Area% of colocalization (1->2): " + d2s(area_1_AND_2/area_1*100,1) + "%");

}