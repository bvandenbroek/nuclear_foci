/*
 * Macro to merge results files.
 * Bram van den Broek, Netherlands Cancer Institute, 2013-2019
 * 
 */

#@ String(label = "String that ends the filenames", value="_ch1.tsv") string
#@ File (label = "Results folder", style = "directory") dir
dir = dir+File.separator;

print("\\Clear");
run("Clear Results");
run("Input/Output...", "save_column");

//Find all result files in directory
//dir = getDirectory("Choose a directory containing the experiment.");
showStatus("reading directory");
file_list = getFileList(dir); //get filenames of directory
showStatus("");

//make a list of images with 'extension' as extension.
j=0;
image_list=newArray(file_list.length);
for(i=0; i<file_list.length; i++){
	if (endsWith(file_list[i],string)) {
		image_list[j] = file_list[i];
		j++;
	}
}
image_list = Array.trim(image_list, j);	//Trimming the array of images
if(image_list.length==0) exit("No results files found");
else print("Experiment contains "+image_list.length+" result files:");

//Get column headers from the first file in the list
text_file = File.openAsString(dir+image_list[0]);
lines = split(text_file,"\n");
headers = split(lines[0],"\t");

//import all results into the results window and save
j=1;
m=0;
for(i=0; i<image_list.length; i++) {
	print("appending file: "+image_list[i]);
	text_file = File.openAsString(dir+image_list[i]);
	lines = split(text_file,"\n");
	for(j=1;j<lines.length;j++) {
		values = split(lines[j],"\t");
//		if(j==1) setResult("filename",m,substring(image_list[i], 0, lengthOf(image_list[i])-12));	//Only print file name at the first occurrance
//		else setResult("filename",m,"");
		setResult("filename",m,substring(image_list[i], 0, lengthOf(image_list[i])-lengthOf(string)));
		for(k=0;k<headers.length;k++) {
			setResult(headers[k],m,parseFloat(values[k]));
		}
		m++;
	}
	updateResults();
}
selectWindow("Results");
saveAs("text", dir+"results_appended.tsv");
print("\nAppended Results File saved as "+dir+"results_appended.tsv");
