// Morphology analysis of pericytes after whole-cell patch clamp  
// 
// This script analyzes the morphology of pericytes injected with a fluorescent dye after whole-cell patch clamp. 
// The script asks the user to draw a region of interest around the patching pipette to exclude it from analysis.
// It also asks the user to draw another ROI around the background resulting from the presence of the pipette 
// named pipette shadow in the script. The user also provides manually an estimate of the maximum intensity
// of the pipette shadow. This ROI drawn by the user is used as a mask to reduce the background caused by the pipette.
// Then a series of filtering, and automated thresholding are applied. The user is asked to confirm the thresholding
// or adjust it if not satisfactory. 
// Then object filtering based on size is applied to extract a single cell as a binary image. Then skeleton analysis 
// is run on a MIP image.
// The code output is: 
// 		6 Tiff files: MIP of the original image, morphology in 2D, morphology in 3D, and 3 images for skeleton analysis      
//		3 CSV files: one for morphology parameters, one for skeleton parameters, one for the manually inserted values.
//		2 ROIs: for the pipette and its shadow, saved for reference and reproducibility. 
//		 
// PREREQUISITES: MorphoLibJ (PLEASE INSTALL THIS LIBRARY BEFORE RUNNING)
// 
// Author: Majed Kikhia, majed.kikhia@charite.de
// June 2024
// License: BSD3
// 
// Copyright 2024 Majed Kikhia, Charité – Universitätsmedizin Berlin, Department of Experimental Neurology  
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//

// Before starting, please set an output folder. Then open an image and run the code:  

// Setting an output folder:

output_folder = "/.../.../.../output"

// configure that binary image are black in background, objects are white
setOption("BlackBackground", true);

// Set foreground color as black and background color as wight
setForegroundColor(0, 0, 0);
setBackgroundColor(255, 255, 255);

// Removing the images from the second channel which has faint signal
getDimensions(img_width, img_height, img_channels, img_slices, img_frames);

setSlice(1);
slice1_max = getValue("Max");

setSlice(2);
slice2_max = getValue("Max");

if (slice1_max > slice2_max) {
	first_slice = 2;
} else { 
	first_slice = 1;
}

getDimensions(img_width, img_height, img_channels, img_slices, img_frames); // Get the new image dimensions after removing the second channel to use them in future steps

// Set the voxel size based on the metadata text file
Stack.setXUnit("micron");
run("Properties...", "channels=1 slices=" + img_slices + " frames=1 pixel_width=0.166667 pixel_height=0.166667 voxel_depth=0.3");

// Update image name, dimensions, and save voxel size
original_img = getTitle();
getDimensions(img_width, img_height, img_channels, img_slices, img_frames);
getVoxelSize(width , height , depth , unit);

// Input the end of the pipette in the Z dimension 
Dialog.createNonBlocking("Please check at which Z plane the pipette ends");
Dialog.addNumber("The pipette ends at Z plane", 70);
Dialog.show();

pipette_z_level =Dialog.getNumber();

// Make a maximum intensity projection and ask the user to draw a polygon around the pipette to exclude it 
run("Z Project...", "projection=[Max Intensity]");
saveAs("Tiff", output_folder + "/" + original_img + "_MIP.tif");

setTool(2);
waitForUser("Draw a region areound the pipette. Be careful around the cell edges!");

if (selectionType==-1) exit("you must draw a region around the pipette first");
// Restore the selection on the original image
selectWindow(original_img);
run("Restore Selection");

// Save the pipette selection to the ROI manager 
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "pipette");
roiManager("Save", output_folder + "/" + original_img + "_pipette_ROI.roi");
roiManager("Select", 0);
roiManager("Deselect");
roiManager("Delete");

// Delete the content of the of selected area until the end of the pipette in Z 
selectWindow(original_img);
for (n=1; n < pipette_z_level + 1; n++) {
	setSlice(n);
	fill();
}

run("Select None");

// Ask the user to set manually a threshold based on the intesity of the pipette shadow

Dialog.createNonBlocking("Please enter The maximum intensity of pipette shadow");
Dialog.addNumber("Maximum intensity of pipette shadow", 500);
Dialog.show();

pipette_shadow =Dialog.getNumber();

setColor(pipette_shadow);

// Creating a mask to reduce the shaddow of the pipette

setTool(2);
waitForUser("Draw a region areound the pipette shadow or the part that you want to substract\n\nand click OK afterward. You can choose to spare or reduce brightness of cell segments with this mask!");

if (selectionType==-1) exit("you must draw a region around the pipette shadow first");

newImage("pipette_mask", "16-bit black", img_width, img_height, img_slices);

run("Restore Selection");

for (n=1; n < img_slices + 1; n++) {
	setSlice(n);
	fill();
}

// Save the pipette shadow selection to the ROI manager
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "pipette_shadow");
roiManager("Save", output_folder + "/" + original_img + "_pipette_shadow_ROI.roi");
roiManager("Select", 0);
roiManager("Deselect");
roiManager("Delete");

run("Select None");

// Substract the mask from the orginal image
imageCalculator("Subtract stack", original_img,"pipette_mask");
run("Select None");

// Then you have the image to work with

// Enhance and the filter image. Take care of saturation levels. Probably 0.03 then threshold with RenyiEntropy
selectWindow(original_img);
run("Gaussian Blur...", "sigma=1 stack");
run("Subtract Background...", "rolling=25 stack");
run("Enhance Contrast...", "saturated=0.08 normalize process_all use");
run("Gaussian Blur...", "sigma=1 stack");

run("8-bit");
enhanced_img = getTitle();

// Thresholding and creating a binary image
run("Threshold...");
setAutoThreshold("Li dark stack"); // Other alogorithms suc as MaxEntropy, and RenyiEntropy can be tried 
setOption("BlackBackground", true);

waitForUser("Are you satisfied with the thresholding? If not, please adjust minimally and press OK");
getThreshold(lower, upper);

// Creating a table to save the manaually entered inputs 
newTableName = "Results table with new name";

Table.create(original_img + "_user_input");
Table.set("Filename", 0, original_img);
Table.set("Pipette ends at Z level", 0, pipette_z_level);
Table.set("Value of pipette shadow mask", 0, pipette_shadow);
Table.set("The threshold applied for segmentation", 0, lower);
Table.save(output_folder + "/" + original_img + "_user_input.csv");

// Applying threshold
run("Convert to Mask", "method=Li background=Dark black create");
binary_img = getTitle(); // This is a binary image that we will work with

// Run connected component analysis and keep the largest object as the single cell  
run("Connected Components Labeling", "connectivity=6 type=[16 bits]");

run("Keep Largest Label");

setVoxelSize(width, height, depth, unit);

// Analyse the morphology of the isolated single cell
run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity euler_number equivalent_ellipsoid ellipsoid_elongations max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
single_cell = getTitle();

// Save the numeric morphology results as csv file
saveAs("Results", output_folder + "/" + original_img + "_morph.csv");
close(original_img + "_morph.csv");

// Save the binary z-stack
run("Green");
saveAs("Tiff", output_folder + "/" + original_img + "_binary_3D.tif");
single_cell = getTitle();

run("Duplicate...", "duplicate");

setVoxelSize(width, height, depth, unit);

// Skeleton analysis
run("Z Project...", "projection=[Max Intensity]");
run("Skeletonize");

// Save the skeleton analysis output
saveAs("Tiff", output_folder + "/" + original_img + "_binary_skeleton.tif");

run("Analyze Skeleton (2D/3D)", "prune=none calculate");
saveAs("Results", output_folder + "/" + original_img + "_skeleton.csv");	
close("Results");

selectWindow("Tagged skeleton");
saveAs("Tiff", output_folder + "/" + original_img + "_tagged_skeleton.tif");

selectWindow("Longest shortest paths");
saveAs("Tiff", output_folder + "/" + original_img + "_longest_shortest_path.tif");

selectWindow(single_cell);

// Save the STDIP (Standard deviation intensity projection)
setVoxelSize(width, height, depth, unit);
run("Region Boundaries Labeling");
run("Z Project...", "projection=[Standard Deviation]");
run("Green");
saveAs("Tiff", output_folder + "/" + original_img +  "_STDIP.tiff");

// close all images
close("*");
close(original_img + "_user_input");

// empty the results table
run("Clear Results");
