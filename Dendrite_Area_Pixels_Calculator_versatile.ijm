//Set directory to save results in table to use for analysis
#@ File (style="directory") imageFolder;
dir = File.getDefaultDir;
dir = replace(dir,"\\","/");

//Get pixel dimensions
setTool("hand");
waitForUser("Get pixel dimensions","Click on the non-annotated input image to get the pixel dimensions\nthen press OK");
getPixelSize(unit, pixelWidth, pixelHeight);
getVoxelSize(width, height, depth, unit);
pixel_area = pixelWidth*pixelHeight;

setTool("multipoint");
waitForUser("Select object and background","click on what you want to measure the area and border length of in the annotated image, then on the background\nand then press OK");
getSelectionCoordinates(xpoints, ypoints);
foreground = getValue(xpoints[0], ypoints[0]);
background = getValue(xpoints[1], ypoints[1]);
setTool("hand");

waitForUser("Starting frame","Put the annotated image at the slice from where the process should start\nthen press OK");

//Get image dimensions
w = getWidth();
h = getHeight();

//Set the currently chosen slice number
slice = getSliceNumber();

//Set headers
print("Slice","Slice_depth_in_microns","Area_in_pixels","Area_in_microns^2","Edge_pixels","edge_length_in_microns");

//Set slice for slice collumn
slice_number = slice;

//Set parameters for circularity normalisation factor retrieval
pcn1 = 7.22032E18;
pcn2 = 1.17768;
pcn3 = 7.57994E-16;
pcn4 = 1.26197 ;

//Set arrays for plots
slice_number_array = newArray();
slice_depth_array = newArray();
area_pixels_array = newArray();
area_micron_array = newArray();
edge_pixel_array = newArray();
edge_length_array = newArray(); 
circularity_array = newArray();
circularity_1_array = newArray();
fourier_frequency_array = newArray();

//Loop for current and next slides until end slide of image
fourier_fequency = 0;
for (i=0;i<nSlices-slice+1;i++){
	//Set pixel counter for Area and edge
	area_pixels = 0;
	edge_pixel = 0;
	edge_length = 0;
	//Loop for every pixel in current slice
	for (x=0;x<w;x++){
		for (y=0;y<h;y++){
			if (getPixel(x,y) == foreground){
				area_pixels += 1;
				if (getPixel(x+1,y) == background){ //Prime
				edge_pixel += 1;
				edge_length += pixelHeight;
					if (getPixel(x-1,y) == background){
					edge_length += pixelHeight;
						if (getPixel(x,y+1) == background){
						edge_length += pixelWidth;
							if (getPixel(x,y-1) == background){
							edge_length += pixelWidth;
							}
						}
					}
					else if (getPixel(x,y+1) == background){
					edge_length += pixelWidth;
	                    if (getPixel(x,y-1) == background){
						edge_length += pixelWidth;
						}
					}
					else if (getPixel(x,y-1) == background){
					edge_length += pixelWidth;
					}
				}
				else if (getPixel(x-1,y) == background){ //Prime
				edge_pixel += 1;
				edge_length += pixelWidth;
					if (getPixel(x,y+1) == background){
					edge_length += pixelWidth;
						if (getPixel(x,y-1) == background){
						edge_length += pixelWidth; 
						}
					}
					else if (getPixel(x,y-1) == background){
					edge_length += pixelWidth;
					}
				}
				else if (getPixel(x,y+1) == background){ //Prime
				edge_pixel += 1;
				edge_length += pixelWidth;
					if (getPixel(x,y-1) == background){
					edge_length += pixelWidth;
				    }
				}
				else if (getPixel(x,y-1) == background){ //Prime
				edge_pixel += 1;
				edge_length += pixelWidth;
				}
			}
		}
	}
	//Retrieve slice depth
	slice_depth = (slice_number-1)*depth; 
	//Retrieve area pixel count, area in microns^2 and edge pixel count
	area_micron = area_pixels*pixel_area;
	print(slice_number,slice_depth,area_pixels,area_micron,edge_pixel,edge_length);
	//Retrieve compactness normalisation factor
    cnf = pcn4+(pcn1-pcn4)/(1+Math.pow(area_pixels/pcn3,pcn2)); 
	//Retrieve compactness
	circularity = ((4*3.14159265359*area_pixels)/(edge_pixel*edge_pixel))/cnf;
	//Retrieve arrays
	slice_number_array = Array.concat(slice_number_array,slice_number);
	area_pixels_array = Array.concat(area_pixels_array,area_pixels);
	edge_pixel_array = Array.concat(edge_pixel_array,edge_pixel);
	slice_depth_array = Array.concat(slice_depth_array,slice_depth);
	area_micron_array = Array.concat(area_micron_array,area_micron);
	edge_length_array = Array.concat(edge_length_array,edge_length);
	circularity_array = Array.concat(circularity_array,circularity);
	circularity_1_array = Array.concat(circularity_1_array,1);
	fourier_frequency_array = Array.concat(fourier_frequency_array,fourier_fequency);
	//Go to next slice and update slice number 
	run("Next Slice [>]");
	slice_number += 1;
	fourier_fequency += 1;
}
waitForUser("Progress","Done");

//Get the amplitudes of the frequencies of the dendrite area
area_FT = Array.fourier(area_pixels_array,"Hann");
Array.print(area_FT);

//Print plots
dcomp = getBoolean("Do you want to display the compactness?");
if (dcomp == 1){ 
Plot.create("Dendrite compactness", "distance (microns)", "compactness (Polsby-Popper)");
Plot.add("line", slice_depth_array, circularity_array);
Plot.setLimits((slice-1)*depth,(nSlices-1)*depth, 0, 1);
}

dae = getBoolean("Do you want to display the area in microns^2 and periphery in microns?");
if (dae == 1){ 
Plot.create("Dendrite cross-sectional area (microns)", "slice", "microns");
Plot.setColor("red");

Plot.add("line", slice_depth_array, edge_length_array);
Plot.setLimits((slice-1)*depth,(nSlices-1)*depth, 0, 25);
}

dae = getBoolean("Do you want to display the area and periphery in pixels?");
if (dae == 1){ 
Plot.create("Dendrite cross-sectional area (pixels)", "slice", "pixels");
Plot.add("line", slice_number_array, area_pixels_array);
Plot.setColor("red"); 
Plot.add("line", slice_number_array, edge_pixel_array);
Plot.setLimits(slice, nSlices, 0, 5000);
}

fo = getBoolean("Do you want to display the fourier spectrum of the area (pixels)?");
if (fo == 1){ 
Plot.create("Frequency spectrum of the dendrite area (pixels) as a function of slices", "frequency", "amplitude", fourier_frequency_array, area_FT);
}

save_option = getBoolean("Want to save results?");
if (save_option == 1){
//Make a table containing the arrays
Table.create("Area_Perifery_Compactness_Fourier");
Table.setColumn("depth_slices",slice_number_array);
Table.setColumn("area_pixels",area_pixels_array);
Table.setColumn("edge_pixels",edge_pixel_array);
Table.setColumn("depth_microns",slice_depth_array);
Table.setColumn("area_microns",area_micron_array);
Table.setColumn("edge_microns",edge_length_array);
Table.setColumn("compactness_Poslby-Popper",circularity_array);
Table.setColumn("Frequency",fourier_frequency_array);
Table.setColumn("Fourier",area_FT);
Table.save(dir+"Area_Perifery_Compactness_Fourier"+".csv");
}

//This is to test whether the edge pixel counter works on a part only containig 8 white pixels
edge_pixel = 0;
for (x=169;x<w;x++){
	for (y=105;y<h;y++){
		if (getPixel(x,y) == 255){
			if (getPixel(x+1,y) == background){
			edge_pixel += 1;
			}
			else if (getPixel(x-1,y) == background){
			edge_pixel += 1;
			}
			else if (getPixel(x,y+1) == background){
			edge_pixel += 1;
			}
			else if (getPixel(x,y-1) == background){
			edge_pixel += 1;
			}
		}	
	} 
}
print(edge_pixel);
//This is to determine for a region containing 8 pixels with 10 borders with black pixels what the border length is.
edge_pixel = 0;
edge_length = 0;
for (x=169;x<w;x++){
	for (y=105;y<h;y++){
		if (getPixel(x,y) == 255){
		area_pixels += 1;
			if (getPixel(x+1,y) == background){ //Prime
			edge_pixel += 1;
			edge_length += pixelHeight;
				if (getPixel(x-1,y) == background){
				edge_length += pixelHeight;
					if (getPixel(x,y+1) == background){
					edge_length += pixelWidth;
						if (getPixel(x,y-1) == background){
						edge_length += pixelWidth;
						}
					}
				}
				else if (getPixel(x,y+1) == background){
				edge_length += pixelWidth;
	            	if (getPixel(x,y-1) == background){
					edge_length += pixelWidth;
					}
				}
				else if (getPixel(x,y-1) == background){
				edge_length += pixelWidth;
				}
			}
			else if (getPixel(x-1,y) == background){ //Prime
			edge_pixel += 1;
			edge_length += pixelWidth;
				if (getPixel(x,y+1) == background){
				edge_length += pixelWidth;
					if (getPixel(x,y-1) == background){
					edge_length += pixelWidth; 
					}
				}
				else if (getPixel(x,y-1) == background){
				edge_length += pixelWidth;
				}
			}
			else if (getPixel(x,y+1) == background){ //Prime
			edge_pixel += 1;
			edge_length += pixelWidth;
				if (getPixel(x,y-1) == background){
				edge_length += pixelWidth;
				}
			}
			else if (getPixel(x,y-1) == background){ //Prime
			edge_pixel += 1;
			edge_length += pixelWidth;
			}
		}
	}
}
print(edge_length);
print(depth);