

// @File(label = "Input directory", style = "directory") input
a=0
list = getFileList(input);
list = Array.sort(list);


for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], ".tif")){
    //setTool("line");
    filename = input +"/"+ list[i];
    filenameMod=filename.replace("\\", "/");
    print(filenameMod);

    name=split(filenameMod,".");
    filenameNew= name[0]+"new.tif";  //new file which will be saved in the same folder of the original image
    
    open(filenameMod);
    name=split(list[i],".");

//DIFFERENT DENOISING STEPS , uncomment the technique you want to use before running
    
//first
	//run("Mean...", "radius=1.5");
	//run("Subtract Background...", "rolling=10 sliding");

//second
	run("Despeckle");
	run("Top Hat...", "radius=7 ");	

//third 
	//run("Subtract Background...", "rolling=15 sliding");
	//run("Top Hat...", "radius=10");
	//run("Remove Outliers...", "radius=2 threshold=50 which=Bright");

	saveAs("Tiff", filenameNew);}
	close();}
	