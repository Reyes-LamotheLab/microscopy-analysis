// @File(label = "Input directory", style = "directory") input
a=1
list = getFileList(input);
list = Array.sort(list);


for (i = 0; i < list.length; i++) {
	if (i%2==0 && endsWith(list[i], ".tif")) {
		print(i);
		filename1 = input +"/"+ list[i];
	    //filenameMod1=filename1.replace("\\", "/");
	    //print(filename1);

		open(filename1);
		//open(filenameMod2);
	
		run("Detect Particles", "ch1i ch1l ch1a=3 ch1s=15 rois=Ovals add=Nothing summary=Append");

		selectWindow("Results");
		
		//fileresults =input +"/"+list[i]+"_Results";
	    //filenameMod=fileresults.replace("\\", "/");
		saveAs("Results", input +"/"+list[i]+"_Results");
	    close();
	}

	
	else {
		print(i);
		filename1 = input +"/"+ list[i];
	    //filenameMod1=filename1.replace("\\", "/");
	    print(filename1);

		open(filename1);

	
		run("Detect Particles", "ch1i ch1l ch1a=3 ch1s=15 rois=Ovals add=Nothing summary=Append");

		selectWindow("Results");
		
		//fileresults = input +"/"+list[i]+"_Results2";
	    //filenameMod=fileresults.replace("\\", "/");
		saveAs("Results", input +"/"+list[i]+"_Results");
	    close();

	}
	
	}
