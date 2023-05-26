// @File(label = "Input directory", style = "directory") input

/* INSTRUCTIONS:
 
 - If you have two channels where the odd numbers (starting at "1") are the first channel change the parameters
 in the first "if" statement and the even numbers in the "else" statement
 - If you only have one channel just comment out the "if" line and all that is IN the "else".
 - IF YOU HAVE WINDOWS UNCOMMENT THE lines: " //filename1=filename1.replace("\\", "/"); " 
	   
*/

a=1
list = getFileList(input);
list = Array.sort(list);
newFolder = input + "/particles_result";
File.makeDirectory(newFolder);

for (i = 0; i < list.length; i++) {
	

	//if (i%2==0 && endsWith(list[i], ".tif")) {
		print(i);
		filename1 = input +"/"+ list[i];
		
	
	
	    //filename1=filename1.replace("\\", "/");
	   
		open(filename1);
	
		run("Detect Particles", "ch1i ch1l ch1a=2 ch1s=15 rois=Ovals add=Nothing summary=Append");
	
		selectWindow("Results");
		
	    //input=input.replace("\\", "/");  
	    
		saveAs("Results", newFolder +"/"+list[i]+"_Results");
	    close();
	//}

/*
	else {
		print(i);
		filename1 = input +"/"+ list[i];
		
	    //filename1=filename1.replace("\\", "/");
	   
	    print(filename1);

		open(filename1);

		run("Detect Particles", "ch1i ch1l ch1a=3 ch1s=15 rois=Ovals add=Nothing summary=Append");

		selectWindow("Results");
		
	    //input=inpur.replace("\\", "/");
	    
		saveAs("Results", input +"/"+list[i]+"_Results");
	    close(); 

	} */
	
	}
