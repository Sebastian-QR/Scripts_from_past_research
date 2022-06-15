//Ask the user to select the source directory - tiff stacks from batch conversion
dir1 = getDirectory("Select NEURONS folder");

//Ask the user to select the destination directory - where to save images  and results
dir2 = getDirectory("Select TRACES folder");

//Ask the user to select the destination directory - where to save images  and results
dir3 = getDirectory("Select OUTPUT folder");


list1 = getFileList(dir1); //neurons list
Array.sort(list1);
//Array.print(list1);

list2 = getFileList(dir2); //traces list
Array.sort(list2);
//Array.print(list2);

for (i=0; i<list1.length; i++) {
	filename = dir1 + list1[i];
	tracename = dir2 + list2[i];
	outputfolder = dir3;
	open(filename);
	nameStore = getTitle();
	newname = substring(nameStore, 0, lengthOf(nameStore)-10);
run("Sholl Analysis (Tracings)...", "traces/(e)swc=["+tracename+"] image=[] load center=[Start of main path] radius=15 enclosing=1 #_primary=[1] infer linear polynomial=[Best fitting degree] most semi-log normalizer=Area/Volume save directory=["+outputfolder+"] do");
selectWindow("Sholl Results");
run("Close");
run("Close All");
}

//dir2+"+newname+".traces

//close(newname+" C=0");
//selectWindow(newname+" C=1");
//run("Invert", "stack");
//run("Simple Neurite Tracer", "look_for_previously_traced");
//waitForUser("Run Sholl and Click 'OK'");
//selectWindow("Tracings for "+newname+" C=1_Sholl-Profiles");
//saveAs("Results", dir2+newname+" - C=1_Sholl-Profiles.csv");
//selectWindow("Sholl Results");
//saveAs("Results", dir2+newname+" Meta Results.csv");
//run("Close All");
//}