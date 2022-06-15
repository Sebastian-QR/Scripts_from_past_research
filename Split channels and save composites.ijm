//Ask the user to select the source directory
dir = getDirectory("Select input directory");

//Ask the user to select the destination directory - where to save images  and results
dir2 = getDirectory("Select output");

//generate a list of the files in the selected directory
//open split to single channel images and save
list = getFileList(dir);
Array.sort(list);
for (i=0; i<list.length; i++) {
	filename = dir + list[i];
	open(filename);
	nameStore = getTitle();
newname = substring(nameStore, 0, lengthOf(nameStore)-5);
//Create and save composite of analised image
run("Merge Channels...", "c1=["+newname+"- C=2] c2=["+newname+"- C=1] c3=["+newname+"- C=0] create keep");	
selectWindow("Composite");
saveAs("Jpeg", dir2+newname+"(composite).jpg");
selectWindow(newname+"- C=0");
saveAs("Jpeg", dir2+newname+"(C0).jpg");
selectWindow(newname+"- C=1");
saveAs("Jpeg", dir2+newname+"(C1).jpg");
selectWindow(newname+"- C=2");
saveAs("Jpeg", dir2+newname+"(C2).jpg");
//Close composite and channel 0
close("Composite");
close(newname+"- C=0");
close(newname+"- C=1");
close(newname+"- C=2");
}