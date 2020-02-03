
ID = getImageID();
selectImage(ID);

//split channels
selectImage(ID);
if(channels==1) {
	print("one channel only");
} 
if(channels>=2) {
	run("Split Channels");
}




filename=getInfo("image.filename");
dirname=getInfo("image.directory");

name = getTitle;

print("ID="+ID)
print("filename="+filename)
print("dirname="+dirname)

run("Split Channels");
selectWindow(filename+" (red)");
close();
selectWindow(filename+" (blue)");
close();
selectWindow(filename+" (green)");
run("Gaussian Blur...", "sigma=20");
setAutoThreshold("Otsu");
setOption("BlackBackground", false);
run("Convert to Mask");
//run("Gray Morphology", "radius=20 type=circle operator=close");
//run("Gray Morphology", "radius=20 type=circle operator=open");
