#! /usr/bin/Rscript

# R function which gives an ordered list of indices on the boundary of points in the input
suppressMessages(library('alphahull'))
#~ library('alphahull')

args <- commandArgs(TRUE)
filename <- args[1] #file containing input point-list
alpha <- as.double(args[2]) #value of alpha

# read the file and prepare the input vector
x=read.table(filename)
# CAUTION : the first three points cannot be collinear!!!
plot(x)
x=unique(x)

ahull2py=function(x,alpha){
	ah = ahull(x, alpha = alpha)
	contour_arcs = ah$arcs
	contour_point_indices = contour_arcs[,7]
	return(contour_point_indices)
}

y=ahull2py(x,alpha)

xx = array(0,c(length(y),2))
for(i in 1:length(y)){
	xx[i,1]=x[y[i],1]
	xx[i,2]=x[y[i],2]
}

# function to write contours into Free-D format curve-files
write.curve2freed = function(leaf,filename){
	d=dim(leaf)
	curv_length = d[1]
	dimension = d[2]
	write('# 0',file=filename, append=FALSE)
	write(paste('# ',dimension,sep=""),file=filename, append=TRUE)
	write(paste('# ',curv_length,sep=""),file=filename, append=TRUE)
	for (j in 1:curv_length) {
		write(paste(leaf[j,1],leaf[j,2],sep=" "),file=filename, append=TRUE)
	}
}

write.curve2freed(xx,filename)

