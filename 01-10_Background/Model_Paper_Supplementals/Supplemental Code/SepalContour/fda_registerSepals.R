#! /usr/bin/Rscript

## 
# Prerequisite :
# install.packages("fda")
#
#
# Usage :
# fda_registerSepals.R  Npoints
# 
# structure of leafsX.data :
# label, length, nb. of primery tooth, tipindex, x_coordinates
# 
#
#

suppressMessages(library('fda'))

args <- commandArgs(TRUE)
Npoints <- as.integer(args[1])
Rejlist <- toString(args[2])
print(Npoints)
print(Rejlist)
#regtype <- args[2] 
#X_file_name_in <- args[1]
#Y_file_name_in <- args[2]
#x <- as.double(args[2])

# conversion of Rejlist from string to array of integers
if (!Rejlist=="NA"){
	Rejlist <- strsplit(Rejlist,"-")
	Rejlist <- Rejlist[[1]][2:length(Rejlist[[1]])]
	Rejlist <- as.integer(Rejlist)
} else { Rejlist <- c() }
# the number of points defining the curves
len=2*Npoints

# the number of basis-functions in fd objects
nbasis=len

# number of landmarks per half leaf :
lm_classnumber = c(0)
lm_paramvalues = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 8, 9, 9.5, 9.8, 9.9, 10)
#lm_paramvalues = c(0, 0.2, 0.4, 0.6, 1, 2, 5, 8, 9, 9.4, 9.6, 9.8, 10) # worst parameters
Nclasses = length(lm_classnumber)
maxNlandmarks = length(lm_paramvalues)
tipvalue = (maxNlandmarks-1)/2 +1


read.leafdata = function(fileX,fileY) {
	tableX=read.table(fileX)
	tableY=read.table(fileY)
	dd=dim(tableX)
	leafdata=array(0,c(dd[1],dd[2],2))	
	for(i in 1:dd[1]){
		for(j in 1:dd[2]){
			leafdata[i,j,1]=tableX[i,j]
			leafdata[i,j,2]=tableY[i,j]
		}
	}
	return(leafdata)
}


read.leaftable = function(filename) {
	whole_table=read.table(filename)
	dd=dim(whole_table)
	#print(dd)
	d2 = dd[2]/2
	leafdata=array(0,c(dd[1],d2,2))
	print(dim(leafdata))
	for(i in 1:dd[1]){
		for(j in 1:d2){
			leafdata[i,j,1]=whole_table[i,j]
			leafdata[i,j,2]=whole_table[i,d2+j]
		}
	}
	return(leafdata)
}

order.leaftable = function(mt) {
	# this function order the matrix mt with respect to the first line
	# (with respect to [1,]) - here the leaflength
	mtt = t(mt)
	mtto = mtt[order(mtt[,1]),] 
	mto = t(mtto)
	return(mto)
}

plot.a_leaf = function(leafsLM,i){
	plot(leafsLM[,i,1],leafsLM[,i,2],type='l')
	return(1)
}


# eliminate an individual because of a problem :
eliminate.leaf = function(leafs,i){
	d= dim(leafs)
	new_leafs = array(0,c(d[1],(d[2]-1),d[3]))
	for (j in 1:(i-1)) {
		new_leafs[,j,] = leafs[,j,]
	}
	for (j in i:(d[2]-1)) {
		new_leafs[,j,] = leafs[,j+1,]
	}
	return(new_leafs)
}

# eliminate an individual because of a problem :
eliminate.leaflength = function(leafl,i){
	len= length(leafl)
	new_leafl = array(0,(len-1))
	for (j in 1:(i-1)) {
		new_leafl[j] = leafl[j]
	}
	for (j in i:(len-1)) {
		new_leafl[j] = leafl[j+1]
	}
	return(new_leafl)
}

# eliminate a list of leafs :
eliminate.leafs = function(leafs,Rejlist){
	#d= dim(leafs)
	#new_leafs = array(0,c(d[1],(d[2]-length(Rejlist)),d[3]))
	new_leafs = leafs[,-Rejlist,]
	return(new_leafs)
}

# eliminate a list of leafs :
eliminate.leaflengths = function(leafl,Rejlist){
	#len= length(leafl)
	#new_leafl = array(0,(len-length(Rejlist)))
	new_leafl = leafl[-Rejlist]
	return(new_leafl)
}


superpose.leafs =function(leafsLM, regtype){
	s= dim(leafsLM)
	N_all=s[2]
	Npoints=s[1]
	if (regtype=='cm') {
		for (i in 1:N_all){
			x=leafsLM[,i,1]
			cx=mean(x)
			leafsLM[,i,1]=x-cx
		}
	}
	if (regtype=='petiole') {
		print('registration on petiole')
		for (i in 1:N_all){
			x=leafsLM[,i,1]
			cx=(x[1]+x[Npoints])/2
			leafsLM[,i,1]=x-cx
		}
	}
	# otherwise we leave it superimposed on the tip.
	return(leafsLM)
}


directory = paste('registration_on_common_landmarks',Npoints,sep="")
directory_in = paste('ContoursData_sample',Npoints,'/',sep="")

###############################################################################
# read the data so that we its dimension and construct the final datatable
###############################################################################
	

X_file_name_in = paste(directory_in,'X_contours0.data',sep="")
Y_file_name_in = paste(directory_in,'Y_contours0.data',sep="")
if (file.exists(X_file_name_in)) {
	print('file exists')
	leafsX=read.table(X_file_name_in)
	leafsY=read.table(Y_file_name_in)
	d=dim(leafsX)
	N_leafs = d[1] # number of curves read
} else {
	N_leafs = 0
}



print(paste('Number of leafs in this genotype: ',N_leafs,sep=""))

# ================================
# initialise the table for the registered data
# ================================

# insert the leaflength in the first column 
# in order to be able to reorder the table with respect to it

leaflength_all = array(1,N_leafs)
leafs_X_all=array(0,c(len+1,N_leafs))
leafs_Y_all=array(0,c(len+1,N_leafs))
leafnames_all = array('',N_leafs)

###############################################################################
# read and register the data in each class of leafs
# ( classes are defined with respect to the number of landmarks per half-leaf)
###############################################################################
ind = 1
i0 = 1
print(paste('Data with 0 landmarks per half leaf',sep=""))
print('--------------------------')
X_file_name_in = paste(directory_in,'X_contours0.data',sep="")
Y_file_name_in = paste(directory_in,'Y_contours0.data',sep="")
if (file.exists(X_file_name_in)) {
	# ================================
	# read the data :
	# ================================
	leafsX=read.table(X_file_name_in)
	leafsY=read.table(Y_file_name_in)
	# transforme them into data.frame
	lX=data.frame(leafsX)
	lY=data.frame(leafsY)
	# ================================
	# extract information :
	# ================================
	d=dim(leafsX)
	N = d[1] # number of curves read
	#nb = 3 + Nlandmarks # number of info fields, other than coordinates 
	nb <- 4
	#len=d[2]-nb # length of point_list defining the curve
	print(paste('Number of leafs: ', N,sep=""))
	leafnames=lX[,1]
	leaflength=lX[,2]
	N_primary = lX[,3]
	landmark_index=lX[,4]
	maxleaflength=max(leaflength)
	minleaflength=min(leaflength)
	print(paste('Minimal and maximal leaflength: ', minleaflength,' - ',maxleaflength,sep=""))
	# ================================
	# construct the data table :
	# ================================
	leafsample=array(0,c(len,d[1],2))
	print(dim(leafsample))
	print(dim(lX))
	print(len)
	for(i in 1:d[1]){
		for(j in 1:len){
			leafsample[j,i,1]=lX[i,j+nb]
			leafsample[j,i,2]=lY[i,j+nb]
			leafsample[j,i,1]=lX[i,j+nb]
			leafsample[j,i,2]=lY[i,j+nb]
		}
	}
	#print(regtype)
	#leafsample=superpose.leafs(leafsample,regtype)
	coordinates=list("X","Y")
	replicates=vector("list",d[1])
	for(i in 1:d[1]){
		replicates[[i]]=lX[i,1]
	}
	# deletion of the not ok leafs
	if (length(Rejlist)!=0){
		leafnames <- leafnames[-Rejlist]
		leafsample <- eliminate.leafs(leafsample, Rejlist)
		leaflength <- eliminate.leaflengths(leaflength, Rejlist)
		N_leafs <- N_leafs-length(Rejlist)
	}
	if (N<1) {#(N>1) {
		# ================================
		# construct fd objects :
		# ================================
		t0=lm_paramvalues[1]
		t1=lm_paramvalues[maxNlandmarks]
		fdatime = seq(t0,t1,len=len)
		fdabasis = create.bspline.basis(c(t0,t1),norder=4,nbasis)
		# information for smoothing :
		#fdaPar = fdPar(fdabasis,int2Lfd(2),1e-8)
		fdaPar = fdPar(fdabasis,0,1e-8)
		wfdPar = fdPar(fdabasis,1,1e-8)
		smoothing=smooth.basis(fdatime,leafsample,fdaPar)
		# take the fd object returned by smoothing
		leafsfd = smoothing$fd
		
		# the parameter-values for the landmarks
		landmark_time=landmark_index*(t1-t0)/len
		print("landmark-time")
		print(landmark_time)
		reflandmarks = mean(landmark_time)
		print("mean landmaks")
		print(reflandmarks)
		reflandmarks=lm_paramvalues[tipvalue]
		print("reflandmaks")
		print(reflandmarks)
		# register on landmarks	
		land_regist=landmarkreg(leafsfd,landmark_time,mean(landmark_time),wfdPar,monwrd=FALSE)
		reg_leafsfd=land_regist$regfd
			# evaluate and write the registered data into a file
		leafsLM=eval.fd(fdatime,reg_leafsfd)
	} else {
		# ================================
		# construct fd objects :
		# ================================
		t0=lm_paramvalues[1]
		t1=lm_paramvalues[maxNlandmarks]
		fdatime = seq(t0,t1,len=len)
		fdabasis = create.bspline.basis(c(t0,t1),norder=4,nbasis)
		# information for smoothing :
		#fdaPar = fdPar(fdabasis,int2Lfd(2),1e-8)
		fdaPar = fdPar(fdabasis,0,1e-8)
		wfdPar = fdPar(fdabasis,1,1e-8)
		smoothing=smooth.basis(fdatime,leafsample,fdaPar)
		# take the fd object returned by smoothing
		leafsfd = smoothing$fd
		leafsLM=leafsample
		reg_leafsfd = leafsfd
	}
	write.table(leafsLM,file=paste(directory,"/landmark_registered_leafs_0.data",sep=""))
	
	leaflength_all = array(1,N_leafs)
	leafs_X_all=array(0,c(len+1,N_leafs))
	leafs_Y_all=array(0,c(len+1,N_leafs))
	leafnames_all = array('',N_leafs)
	leafnames_all[1:N_leafs] = as.character(leafnames)
	leaflength_all[1:N_leafs] = leaflength
	leafs_X_all[2:(len+1),1:N_leafs] = leafsLM[,,1]
	leafs_Y_all[2:(len+1),1:N_leafs] = leafsLM[,,2]
	
	# represent the registered curves
	png(paste(directory,'/TheDataXY_lm_reg_0.png',sep=""))
	par(mfrow=c(2,1))
	plot(reg_leafsfd)
	dev.off()
} else {
	print(paste('Number of leafs: ', 0,sep=""))
}


###############################################################################
# si les fichiers existent déjà, on les lit :
###############################################################################



print('Sorti de la boucle, on trie')

leafs_X_all[1,] = leaflength_all
leafs_Y_all[1,] = leaflength_all

# order X coordinates
#--------------------
#leafs_X_ordered=order.leaftable(leafs_X_all)

# order Y coordinates
#--------------------
#leafs_Y_ordered=order.leaftable(leafs_Y_all)

# order leafnames
#--------------------
lntable <- data.frame(leafl=leaflength_all, leafn=leafnames_all)
#lntable[with(lntable, order(leafl)),]
#print(leafnames_all)
#print(leaflength_all)

# construct the whole data structure with ordered leaves
# --------------------------------------------------

#leafsLM=array(0,c(len,N_leafs,2))
#leafsLM[,,1] = leafs_X_ordered[2:(len+1),]
#leafsLM[,,2] = leafs_Y_ordered[2:(len+1),]



#print('after ordering')
#leaflength = leafs_X_ordered[1,]
leaflength = leafs_X_all[1,]

# function to write leaf contours into Free-D format curve-files
write.leafcurves2freed = function(leafs,namelist,directory_name,suffix){
	#dir.create(directory_name)
	d=dim(leafs)
	curv_nb = d[2]
	curv_length = d[1]
	dimension = d[3]
	for (i in 1:curv_nb){
		filename = paste(directory_name,'/',namelist[i],suffix,sep="")
		write('# 0',file=filename, append=FALSE)
		write(paste('# ',dimension,sep=""),file=filename, append=TRUE)
		write(paste('# ',curv_length,sep=""),file=filename, append=TRUE)
		for (j in 1:curv_length) {
			write(paste(leafs[j,i,1],leafs[j,i,2],sep=" "),file=filename, append=TRUE)
		}
	}
}

print('On ecrit les résultats')
curvenamelist = lntable[,2]
print(curvenamelist)
suffix=paste('_Registered',str(2*Npoints),'.cv',sep="")
write.leafcurves2freed(leafsLM,curvenamelist,'registered_curves',suffix)


# store :
#---------------
write.table(leafsLM,file=paste(directory,"/landmark_registered_leafs_all.data",sep=""))
write.table(leaflength,file=paste(directory,"/landmark_registered_leafs_all-leaflength.data",sep=""))

print('On normalise')
leafsLM_n = leafsLM
for (i in 1:N_leafs){
	leafsLM_n[,i,] = leafsLM[,i,]/leaflength[i]
}

print('On ecrit les résultats')
# store :
#---------------
write.table(leafsLM_n,file=paste(directory,"/landmark_registered_leafs_all_normalised.data",sep=""))



