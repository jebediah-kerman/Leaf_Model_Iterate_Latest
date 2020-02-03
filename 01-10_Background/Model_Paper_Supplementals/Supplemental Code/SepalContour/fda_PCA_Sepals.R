#! /usr/bin/Rscript

## 
# Prerequisite :
# install.packages("fda")
#
#
# Usage :
# fda_PCA_Sepals.R  Npoints #gen1 gen2 gen3 gen4
# 
# structure of leafsX.data :
# label, length, nb. of primery tooth, tipindex, x_coordinates
# 
#
#
                                                           

args <- commandArgs(TRUE)

Npoints=400
Npoints <- as.integer(args[1])
len=2*Npoints

# here is the list of samples to compare :
#genes = c("col-0")#c("mbd")#c("col-0")#c("ws4")#c("N6549xmbd")#c("bot1.7xmbd","mbd","col-0","ws4","N6549xmbd")
#c("test")#c("aba-pom2.4")#c("aba-pom2.4","12") #c("1", "3", "5","7","9","11","JL")
#types <- c('aba')

genes <- toString(args[2])
scale_type <- toString(args[3])

# genes & types formatage treatment	
genes <- strsplit(genes,"/")
genes <- genes[[1]][1:length(genes[[1]])]


N_gentypes = length(genes) #nb de genotypes
N_types <- 1
N_subclasses <- length(genes)

# number of classes with respect to the leaflength # we make no difference between lengths
N_classes=1

Npoints=400
Npoints <- as.integer(args[1])
len=2*Npoints


all_colors = c("black", "red", "blue", "green", "cyan", "magenta", "orange","purple")
#all_colors = c( "red", "blue", "green", "cyan", "magenta", "orange")
gen_colors = all_colors[1:N_gentypes]


suppressMessages(library('fda'))

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


read.lltable = function(filename) {
	whole_table=read.table(filename)
	dd=dim(whole_table)
	leafdata=array(0,c(dd[1]))
	#print(dim(leafdata))
	for(i in 1:dd[1]){
			leafdata[i]=whole_table[i,1]
	}
	return(leafdata)
}


directory_in = paste('/registration_on_common_landmarks',Npoints,sep="")

directory=paste("../ImagesPCA_",scale_type,"/ImagesPCA",sep='')

for (i in genes) {
	print(i)
	directory=paste(directory,i,sep="_")
	}
dir.create(directory)

determine_directory = function(gen) {
	return(paste(gen,directory_in,sep=""))
}


determine_datafile = function(gen) {
	#return(paste(gen,directory_in,"/registered_landmarks_normalised_leafs_and_mean_basis200.data",sep=""))
	#return(paste(gen,directory_in,"/landmark_registered_normalised_leafs.data",sep=""))
	return(paste(gen,directory_in,"/landmark_registered_leafs_all_normalised.data",sep=""))
}

# non ecrit dans l'etape d'avant
determine_mdatafile = function(gen) {
	#return(paste(gen,directory_in,"/registered_landmarks_normalised_leafs_and_mean_basis200.data",sep=""))
	#return(paste(gen,directory_in,"/registered_landmarks_means_normed.data",sep=""))
	return(paste(gen,directory_in,"/registered_landmarks_means.data",sep=""))
}

determine_leaflengthfile = function(gen) {
	return(paste(gen,directory_in,"/landmark_registered_leafs_all-leaflength.data",sep=""))
}

determine_mleaflengthfile = function(gen) {
	return(paste(gen,directory_in,"/registered_landmarks_class_leaflength.data",sep=""))
}

Nsubclasses <-array(0,c(N_subclasses))
Ngenes = array(0,c(N_gentypes))


for(i in 1:N_gentypes){
	temp <- 0
	file_ll = determine_leaflengthfile(genes[i])
	llgen=read.lltable(file_ll)
	Nsubclasses[(i-1) + t]=length(llgen)
	temp <- temp + length(llgen)
	}
	Ngenes[i] <- temp
}


# pop the zeros of a vector :
pop_zero=function(a){
	b=a
	for (i in length(a):1){
		if (a[i]==0) {b=b[-i]}
	}
	return(b)
}



N_all=sum(Ngenes)	#nb total d'images

leafsample = array(0,c(len,N_all,2))
leaflength = array(0,c(N_all))

# donnees a partir d'un fichier non ecrit
mleafs = array(0,c(N_gentypes,len,N_classes,2))		
#mll = array(0,c(N_gentypes,N_classes))


indstart=1
for(i in 1:N_gentypes){
	print(i)
	dirgen = determine_directory(genes[i])
	filegen = determine_datafile(genes[i])
	file_ll = determine_leaflengthfile(genes[i])
	#mfile = determine_mdatafile(genes[i])
	#mfile_ll = determine_mleaflengthfile(genes[i])
	indstop=indstart+Nsubclasses[(i-1)+t]-1
	leafsample[,indstart:indstop,]=read.leaftable(filegen)			
	leaflength[indstart:indstop]=read.lltable(file_ll)
	#mleafs[i,,,]=read.leaftable(mfile)		# donnees a partir d'un fichier non ecrit
	#mll[i,]=read.lltable(mfile_ll)
	indstart=indstop+1

	}

# ================================
# represent mean evolutions comparatively :
# ================================

# 1. compare mutants for each class
# ----------------------------------------------------

# ! mleafs n existe pas ici : liste de landmarks
#for (i in 1:N_classes){	
#	plot(0,0, main=paste("Mean contours, class no. ",i, ", meanlength=",round(mll[1,i]),"um", sep=""),
#		xlab="X", ylab="Y",
#		xlim=c( 0,1), ylim=c(-0.5,0.5))
#	for (j in 1:N_gentypes){
#		lines(mleafs[j,,i,1],mleafs[j,,i,2],lwd=2, col=gen_colors[j])#,col="black")
#	}
#	legend("topright", inset=.05, genes, fill=gen_colors)
#	grid()
#	dev.copy(png,paste(directory,"/means_contours_",i,".png",sep=""))
#	dev.off()
#}


# ================================
# assemble the data construct the "leafs" data table :
# ================================

# constructions of fd objects :
# ------------------------------
t0=0
t1=10

fdatime = seq(t0,t1,len=len)
fdabasis = create.bspline.basis(c(t0,t1),norder=4,nbasis=200)

# information for smoothing :
# ---------------------------
fdaPar = fdPar(fdabasis,int2Lfd(2),1e-8)

shift.leafsample = function(leafsLM, a){
	# shifts the leaf leafsLM[,i,] by the amount a[i-1] on the x axis
	# the first leaf stays in place
	d = dim(leafsLM)
	leafsLM1=leafsLM
	N_all = d[2]
	for (i in 2:N_all){
		leafsLM1[,i,1]=leafsLM[,i,1]+a[i-1]
		}
	return(leafsLM1)
}

total.variability = function(a) {
	leafsample1=shift.leafsample(leafsample,a)
	# construct the fd object correspondig to the sample
	leafsfd=smooth.basis(fdatime,leafsample1,fdaPar)$fd
	# Bivariate PCA
	#N_harmonics=5
	fdapca = pca.fd(leafsfd)#,nharm=N_harmonics)
	lambda=fdapca$varprop
	cum_lambda=cumsum(lambda)
	mu = fdapca$values
	smu=sum(mu)
	return(smu)
}


# ---------------------


# constructing the fd object
leafsfd=smooth.basis(fdatime,leafsample,fdaPar)$fd

# ================================
# Bivariate PCA
# ================================

N_harmonics=6

fdapca = pca.fd(leafsfd,nharm=N_harmonics)
lambda=fdapca$varprop
cum_lambda=cumsum(lambda)
mu = fdapca$values
smu=sum(mu)
 
png(paste(directory,"/lambda.png",sep=""))
plot(lambda,type='o', pch=20,xlab="",ylab="")
title(main=paste('Lambda (Sum of eigenvalues =',smu,')',sep=""))
grid()
#dev.copy(png,paste(directory,"/lambda.png",sep=""))
dev.off()

png(paste(directory,"/lambda_cumulated.png",sep=""))
plot(cum_lambda,type='o', pch=20,xlab="",ylab="")
title(main='Cumulated lambda')
grid()
dev.off()


# the mean contour
# -------------------
meanfd=fdapca$meanfd
# plot the mean contour

# evaluated on the same timeintervall as the data :
meanvec=eval.fd(fdatime, meanfd)
# plot the mean vector
png(paste(directory,"/mean_sepal.png",sep=""))
plot(meanvec[,,1],meanvec[,,2],pch=20,asp=1,xlab="",ylab="")
grid(col='black')
mean_vec=meanvec[,1,]
dev.off()

print.contour=function(mean_vec,filename){
	#plot(meanvec[,,1],meanvec[,,2],pch=20,xaxt='n', yaxt='n',ann=FALSE)
	png(filename)
	plot(mean_vec[,1],mean_vec[,2],pch=20,axes=FALSE,ann=FALSE,asp=1)
	#plot(mean_vec[,1],mean_vec[,2],asp=1,pch=20,xaxt='n', yaxt='n',ann=FALSE)
	#grid(col='black')
	dev.off()
}

print.contour(mean_vec,paste(directory,"/mean_sepal_noaxis.png",sep=""))

##############################################
# a matrix of scores on the principal components
##############################################

scores=fdapca$scores

# the total score is the sum of x and y scores
scores_xy=scores[,,1]+scores[,,2]

##############################################
# interpret the different principal components
##############################################
harmonicsvec=eval.fd(fdatime, fdapca$harmonics)

# check the length of principal vectors
PC_vec_module = array(0,N_harmonics)
for (i in 1:N_harmonics){
	PC_fd = fdapca$harmonics[i]
	PC_vec_module[i] = sqrt(inprod(PC_fd[1],PC_fd[1])+inprod(PC_fd[2],PC_fd[2]))
}


# a function which projects a curve on a principal axis ?
project.on.PC = function(vec,pcindex){
	vec_centered=vec-meanvec[,1,]
	leaf_fd=(smooth.basis(fdatime,vec_centered,fdaPar))$fd
	PC_fd=fdapca$harmonics[pcindex]
	score = diag(inprod(leaf_fd,PC_fd))
	return(score)
}

# project the mean leaf on each PC axis
meanscore=array(0, c(N_harmonics))
meanleaf=array(0,c(len,2))
meanleaf[,]=meanvec[,1,]

for (i in 1:N_harmonics){ 
	meanscore[i] = sum(project.on.PC(meanleaf,i))
	}

# project the mean curves on PC axes :
mscore_xy = array(0,c(N_gentypes,N_classes,N_harmonics))
for (i in 1:N_gentypes){
	for (j in 1:N_classes){
		for (k in 1:N_harmonics){
			mscore_xy[i,j,k] = sum(project.on.PC(mleafs[i,,j,],k))-meanscore[k]
		}
	}
}

PC_vec=harmonicsvec[,1,]
mPC1_vecp=mean_vec+1*PC_vec
score_PC1 = sum(project.on.PC(mPC1_vecp,1))


# Plotting components as perturbations of the mean
# --------------------------------------------

plot.PC = function(pcindex) {	
	PC_vec=2e-1*harmonicsvec[,pcindex,]
	mPC1_vecp=mean_vec+1*PC_vec
	mPC1_vecm=mean_vec-1*PC_vec
	mPC2_vecp=mean_vec+0.5*PC_vec
	mPC2_vecm=mean_vec-0.5*PC_vec

	max_x=max(max(mPC1_vecp[,1]),max(mPC1_vecm[,1]))
	min_x=min(min(mPC1_vecp[,1]),min(mPC1_vecm[,1]))

	max_y=max(max(mPC1_vecp[,2]),max(mPC1_vecm[,2]))
	min_y=min(min(mPC1_vecp[,2]),min(mPC1_vecm[,2]))

	xlim <- c(min_x,max_x)
	ylim <- c(min_y,max_y)

	png(paste(directory,"/PC",pcindex,".png",sep=""))
	plot(mean_vec[,1],mean_vec[,2],pch=10, type='l',xlim=xlim,ylim=ylim,asp=1,xlab="",ylab="")
	title(main=paste("PC",pcindex,"  -  variability percentage ", round(lambda[pcindex]*100,1),"%",sep=""))
	lines(mPC1_vecp[,1],mPC1_vecp[,2],col='green')
	lines(mPC1_vecm[,1],mPC1_vecm[,2],col='red')
	lines(mPC2_vecp[,1],mPC2_vecp[,2],lty=2,col='green')
	lines(mPC2_vecm[,1],mPC2_vecm[,2],lty=2,col='red')
	grid()
	dev.off()
	return(1)
}

for (i in 1:N_harmonics){
	plot.PC(i)
}

print.PCicons=function(pcindex){
	PC_vec=2e-1*harmonicsvec[,pcindex,]
	alpha = 1
	mPC1_vecp=mean_vec+alpha*PC_vec
	mPC1_vecm=mean_vec-alpha*PC_vec
	print.contour(mPC1_vecp,paste(directory,"/PC",pcindex,"p_icon.png",sep=""))
	print.contour(mPC1_vecm,paste(directory,"/PC",pcindex,"m_icon.png",sep=""))
}


for (i in 1:N_harmonics){
	print.PCicons(i)
}


colorcode <- c()
for (ind in 1:N_gentypes){
	colorcode <- c(colorcode, rep(gen_colors[ind], Ngenes[ind]))
}
#shapecode <- c()
#for (ind in 1:N_subclasses){
#	t <- ind%%N_types
#	if (t==0){t <- N_types}
#	shapecode <- c(shapecode, rep(typ_shapes[t], Nsubclasses[ind]))
#}
 
plot.PCplain = function(pcindex1,pcindex2) {    
	png(paste(directory,"/scores_PC",pcindex1,"-",pcindex2,".png",sep=""))
    plot(scores_xy[,pcindex1],scores_xy[,pcindex2], col=colorcode,xlab=paste("PC",pcindex1,sep=""),
    ylab=paste("PC",pcindex2,sep=""))
    title(main=paste("PC",pcindex1,"  -  PC", pcindex2," plane",sep=""))
    legend("topright", inset=.05,  genes ,fill=gen_colors)
    grid()
    dev.off()
    return(1)
}
 
plot.PCplain(1,2) 
