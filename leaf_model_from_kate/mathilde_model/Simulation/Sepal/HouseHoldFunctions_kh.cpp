//////////////////////////////////////////////////
//////////////////////////////////////////////////
//                                              //
//                 Freefem++ script             //
//                                              //
//  author: Mathilde Dumond                     //
//                                              //
//  Supplemental Information of the paper:      //
//  Variable cell growth yields reproducible    //
//  organ development through spatiotemporal    //
//  averaging                                   //
//                                              //
//  Hong et al.                                 //
//  Dev Cell 2016                               //
//                                              //
//////////////////////////////////////////////////
//////////////////////////////////////////////////

// File with all the house hold functions

// EndSim: determines if this step is the final step
// GetFolder: determines which folder to put the plot in
// FindyLim: computes ymin and ymax for a given mesh
// GetElastFromMecaParam: picks an elastic value in the given normal distribution
// GetElastVert: computes the new elastic value of a vertex based on the current elastic value and the renewal rate

// Function determining if the step is the last step
func bool EndSim(bool ROSFront, real ROSheight, real areacurr, real MaxArea){
	bool EndSimul = 0;
	if (ROSFront){
		if (ROSheight < 0.){
			EndSimul=1;
		}
	} else {
		if (areacurr > MaxArea){
			EndSimul=1;
		}
	}
	return EndSimul;
}

// Function determining which folder to put the plot in
func string GetFolder(string folder){
	string folderOut = folder;
	if ((folder=="AR") | (folder=="BR") |
	(folder=="IR") | (folder=="GR") | (folder=="ER") |
	(folder=="FR")) {
		folderOut = "RotMat";
	}
	else if ((folder=="A1") | (folder=="A2") | (folder=="B12") |
	(folder=="C3")) {
		folderOut = "Matrix";
	}
	return folderOut;
}

// Computes ymax, ymin
func real[int] FindyLim(mesh sepal, int nbvertices){
	real[int] Values(4); //imin, ymin, imax, ymax
	Values[1] = 0;
	Values[3] = 0;
	for (int i=0;i<nbvertices;i++){
		if (sepal(i).y>Values[3]){
			Values[3] = sepal(i).y;
			Values[2] = i;
		}
		if (sepal(i).y<Values[1]){
			Values[1] = sepal(i).y;
			Values[0] = i;
		}
	}
	return Values;
}

// Get new Elasticity value
func real GetElastFromMecaParam(real ElastMean, real ElastSd, real MinElast){
	real U = randreal1();
	real V = randreal1();
	real NewElast = ElastMean + sqrt(-2*log(U))*cos(2*pi*V)*ElastSd; // Box-Muller transform
	while (NewElast<MinElast){
		U = randreal1();
		V = randreal1();
		NewElast = ElastMean + sqrt(-2*log(U))*cos(2*pi*V)*ElastSd; // Box-Muller transform
	}
	return NewElast;
}

// Get new Elasticity value
// same as the function before but if no MinElast argument is given, then no while loop
// (to be used in a function where there is the while loop)
func real GetElastFromMecaParam(real ElastMean, real ElastSd){
	real U = randreal1();
	real V = randreal1();
	real NewElast = ElastMean + sqrt(-2*log(U))*cos(2*pi*V)*ElastSd; // Box-Muller transform
	return NewElast;
}

// Computes New ElastVertices with or without memory (ElastCoefTimeVar)
// for each vertex

// assign new Elasticity w effects from Growth Arrest front varied
func real GetElastFromGrowthFront(real ElastMean, real ElastCoefTimeVar, real CurrElastMean, real ElastSd, real CurrElastSd, real CurrElastV, real MinElast){
	real AimElastMean = ElastMean - (1.-ElastCoefTimeVar)*CurrElastMean;
	real AimElastSd = sqrt(ElastSd^2. - (1.-ElastCoefTimeVar)^2.*CurrElastSd^2.);
	real AimElast = GetElastFromMecaParam(AimElastMean, AimElastSd);
	real NewElast = CurrElastV*(1.-ElastCoefTimeVar) + ElastCoefTimeVar*AimElast;
	// compute the NewElast
	while (NewElast<MinElast){
		AimElast = GetElastFromMecaParam(AimElastMean, AimElastSd);
		NewElast = CurrElastV*(1.-ElastCoefTimeVar) + ElastCoefTimeVar*AimElast;
	}
	return NewElast;
}

func real CalcGrowthFactorConc(real area, real MolesGrowthFactor){
	return MolesGrowthFactor / area;
}

func real GetElastVert(int i, real ElastMean, real ElastSd, real CurrElastMean, real CurrElastSd, real MinElast, real CurrElastV, real ElastCoefTimeVar, mesh sepal){
	// get AimElast
	// first, get the parameters of the aim distrib
	real NewElast;
	if (ElastCoefTimeVar == 1.){
		// no memory
		NewElast = GetElastFromMecaParam(ElastMean, ElastSd, MinElast);
	} else if (ElastCoefTimeVar == 0.){
		// full memory
		NewElast = CurrElastV;
	} else {
		NewElast = GetElastFromGrowthFront(ElastMean, ElastCoefTimeVar, CurrElastMean, ElastSd, CurrElastSd, CurrElastV, MinElast);
	}
	return NewElast;
}
