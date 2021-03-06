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

// Note by Shuyao on 01/16/20:
// The parameters are adjusted to reflect the ones the authors used in Figure 7C, which is a WT-like model

// define the number of temporal steps
int step, maxstep=1000;
//size of the triangles
real TriangleSize = 1./4.; 		// To avoid segmentation fault
//real TriangleSize = 1./5.; 	// all figures except fig 2H and fig 7D
//real TriangleSize = 1./1.5; 	// fig 2H
//real TriangleSize = 1./3.5; 	// fig 7D
radius=1;	// initial radius of leaf primordium

// define the mechanical parameters
real p=5.e5;		// Turgor pressure
real nu = 0.45;  	// Poisson's ratio
// parameter for mechanical properties
real rho = .5;		// Non-dimensional modulus for the elasiticity matrix to be well-defined

// define the parameters of the gaussian distribution of the elasticity
real ElastMean = 3270000;
// real ElastSd = 0; 	// fig 2A; fig 7B
real ElastSd = 2700000; // all figures except fig 2A and fig 7B
real MinElast = 100000;
// how fast is the change from the
// current value to the new random
// value ElastCoelTimeVar = 1 = Full resampling, = 0 = No resampling
 real ElastCoefTimeVar = 1.; 		// fig 2D
// real ElastCoefTimeVar = 0.; 		// fig 2A; 2F; 2H; fig 7B
// real ElastCoefTimeVar = 0.1; 		// fig S2C; S2F-G; fig 7C; 7D
// real ElastCoefTimeVar = 0.9; 	// fig S2C
// real ElastCoefTimeVar = 0.003; 	// fig S2C

/////////////////////////Growth factor related
 int prefaCurve = 4;	// Which response curve should be used? 1 for tanh, 2 for linear, 3 for step, 4 for hill (default)
 real DegRate =.1; 		// Degradation rate. It should be  0 < DegRate << 1
 real D = .1; 			// Diffusion constant. It should be 0 < D < 1
 real J = .01; 			// Inward flux at the base. It should be 0 < J < 1
 real DensityI = 0; 	// Initial density
 real RelEl= 10.; 		// The ratio beween elasticity for very low and very high density 1 < RelEl
 real RelElFactor = log(RelEl)/2;
 real Rhz = .03; 		// The density for the transition from high to low density regime.
 real dRho = 3.; 		// sensitivity of the elasticity with respect to the density at the transition from low to high density regime.
 if(prefaCurve == 4){
 	if(dRho > 10){		// Hill function can be substituted by Step function when the power > 10
 		prefaCurve = 3;
 	}
 }


 // define the anisotropy parameters
 real Anis = .2;
 real Theta = 0; // orientation of the anisotropy

 // Growth Front Arrest:
 bool frontArrest = 1; // boolean: is the simulation stopping because of the growth front arrest or another factor (MaxArea)
 //real frontArrHeightIni = 6.; // so that the leaves are better developed
 // real frontArrHeightIni = 3.; // all figures except fig 7B(ii)and fig 7D
 // real frontArrHeightIni = 2.7; // fig 7B(ii); 7D
 // real frontArrHeightIniSD = 0.; // all figures except fig 7B; fig 7C and fig 7D
 //real frontArrHeightIniSD = 0.05; // fig 7B(i)
//real frontArrHeightIniSD = 0.5; // fig 7B(ii)
real frontArrHeightIniSD = 0.08; // fig 7C
//real frontArrHeightIniSD = 0.15; // fig 7D
// Adding Variation
real fAa = randreal1();
real fAb = randreal1();
real fAheight = max(frontArrHeightIni + sqrt(-2*log(fAa))*cos(2*pi*fAb)*frontArrHeightIniSD, 0.01); // Box-Muller transform
frontArrHeightIni = fAheight;
real fAspeed = 0.1; // speed of the growth front arrest towards the bottom
real fAElastFactor = 2.;	// The multiplicative factor on Elasticity if a point is beyond arrest front

// if the simulation does not end with the growth front arrest
// at which area the simulation ends
real MaxArea = 100;

// save/plot the data or not, and the frequency of saving/plotting
int savePic=1;
int picturestep=5;

// whether to scale the plots to screen, or plot actual sizes (specify a boundary)
bool savePicScaled=0;
func savePicBounds = [[-10.,0.],[10.,20.]];

// For the plots, isovalues to be plotted
real[int] visoElast = [0.e6, 5.e6, 10.e6, 15.e6, 20.e6, 25.e6, 30.e6, 35.e6, 40.e6, 45.e6, 50.e6, 55.e6, 60.e6, 65.e6, 70.e6, 75.e6, 80.e6];

// Evaluation of leaf shape
real targetMatureHeightAFInitHeightRatio = 2.;
real targetHeightWidthRatio = 3.;













