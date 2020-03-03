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
int step, maxstep=200;
// define the minimum height gain between two consecutive steps for the program to terminate
real minYgain = 0.05;
//size of the triangles
real TriangleSize = 1./3.5; 		// To avoid segmentation fault
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
 //real J = .1; 			// Inward flux at the base. It should be 0 < J < 1
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


 /////////////////////////Growth inhibitor related
 int prefbCurve = 4;	// Which response curve should be used? 1 for tanh, 2 for linear, 3 for step, 4 for hill (default)
 real DegRateInh =.1; 	// Degradation rate. It should be  0 < DegRate << 1
 real DInh = .4; 		// Diffusion constant. It should be 0 < D < 1
 //real JInh = .01; 	// Inward flux at the base. It should be 0 < J < 1
 real DensityInhI = 0; 	// Initial density
 real RelElInh = 100.; 	// The ratio beween elasticity for very low and very high density 1 < RelEl
 real RelElFactorInh = log(RelElInh);
 real RhzInh = 2; 		// The density for the transition from high to low density regime.
 real dRhoInh = 6.; 	// sensitivity of the elasticity with respect to the density at the transition from low to high density regime.
 if(prefbCurve == 4){
 	if(dRhoInh > 10){	// Hill function can be substituted by Step function when the power > 10
 		prefbCurve = 3;
 	}
 }


 /////////////////////////Production of the growth inhibitor 2 (from the tip) based on growth factor concentration
 real GI2maxProd = 1.;					// maximum rate of production (at absolute void of the growth inhibitor)
 //real GI2transitionalGFconc = 0.00001; 	// The density of the growth factor for the transition from high to low production of the growth inhibitor
 real GI2transitionalGFconc = 0.04;
 real GI2power = 5.; 					// sensitivity of the growth inhibitor production with respect to the growth factor density.


 // define the anisotropy parameters
 //real Anis = .2;
 real AnisStart = 0.4;
 real AnisEnd = 0.1;
 real Theta = 0; // orientation of the anisotropy

 // Growth Front Arrest:
 bool frontArrest = 0; // boolean: is the simulation stopping because of the growth front arrest or another factor (MaxArea)
 real frontArrHeightIni = 6.; // so that the leaves are better developed
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
real fAElastFactor = 2.;		// The multiplicative factor on Elasticity if a point is beyond arrest front
//real fADiffusionFactor = 0.; // The multiplicative factor on diffusion constant D if a point is beyond arrest front

// if the simulation does not end with the growth front arrest
// at which area the simulation ends
real MaxArea = 1000;

// save/plot the data or not, and the frequency of saving/plotting
int savePic=1;
int picturestep=2;

// whether to scale the plots to screen, or plot actual sizes (specify a boundary)
bool savePicScaled=0;
func savePicBounds = [[-8.,0.],[8.,20.]];

// For the plots, isovalues to be plotted
real[int] visoElast = [0.e7, 1.e7, 2.e7, 3.e7, 4.e7, 5.e7, 6.e7, 7.e7, 8.e7, 9.e7, 10.e7, 1.e10];
real[int] visoGrowthRate = [0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 2.];
real[int] visoGrowthFac = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 10.];
real[int] visoGrowthInh = [0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 100.];
//real[int] visoGrowthInh = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 5., 100.];

// Evaluation of leaf shape
real targetMatureHeightAFInitHeightRatio = 2.;
real targetHeightWidthRatio = 3.;













