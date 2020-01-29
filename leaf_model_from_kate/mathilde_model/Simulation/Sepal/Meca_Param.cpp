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

// define the number of temporal steps
int step, maxstep=500;
//size of the triangles
real TriangleSize = 1./5.; // all figures except fig 2H and fig 7D
//real TriangleSize = 1./1.5; // fig 2H
//real TriangleSize = 1./3.5; // fig 7D
radius=1.1;

// define the mechanical parameters
real p=5.e5;
real nu = 0.45;  // Poisson Coef
// parameter for mechanical properties
real rho = 0.5;

// define the parameters of the gaussian distrib of the elasticity
real ElastMean = 3270000;
// real ElastSd = 0; 2700000// fig 2A; fig 7B
real ElastSd = 0; // all figures except fig 2A and fig 7B
real MinElast = 100000;
// how fast is the change from the
// current value to the new random
// value ElastCoelTimeVar = 1 = Full resampling, = 0 = No resampling
// real ElastCoefTimeVar = 1.; // fig 2D
real ElastCoefTimeVar = 1.; // fig 2A; 2F; 2H; fig 7B
// real ElastCoefTimeVar = 0.1; // fig S2C; S2F-G; fig 7C; 7D
// real ElastCoefTimeVar = 0.9; // fig S2C
// real ElastCoefTimeVar = 0.003; // fig S2C


// define the anisotropy parameters
real Anis = 0.2;
real Theta = 0; // orientation of the anisotropy

// Growth Front Arrest:
bool frontArrest = 1; // boolean: is the simulation stopping because of the growth front arrest or another factor (MaxArea)
real frontArrHeightIni = 3.; // all figures except fig 7B(ii)and fig 7D
// real frontArrHeightIni = 2.7; // fig 7B(ii); 7D
real frontArrHeightIniSD = 0.; // all figures except fig 7B; fig 7C and fig 7D
//real frontArrHeightIniSD = 0.05; // fig 7B(i)
//real frontArrHeightIniSD = 0.5; // fig 7B(ii)
//real frontArrHeightIniSD = 0.08; // fig 7C
//real frontArrHeightIniSD = 0.15; // fig 7D
// Adding Variation
real fAa = randreal1();
real fAb = randreal1();
real fAheight = max(frontArrHeightIni + sqrt(-2*log(fAa))*cos(2*pi*fAb)*frontArrHeightIniSD, 0.01); // Box-Muller transform
real fAspeed = 0.05; // speed of the growth front arrest towards the bottom

// if the simulation does not end with the growth front arrest
// at which area the simulation ends
real MaxArea = 200.0;

// save/plot the data or not, and the frequency of saving/plotting
int savePic=1;
int picturestep=30;

//For "Stupid GF"
real MaxOrganHeight = sqrt(MaxArea);

// For "Concentration GF"
real MolesGrowthFactor = 500;
real alpha = -1;

// For Spatially dependent Growth Factor
//
real DegradationRate_tau = 20;
real DiffusionConstant_D = 60;
real BaseFlux_J = 400;
real InitialDensity_p_i = 30;
