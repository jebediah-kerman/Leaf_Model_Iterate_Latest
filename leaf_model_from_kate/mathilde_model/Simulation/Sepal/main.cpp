
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

int simnumber=2;
int seed=2;
randinit(seed);

string f,numb;
// colors for the plot
real[int] colors(120);
colors=(0.e6:2e6:50e6);

// define the geometric variables and parameters
real radius, plotsize;

// all mechanical properties for the simulation - only file to be modified
include "Meca_Param.cpp"

// initial shape
plotsize=6*radius;
step=0;
real ymin=0;
real ymax=radius;

// compute the mesh
//definition of the domain for the problem 
border out(t=0,pi){x=radius*cos(t);y=radius*sin(t);label=1;} // half circle
border bot(t=-radius,radius){x=t;y=0;label=4;}
border box1(t=-plotsize,plotsize){x=t;y=-plotsize/20;label=3;} // the box for the plot
border box2(t=-plotsize,plotsize){x=t;y=2*plotsize;label=3;}
border box3(t=-plotsize/20,2*plotsize){x=plotsize;y=t;label=3;}
border box4(t=-plotsize/20,2*plotsize){x=-plotsize;y=t;label=3;}
mesh sepal=buildmesh(out(63)+bot(40)); // the sepal mesh
mesh box=buildmesh(box1(1)+box3(1)+box2(-1)+box4(-1)); //the box for the plot
// choose the finite elements
//define the finite element space
fespace fsepal(sepal,P2); // P2 elements for the sepal
fespace fsepal0(sepal,P0);
fespace fsepal1(sepal,P1);

include "HouseHoldFunctions_kh.cpp"
real area;
int nbvertices = sepal.nv;
real[int] ElastVertices(nbvertices);
real[int] DensityVertices;

bool fAActivated;
bool EndSimul;
real CurrElastV;
real CurrElastMean;
real CurrElastSd;

// displacement (ux,uy), test function (wx,wy)
// define the field variables for the finite element space
fsepal ux=0,uy=0,wx,wy,ur;
fsepal ro = d
// a variable used just to plot combinations of the others
fsepal variableplot=0;
//fbox backgroundcolor=4;

area= int2d(sepal)(1.);

// initialization: randomization of elasticity
for (int i=0;i<nbvertices;i++){
	ElastVertices[i] = GetElastFromMecaParam(ElastMean, ElastSd, MinElast);
}

real GrowthFactorConc = CalcGrowthFactorConc(area, MolesGrowthFactor);

// conversion from "vertice-based" data to "x,y-based" data
include "ConvVertice-xy.cpp"

// definition of the Hooke's matrix for orthotropic material in 2D [A1, B, 0; B, A2, 0; 0, 0, C]

/**
"Hooke's Matrix Classic"
Mathilde's original code, E only based on field of elasticity
func A1 = max(Elastxyh(x,y)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1+Anis/2.),0.);
fsepal1 A1h = A1;
func A2 = max(Elastxyh(x,y)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1-Anis/2.),0.);
fsepal1 A2h = A2;
func B12 = rho*sqrt(A1*A2);
fsepal1 B12h = B12;
func C3 = Elastxyh(x,y)/(2.*(1.+nu));
fsepal1 C3h = C3;
**/

/**
"Stupid Growth Factor"
Elasticity increases linearly with position relative to total final height of sepal 

func A1 = max((Elastxyh(x,y)+ElastMean*y/10)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1+Anis/2.),0.);
fsepal1 A1h = A1;
func A2 = max((Elastxyh(x,y)+ElastMean*y/10)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1-Anis/2.),0.);
fsepal1 A2h = A2;
func B12 = rho*sqrt(A1*A2);
fsepal1 B12h = B12;
func C3 = (Elastxyh(x,y)+ElastMean*y/10)/(2.*(1.+nu));
fsepal1 C3h = C3;

**/

/**
"DelConcentration Growth Factor"
func A1 = max(max(Elastxyh(x,y)+alpha*GrowthFactorConc, 0.)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1+Anis/2.),0.);
fsepal1 A1h = A1;
func A2 = max(max(Elastxyh(x,y)+alpha*GrowthFactorConc, 0.)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1-Anis/2.),0.);
fsepal1 A2h = A2;
func B12 = rho*sqrt(A1*A2);
fsepal1 B12h = B12;
func C3 = max(Elastxyh(x,y)+alpha*GrowthFactorConc, 0.)/(2.*(1.+nu));
fsepal1 C3h = C3;
**/

func A1 = max(max(Elastxyh(x,y)+alpha*GrowthFactorConc, 0.)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1+Anis/2.),0.);
fsepal1 A1h = A1;
func A2 = max(max(Elastxyh(x,y)+alpha*GrowthFactorConc, 0.)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1-Anis/2.),0.);
fsepal1 A2h = A2;
func B12 = rho*sqrt(A1*A2);
fsepal1 B12h = B12;
func C3 = max(Elastxyh(x,y)+alpha*GrowthFactorConc, 0.)/(2.*(1.+nu));
fsepal1 C3h = C3;


// !! This matrix is defined in the eigen direction of the material, and if left as it is, allows
// anisotropy only in the x-y orientations
// -> need to rotate this matrix to be general - and allow the anisotropy to not be horizontal or vertical

// computation of the new Hooke's matrix after rotation [A B G; B E F; G F I]
include "RotMat_Comput.cpp"


// define the variationnal equation to solve - developped
// functions describing the elements in the mesh and on the boundary ("on") based on weak form of equation governing deformation due to
// turgor pressure
problem elasticity([ux,uy],[wx,wy])=
	int2d(sepal)((
	ARh*dx(ux)*dx(wx) + ERh*dy(uy)*dy(wy)
	+ BRh*(dx(ux)*dy(wy) + dy(uy)*dx(wx))
	+ 0.5*IRh*(dx(uy)+dy(ux))*(dx(wy)+dy(wx))
	+ GRh*(0.5*(dx(uy)+dy(ux))*dx(wx)+dx(ux)*(dx(wy)+dy(wx)))
	+ FRh*(0.5*(dx(uy)+dy(ux))*dy(wy)+dy(uy)*(dx(wy)+dy(wx)))))
	// N is the normal to the surface
	-int1d(sepal,1)(p*(wx*N.x+wy*N.y))
	+on(4,ux=0, uy=0); // boundary condition to fix the bottom of the structure

elasticity;

real[int] tmp(ux[].n);
// compute the deformed sepal
sepal=movemesh(sepal,[x+ux,y+uy]);
include "VarDispl.cpp"

if (savePic==1 && step%picturestep==0){
	include "SavePlots.cpp"
}
// temporal loop
include "SaveData.cpp"
step +=1;
fAActivated = 0;
EndSimul = EndSim(frontArrest, fAheight, area, MaxArea);
while (step<maxstep && !(EndSimul)){
	cout << "step : " << step << endl;
	cout << "area: " << area << endl;

	ElastVertices.resize(sepal.nv);
	// Elasticity distribution computation
	CurrElastMean = int2d(sepal)(Elastxyh)/int2d(sepal)(1.);
	CurrElastSd = sqrt(int2d(sepal)((Elastxyh-CurrElastMean)^2.)/int2d(sepal)(1.));
	// cout << "((Elastxyh-CurrElastMean)^2.)" << ((Elastxyh-CurrElastMean)^2.) << endl;
	for (int i=0;i<sepal.nv;i++){
		// Elast update
		x = sepal(i).x;
		y = sepal(i).y;
//		CurrElastV = Elastxyh(sepal(i).x, sepal(i).y) *(1+GradientFactor/(1+y));
		CurrElastV = Elastxyh(sepal(i).x, sepal(i).y);	

		ElastVertices(i) = GetElastVert(i, ElastMean, ElastSd, CurrElastMean, CurrElastSd, MinElast, CurrElastV, ElastCoefTimeVar, sepal);
	}

  // Update of the xyh data
	nbTr = sepal(x,y).nuTriangle;


  V0x = (sepal[nbTr(x,y)][0]).x;
  V0y = (sepal[nbTr(x,y)][0]).y;
  V1x = (sepal[nbTr(x,y)][1]).x;
  V1y = (sepal[nbTr(x,y)][1]).y;
  V2x = (sepal[nbTr(x,y)][2]).x;
  V2y = (sepal[nbTr(x,y)][2]).y;
  Elastxyh = Elastxy;
	A1h = A1; A2h = A2; B12h = B12; C3h = C3;
	ARh = AR; BRh = BR; ERh = ER; FRh = FR; GRh = GR; IRh = IR;

	// compute the displacement
	elasticity;


	// displace mesh and variables
	real[int] tmp(ux[].n);
	sepal=movemesh(sepal,[x+ux,y+uy]);
	include "VarDispl.cpp"
	area= int2d(sepal)(1.);
	real[int] Values = FindyLim(sepal, sepal.nv);
	ymin = Values[1];
	ymax = Values[3];

	//recalculate GrowthFactorConc based on new area
	GrowthFactorConc = CalcGrowthFactorConc(area, MolesGrowthFactor);
	// Computation of where is the front arrest limit now
	if (frontArrest){
		if (fAActivated == 1){
			fAheight = fAheight - fAspeed;
		}
		cout << "height front arrest limit " << fAheight << endl;
		if (ymax > frontArrHeightIni){
			fAActivated = 1;
		}
	}
		// save the picture
	if (savePic==1 && step%picturestep==0){
		include "SavePlots.cpp"
	}

	include "SaveData.cpp"
	// Is this step the last step?
	EndSimul = EndSim(frontArrest, fAheight, area, MaxArea);

	step++;
} // end of the loop over the time
