int simnumber=17;
int seed=17;
real Rhz=0.05;
real junk2=1.;
real junk3=1.;

// Initialize random number generation using seed, which is D from Repeat.sh
randinit(seed);

string f,numb;
// colors for the plot
real[int] colors(120);
colors=(0.e6:2e6:50e6);

// define the geometric variables and parameters
real radius, plotsize;

// all mechanical properties for the simulation - only file to be modified
include "Meca_Param.cpp"
include "HouseHoldFunctions.cpp"

// initial shape
plotsize=6*radius;
step=0;
real ymin=0;
real ymax=radius;
real ygain=ymax;
real Anis=getAnis(ymax);

// compute the mesh
border out(t=0,pi){x=radius*cos(t);y=radius*sin(t);label=1;} // half circle
border bot(t=-radius,radius){x=t;y=0;label=4;}
border box1(t=-plotsize,plotsize){x=t;y=-plotsize/20;label=3;} // the box for the plot
border box2(t=-plotsize,plotsize){x=t;y=2*plotsize;label=3;}
border box3(t=-plotsize/20,2*plotsize){x=plotsize;y=t;label=3;}
border box4(t=-plotsize/20,2*plotsize){x=-plotsize;y=t;label=3;}
mesh sepal=buildmesh(out(63)+bot(40)); // the sepal mesh
mesh box=buildmesh(box1(1)+box3(1)+box2(-1)+box4(-1)); //the box for the plot
// choose the finite elements
fespace fsepal(sepal,P2); // P2 elements for the sepal
fespace fsepal0(sepal,P0);
fespace fsepal1(sepal,P1);


real area;
int nbvertices = sepal.nv;
real[int] ElastVertices(nbvertices);
real[int] DensityVertices(nbvertices);// added by AF
real[int] DensityInhVertices(nbvertices);// added by Shuyao
bool fAActivated;
bool EndSimul;
real CurrElastV;
real CurrElastMean;
real CurrElastSd;
real CurrDensityV;
real CurrDensityInhV;
// displacement (ux,uy), test function (wx,wy)
fsepal ux=0,uy=0,wx,wy,ur;
fsepal Rho=DensityI,RhoInh=DensityInhI,w,wInh;
// a variable used just to plot combinations of the others
fsepal variableplot=0;
//fbox backgroundcolor=4;

area= int2d(sepal)(1.);

// initialization: randomization of elasticity
for (int i=0;i<nbvertices;i++){
	//ElastVertices[i] = prefa(DensityI, RelElFactor, Rhz, dRho, prefaCurve)*prefb(DensityInhI, RelElFactorInh, RhzInh, dRhoInh, prefbCurve)*GetElastFromMecaParam(ElastMean, ElastSd, MinElast);
	ElastVertices[i] = prefab(DensityI, 0., RelElFactor, Rhz, dRho, prefaCurve)*GetElastFromMecaParam(ElastMean, ElastSd, MinElast);
	// ElastVertices[i] = GetElastFromMecaParam(ElastMean, ElastSd, MinElast);
	DensityVertices[i] = DensityI;
	DensityInhVertices[i] = DensityInhI;
}

// conversion from "vertice-based" data to "x,y-based" data
include "ConvVertice-xy.cpp"

// definition of the Hooke's matrix for orthotropic material in 2D [A1, B, 0; B, A2, 0; 0, 0, C]
func A1 = max(Elastxyh(x,y)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1+Anis/2.),0.);

fsepal1 A1h = A1;
func A2 = max(Elastxyh(x,y)*(1.-nu)/((1.+nu)*(1.-2.*nu))*(1-Anis/2.),0.);

fsepal1 A2h = A2;
func B12 = rho*sqrt(A1*A2);
fsepal1 B12h = B12;
func C3 = Elastxyh(x,y)/(2.*(1.+nu));
fsepal1 C3h = C3;
// !! This matrix is defined in the eigen direction of the material, and if left as it is, allows
// anisotropy only in the x-y orientations
// -> need to rotate this matrix to be general - and allow the anisotropy to not be horizontal or vertical

// computation of the new Hooke's matrix after rotation [A B G; B E F; G F I]
include "RotMat_Comput.cpp"



// define the variationnal equation to solve - developped
problem elasticity([ux,uy],[wx,wy])=
	int2d(sepal)((
	ARh*dx(ux)*dx(wx) + ERh*dy(uy)*dy(wy)
	+ BRh*(dx(ux)*dy(wy) + dy(uy)*dx(wx))
	+ 0.5*IRh*(dx(uy)+dy(ux))*(dx(wy)+dy(wx))
	+ GRh*(0.5*(dx(uy)+dy(ux))*dx(wx)+dx(ux)*(dx(wy)+dy(wx)))
	+ FRh*(0.5*(dx(uy)+dy(ux))*dy(wy)+dy(uy)*(dx(wy)+dy(wx)))))
	-int1d(sepal,1)(p*(wx*N.x+wy*N.y))
	+on(4,ux=0, uy=0); // boundary condition to fix the bottom of the structure
//AF problem writen by AF dillution to be added

problem density(Rho,w)=
	int2d(sepal)(
		Rho*w
		+ D*(dx(Rho) * dx(w)+ dy(Rho) * dy(w))
	)
	+ int2d(sepal)(
		Rho*(dx(ux)+dy(uy))*w  // Ne marche pas?
	)
//    + int1d(sepal,1)(0*w)
	- int1d(sepal,4)(
		J*w
	)
	- int2d(sepal)(
		(1-DegRate)*Densityxyh*w
	);

// densityInh problem duplicated from above, by Shuyao
// In this version, the growth inhibitor is one being produced whenever growth factor concentration is low, mimicing miR396
problem densityInh(RhoInh,wInh)=
	int2d(sepal)(
		RhoInh*wInh
		+ DInh*(dx(RhoInh) * dx(wInh)+ dy(RhoInh) * dy(wInh))
	)
	+ int2d(sepal)(
		RhoInh*(dx(ux)+dy(uy))*wInh  // Ne marche pas?
	)
//    + int1d(sepal,1)(0*w)
	//- int1d(sepal,4)(
	//	JInh*wInh
	//)
	- int2d(sepal)(
		getGI2fromGF(Rho)*wInh
	)
	- int2d(sepal)(
		(1-DegRateInh)*DensityInhxyh*wInh
	);

elasticity;
density;
densityInh;

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
EndSimul = EndSim(frontArrest, fAheight, area, MaxArea, ygain);
while (step<maxstep && !(EndSimul)){
	cout << "step : " << step << endl;
	cout << "area: " << area << endl;
	cout << "height: " << ymax << endl;
	cout << "ygain: " << ygain << endl;


	ElastVertices.resize(sepal.nv);
	DensityVertices.resize(sepal.nv);
	DensityInhVertices.resize(sepal.nv);
	// Elasticity distribution computation
	CurrElastMean = int2d(sepal)(Elastxyh)/int2d(sepal)(1.);
	CurrElastSd = int2d(sepal)((Elastxyh-CurrElastMean)^2.)/int2d(sepal)(1.);

	for (int i=0;i<sepal.nv;i++){
		
		// Elast update
		CurrElastV = Elastxyh(sepal(i).x, sepal(i).y);
		CurrDensityV = Rho(sepal(i).x, sepal(i).y);
		CurrDensityInhV = RhoInh(sepal(i).x, sepal(i).y);
		//ElastVertices(i) = prefa(CurrDensityV, RelElFactor, Rhz, dRho, prefaCurve)*prefb(CurrDensityInhV, RelElFactorInh, RhzInh, dRhoInh, prefbCurve)*GetElastVert(i, ElastMean, ElastSd, CurrElastMean, CurrElastSd, MinElast, CurrElastV, ElastCoefTimeVar, sepal);
		ElastVertices(i) = prefab(CurrDensityV, CurrDensityInhV, RelElFactor, Rhz, dRho, prefaCurve)*GetElastVert(i, ElastMean, ElastSd, CurrElastMean, CurrElastSd, MinElast, CurrElastV, ElastCoefTimeVar, sepal);
		DensityVertices(i) = CurrDensityV;
		DensityInhVertices(i) = CurrDensityInhV;

		// If a vertice is beyond the arrest front
		if(sepal(i).y > fAheight){
			ElastVertices(i) = ElastVertices(i) * fAElastFactor;	// Increase its elasticity
		}
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
	Densityxyh = Densityxy;
	DensityInhxyh = DensityInhxy;

	// compute the displacement
	elasticity;
	density;
	densityInh;

	// displace mesh and variables
	real[int] tmp(ux[].n);
	sepal=movemesh(sepal,[x+ux,y+uy]);
	include "VarDispl.cpp"
	area= int2d(sepal)(1.);
	real[int] Values = FindyLim(sepal, sepal.nv);
	ymin = Values[1];
	ygain = Values[3]-ymax;
	ymax = Values[3];
	Anis=getAnis(ymax);


	// Computation of where is the front arrest limit now
	if (frontArrest){
		if (fAActivated == 1){
			fAheight = fAheight - fAspeed;
		}
		cout << "height front arrest limit " << fAheight << endl;
		if (fAActivated == 0){
			if (ymax > frontArrHeightIni){
				fAActivated = 1;
			}
		}
	}
		// save the picture
	if (savePic==1 && step%picturestep==1){
		include "SavePlots.cpp"
	}

	include "SaveData.cpp"
	// Is this step the last step?
	EndSimul = EndSim(frontArrest, fAheight, area, MaxArea, ygain);

	step++;
} // end of the loop over the time



// Calculate the criteria to prepare for output
real[int] Values = FindxLim(sepal, sepal.nv);
real CritAFInit = pow(ymax / frontArrHeightIni - targetMatureHeightAFInitHeightRatio, 2);
real CritHWRatio = pow(ymax / (Values[3]-Values[1]) - targetHeightWidthRatio, 2);


//Remplissage du fichier pour repeter meme s'il y a des erreurs
string ferr;
ferr = "../../Output_Summary/good_output.txt";
ofstream fferr(ferr, append);











fferr << "17\t0.05\t1.\t1.\t" << CritAFInit << "\t" << CritHWRatio << endl;
