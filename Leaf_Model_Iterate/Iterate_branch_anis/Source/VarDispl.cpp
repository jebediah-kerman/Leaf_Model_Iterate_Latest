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

// displace the variables
tmp=ux[]; ux=0; ux[]=tmp;
tmp=uy[]; uy=0; uy[]=tmp;
tmp=Rho[]; Rho=0; Rho[]=tmp;
tmp=nbTr[]; nbTr=0; nbTr[]=tmp;
tmp=V0x[]; V0x=0; V0x[]=tmp;
tmp=V0y[]; V0y=0; V0y[]=tmp;
tmp=V1x[]; V1x=0; V1x[]=tmp;
tmp=V1y[]; V1y=0; V1y[]=tmp;
tmp=V2x[]; V2x=0; V2x[]=tmp;
tmp=V2y[]; V2y=0; V2y[]=tmp;
tmp=Densityxyh[]; Densityxyh=0; Densityxyh[]=tmp;
tmp=Elastxyh[]; Elastxyh=0; Elastxyh[]=tmp;
tmp=A1h[]; A1h=0; A1h[]=tmp;
tmp=A2h[]; A2h=0; A2h[]=tmp;
tmp=B12h[]; B12h=0; B12h[]=tmp;
tmp=C3h[]; C3h=0; C3h[]=tmp;
tmp=ARh[];ARh=0; ARh[]=tmp;
tmp=BRh[];BRh=0; BRh[]=tmp;
tmp=ERh[];ERh=0; ERh[]=tmp;
tmp=FRh[];FRh=0; FRh[]=tmp;
tmp=GRh[];GRh=0; GRh[]=tmp;
tmp=IRh[];IRh=0; IRh[]=tmp;
// re-mesh the sepal
//plot(sepal);
sepal=adaptmesh(sepal,TriangleSize,IsMetric=1,nbvx=10000,verbosity=1);
//plot(sepal);
sepal=adaptmesh(sepal,TriangleSize,IsMetric=1,nbvx=10000,verbosity=1);
//plot(sepal);
sepal=adaptmesh(sepal,TriangleSize,IsMetric=1,nbvx=10000,verbosity=1);
//plot(sepal);
// delete the variables on the old mesh
ux=ux;uy=uy;
Rho=Rho;
nbTr=nbTr;V0x=V0x;V0y=V0y;V1x=V1x;V1y=V1y;V2x=V2x;V2y=V2y;
Elastxyh=Elastxyh;
Densityxyh=Densityxyh;
A1h=A1();A2h=A2();B12h=B12;C3h=C3;
ARh=AR;BRh=BR;ERh=ER;FRh=FR;GRh=GR;IRh=IR;
