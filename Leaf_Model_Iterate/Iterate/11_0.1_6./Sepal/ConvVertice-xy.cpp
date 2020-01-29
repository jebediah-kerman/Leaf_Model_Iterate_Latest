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

// conversion from "vertice-based" data to "x,y-based" data
// sepal[sepal(x,y).nuTriangle][1] : 1st (2nd) vertex of the triangle containing (x,y)
// meanE(x,y) = area(PAB)/area(ABC)*mE at C + area(PCB)/area(ABC)*mE at A + area(PAC)/area(ABC)*mE at B
// P (x,y), ABC vertices of the triangle
// calcul surface : 1/2*((xB-xA)(yC-yA)-(xC-xA)(yB-yA))
// names: A0 is the triangle NOT containing A
fsepal0 nbTr = sepal(x,y).nuTriangle;
fsepal0 V0x = (sepal[nbTr(x,y)][0]).x;
fsepal0 V0y = (sepal[nbTr(x,y)][0]).y;
fsepal0 V1x = (sepal[nbTr(x,y)][1]).x;
fsepal0 V1y = (sepal[nbTr(x,y)][1]).y;
fsepal0 V2x = (sepal[nbTr(x,y)][2]).x;
fsepal0 V2y = (sepal[nbTr(x,y)][2]).y;
func Ai0 = (abs((V2x(x,y)-x)*(V1y(x,y)-y)-(V1x(x,y)-x)*(V2y(x,y)-y)))/2.;
func Ai1 = (abs((V0x(x,y)-x)*(V2y(x,y)-y)-(V2x(x,y)-x)*(V0y(x,y)-y)))/2.;
func Ai2 = (abs((V0x(x,y)-x)*(V1y(x,y)-y)-(V1x(x,y)-x)*(V0y(x,y)-y)))/2.;
func Elastxy = (ElastVertices[sepal[nbTr(x,y)][0]]*Ai0(x,y)
		  +ElastVertices[sepal[nbTr(x,y)][1]]*Ai1(x,y)
		  +ElastVertices[sepal[nbTr(x,y)][2]]*Ai2(x,y))
                /(sepal[nbTr(x,y)].area);
func Densityxy = (DensityVertices[sepal[nbTr(x,y)][0]]*Ai0(x,y)
		  +DensityVertices[sepal[nbTr(x,y)][1]]*Ai1(x,y)
		  +DensityVertices[sepal[nbTr(x,y)][2]]*Ai2(x,y))
                /(sepal[nbTr(x,y)].area);
fsepal1 Elastxyh = Elastxy;
fsepal1 Densityxyh = Densityxy;
