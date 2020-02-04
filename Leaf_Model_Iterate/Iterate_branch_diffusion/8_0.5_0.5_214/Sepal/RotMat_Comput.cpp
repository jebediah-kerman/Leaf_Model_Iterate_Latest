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


// Computation of the new matrix after rotation: [A B G; B E F; G F I]
func AR = A1h*(cos(Theta))^4 + 2.*B12h*((cos(Theta))^2)*((sin(Theta))^2) + A2h*(sin(Theta))^4 + C3h*(sin(2*Theta))^2;
fsepal1 ARh = AR;
func BR = (A1h + A2h + 6.*B12h - 4.*C3h - (A1h + A2h - 2.*B12h - 4.*C3h)*cos(4*Theta))/8.;
fsepal1 BRh = BR;
func ER = A2h*(cos(Theta))^4 + 2.*B12h*((cos(Theta))^2)*((sin(Theta))^2) + A1h*(sin(Theta))^4 + C3h*(sin(2.*Theta))^2;
fsepal1 ERh = ER;
func FR = ((-A1h + A2h + (A1h + A2h - 2.*B12h - 4.*C3h)*cos(2.*Theta))*sin(-2.*Theta))/4.;
fsepal1 FRh = FR;
func GR = -((A1h - A2h + (A1h + A2h - 2.*B12h - 4.*C3h)*cos(2.*Theta))*sin(-2.*Theta))/4.;
fsepal1 GRh = GR;
func IR = (A1h + A2h - 2.*B12h + 4.*C3h - (A1h + A2h - 2.*B12h - 4.*C3h)*cos(4.*Theta))/8.;
fsepal1 IRh = IR;
