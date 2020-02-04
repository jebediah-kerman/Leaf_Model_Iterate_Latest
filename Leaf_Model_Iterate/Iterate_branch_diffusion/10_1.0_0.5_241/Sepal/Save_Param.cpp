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


int nbOutput=13;
fsepal1[int] listOutput(nbOutput);
listOutput[0] = Elastxyh;
listOutput[1] = ARh;
listOutput[2] = BRh;
listOutput[3] = IRh;
listOutput[4] = GRh;
listOutput[5] = FRh;
listOutput[6] = A1h;
listOutput[7] = A2h;
listOutput[8] = B12h;
listOutput[9] = C3h;
listOutput[10] = Densityxyh;
listOutput[11] = Dc(y);
listOutput[12] = dx(ux)+dy(uy);


int nbOutputV=1;
fsepal[int] listOutputV(nbOutputV*2);
listOutputV[0] = ux; // fsepal
listOutputV[1] = uy; // fsepal

string[int] listNamesOutput(nbOutput);
listNamesOutput[0] = "Elast";
listNamesOutput[1] = "AR";
listNamesOutput[2] = "BR";
listNamesOutput[3] = "IR";
listNamesOutput[4] = "GR";
listNamesOutput[5] = "FR";
listNamesOutput[6] = "A1";
listNamesOutput[7] = "A2";
listNamesOutput[8] = "B12";
listNamesOutput[9] = "C3";
listNamesOutput[10] = "Density";
listNamesOutput[11] = "DiffusionConst";
listNamesOutput[12] = "GrowthRate";

string[int] listNamesOutputV(nbOutputV);
listNamesOutputV[0] = "Displ";
