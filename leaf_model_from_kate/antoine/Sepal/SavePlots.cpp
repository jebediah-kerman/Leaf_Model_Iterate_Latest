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


include "Save_Param.cpp"

string numb;
string f;
if (step<10){ numb="\\0000";}
if (step<100 && step>9){ numb="\\000";}
if (step<1000 && step>99){ numb="\\00";}
if (step<10000 && step>999){ numb="\\0";}
string folder;

// cout << std::chrono::system_clock::now() << endl;

// first, saving of the mesh
f="../Plot/Mesh/"+simnumber+"_Mesh"+numb+step+".eps";
plot(sepal,ps=f);
// then, saving of the scalar plots
for (int i=0;i<nbOutput;i++){
  //cout << i << listNamesOutput[i] << endl;
  folder = listNamesOutput[i];
  folder = GetFolder(folder);
  f="../Plot/"+folder+"/"+simnumber+"_"+listNamesOutput[i]+numb+step+".eps";
  plot(listOutput[i],ps=f,fill=1,grey=1,value=1);
 }
// then, saving the vectorial plots
for (int i=0;i<nbOutputV;i++){
  folder = listNamesOutputV[i];
  folder = GetFolder(folder);
  f="../Plot/"+folder+"/"+simnumber+"_"+listNamesOutputV[i]+numb+step+".eps";
  plot([listOutputV[2*i],listOutputV[2*i+1]],ps=f,value=1);
 }
