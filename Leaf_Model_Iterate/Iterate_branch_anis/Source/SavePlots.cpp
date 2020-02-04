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
if (step<10){ numb="0000";}
if (step<100 && step>9){ numb="000";}
if (step<1000 && step>99){ numb="00";}
if (step<10000 && step>999){ numb="0";}
string folder;

// cout << std::chrono::system_clock::now() << endl;

// first, saving of the mesh
f="../Plot/Mesh/"+simnumber+"_Mesh"+numb+step+".eps";
if (savePicScaled){
	plot(sepal, ps=f);
} else{
	plot(sepal, ps=f, bb=savePicBounds);
}

// then, saving of the scalar plots
// Elasticity
  folder = listNamesOutput[0];
  folder = GetFolder(folder);
  f="../Plot/"+folder+"/"+simnumber+"_"+listNamesOutput[0]+numb+step+".eps";
  if (savePicScaled){
    plot(listOutput[0],ps=f,fill=1,grey=1,value=1, viso=visoElast);
  } else{
    plot(listOutput[0],ps=f,fill=1,grey=1,value=1, viso=visoElast, bb=savePicBounds);
  }


// Others
for (int i=1;i<nbOutput;i++){
  //cout << i << listNamesOutput[i] << endl;
  folder = listNamesOutput[i];
  folder = GetFolder(folder);
  f="../Plot/"+folder+"/"+simnumber+"_"+listNamesOutput[i]+numb+step+".eps";
  if (savePicScaled){
    plot(listOutput[i],ps=f,fill=1,grey=1,value=1);
  } else{
    plot(listOutput[i],ps=f,fill=1,grey=1,value=1, bb=savePicBounds);
  }

 }
// then, saving the vectorial plots
for (int i=0;i<nbOutputV;i++){
  folder = listNamesOutputV[i];
  folder = GetFolder(folder);
  f="../Plot/"+folder+"/"+simnumber+"_"+listNamesOutputV[i]+numb+step+".eps";
  if (savePicScaled) {
    plot([listOutputV[2*i],listOutputV[2*i+1]],ps=f,value=1);
  } else{
    plot([listOutputV[2*i],listOutputV[2*i+1]],ps=f,value=1, bb=savePicBounds);
  }

 }



