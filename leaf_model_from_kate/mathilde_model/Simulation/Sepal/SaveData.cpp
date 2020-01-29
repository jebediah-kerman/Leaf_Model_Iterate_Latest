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


if (step<10){ numb="__0000";}
if (step<100 && step>9){ numb="__000";}
if (step<1000 && step>99){ numb="__00";}
if (step<10000 && step>999){ numb="__0";}

f = "../Data/" + simnumber+numb+step+"__data.txt";
ofstream fff(f);
fff << "Index_Vertex nbTriangle nbInTr x y ux uy " <<
"Elast " <<
"AR BR IR GR ER FR A1 A2 B12 C3 " <<
"Defxx Defxy Defyx Defyy " << endl;
for (int i=0;i<sepal.nt;i++){
  for (int j=0;j<3;j++){
    fff << sepal[i][j] << " " <<
    i << " " <<
    j << " " <<
    sepal[i][j].x << " " <<
    sepal[i][j].y << " " <<
    ux(sepal[i][j].x,sepal[i][j].y) << " " <<
    uy(sepal[i][j].x,sepal[i][j].y) << " " <<
    Elastxyh(sepal[i][j].x,sepal[i][j].y) << " " <<
    ARh(sepal[i][j].x,sepal[i][j].y) << " " <<
    BRh(sepal[i][j].x,sepal[i][j].y) << " " <<
    IRh(sepal[i][j].x,sepal[i][j].y) << " " <<
    GRh(sepal[i][j].x,sepal[i][j].y) << " " <<
    ERh(sepal[i][j].x,sepal[i][j].y) << " " <<
    FRh(sepal[i][j].x,sepal[i][j].y) << " " <<
    A1h(sepal[i][j].x,sepal[i][j].y) << " " <<
    A2h(sepal[i][j].x,sepal[i][j].y) << " " <<
    B12h(sepal[i][j].x,sepal[i][j].y) << " " <<
    C3h(sepal[i][j].x,sepal[i][j].y) << " " <<
    dx(ux)(sepal[i][j].x,sepal[i][j].y) <<" " <<
    dx(uy)(sepal[i][j].x,sepal[i][j].y) <<" " <<
    dy(ux)(sepal[i][j].x,sepal[i][j].y) <<" " <<
    dy(uy)(sepal[i][j].x,sepal[i][j].y) << endl;

  }
 }
