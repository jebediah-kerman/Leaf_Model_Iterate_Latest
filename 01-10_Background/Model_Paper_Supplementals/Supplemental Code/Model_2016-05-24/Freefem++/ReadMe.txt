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

Meca_Param.cpp:
Only file to be modified.
Contains all the parameters of the simulations that can be tuned.

main.cpp:
Main file. The one to be excecuted in Freefem++ (calls all the other files).


////////


ConvVertice-xy.cpp:
Does the convertion from vertex-based data to x-y based data.

HouseHoldFunctions.cpp:
Contains functions for the simulations.
// EndSim: determines if this step is the final step
// GetFolder: determines which folder to put the plot in
// FindyLim: computes ymin and ymax for a given mesh
// GetElastFromMecaParam: picks an elastic value in the given normal distribution
// GetElastVert: computes the new elastic value of a vertex based on the current
elastic value and the renewal rate

RotMat_Comput.cpp:
Computes the rotation of a 2x2x2x2 matrix C_{ijkl} (the Hooke's matrix, that can
be gotten from [A1, B, 0; B, A2, 0; 0, 0, C]).
The new matrix is:
C_{abcd} = R^{-1}_{ai} R_{jb} R_{kc} R^{-1}_{dl} C_{ijkl}
where R = [cos(Theta), -sin(Theta); sin(Theta), cos(Theta)] is the rotation matrix

Save_Param.cpp:
List of the parameters saved (plots done by Freefem++).

SaveData.cpp:
Save of the numerical data.

SavePlots.cpp:
Save of the plots specified in Save_Params.cpp

VarDispl.cpp:
Deals with the displacement of the data defined on the mesh when the mesh is
modified - when the structure grows
