(* ::Package:: *)

(* ::Title:: *)
(*This package holds several ofter used stuff for FEM analysis.*)


BeginPackage["FEMUtils`"];


(* ::Section:: *)
(*Definition of Public Functions and Parameters.*)


fieldFunction::usage="TODO";
discreteValues::usage="TODO";
nodalCoordiantes3D::usage="TODO";
xi::usage="TODO";
xi1::usage="TODO";
xi2::usage="TODO";
xi3::usage="TODO";
posDof::usage="TODO";
tanDof::usage="TODO";
lenDof::usage="TODO";


(* ::Section:: *)
(*Private package contents.*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Standard Lagrange shape functions in 1D.*)


(* Lagrange shape functions *)
shapeFunctionsLagrange1D[order_] := Module[{xiPoints},
  xiPoints = Flatten[{-1, 1, Drop[Drop[Table[xi, {xi, -1, 1, 2/order}], 1], -1]}];
  Table[  	
   Simplify[
    InterpolatingPolynomial[Transpose[{xiPoints, 
       IdentityMatrix[order + 1][[i]]}], xi]]
   , {i, order + 1}]
  ]


(* ::Subsection:: *)
(*Shape functions for standard elements.*)


(*Default function (throws error).*)
shapeFunctions[var_:None] := (Print[StringForm["Shape function `` is not implemented!", var]];Abort[])
discreteValues[var_:None] := (Print[StringForm["Discrete values for `` is not implemented!", var]];Abort[])
discreteValuesOrder[None] := (Print[StringForm["Discrete values order for `` is not implemented!", None]];Abort[])

(*Calculate the field function with the nodal values.*)
fieldFunction[var_] := shapeFunctions[var] . discreteValues[var][[ discreteValuesOrder[var] ]];

(*If discreteValues is not given use the number of nodes.*)
discreteValues[var_] := Table[posDof[i], {i, Length[shapeFunctions[var]]}];

(*If discreteValuesOrder is not given use the ordering in discreteValues.*)
discreteValuesOrder[var_] := Table[i, {i, Length[discreteValues[var]]}];

(*Dimension of the element.*)
elementDimension[var_] := Module[{field, dim},
    field = fieldFunction[var];
    dim = 0;
    If[!FreeQ[field, xi], dim = 1;];
    If[FreeQ[field, xi1],(Print[StringForm["Error in elementDimension!", None]];Abort[])];
    If[!FreeQ[field, xi2], dim = 2;];
    If[!FreeQ[field, xi3], dim = 3;];
    If[dim==0,(Print[StringForm["Error in elementDimension!", None]];Abort[])];
    dim
    ]

(*Nodal coordinates for 3D elements.*)
nodalCoordiantes3D[var_] := Module[{shapeFunct, solCoordinates},
    shapeFunct=D[fieldFunction[var],{discreteValues[var]}];
    solCoordinates=Table[Solve[Table[If[i==iN,1,0]==shapeFunct[[i]],{i,Length[shapeFunct]}],{xi1,xi2,xi3}],{iN,Length[shapeFunct]}];
    Flatten[{xi1,xi2,xi3}/.solCoordinates,1]
]


(* ::Subsubsection:: *)
(*1D elements.*)


(*Line with linear shape functions.*)
shapeFunctions["line2"] := shapeFunctionsLagrange1D[1];

(*Line with quadratic shape functions.*)
shapeFunctions["line3"] := shapeFunctionsLagrange1D[2];

(*Line with cubic shape functions.*)
shapeFunctions["line4"] := shapeFunctionsLagrange1D[3];

(*Hermite polynomials of 3rd order.*)
discreteValues["hermite3"] := {posDof[1], tanDof[1], posDof[2], tanDof[2]};
shapeFunctions["hermite3"] := {
  1/4 (2 + xi) (1 - xi)^2, 1/4 (1 + xi) (1 - xi)^2 lenDof/2,
  1/4 (2 - xi) (1 + xi)^2, -(1/4) (1 - xi) (1 + xi)^2 lenDof/2
  };


(* ::Section:: *)
(*2D elements.*)


(*Bilinear quad element.*)
discreteValuesOrder["quad4"] := {1, 4, 2, 3};
shapeFunctions["quad4"] := Flatten[ Transpose[{shapeFunctions["line2"]} /. xi -> xi1] . ({shapeFunctions["line2"]} /. xi -> xi2) ];

(*Quadratic rectangular element.*)
discreteValuesOrder["quad9"] := {1, 4, 8, 2, 3, 6, 5, 7, 9};
shapeFunctions["quad9"] := Flatten[ Transpose[{shapeFunctions["line3"]} /. xi -> xi1] . ({shapeFunctions["line3"]} /. xi -> xi2) ];

(*Seredipity element.*)
shapeFunctions["quad8"] := {1/4 ((1 - xi1) (1 - xi2) - (1 - xi1^2) (1 - xi2) - (1 - xi1) (1 - xi2^2)),
    1/4 ((1 + xi1) (1 - xi2) - (1 - xi1^2) (1 - xi2) - (1 + xi1) (1 - xi2^2)),
    1/4 ((1 + xi1) (1 + xi2) - (1 - xi1^2) (1 + xi2) - (1 + xi1) (1 - xi2^2)),
    1/4 ((1 - xi1) (1 + xi2) - (1 - xi1^2) (1 + xi2) - (1 - xi1) (1 - xi2^2)),
    1/2 (1 - xi1^2) (1 - xi2), 1/2 (1 + xi1) (1 - xi2^2), 1/2 (1 - xi1^2) (1 + xi2),
    1/2 (1 - xi1) (1 - xi2^2)};


(* ::Section:: *)
(*3D elements.*)


(*Hex8 element.*)
shapeFunctions["hex8"] := {1/8 (1 - xi1) (1 - xi2) (1 - xi3), 1/8 (1 + xi1) (1 - xi2) (1 - xi3),
  1/8 (1 + xi1) (1 + xi2) (1 - xi3), 
 1/8 (1 - xi1) (1 + xi2) (1 - xi3), 1/8 (1 - xi1) (1 - xi2) (1 + xi3),
  1/8 (1 + xi1) (1 - xi2) (1 + xi3), 
 1/8 (1 + xi1) (1 + xi2) (1 + xi3), 1/8 (1 - xi1) (1 + xi2) (1 + xi3)};


(*Hex8 element.*)
shapeFunctions["hex20"] := {1/8 (1 - xi1) (1 - xi2) (1 - xi3) (-2. - xi1 - xi2 - xi3), 
 1/8 (1 + xi1) (1 - xi2) (1 - xi3) (-2. + xi1 - xi2 - xi3), 
 1/8 (1 + xi1) (1 + xi2) (1 - xi3) (-2. + xi1 + xi2 - xi3), 
 1/8 (1 - xi1) (1 + xi2) (1 - xi3) (-2. - xi1 + xi2 - xi3), 
 1/8 (1 - xi1) (1 - xi2) (1 + xi3) (-2. - xi1 - xi2 + xi3), 
 1/8 (1 + xi1) (1 - xi2) (1 + xi3) (-2. + xi1 - xi2 + xi3), 
 1/8 (1 + xi1) (1 + xi2) (1 + xi3) (-2. + xi1 + xi2 + xi3), 
 1/8 (1 - xi1) (1 + xi2) (1 + xi3) (-2. - xi1 + xi2 + xi3), 
 1/4 (1 - xi1^2) (1 - xi2) (1 - xi3), 
 1/4 (1 + xi1) (1 - xi2^2) (1 - xi3), 
 1/4 (1 - xi1^2) (1 + xi2) (1 - xi3), 
 1/4 (1 - xi1) (1 - xi2^2) (1 - xi3), 
 1/4 (1 - xi1) (1 - xi2) (1 - xi3^2), 
 1/4 (1 + xi1) (1 - xi2) (1 - xi3^2), 
 1/4 (1 + xi1) (1 + xi2) (1 - xi3^2), 
 1/4 (1 - xi1) (1 + xi2) (1 - xi3^2), 
 1/4 (1 - xi1^2) (1 - xi2) (1 + xi3), 
 1/4 (1 + xi1) (1 - xi2^2) (1 + xi3), 
 1/4 (1 - xi1^2) (1 + xi2) (1 + xi3), 
 1/4 (1 - xi1) (1 - xi2^2) (1 + xi3)};


(*Hex27 element.*)
shapeFunctions["hex27"] := {1/8 (-1 + xi1) xi1 (-1 + xi2) xi2 (-1 + xi3) xi3, 
 1/8 xi1 (1 + xi1) (-1 + xi2) xi2 (-1 + xi3) xi3, 
 1/8 xi1 (1 + xi1) xi2 (1 + xi2) (-1 + xi3) xi3, 
 1/8 (-1 + xi1) xi1 xi2 (1 + xi2) (-1 + xi3) xi3, 
 1/8 (-1 + xi1) xi1 (-1 + xi2) xi2 xi3 (1 + xi3), 
 1/8 xi1 (1 + xi1) (-1 + xi2) xi2 xi3 (1 + xi3), 
 1/8 xi1 (1 + xi1) xi2 (1 + xi2) xi3 (1 + xi3), 
 1/8 (-1 + xi1) xi1 xi2 (1 + xi2) xi3 (1 + xi3), 
 1/4 (1 - xi1^2) (-1 + xi2) xi2 (-1 + xi3) xi3, 
 1/4 xi1 (1 + xi1) (1 - xi2^2) (-1 + xi3) xi3, 
 1/4 (1 - xi1^2) xi2 (1 + xi2) (-1 + xi3) xi3, 
 1/4 (-1 + xi1) xi1 (1 - xi2^2) (-1 + xi3) xi3, 
 1/4 (-1 + xi1) xi1 (-1 + xi2) xi2 (1 - xi3^2), 
 1/4 xi1 (1 + xi1) (-1 + xi2) xi2 (1 - xi3^2), 
 1/4 xi1 (1 + xi1) xi2 (1 + xi2) (1 - xi3^2), 
 1/4 (-1 + xi1) xi1 xi2 (1 + xi2) (1 - xi3^2), 
 1/4 (1 - xi1^2) (-1 + xi2) xi2 xi3 (1 + xi3), 
 1/4 xi1 (1 + xi1) (1 - xi2^2) xi3 (1 + xi3), 
 1/4 (1 - xi1^2) xi2 (1 + xi2) xi3 (1 + xi3), 
 1/4 (-1 + xi1) xi1 (1 - xi2^2) xi3 (1 + xi3), 
 1/2 (1 - xi1^2) (1 - xi2^2) (-1 + xi3) xi3, 
 1/2 (1 - xi1^2) (-1 + xi2) xi2 (1 - xi3^2), 
 1/2 xi1 (1 + xi1) (1 - xi2^2) (1 - xi3^2), 
 1/2 (1 - xi1^2) xi2 (1 + xi2) (1 - xi3^2), 
 1/2 (-1 + xi1) xi1 (1 - xi2^2) (1 - xi3^2), 
 1/2 (1 - xi1^2) (1 - xi2^2) xi3 (1 + xi3), (1 - xi1^2) (1 - 
    xi2^2) (1 - xi3^2)};


(*Tet4 element.*)
shapeFunctions["tet4"] := {1 - xi1 - xi2 - xi3, xi1, xi2, xi3};


(*Tet10 element.*)
shapeFunctions["tet10"] := {(-1 + 2 (1 - xi1 - xi2 - xi3)) (1 - xi1 - xi2 - xi3), 
 xi1 (-1 + 2 xi1), xi2 (-1 + 2 xi2), xi3 (-1 + 2 xi3), 
 4 xi1 (1 - xi1 - xi2 - xi3), 4 xi1 xi2, 4 xi2 (1 - xi1 - xi2 - xi3), 
 4 (1 - xi1 - xi2 - xi3) xi3, 4 xi1 xi3, 4 xi2 xi3};


(*End private.*)
End[];

(*End the package.*)
EndPackage[];
