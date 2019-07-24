(* ::Package:: *)

(* ::Title:: *)
(*This package holds several ofter used stuff for FEM analysis.*)


BeginPackage["FEMUtils`"];


(* ::Section:: *)
(*Definition of Public Functions and Parameters.*)


fieldFunction::usage = "TODO";
shapeFunctions::usage = "TODO";
discreteValues::usage = "TODO";
discreteValuesElement::usage = "TODO";
discreteValuesNode::usage = "TODO";
nodalCoordiantes3D::usage = "TODO";
xi::usage = "TODO";
xi1::usage = "TODO";
xi2::usage = "TODO";
xi3::usage = "TODO";
posDof::usage = "TODO";
tanDof::usage = "TODO";
lenDof::usage = "TODO";


(* ::Section:: *)
(*Private package contents.*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Standard Lagrange shape functions in 1D.*)


(* Lagrange shape functions *)
shapeFunctionsInternalLagrange1D[order_] :=
    Module[ {xiPoints},
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
shapeFunctionsInternal[var_:None] :=
    (Print[StringForm["Shape function `` is not implemented!", var]];
     Abort[])
discreteValues[var_:None] :=
    (Print[StringForm["Discrete values for `` is not implemented!", var]];
     Abort[])
discreteValuesElement[var_:None] :=
    (Print[StringForm["Discrete values for `` is not implemented!", var]];
     Abort[])
discreteValuesNode[var_:None] :=
    (Print[StringForm["Discrete values for `` is not implemented!", var]];
     Abort[])
discreteValuesOrder[None] :=
    (Print[StringForm["Discrete values order for `` is not implemented!", None]];
     Abort[])

(*Function to get the actual shape functions*)
shapeFunctions[var_] :=
    D[fieldFunction[var], {discreteValues[var]}];

(*Calculate the field function with the nodal values.*)
fieldFunction[var_] :=
    shapeFunctionsInternal[var] . discreteValues[var][[ discreteValuesOrder[var] ]];

(*Get the field function and discrete values for multidimensional fields.*)
replaceDiscreteValuesDim[nDim_Integer] := {varName_[i_] :> varName[dim, i]};
fieldFunction[var_, nDim_Integer] :=
    Module[ {},
        Table[
            fieldFunction[var] /. replaceDiscreteValuesDim[nDim]
            ,{dim, nDim}]
    ]
discreteValues[var_, nDim_Integer] := Table[
	discreteValues[var] /. replaceDiscreteValuesDim[dim]
	,{dim, nDim}];

(*If discreteValues is not given use the number of nodes.*)
discreteValues[var_] :=
    Table[posDof[i], {i, Length[shapeFunctionsInternal[var]]}];

(*If discreteValuesOrder is not given use the ordering in discreteValues.*)
discreteValuesOrder[var_] :=
    Table[i, {i, Length[discreteValues[var]]}];

(*Dimension of the element.*)
elementDimension[var_] :=
    Module[ {field, dim},
        field = fieldFunction[var];
        dim = 0;
        If[ !FreeQ[field, xi],
            dim = 1;
        ];
        If[ FreeQ[field, xi1],
            (Print[StringForm["Error in elementDimension!", None]];
             Abort[])
        ];
        If[ !FreeQ[field, xi2],
            dim = 2;
        ];
        If[ !FreeQ[field, xi3],
            dim = 3;
        ];
        If[ dim==0,
            (Print[StringForm["Error in elementDimension!", None]];
             Abort[])
        ];
        dim
    ]

(*Nodal coordinates for 3D elements.*)
nodalCoordiantes3D[var_] :=
    Module[ {shapeFunct, solCoordinates},
        shapeFunct = D[fieldFunction[var],{discreteValues[var]}];
        solCoordinates = Table[Solve[Table[If[ i==iN,
                                               1,
                                               0
                                           ]==shapeFunct[[i]],{i,Length[shapeFunct]}],{xi1,xi2,xi3}],{iN,Length[shapeFunct]}];
        Flatten[{xi1,xi2,xi3}/.solCoordinates,1]
    ]


(* ::Subsubsection:: *)
(*1D elements.*)


(*Line with linear shape functions.*)
discreteValuesElement["line2"] := {};
discreteValuesNode["line2"] := {1};
shapeFunctionsInternal["line2"] :=
    shapeFunctionsInternalLagrange1D[1];

(*Line with quadratic shape functions.*)
discreteValuesElement["line3"] := {3};
discreteValuesNode["line3"] := {1};
shapeFunctionsInternal["line3"] :=
    shapeFunctionsInternalLagrange1D[2];

(*Line with cubic shape functions.*)
discreteValuesElement["line4"] := {3, 4};
discreteValuesNode["line4"] := {1};
shapeFunctionsInternal["line4"] :=
    shapeFunctionsInternalLagrange1D[3];

(*Hermite polynomials of 3rd order.*)
discreteValues["hermite3"] :=
    {posDof[1], tanDof[1], posDof[2], tanDof[2]};
discreteValuesElement["hermite3"] := {};
discreteValuesNode["hermite3"] := {1, 2};
shapeFunctionsInternal["hermite3"] :=
    {
    1/4 (2 + xi) (1 - xi)^2, 1/4 (1 + xi) (1 - xi)^2 lenDof/2,
    1/4 (2 - xi) (1 + xi)^2, -(1/4) (1 - xi) (1 + xi)^2 lenDof/2
    };


(* ::Section:: *)
(*2D elements.*)


(*Bilinear quad element.*)
discreteValuesOrder["quad4"] :=
    {1, 4, 2, 3};
shapeFunctionsInternal["quad4"] :=
    Flatten[ Transpose[{shapeFunctionsInternal["line2"]} /. xi -> xi1] . ({shapeFunctionsInternal["line2"]} /. xi -> xi2) ];

(*Quadratic rectangular element.*)
discreteValuesOrder["quad9"] :=
    {1, 4, 8, 2, 3, 6, 5, 7, 9};
shapeFunctionsInternal["quad9"] :=
    Flatten[ Transpose[{shapeFunctionsInternal["line3"]} /. xi -> xi1] . ({shapeFunctionsInternal["line3"]} /. xi -> xi2) ];

(*Seredipity element.*)
shapeFunctionsInternal["quad8"] :=
    {1/4 ((1 - xi1) (1 - xi2) - (1 - xi1^2) (1 - xi2) - (1 - xi1) (1 - xi2^2)),
    1/4 ((1 + xi1) (1 - xi2) - (1 - xi1^2) (1 - xi2) - (1 + xi1) (1 - xi2^2)),
    1/4 ((1 + xi1) (1 + xi2) - (1 - xi1^2) (1 + xi2) - (1 + xi1) (1 - xi2^2)),
    1/4 ((1 - xi1) (1 + xi2) - (1 - xi1^2) (1 + xi2) - (1 - xi1) (1 - xi2^2)),
    1/2 (1 - xi1^2) (1 - xi2), 1/2 (1 + xi1) (1 - xi2^2), 1/2 (1 - xi1^2) (1 + xi2),
    1/2 (1 - xi1) (1 - xi2^2)};


(* ::Section:: *)
(*3D elements.*)


(*Hex8 element.*)
shapeFunctionsInternal["hex8"] :=
    {1/8 (1 - xi1) (1 - xi2) (1 - xi3), 1/8 (1 + xi1) (1 - xi2) (1 - xi3),
    1/8 (1 + xi1) (1 + xi2) (1 - xi3), 
    1/8 (1 - xi1) (1 + xi2) (1 - xi3), 1/8 (1 - xi1) (1 - xi2) (1 + xi3),
    1/8 (1 + xi1) (1 - xi2) (1 + xi3), 
    1/8 (1 + xi1) (1 + xi2) (1 + xi3), 1/8 (1 - xi1) (1 + xi2) (1 + xi3)};


(*Hex8 element.*)
shapeFunctionsInternal["hex20"] :=
    {1/8 (1 - xi1) (1 - xi2) (1 - xi3) (-2. - xi1 - xi2 - xi3), 
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
shapeFunctionsInternal["hex27"] :=
    {1/8 (-1 + xi1) xi1 (-1 + xi2) xi2 (-1 + xi3) xi3, 
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
shapeFunctionsInternal["tet4"] :=
    {1 - xi1 - xi2 - xi3, xi1, xi2, xi3};


(*Tet10 element.*)
shapeFunctionsInternal["tet10"] :=
    {(-1 + 2 (1 - xi1 - xi2 - xi3)) (1 - xi1 - xi2 - xi3), 
    xi1 (-1 + 2 xi1), xi2 (-1 + 2 xi2), xi3 (-1 + 2 xi3), 
    4 xi1 (1 - xi1 - xi2 - xi3), 4 xi1 xi2, 4 xi2 (1 - xi1 - xi2 - xi3), 
    4 (1 - xi1 - xi2 - xi3) xi3, 4 xi1 xi3, 4 xi2 xi3};


(*End private.*)
End[];

(*End the package.*)
EndPackage[];
