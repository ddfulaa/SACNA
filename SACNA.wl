(* ::Package:: *)

(* ::Input::Initialization:: *)
(*Semi Algebraic Chemical Network Analyzer Package SACNA *)
BeginPackage["SACNA`"]
Unprotect@@Names["SACNA`*"];
ClearAll@@Names["SACNA`*"];

ClausuraDual::usage="ClausuraDual[Reacciones]"
GetRates::usage="GetRates[Reacciones]"
CheckReactions::usage="CheckReactions[r]"
ExtraerEspecies::usage="ExtraerEspecies[reaccion]"
GetReactionMatrices::usage="GetReactionMatrices[Reacciones, Especies]"
OrList::usage="OrList[x_List]"
AndList::usage="AndList[x_List]"
GetSSEqns::usage="GetSSEqns[Reacciones, rates]"
GetRacemicConditions::usage="GetRacemicConditions[Reacciones]"
GetRacemicSubstitutionList::usage="GetRacemicSubstitutionList[Reacciones]"
GetLSpecies::usage="GetLSpecies[Reacciones]"
GetDSpecies::usage="GetDSpecies[Reacciones]"
GetDiffusionParameters::usage="GetDiffusionParameters[Reacciones]"
GetAllDiffusionParameters::usage="GetAllDiffusionParameters[Reacciones]"
GetDiffusionExpressions::usage="GetDiffusionExpressions[Reacciones]"
GetODEs::usage="GetODEs[Reacciones, rates]"
GetPDEs::usage="GetPDEs[Reacciones, rates]"
GetJacobianMatrix::usage="GetJacobianMatrix[Reacciones, rates]"
DeletePositiveCoefficientsPolynomialOR::usage="DeletePositiveCoefficientsPolynomialOR[polynomial]"
DeletePositiveCoefficientsPolynomialAND::usage="DeletePositiveCoefficientsPolynomialAND[polynomial]"
GetRouthKurwitzTable::usage="GetRouthKurwitzTable[Matrix,StabilityChoice]"
RunSemiAlgebraicAnalysis::usage="RunSemiAlgebraicAnalysis[Reacciones, rates, tiempo]"
ReactionSystemSimulator::usage="ReactionSystemSimulator[Reacciones,Rates,ratesInitialValues,concentrationsInitialValues,t,ta,tb]"
ReactionDiffusionSystemSimulator::usage="ReactionDiffusionSystemSimulator[Lista de Reacciones, Lista de Razones correspondientes a cada reacci\[OAcute]n,Valores iniciales de las razones, Valores iniciales de las concentraciones,Valores de difusi\[OAcute]n,Variable Espacial(Escribir x),Variable Temporal(esribir t), Tiempo Inicial, Tiempo Final]"

Begin["`Private`"]
(*ClausuraDual requiere de una lista de reacciones. 
Devuelve una lista de reacciones y sus reacciones duales correspondientes.*)

ClausuraDual[Reacciones_]:= Module[{ReaccionesDualesTE,ReaccionesTE,Dual},
ReaccionesDualesTE=ToExpression[StringReplace[Reacciones,{"L"->"D", "D"-> "L" }]];
ReaccionesTE=ToExpression[Reacciones];
Dual=StringSplit[StringReplace[ToString[DeleteDuplicates[Join[ReaccionesTE,ReaccionesDualesTE]]],{" "->"", "{"->"", "}"->""}],","];
Return[Dual]
];

(* GetRates necesita una lista de reacciones, cada una de las cuales debe tener su dual correspondiente. 
Devuelve la lista de constantes para cada reacci\[OAcute]n, teniendo en cuenta que una reacci\[OAcute]n y su dual tienen la misma constante. *)
GetRates[Reacciones_] := Module[{DualPositions,rates},DualPositions=DeleteDuplicates[Table[Sort[DeleteCases[{i,FirstPosition[ToExpression[Reacciones],ToExpression[StringReplace[Reacciones[[i]],{"L"->"D","D"->"L"}]]][[1]]},"NotFound"]] , {i,Length[Reacciones]}]];
rates=Table["", {i,Length[Reacciones]}];
For[i=1,i<=Length[DualPositions],i++,
rates[[DualPositions[[i,1]]]]= StringJoin["k",ToString[i]];
rates[[DualPositions[[i,2]]]]= StringJoin["k",ToString[i]];];
rates=ToExpression[rates];
Return[rates]
];

(*CheckReactions revisa si todas las reacciones cumplen con la sintaxis correcta. 
Devuelve una lista de valores booleanos. Si todos los valores son True, no hay problemas con las redes ingresadas *)
CheckReactions[r_] := Module[{}, StringMatchQ[r, RegularExpression["((\\d*)([LDZ]\\d*)\\+*)*->((\\d*)([LDZN]\\d*)\\+*)*"]]];

(*ExtraerEspecies devuelve la lista ordenada con todas las especies incluidas en las reacciones ingresadas*)
ExtraerEspecies[reaccion_] := Module[{},Sort[Union @@ StringCases[reaccion, RegularExpression["[LDZ]\\d+"]]]];(*GetReactionMatrices Obtiene las matrices c y d del sistema formado por las reacciones.*)
GetReactionMatrices[Reacciones_, Especies_]:= Module[{ExtraerC,ExtraerD,TablaC,TablaD,c,d},
ExtraerC[reaccion_] := StringCases[StringSplit[reaccion, "->"][[1]], RegularExpression["(\\d*)([LDZ]\\d+)"] -> {"$1", "$2"}];
ExtraerD[reaccion_] := StringCases[StringSplit[reaccion, "->"][[2]], RegularExpression["(\\d*)([LDZN]\\d+)"] -> {"$1", "$2"}];
TablaC = Table[ExtraerC@Reacciones[[i]], {i, Length[Reacciones]}];
TablaD = Table[ExtraerD@Reacciones[[i]], {i, Length[Reacciones]}];
construirFilas[especie_, fila_] := Table[Total[Table[If[especie[[i]] == fila[[j, 2]], If[fila[[j, 1]] == "", 1, ToExpression[fila[[j, 1]]]], 0], {i, Length[especie]}, {j, Length[fila]}][[k]]], {k, Length[especie]}];
c = Table[construirFilas[Especies, TablaC[[l]]], {l, 1, Length[Reacciones]}];
d = Table[construirFilas[Especies, TablaD[[l]]], {l, 1, Length[Reacciones]}] ;
Return[{c,d}]
];
(*Implementaci\[OAcute]n de Or para elementos de una lista*)
OrList [x_List] := Module[{}, Or @@ x];

(*Implementaci\[OAcute]n de Or para elementos de una lista*)
AndList [x_List] := Module[{}, And @@ x];

(*GetSSEqns: Obtiene una lista con las ecuaciones del estado estacionario.*)
GetSSEqns[Reacciones_,rates_] := Module[{ccs,eqns,c,d,Especies},
(*rates=Quiet[GetRates[Reacciones]];*)
Especies = ExtraerEspecies[Reacciones];
ccs=ToExpression[Especies];
{c,d}=GetReactionMatrices[Reacciones,Especies];
eqns = Table[(Sum[rates[[j]]*(d[[j, i]] - c[[j, i]])*(Product[ccs[[l]]^c[[j, l]], {l, 1, Length[ccs]}]), {j, 1, Length[rates]}]), {i, 1, Length[ccs]}]; 
Return[eqns]
];

(*GetRacemicConditions: Obtiene las condiciones rac\[EAcute]micas *)
GetRacemicConditions[Reacciones_] := Module[{Especies},
Especies = ExtraerEspecies[Reacciones];
rc=AndList[ToExpression[Union @@ Union @@ Table[StringCases[Especies[[i]], RegularExpression["L(\\d+)"] -> {"L$1 == D$1"}], {i, Length[Especies]}]]];
Return[rc]
];

(*GetRacemicSubstitutionList: Obtiene una lista con sustituciones de las condiciones rac\[EAcute]micas, para simplificar t\[EAcute]rminos*)
GetRacemicSubstitutionList[Reacciones_] := Module[{Especies,racemicSubstitution},
Especies = ExtraerEspecies[Reacciones];
racemicSubstitution = ToExpression[Union @@ Union @@ Table[StringCases[Especies[[i]], RegularExpression["L(\\d+)"] -> {"L$1 -> D$1"}], {i, Length[Especies]}]];
Return[racemicSubstitution]
];

(*GetLSpecies: Obtiene una lista con todas las especies de tipo L*)
GetLSpecies[Reacciones_] := Module[{Especies,LSpecies},
Especies = ExtraerEspecies[Reacciones];
LSpecies = ToExpression[Union @@ Union @@ Table[StringCases[Especies[[i]], RegularExpression["L(\\d+)"] -> {"L$1"}], {i, Length[Especies]}]];
Return[LSpecies]
];

(*GetDSpecies: Obtiene una lista con todas las especies de tipo D*)
GetDSpecies[Reacciones_] := Module[{Especies,DSpecies},
Especies = ExtraerEspecies[Reacciones];
DSpecies = ToExpression[Union @@ Union @@ Table[StringCases[Especies[[i]], RegularExpression["D(\\d+)"] -> {"D$1"}], {i, Length[Especies]}]];
Return[DSpecies]
];

(*GetDiffusionParameters: Obtiene una lista con los par\[AAcute]metros para el modelo de reacci\[OAcute]n-difusi\[OAcute]n*)
GetDiffusionParameters[Reacciones_] := Module[{Especies,DiffusionParameters},
Especies = ExtraerEspecies[Reacciones];
DiffusionParameters = ToExpression[Union @@ Union @@ Table[StringCases[Especies[[i]], RegularExpression["L(\\d+)"] -> {"\[Delta]$1"}], {i, Length[Especies]}]];
Return[DiffusionParameters]
];
(*GetAllDiffusionParameters: Obtiene una lista con los par\[AAcute]metros para el modelo de reacci\[OAcute]n-difusi\[OAcute]n, pero les agrega el nombre de las especies, esto para el simulador gr\[AAcute]fico*)
GetAllDiffusionParameters[Reacciones_] := Module[{Especies,DiffusionParameters},
Especies = ExtraerEspecies[Reacciones];
DiffusionParameters = ToExpression[Union @@ Union @@ Table[StringCases[Especies[[i]], RegularExpression["([LDZ]\\d+)"] -> {"\[Delta]$1"}], {i, Length[Especies]}]];
Return[DiffusionParameters]
];
(*GetDiffusionExpressions: Obtiene una lista con las derivadas espaciales para el caso reacci\[OAcute]ni difusi\[OAcute]n*)
GetDiffusionExpressions[Reacciones_] := Module[{Especies,DiffusionParameters},
Especies = ExtraerEspecies[Reacciones];
DiffusionParameters = ToExpression[Union @@ Union @@ Table[StringCases[Especies[[i]], RegularExpression["([LDZ]\\d+)"] -> {"\[Delta]$1 D[D[$1[t,x],x],x]"}], {i, Length[Especies]}]];
Return[DiffusionParameters]
];
(*GetODEqns: Obtiene una lista con las ecuaciones diferenciales del sistema en el caso reactivo.
sndDerivatives=GetDiffusionExpressions[Reacciones];
difParams=GetAllDiffusionParameters[Reacciones];
*)
GetODEs[Reacciones_,rates_,t_] := Module[{ccs,eqns,c,d,Especies,DerivadasDeEspecies},
(*rates=Quiet[GetRates[Reacciones]];*)
Especies = ExtraerEspecies[Reacciones];
ccs=ToExpression[Especies];
DerivadasDeEspecies=Table[ToExpression[StringJoin["D[",Especies[[i]],"[t],t]"]],{i,Length[Especies]}];
{c,d}=GetReactionMatrices[Reacciones,Especies];
eqns = Table[(Sum[rates[[j]]*(d[[j, i]] - c[[j, i]])*(Product[ccs[[l]][t]^c[[j, l]], {l, 1, Length[ccs]}]), {j, 1, Length[rates]}]) == DerivadasDeEspecies[[i]] , {i, 1, Length[ccs]}]; 
Return[eqns]
];
(*GetPDEqns: Obtiene una lista con las ecuaciones diferenciales parciales del sistema en el caso reacci\[OAcute]n difusi\[OAcute]n.
*)
GetPDEs[Reacciones_,rates_,t_,x_] := Module[{ccs,eqns,c,d,Especies,DerivadasDeEspecies,sndDerivatives,difParams},
(*rates=Quiet[GetRates[Reacciones]];*)
Especies = ExtraerEspecies[Reacciones];
ccs=ToExpression[Especies];
DerivadasDeEspecies=Table[ToExpression[StringJoin["D[",Especies[[i]],"[t,x],t]"]],{i,Length[Especies]}];
sndDerivatives=GetDiffusionExpressions[Reacciones];
difParams=GetAllDiffusionParameters[Reacciones];
{c,d}=GetReactionMatrices[Reacciones,Especies];
eqns = Table[(Sum[rates[[j]]*(d[[j, i]] - c[[j, i]])*(Product[ccs[[l]][t,x]^c[[j, l]], {l, 1, Length[ccs]}]), {j, 1, Length[rates]}]) +difParams[[i]]*sndDerivatives[[i]]== DerivadasDeEspecies[[i]] , {i, 1, Length[ccs]}]; 
Return[eqns]
];

(*GetJacobianMatrix: Obtiene la matriz jacobiana del sistema*)
GetJacobianMatrix[Reacciones_,rates_] := Module[{}, Grad[GetSSEqns[Reacciones,rates],ToExpression[ExtraerEspecies[Reacciones]]]];
(*Convierte un polinomio con coeficientes positivos en una expresion False, esto porque se quiere obtener la expresion poli< 0, lo cual es Falso si todos sus coeficientes son positivos*)
DeletePositiveCoefficientsPolynomialOR[poli_] := Module[{},If[Or@@Thread[ToExpression[Join@@StringCases[ToString[CoefficientRules[poli]],RegularExpression["->\\s*(-*\\d+)\\s*,*"]->{"$1"}]]<0],poli <0,False]];
(*Convierte un polinomio con coeficientes positivos en una expresion True, esto porque se quiere obtener la expresion poli\[GreaterEqual] 0, lo cual es verdad si todos sus coeficientes son positivos*)
DeletePositiveCoefficientsPolynomialAND[poli_] := Module[{},If[And@@Thread[ToExpression[Join@@StringCases[ToString[CoefficientRules[poli]],RegularExpression["->\\s*(-*\\d+)\\s*,*"]->{"$1"}]]>= 0],True,poli >= 0]];
(*GetRouthKurwitzTable: Obtiene la primer columna de la tabla de Routh Kurwitz. Si el parametro Choice es True, regresa la tabla con condiciones
para revisar si la matriz es estable. En caso contrario las condiciones indicadas ser\[AAcute]n para revisar si la matriz es inestable.*)
GetRouthKurwitzTable[Matrix_, Choice_] := Module[{positive, polinomio,\[Lambda],ListaDeCoeficientes,n,routhHurwitzFunction,RKconditionsTable},
(*positive: Convierte un polinomio caracter\[IAcute]stico a uno m\[OAcute]nico. Es decir, ajusta para que el coeficiente del t\[EAcute]rmino de mayor grado sea positivo.*)
positive[poly_, var_] := poly * (-1)^(Exponent[poly,var]);
polinomio = CharacteristicPolynomial[Matrix, \[Lambda]];
polinomio = positive[polinomio, \[Lambda]];
ListaDeCoeficientes = Reverse[CoefficientList[polinomio, \[Lambda]]];
routhHurwitzFunction[1, i_, lista_] := routhHurwitzFunction[1, i, lista] = TakeDrop[lista, {1, -1, 2}][[1]][[i]];
routhHurwitzFunction[2, i_, lista_] := routhHurwitzFunction[2, i, lista] = TakeDrop[lista, {1, -1, 2}][[2]][[i]];
routhHurwitzFunction[i_, j_, lista_] := If[ j < Length[lista]/2, routhHurwitzFunction[i, j, lista] = routhHurwitzFunction[i - 2, j + 1, lista] - routhHurwitzFunction[i - 2, 1, lista]*routhHurwitzFunction[i - 1, j + 1, lista]/routhHurwitzFunction[i - 1, 1, lista], 0] ;
n = Length[ListaDeCoeficientes];
If[Mod[Length[ListaDeCoeficientes], 2] > 0, ListaDeCoeficientes = PadRight[ListaDeCoeficientes, Length[ListaDeCoeficientes] + 1]];
 ListaDeCoeficientes = PadRight[ListaDeCoeficientes, Length[ListaDeCoeficientes] + 2];
RKconditionsTable = If[Choice==False,Table[DeletePositiveCoefficientsPolynomialOR[Numerator[Together[routhHurwitzFunction[i, 1, ListaDeCoeficientes]]]] , {i, 2, n}] ,Table[DeletePositiveCoefficientsPolynomialAND[Numerator[Together[routhHurwitzFunction[i, 1, ListaDeCoeficientes]]]] , {i, 2, n}]];
Return[RKconditionsTable];
];

RunSemiAlgebraicAnalysis[reacciones_, Rates_, tiempo_] := Module[{Reacciones,rates,Especies,SSEqns,racemicSubstitution,RacemicConditions,EspeciesL,ParametrosDeDifusion,Js,Jsracemic,rcond,MatrizA,MatrizB,MatrizDiagonal,RKStablilityConditions1,RKStablilityConditions2,RKUnstablilityConditionsR,RKUnstablilityConditionsRD,Condiciones,CondicionesSemialgebraicas,CondicionesSemialgebraicas2,CondicionesSemialgebraicasReaccion,CondicionesSemialgebraicasReaccionDifusion,ejemplo},
Print["Reacciones ingresadas"];
Print[reacciones];
Print["Clausura Quiral de las reacciones ingresadas"];
Reacciones=Quiet[ClausuraDual[reacciones]];
Print[Reacciones];
Print["Revisi\[OAcute]n de la sintaxis de las reacciones ingresadas"];
Print[CheckReactions[Reacciones]];
Print["Constantes de las reacciones"];
Print[If[Length[Rates]==0,GetRates[Reacciones], Rates]];
rates=If[Length[Rates]==0,GetRates[Reacciones], Rates];

Print[rates];
Print["Especies"];
Especies = ExtraerEspecies[Reacciones];
Print[Especies];
Print["Ecuaciones del sistema"];
SSEqns=GetSSEqns[Reacciones,rates];
Print[SSEqns==0];
Print["Sustituciones rac\[EAcute]micas"];
racemicSubstitution = GetRacemicSubstitutionList[Reacciones];
Print[racemicSubstitution];
Print[SSEqns == 0 /. racemicSubstitution];
Print["Condiciones rac\[EAcute]micas"];
RacemicConditions = GetRacemicConditions[Reacciones];
Print[RacemicConditions];
EspeciesL = GetLSpecies[Reacciones];
ParametrosDeDifusion = GetDiffusionParameters[Reacciones];
Print["Matriz Jacobiana del sistema"];
Js = GetJacobianMatrix[Reacciones,rates];
Print[Js //MatrixForm];
Jsracemic = Js /. racemicSubstitution;Print["Matriz Jacobiana del sistema con sustituci\[OAcute]n de las condiciones rac\[EAcute]micas"];
Print[Jsracemic //MatrixForm];
rcond = Length[racemicSubstitution];
MatrizA = Jsracemic[[1 ;; rcond, 1 ;; rcond]];
Print["Matriz A"];
Print[MatrizA //MatrixForm];
MatrizB = Jsracemic[[1 ;; rcond, rcond + 1 ;; 2 rcond]];
Print["Matriz B"];
Print[MatrizB //MatrixForm];
Print["Matriz A+B"];
Print[MatrizA+MatrizB //MatrixForm];
Print["Matriz A-B"];
Print[MatrizA-MatrizB //MatrixForm];
MatrizDiagonal = DiagonalMatrix[ParametrosDeDifusion];
Print["Matriz A-B-\[CapitalDelta]"];
Print[MatrizA-MatrizB-MatrizDiagonal //MatrixForm];(*Routh Kurwitz Conditions for stability in the Reaction-diffussion case*)
RKStablilityConditions1=GetRouthKurwitzTable[MatrizA+MatrizB,True];
Print["Condiciones para la estabilidad en el caso de reacci\[OAcute]n difusi\[OAcute]n para la matriz A+B"];
Print[RKStablilityConditions1];
(*Routh Kurwitz Conditions for stability in the Reaction-diffussion case*)
RKStablilityConditions2=GetRouthKurwitzTable[MatrizA-MatrizB,True];
Print["Condiciones para la estabilidad en el caso de reacci\[OAcute]n difusi\[OAcute]n para la matriz A-B"];
Print[RKStablilityConditions2];
(*Routh Kurwitz Conditions for unstability in the Reaction-diffussion case*)
RKUnstablilityConditionsRD=GetRouthKurwitzTable[MatrizA-MatrizB-MatrizDiagonal,False];
Print["Condiciones para la inestabilidad en el caso de reacci\[OAcute]n difusi\[OAcute]n para la matriz A-B-\[CapitalDelta]"];
Print[RKUnstablilityConditionsRD];
(*Routh Kurwitz Conditions for the Reaction case*)
Print["Condiciones para la inestabilidad en el caso de reacci\[OAcute]n difusi\[OAcute]n para la matriz A-B"];
RKUnstablilityConditionsR=GetRouthKurwitzTable[MatrizA-MatrizB,False];
Print[RKUnstablilityConditionsR];
CondicionesSemialgebraicasReaccion= (OrList[RKUnstablilityConditionsR] );

CondicionesSemialgebraicasReaccionDifusion = (OrList[RKUnstablilityConditionsRD] ) && AndList[RKStablilityConditions1] && AndList[RKStablilityConditions2] &&  AndList[Thread[ParametrosDeDifusion>= 0]];
Print["Condiciones Semialgebraicas del sistema reacci\[OAcute]n-difusi\[OAcute]n"];
(*Condiciones= ToExpression[Especies] >=   0 && ToExpression[Especies] != 0  && SSEqns == 0 && Union[rates] > 0 &&RacemicConditions;*)
Condiciones= AndList[Thread[ToExpression[Especies] >=   0]] && AndList[Thread[SSEqns == 0]] &&AndList[Thread[Union[rates] > 0]] &&RacemicConditions;CondicionesSemialgebraicas =Simplify[Condiciones &&CondicionesSemialgebraicasReaccionDifusion /. racemicSubstitution];
Print[CondicionesSemialgebraicas];
Print["Soluci\[OAcute]n del sistema reacci\[OAcute]n-difusi\[OAcute]n"];
Print[TimeConstrained[GenericCylindricalDecomposition[CondicionesSemialgebraicas, Union[ParametrosDeDifusion,ToExpression[Especies], rates]] ,tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"]];
Print["Por ejemplo"];
Print[TimeConstrained[FindInstance[CondicionesSemialgebraicas, Union[ParametrosDeDifusion,ToExpression[Especies], rates], Reals] ,tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"]];
Print["Condiciones semialgebraicas del sistema reacci\[OAcute]n"];
CondicionesSemialgebraicas2 =Simplify[Condiciones &&CondicionesSemialgebraicasReaccion /. racemicSubstitution];
Print[CondicionesSemialgebraicas2];
Print["Soluci\[OAcute]n del sistema reacci\[OAcute]n"];
(*Reduce[CondicionesSemialgebraicas2, Union[ToExpression[Especies], rates], Reals]*)
(*Print[Reduce[CondicionesSemialgebraicas2, Union[ToExpression[Especies], rates]]]*)

(*Reduce[CondicionesSemialgebraicas2, Union[ToExpression[Especies], rates]]*)
Print[TimeConstrained[Reduce[CondicionesSemialgebraicas2, Union[ToExpression[Especies]/.racemicSubstitution, rates]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"]];
ejemplo=TimeConstrained[First[FindInstance[CondicionesSemialgebraicas2, Union[ToExpression[Especies]/.racemicSubstitution, rates], Reals]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"];
Print["Por ejemplo"];
Print[ejemplo];
Return[True]
];

ReactionSystemSimulator[reacciones_,Rates_,ratesInitialValues_,concentrationsInitialValues_,t_,ta_,tb_] := Module[{Especies,Reacciones,rates,err,rangoTiempo,DEqns,initialConditionsExpr,EspeciesL,Sols,Soluciones,EspeciesLStrings,GraficasSolucionesQuirales,ccs,DerivadasDeEspecies,c,d,eqns},Reacciones=Quiet[ClausuraDual[reacciones]];Especies=ExtraerEspecies[Reacciones];EspeciesL=GetLSpecies[Reacciones];rates=If[Length[Rates]==0,GetRates[Reacciones],Rates];
ccs=ToExpression[Especies];
DerivadasDeEspecies=Table[ToExpression[StringJoin["D[",Especies[[i]],"[t],t]"]],{i,Length[Especies]}];
{c,d}=GetReactionMatrices[Reacciones,Especies];eqns = Table[(Sum[rates[[j]]*(d[[j, i]] - c[[j, i]])*(Product[ccs[[l]][t]^c[[j, l]], {l, 1, Length[ccs]}]), {j, 1, Length[rates]}]) == DerivadasDeEspecies[[i]] , {i, 1, Length[ccs]}]; 

Print[DeleteDuplicates[rates]];Print[Especies];
Print[Thread[DeleteDuplicates[rates]->ratesInitialValues]];
rangoTiempo={t,ta,tb};
(*Se obtiene el sistema de ecuaciones diferenciales*)DEqns=eqns /.Thread[DeleteDuplicates[rates]->ratesInitialValues];Print["Ecuaciones diferenciales"];Print[DEqns];(*Obtiene una lista con los valores de las concentraciones iniciales*)initialConditionsExpr=Table[ToExpression[StringJoin[Especies[[i]],"[0]"]]==concentrationsInitialValues[[i]],{i,Length[Especies]}];Print["Condiciones Iniciales"];Print[initialConditionsExpr];
(*Obtiene la primer soluci\[OAcute]n del sistema de EDO num\[EAcute]ricamente*)Sols=NDSolve[Join[DEqns,initialConditionsExpr],ToExpression[Especies],rangoTiempo];Print["Species Concentrations Graphic"];EspeciesLStrings=Table[ToString[EspeciesL[[i]]],{i,Length[EspeciesL]}];Soluciones=Table[ToExpression[StringJoin["{",EspeciesLStrings[[i]],"[t],",StringReplace[EspeciesLStrings[[i]],"L"->"D"],"[t]}"]],{i,Length[EspeciesL]}];GraficasSolucionesQuirales=Table[Plot[{Soluciones[[i,1]]/.Sols, Soluciones[[i,2]]/.Sols},rangoTiempo,PlotStyle->{{Red,Thick},{Blue,Thick}}, PlotLegends->{EspeciesLStrings[[i]],StringReplace[EspeciesLStrings[[i]], "L"->"D"]},Frame->True,FrameLabel->{"Time","Concentrations"},LabelStyle->Directive[Black,Bold, Medium], PlotRange->All],{i,Length[EspeciesL]}];For[i=1,i<=Length[GraficasSolucionesQuirales],i++, Print[GraficasSolucionesQuirales[[i]]]];Return[GraficasSolucionesQuirales];
];

ReactionDiffusionSystemSimulator[reacciones_,Rates_,ratesInitialValues_,concentrationsInitialValues_,diffusionValues_,x_,t_,tinicial_, tfinal_]:=Module[{Especies,EspeciesL,Reacciones,rates,rangoTiempo,rangoEspacio,diffusionConstants,DEqns,initialTemporalConditions,Sols,EspeciesLStrings,Soluciones,Graficas3D},
Reacciones=Quiet[ClausuraDual[reacciones]];Especies=ExtraerEspecies[Reacciones];EspeciesL=GetLSpecies[Reacciones];rates=If[Length[Rates]==0,GetRates[Reacciones],Rates];
Print[DeleteDuplicates[rates]];
Print[Thread[DeleteDuplicates[rates]->ratesInitialValues]];
Print[Especies];
rangoTiempo={t,tinicial,tfinal};
rangoEspacio={x,0,2*Pi};
diffusionConstants=GetAllDiffusionParameters[Reacciones];
Print[diffusionConstants];
Print[Thread[diffusionConstants->diffusionValues]];
(*Se obtiene el sistema de ecuaciones diferenciales*)
DEqns=GetPDEs[Reacciones,rates,t,x] /.Thread[DeleteDuplicates[rates]->ratesInitialValues]/. Thread[diffusionConstants->diffusionValues];
Print["Ecuaciones diferenciales parciales"];
Print[DEqns];
(*Obtiene una lista con los valores de las concentraciones iniciales*)
initialTemporalConditions=Table[ToExpression[StringJoin[Especies[[i]],"[0,x]"]]==concentrationsInitialValues[[i]],{i,Length[Especies]}];
Print["Condiciones iniciales temporales"];
Print[initialTemporalConditions];

(*Cuando no se especifican las dem\[AAcute]s condiciones, mathematica supone que son las condiciones de Neumann \[Equal] 0, buscar Neumann Zero en la documentaci\[OAcute]n https://reference.wolfram.com/language/FEMDocumentation/tutorial/SolvingPDEwithFEM.html*)
Sols=NDSolve[Join[DEqns,initialTemporalConditions],ToExpression[Especies],rangoTiempo,rangoEspacio,Method->"FiniteElement"];
EspeciesLStrings=Table[ToString[EspeciesL[[i]]],{i,Length[EspeciesL]}];
Soluciones=Table[ToExpression[StringJoin["{",EspeciesLStrings[[i]],"[t,x],",StringReplace[EspeciesLStrings[[i]],"L"->"D"],"[t,x]}"]],{i,Length[EspeciesL]}];
Graficas3D=Table[Plot3D[{Soluciones[[i,1]]/.Sols, Soluciones[[i,2]]/.Sols},rangoTiempo,rangoEspacio,AxesLabel->Automatic,PlotRange->All,PlotLegends->{EspeciesLStrings[[i]],StringReplace[EspeciesLStrings[[i]], "L"->"D"]}],{i,Length[EspeciesL]}];
For[i=1,i<=Length[Graficas3D],i++, Print[Graficas3D[[i]]]];
Return[Graficas3D];
];

End[]
Protect@@Names["SACNA`*"];
EndPackage[]

