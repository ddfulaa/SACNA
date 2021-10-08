(* ::Package:: *)

(* ::Input::Initialization:: *)
(*Semi Algebraic Chemical Network Analyzer Package SACNA *)
BeginPackage["SACNA`"]
Unprotect@@Names["SACNA`*"];
ClearAll@@Names["SACNA`*"];

ClausuraDual::usage="ClausuraDual[List of reactions]"
GetRates::usage="GetRates[List of reactions]"
CheckReactions::usage="CheckReactions[Lis of reactions]"
ExtraerEspecies::usage="ExtraerEspecies[List of reactions]"
GetReactionMatrices::usage="GetReactionMatrices[List of reactions, List of species]"
OrList::usage="OrList[List]"
AndList::usage="AndList[List]"
GetSSEqns::usage="GetSSEqns[List of reactions, List of the respective reaction rates]"
GetRacemicConditions::usage="GetRacemicConditions[List of reactions]"
GetRacemicSubstitutionList::usage="GetRacemicSubstitutionList[List of reactions]"
GetLSpecies::usage="GetLSpecies[List of reactions]"
GetDSpecies::usage="GetDSpecies[List of reactions]"
GetODEs::usage="GetODEs[List of reactions, List of the respective reaction rates]"
GetPDEs::usage="GetPDEs[List of reactions, List of the respective reaction rates]"
GetJacobianMatrix::usage="GetJacobianMatrix[List of reactions, List of the respective reaction rates]"
DeletePositiveCoefficientsPolynomialOR::usage="DeletePositiveCoefficientsPolynomialOR[polynomial]"
DeletePositiveCoefficientsPolynomialAND::usage="DeletePositiveCoefficientsPolynomialAND[polynomial]"
GetRouthKurwitzTable::usage="GetRouthKurwitzTable[Matrix,StabilityChoice]"
GetRouthKurwitzTableV2::usage="GetRouthKurwitzTableV2[Matrix,StabilityChoice]"
RunSemiAlgebraicAnalysis::usage="RunSemiAlgebraicAnalysis[Reacciones, rates, tiempo]"
ReactionSystemSimulator::usage="ReactionSystemSimulatorV2[Reacciones,Rates,CADsolutionInstance: Una muestra de la CAD de la soluci\[OAcute]n,Variable temporal: t (escribir simplemente t),Tiempo inicial: ta,Tiempo final: tb]"
RunSimplifiedAnalysis::usage="RunSimplifiedAnalysis[Lista de Reacciones, Lista de Razones correspondientes a cada reacci\[OAcute]n, tiempo, n\[UAcute]mero de condici\[OAcute]n en la tabla de Routh-Huwrwitz]"
RunSemiAlgebraicAnalysisWithInitialValues::usage="RunSemiAlgebraicAnalysisWithInitialValues[reacciones, Razones, tiempo, Lista de par\[AAcute]metros a sustituir] "
RunSimplifiedAnalysisWithInitialValues::usage="RunSimplifiedAnalysisWithInitialValues[reacciones, Razones, tiempo, n\[UAcute]mero de condici\[OAcute]n en la tabla de Routh-Huwrwitz, Lista de par\[AAcute]metros a sustituir] "
ExportToChemKinLator::usage="ExportToChemKinLator[reacciones_,Rates_,CADsolutionInstance_,ATOL_, delta_, RTOL_, T0_, Tmax_] "

Begin["`Private`"]


(*ClausuraDual requiere de una lista de reacciones. 
Devuelve una lista de reacciones y sus reacciones duales correspondientes.*)
(* Returns the clousure of a list of reactions (reactions and their dual reactions)
*)
ClausuraDual[Reacciones_]:= Module[{ReaccionesDualesTE,ReaccionesTE,Dual},
ReaccionesDualesTE=ToExpression[StringReplace[Reacciones,{"L"->"D", "D"-> "L" }]];
ReaccionesTE=ToExpression[Reacciones];
Dual=StringSplit[StringReplace[ToString[DeleteDuplicates[Join[ReaccionesTE,ReaccionesDualesTE]]],{" "->"", "{"->"", "}"->""}],","];
Return[Dual]
];

(* GetRates necesita una lista de reacciones, cada una de las cuales debe tener su dual correspondiente. 
Devuelve la lista de constantes para cada reacci\[OAcute]n, teniendo en cuenta que una reacci\[OAcute]n y su dual tienen la misma constante. *)
(*
Returns the list of constants for each reaction, considering thaa a reaction and its dual have the same rate constant
*)
GetRates[Reacciones_] := Module[{DualPositions,rates},
DualPositions=DeleteDuplicates[Table[Sort[DeleteCases[{i,FirstPosition[ToExpression[Reacciones],ToExpression[StringReplace[Reacciones[[i]],{"L"->"D","D"->"L"}]]][[1]]},"NotFound"]] , {i,Length[Reacciones]}]];
rates=Table["", {i,Length[Reacciones]}];
For[i=1,i<=Length[DualPositions],i++,
rates[[DualPositions[[i,1]]]]= StringJoin["k",ToString[i]];
rates[[DualPositions[[i,2]]]]= StringJoin["k",ToString[i]];];
rates=ToExpression[rates];
Return[rates]
];

(*CheckReactions revisa si todas las reacciones cumplen con la sintaxis correcta. 
Devuelve una lista de valores booleanos. Si todos los valores son True, no hay problemas con las redes ingresadas *)
(*
Returns a list of boolean values. If all values are True, the input networks are correct.
*)
CheckReactions[r_] := Module[{}, StringMatchQ[r, RegularExpression["((\\d*)([LDZN]\\d*)\\+*)*->((\\d*)([LDZN]\\d*)\\+*)*"]]];

(*ExtraerEspecies devuelve la lista ordenada con todas las especies incluidas en las reacciones ingresadas*)
(**)

ExtraerEspecies[reaccion_] := Module[{},
Sort[Union @@ StringCases[reaccion, RegularExpression["[LDZ]\\d+"]]]
];

(*GetReactionMatrices Obtiene las matrices c y d del sistema formado por las reacciones.*)
GetReactionMatrices[Reacciones_, Especies_]:= Module[{ExtraerC,ExtraerD,TablaC,TablaD,c,d,construirFilas},
ExtraerC[reaccion_] := StringCases[StringSplit[reaccion, "->"][[1]], RegularExpression["(\\d*)([LDZN]\\d+)"] -> {"$1", "$2"}];
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
GetRacemicConditions[Reacciones_] := Module[{Especies,rc},
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

(*GetRouthKurwitzTable: Obtiene la primer columna de la tabla de Routh Kurwitz. Si el parametro Choice es True, regresa la tabla con condiciones
para revisar si la matriz es estable. En caso contrario las condiciones indicadas ser\[AAcute]n para revisar si la matriz es inestable.*)
GetRouthKurwitzTableV2[Matrix_, Choice_] := Module[{positive, polinomio,\[Lambda],ListaDeCoeficientes,n,routhHurwitzFunction,RKconditionsTable},
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
RKconditionsTable = If[Choice==False,Table[Numerator[Together[routhHurwitzFunction[i, 1, ListaDeCoeficientes]]]<0 , {i, 2, n}] ,Table[Numerator[Together[routhHurwitzFunction[i, 1, ListaDeCoeficientes]]]>= 0 , {i, 2, n}]];
Return[RKconditionsTable];
];

(*RunSemiAlgebraicAnalysis*)
RunSemiAlgebraicAnalysis[reacciones_, Rates_, tiempo_] := Module[{Reacciones,rates,Especies,SSEqns,racemicSubstitution,RacemicConditions,EspeciesL,ParametrosDeDifusion,Js,Jsracemic,rcond,MatrizA,MatrizB,MatrizDiagonal,RKStablilityConditions1,RKStablilityConditions2,RKUnstablilityConditionsR,RKUnstablilityConditionsRD,Condiciones,CondicionesSemialgebraicas,CondicionesSemialgebraicas2,CondicionesSemialgebraicasReaccion,CondicionesSemialgebraicasReaccionDifusion,ejemplo,RacemicSSEqns,solucionCAD,condition},
Reacciones=Quiet[ClausuraDual[reacciones]];
rates=If[Length[Rates]==0,GetRates[Reacciones], Rates];
Especies = ExtraerEspecies[Reacciones];
SSEqns=GetSSEqns[Reacciones,rates];
racemicSubstitution = GetRacemicSubstitutionList[Reacciones];
RacemicSSEqns= DeleteDuplicates[SSEqns /.racemicSubstitution];
RacemicConditions = GetRacemicConditions[Reacciones];
EspeciesL = GetLSpecies[Reacciones];
Js = GetJacobianMatrix[Reacciones,rates];
Jsracemic = Js /. racemicSubstitution;
rcond = Length[racemicSubstitution];
MatrizA = Jsracemic[[1 ;; rcond, 1 ;; rcond]];
MatrizB = Jsracemic[[1 ;; rcond, rcond + 1 ;; 2 rcond]];
RKUnstablilityConditionsR=GetRouthKurwitzTableV2[MatrizA-MatrizB,False];
condition=Input[StringJoin["Choose a Routh-Hurwitz contition (Enter a number from 1 to ",ToString[Length[RKUnstablilityConditionsR]],")."]];
CondicionesSemialgebraicasReaccion= RKUnstablilityConditionsR[[condition]];
CondicionesSemialgebraicasReaccionDifusion = (OrList[RKUnstablilityConditionsRD] ) && AndList[RKStablilityConditions1] && AndList[RKStablilityConditions2] &&  AndList[Thread[ParametrosDeDifusion>= 0]];
Condiciones= AndList[Thread[ToExpression[Especies] >=    0]] && AndList[Thread[SSEqns == 0]] &&AndList[Thread[Union[rates] > 0]] &&RacemicConditions;
CondicionesSemialgebraicas2 =Simplify[Condiciones &&CondicionesSemialgebraicasReaccion /. racemicSubstitution];
solucionCAD=TimeConstrained[CylindricalDecomposition[CondicionesSemialgebraicas2, Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"];
ejemplo=If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",TimeConstrained[First[FindInstance[CondicionesSemialgebraicas2,Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]], Reals]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"]];
If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print["Intento de buscar un ejemplo sin CAD"]];
If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print[ejemplo]];solucionCAD=If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",ejemplo,solucionCAD];
Return[solucionCAD];
];

ReactionSystemSimulator[reacciones_,Rates_,CADsolutionInstance_,especiesMaxInitialDifference_,t_,ta_,tb_] := Module[{Especies,Reacciones,rates,rangoTiempo,DEqns,initialConditionsExpr,EspeciesL,Sols,Soluciones,EspeciesLStrings,GraficasSolucionesQuirales,ccs,DerivadasDeEspecies,c,d,eqns,ratesInitialValues,concentrationsInitialValues,auxList,racemicSubstitution},Reacciones=Quiet[ClausuraDual[reacciones]];Especies=ExtraerEspecies[Reacciones];EspeciesL=GetLSpecies[Reacciones];rates=If[Length[Rates]==0,GetRates[Reacciones],Rates];
ccs=ToExpression[Especies];
DerivadasDeEspecies=Table[ToExpression[StringJoin["D[",Especies[[i]],"[t],t]"]],{i,Length[Especies]}];
{c,d}=GetReactionMatrices[Reacciones,Especies];eqns = Table[(Sum[rates[[j]]*(d[[j, i]] - c[[j, i]])*(Product[ccs[[l]][t]^c[[j, l]], {l, 1, Length[ccs]}]), {j, 1, Length[rates]}]) == DerivadasDeEspecies[[i]] , {i, 1, Length[ccs]}]; 
ratesInitialValues=DeleteDuplicates[rates] /. CADsolutionInstance;
racemicSubstitution = GetRacemicSubstitutionList[Reacciones];
auxList=ToExpression[Especies]/. racemicSubstitution /.CADsolutionInstance;
concentrationsInitialValues=auxList+RandomReal[{0,especiesMaxInitialDifference},Length[auxList]];
rangoTiempo={t,ta,tb};
(*Se obtiene el sistema de ecuaciones diferenciales*)
DEqns=eqns /.Thread[DeleteDuplicates[rates]->ratesInitialValues];
(*Obtiene una lista con los valores de las concentraciones iniciales*)
initialConditionsExpr=Table[ToExpression[StringJoin[Especies[[i]],"[0]"]]==concentrationsInitialValues[[i]],{i,Length[Especies]}];
(*Obtiene la primer soluci\[OAcute]n del sistema de EDO num\[EAcute]ricamente*)
Sols=NDSolve[Join[DEqns,initialConditionsExpr],ToExpression[Especies],rangoTiempo,WorkingPrecision->MachinePrecision];
Print["Species Concentrations Graphic"];
EspeciesLStrings=Table[ToString[EspeciesL[[i]]],{i,Length[EspeciesL]}];
Soluciones=Table[ToExpression[StringJoin["{",EspeciesLStrings[[i]],"[t],",StringReplace[EspeciesLStrings[[i]],"L"->"D"],"[t]}"]],{i,Length[EspeciesL]}];
GraficasSolucionesQuirales=Table[Plot[{Soluciones[[i,1]]/.Sols, Soluciones[[i,2]]/.Sols},rangoTiempo,PlotStyle->{{Red,Thick},{Blue,Thick}}, PlotLegends->{EspeciesLStrings[[i]],StringReplace[EspeciesLStrings[[i]], "L"->"D"]},Frame->True,FrameLabel->{"Time","Concentrations"},LabelStyle->Directive[Black,Bold, Medium], PlotRange->All],{i,Length[EspeciesL]}];
For[i=1,i<=Length[GraficasSolucionesQuirales],i++, Print[GraficasSolucionesQuirales[[i]]]];
Return[GraficasSolucionesQuirales];
];
(*RunSimplifiedAnalysis solo toma la primera y la \[UAcute]ltima condici\[OAcute]n de Routh Hurwitz para simplificar los c\[AAcute]lculos y poder encontrar alguna soluci\[OAcute]n al
sistema semialgebraico. Usualmente la soluci\[OAcute]n se encuentra aqu\[IAcute] y estas condiciones son las m\[AAcute]s cortas*)
RunSimplifiedAnalysis[reacciones_, Rates_, tiempo_,conditionNumber_] := Module[{Reacciones,rates,Especies,SSEqns,racemicSubstitution,RacemicConditions,EspeciesL,ParametrosDeDifusion,Js,Jsracemic,rcond,MatrizA,MatrizB,MatrizDiagonal,RKStablilityConditions1,RKStablilityConditions2,RKUnstablilityConditionsR,RKUnstablilityConditionsRD,Condiciones,CondicionesSemialgebraicas,CondicionesSemialgebraicas2,CondicionesSemialgebraicasReaccion,CondicionesSemialgebraicasReaccionDifusion,ejemplo,RacemicSSEqns, solucionCAD},
Reacciones=Quiet[ClausuraDual[reacciones]];
rates=If[Length[Rates]==0,GetRates[Reacciones], Rates];
Especies = ExtraerEspecies[Reacciones];
SSEqns=GetSSEqns[Reacciones,rates];
racemicSubstitution = GetRacemicSubstitutionList[Reacciones];
RacemicSSEqns= DeleteDuplicates[SSEqns /.racemicSubstitution];
RacemicConditions = GetRacemicConditions[Reacciones];
EspeciesL = GetLSpecies[Reacciones];
Js = GetJacobianMatrix[Reacciones,rates];
Jsracemic = Js /. racemicSubstitution;Print["Matriz Jacobiana del sistema con sustituci\[OAcute]n de las condiciones rac\[EAcute]micas"];
rcond = Length[racemicSubstitution];
MatrizA = Jsracemic[[1 ;; rcond, 1 ;; rcond]];
MatrizB = Jsracemic[[1 ;; rcond, rcond + 1 ;; 2 rcond]];
(*Routh Kurwitz Conditions for the Reaction case*)
RKUnstablilityConditionsR=GetRouthKurwitzTableV2[MatrizA-MatrizB,False];
CondicionesSemialgebraicasReaccion= If[conditionNumber==-1,Last[RKUnstablilityConditionsR],RKUnstablilityConditionsR[[conditionNumber]]];
Condiciones= AndList[Thread[ToExpression[Especies] >=    0]] && AndList[Thread[SSEqns == 0]] &&AndList[Thread[Union[rates] > 0]] &&RacemicConditions;
CondicionesSemialgebraicas2 =Simplify[Condiciones &&CondicionesSemialgebraicasReaccion /. racemicSubstitution];
solucionCAD=TimeConstrained[CylindricalDecomposition[CondicionesSemialgebraicas2, Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"];
ejemplo=If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",TimeConstrained[First[FindInstance[CondicionesSemialgebraicas2,Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]], Reals]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"]];
If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print["Intento de buscar un ejemplo sin CAD"]];
If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print[ejemplo]];solucionCAD=If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",ejemplo,solucionCAD];
Return[solucionCAD]
];

(*Esta funci\[OAcute]n toma valores iniciales de algunos par\[AAcute]metros para simplificar el sistema semialgebraico*)
RunSemiAlgebraicAnalysisWithInitialValues[reacciones_, Rates_, tiempo_, initialValuesList_] := Module[{Reacciones,rates,Especies,SSEqns,racemicSubstitution,RacemicConditions,EspeciesL,ParametrosDeDifusion,Js,Jsracemic,rcond,MatrizA,MatrizB,MatrizDiagonal,RKStablilityConditions1,RKStablilityConditions2,RKUnstablilityConditionsR,RKUnstablilityConditionsRD,Condiciones,CondicionesSemialgebraicas,CondicionesSemialgebraicas2,CondicionesSemialgebraicasReaccion,CondicionesSemialgebraicasReaccionDifusion,ejemplo,RacemicSSEqns, solucionCAD,initialValues, condition},
initialValues=Rationalize[initialValuesList];
Reacciones=Quiet[ClausuraDual[reacciones]];
rates=If[Length[Rates]==0,GetRates[Reacciones], Rates];
Especies = ExtraerEspecies[Reacciones];
SSEqns=GetSSEqns[Reacciones,rates];
racemicSubstitution = GetRacemicSubstitutionList[Reacciones];
RacemicSSEqns= DeleteDuplicates[SSEqns /.racemicSubstitution];
RacemicConditions = GetRacemicConditions[Reacciones];
EspeciesL = GetLSpecies[Reacciones];
Js = GetJacobianMatrix[Reacciones,rates];
Jsracemic = Js /. racemicSubstitution;Print["Matriz Jacobiana del sistema con sustituci\[OAcute]n de las condiciones rac\[EAcute]micas"];
rcond = Length[racemicSubstitution];
MatrizA = Jsracemic[[1 ;; rcond, 1 ;; rcond]];
MatrizB = Jsracemic[[1 ;; rcond, rcond + 1 ;; 2 rcond]];

(*Routh Kurwitz Conditions for the Reaction case*)
RKUnstablilityConditionsR=GetRouthKurwitzTable[MatrizA-MatrizB/.initialValues,False];
condition=Input[StringJoin["Choose a Routh-Hurwitz contition (Enter a number from 1 to ",ToString[Length[RKUnstablilityConditionsR]],")."]];
CondicionesSemialgebraicasReaccion= RKUnstablilityConditionsR[[condition]];
Condiciones= AndList[Thread[ToExpression[Especies] >=    0]] && AndList[Thread[SSEqns==0]] &&AndList[Thread[Union[rates] >  0]] &&RacemicConditions;
CondicionesSemialgebraicas2 =Condiciones &&CondicionesSemialgebraicasReaccion /. racemicSubstitution /. initialValues;
solucionCAD=TimeConstrained[CylindricalDecomposition[CondicionesSemialgebraicas2,Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"];
ejemplo=If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",TimeConstrained[First[FindInstance[CondicionesSemialgebraicas2,Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]], Reals]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"]];
If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print["Intento de buscar un ejemplo sin CAD"]];
If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print[ejemplo]];solucionCAD=If[solucionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",ejemplo,solucionCAD];
Return[solucionCAD]
];

(*Esta funci\[OAcute]n toma valores iniciales de algunos par\[AAcute]metros para simplificar el sistema semialgebraico*)
RunSimplifiedAnalysisWithInitialValues[reacciones_, Rates_, tiempo_, conditionNumber_, initialValuesList_] := Module[{Reacciones,rates,Especies,SSEqns,racemicSubstitution,RacemicConditions,EspeciesL,ParametrosDeDifusion,Js,Jsracemic,rcond,MatrizA,MatrizB,MatrizDiagonal,RKStablilityConditions1,RKStablilityConditions2,RKUnstablilityConditionsR,RKUnstablilityConditionsRD,Condiciones,CondicionesSemialgebraicas,CondicionesSemialgebraicas2,CondicionesSemialgebraicasReaccion,CondicionesSemialgebraicasReaccionDifusion,ejemplo,RacemicSSEqns, solutionCAD, initialValues},
initialValues=Rationalize[initialValuesList];
Reacciones=Quiet[ClausuraDual[reacciones]];
rates=If[Length[Rates]==0,GetRates[Reacciones], Rates];
Especies = ExtraerEspecies[Reacciones];
SSEqns=GetSSEqns[Reacciones,rates];
racemicSubstitution = GetRacemicSubstitutionList[Reacciones];
RacemicSSEqns= DeleteDuplicates[SSEqns /.racemicSubstitution];
RacemicConditions = GetRacemicConditions[Reacciones];
EspeciesL = GetLSpecies[Reacciones];
Js = GetJacobianMatrix[Reacciones,rates];
Jsracemic = Js /. racemicSubstitution;Print["Matriz Jacobiana del sistema con sustituci\[OAcute]n de las condiciones rac\[EAcute]micas"];
rcond = Length[racemicSubstitution];
MatrizA = Jsracemic[[1 ;; rcond, 1 ;; rcond]];
MatrizB = Jsracemic[[1 ;; rcond, rcond + 1 ;; 2 rcond]];
(*Routh Kurwitz Conditions for the Reaction case*)
RKUnstablilityConditionsR=GetRouthKurwitzTableV2[MatrizA-MatrizB /.initialValues,False];
CondicionesSemialgebraicasReaccion= If[conditionNumber==-1,Last[RKUnstablilityConditionsR],RKUnstablilityConditionsR[[conditionNumber]]];
Condiciones= AndList[Thread[ToExpression[Especies] >=    0]] && AndList[Thread[RacemicSSEqns ==0]] &&AndList[Thread[Union[rates] > 0]] ;
CondicionesSemialgebraicas2 =Condiciones &&CondicionesSemialgebraicasReaccion /. racemicSubstitution /. initialValues;
solutionCAD=TimeConstrained[CylindricalDecomposition[CondicionesSemialgebraicas2,Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"];
ejemplo=If[solutionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",TimeConstrained[First[FindInstance[CondicionesSemialgebraicas2,Sort[DeleteDuplicates@Cases[CondicionesSemialgebraicas2,_Symbol,Infinity]], Reals]],tiempo,"Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n"]];
If[solutionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print["Intento de buscar un ejemplo sin CAD"]];
If[solutionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",Print[ejemplo]];solutionCAD=If[solutionCAD==="Se alcanz\[OAcute] el tiempo impuesto sin llegar a una soluci\[OAcute]n",ejemplo,solutionCAD];
Return[solutionCAD]
];

ExportToChemKinLator[reacciones_,Rates_,CADsolutionInstance_,ATOL_, RTOL_, T0_, Tmax_, samples_] := Module[{Especies,Reacciones,rates,initialConditionsExpr,EspeciesL,ratesInitialValues,concentrationsList,racemicSubstitution,reactionratesList,jsoncontent,delta},
delta=(Tmax-T0)/samples;
Reacciones=Quiet[ClausuraDual[reacciones]];
Especies=ExtraerEspecies[Reacciones];
EspeciesL=GetLSpecies[Reacciones];rates=If[Length[Rates]==0,GetRates[Reacciones],Rates];
ratesInitialValues=DeleteDuplicates[rates] /. CADsolutionInstance;
racemicSubstitution = GetRacemicSubstitutionList[Reacciones];
concentrationsList=N[Table[Especies[[i]] ->ToExpression[Especies[[i]]], {i,Length[Especies]}] /. racemicSubstitution /. CADsolutionInstance,8];
Reacciones=Thread[StringDelete[Reacciones,"N1"]];
Reacciones=Thread[StringReplace[Reacciones,{"+"->" + ","->"->" -> ", RegularExpression["(\\d+)([DLZ]\\d+)"]->"$1 $2"}]];
Reacciones=Thread[StringReplace[Reacciones,{RegularExpression["^\\s"]->"",RegularExpression["\\s$"]->""}]];
reactionratesList=N[Table[Reacciones[[i]] -> rates[[i]],{i,Length[Reacciones]}] /. CADsolutionInstance,8];
jsoncontent={"_WARNING!_"->"This file was generated using CHEMKINLATOR. Any change to the file could make it unreadable, modify it at your own risk.","network"->{"concentrations"->concentrationsList,"rates"->reactionratesList},"simulation details"->{"ATOL"->ATOL,"Delta"->delta,"RTOL"->RTOL,"T_0"->T0,"Tmax"->Tmax},"species to plot"->{},"tabs"->{{"name"->"Time Series","type"->"timeseries"}},"version"->0.3`};
Export[SystemDialogInput["FileSave","simulacion.simu.json"],jsoncontent,"JSON"];
Return[jsoncontent];
];

End[]
Protect@@Names["SACNA`*"];
EndPackage[]

