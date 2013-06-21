(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Mathematica 8.0' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       157,          7]
NotebookDataLength[    100013,       2042]
NotebookOptionsPosition[     98570,       1994]
NotebookOutlinePosition[     98906,       2009]
CellTagsIndexPosition[     98863,       2006]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"BeginPackage", "[", "\"\<MultifieldFunctions`\>\"", "]"}], 
  "*)"}]], "Input",
 CellChangeTimes->{{3.580777849034735*^9, 3.580777874451135*^9}, 
   3.580778363676653*^9, {3.580778488652233*^9, 3.580778495630071*^9}}],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{
   RowBox[{"H", "::", "usage"}], "=", 
   RowBox[{
    RowBox[{"\"\<Hubble factor\>\"", "\[IndentingNewLine]", 
     RowBox[{"weqnstate", "::", "usage"}]}], "=", 
    RowBox[{
     RowBox[{
     "\"\<Equation of state parameter.  w=p/\[Rho]\>\"", 
      "\[IndentingNewLine]", 
      RowBox[{"\[Epsilon]", "::", "usage"}]}], "=", 
     RowBox[{
      RowBox[{
      "\"\<First slow roll parameter: \
\[Epsilon]=-(dH/dt)/\!\(\*SuperscriptBox[\(H\), \(2\)]\)=(1/2)d\[Phi].d\[Phi]\
\>\"", "\[IndentingNewLine]", 
       RowBox[{"a", "::", "usage"}]}], "=", 
      RowBox[{
       RowBox[{"\"\<Scale factor.\>\"", "\[IndentingNewLine]", 
        RowBox[{"kgeqn", "::", "usage"}]}], "=", 
       "\"\<Background equation of motion.\>\""}]}]}]}]}], "*)"}]], "Input",
 CellChangeTimes->{{3.580778164024135*^9, 3.580778297511432*^9}, {
  3.580778491640447*^9, 3.580778500001092*^9}}],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"Begin", "[", "\"\<`Private`\>\"", "]"}], "*)"}]], "Input",
 CellChangeTimes->{{3.580778307029661*^9, 3.580778317810142*^9}, {
  3.580778503107834*^9, 3.580778504006513*^9}}],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{"Library", " ", "Call"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
    "Needs", "[", 
     "\"\<DifferentialEquations`InterpolatingFunctionAnatomy`\>\"", "]"}], 
    ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
    "Needs", "[", "\"\<DifferentialEquations`NDSolveProblems`\>\"", "]"}], 
    ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
    "Needs", "[", "\"\<DifferentialEquations`NDSolveUtilities`\>\"", "]"}], 
    ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Needs", "[", "\"\<FunctionApproximations`\>\"", "]"}], 
    ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{
    RowBox[{"Potential", " ", "&"}], " ", "deriv"}], "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"m2", "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"10", "^", 
        RowBox[{"(", 
         RowBox[{"-", "10.422895047377294"}], ")"}]}], ",", 
       RowBox[{"10", "^", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "10.422895047377294"}], "+", 
          RowBox[{"Log10", "[", "81.0", "]"}]}], ")"}]}]}], "}"}]}], ";"}], 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"V", "[", "\[Phi]_", "]"}], ":=", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1.", "/", "2"}], ")"}], 
     RowBox[{"m2", ".", 
      RowBox[{"(", 
       RowBox[{"\[Phi]", "*", "\[Phi]"}], ")"}]}]}]}], "   ", 
   RowBox[{"(*", "Scalar", "*)"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"dVd\[Phi]", "[", 
     RowBox[{"\[Phi]_", "?", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"VectorQ", "[", 
         RowBox[{"#", ",", " ", "NumericQ"}], "]"}], " ", "&"}], ")"}]}], 
     "]"}], ":=", 
    RowBox[{"m2", "*", "\[Phi]"}]}], "  ", 
   RowBox[{"(*", "Vector", "*)"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"m2matrix", "[", 
     RowBox[{"\[Phi]_", "?", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"VectorQ", "[", 
         RowBox[{"#", ",", " ", "NumericQ"}], "]"}], " ", "&"}], ")"}]}], 
     "]"}], ":=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"m2", "[", 
         RowBox[{"[", "1", "]"}], "]"}], ",", "0"}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"0", ",", 
        RowBox[{"m2", "[", 
         RowBox[{"[", "2", "]"}], "]"}]}], "}"}]}], "}"}], " ", 
    RowBox[{"(*", 
     RowBox[{
      RowBox[{"Hessian", " ", "matrix", " ", "for", " ", "V"}], ",", " ", 
      RowBox[{
       RowBox[{"ddV", "/", "d\[Phi]_id\[Phi]"}], "_j"}]}], 
     "*)"}]}]}]}]], "Input",
 CellChangeTimes->{{3.58076651484524*^9, 3.580766516841665*^9}, {
   3.580766684887984*^9, 3.58076670950583*^9}, {3.580771325628658*^9, 
   3.58077137533566*^9}, 3.580777734726744*^9, 3.580778534084891*^9}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"\[Phi]0", "=", 
   RowBox[{"{", 
    RowBox[{"10.31001", ",", "12.93651"}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"d\[Phi]0", "=", 
   RowBox[{
    RowBox[{"-", 
     RowBox[{"dVd\[Phi]", "[", "\[Phi]0", "]"}]}], "/", 
    RowBox[{"(", 
     RowBox[{"3", " ", 
      RowBox[{"H2", "[", 
       RowBox[{"\[Phi]0", ",", 
        RowBox[{"{", 
         RowBox[{"0.", ",", "0."}], "}"}]}], "]"}]}], ")"}]}]}], ";", 
  RowBox[{"(*", 
   RowBox[{"slow", " ", "roll", " ", 
    RowBox[{"w", "/", "negligble"}], " ", "KE", " ", "contrib"}], "*)"}], 
  "\[IndentingNewLine]", "\[IndentingNewLine]", 
  RowBox[{"(*", 
   RowBox[{"Set", " ", "other", " ", "ICs"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{"a0", "=", "1."}], ";"}], "\[IndentingNewLine]"}], "Input",
 CellChangeTimes->{{3.580682758555876*^9, 3.580682794439814*^9}, {
   3.580682951205136*^9, 3.58068295706301*^9}, {3.580683230166759*^9, 
   3.580683235082072*^9}, 3.580683732992644*^9, {3.580684287353712*^9, 
   3.580684288003732*^9}, {3.580684378266436*^9, 3.580684412512316*^9}, {
   3.580684858044544*^9, 3.58068486185786*^9}, {3.580759776123397*^9, 
   3.580759782894651*^9}, 3.58076011008255*^9, {3.5807603594566*^9, 
   3.580760360177565*^9}, {3.580762317702638*^9, 3.580762327696071*^9}, {
   3.580762549502399*^9, 3.58076256834792*^9}, {3.580762631476859*^9, 
   3.580762636428616*^9}, {3.580762747121856*^9, 3.58076275941746*^9}, {
   3.580764562793563*^9, 3.580764636488595*^9}, {3.580771765088671*^9, 
   3.580771769073398*^9}, {3.580776595901836*^9, 3.580776599757095*^9}, 
   3.580777591175517*^9, 3.58077854044252*^9}],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{
    RowBox[{"KE", " ", 
     RowBox[{"densities", " ", "--"}]}], "-", " ", "background"}], "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"(*", 
   RowBox[{"All", " ", "are", " ", "scalars"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"KEoverH2", "[", "d\[Phi]_", "]"}], ":=", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1.", "/", "2"}], ")"}], 
     RowBox[{"(", 
      RowBox[{"d\[Phi]", ".", "d\[Phi]"}], ")"}]}]}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"H", "[", 
     RowBox[{"\[Phi]_", ",", "d\[Phi]_"}], "]"}], ":=", 
    RowBox[{"Sqrt", "[", 
     RowBox[{
      RowBox[{"V", "[", "\[Phi]", "]"}], "/", 
      RowBox[{"(", 
       RowBox[{"3.", "-", 
        RowBox[{"KEoverH2", "[", "d\[Phi]", "]"}]}], ")"}]}], "]"}]}], 
   RowBox[{"(*", 
    RowBox[{"G", "=", 
     RowBox[{
      RowBox[{
       SuperscriptBox[
        SubscriptBox["M", "P"], "2"], "/", "8"}], "\[Pi]"}]}], "*)"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"H2", "[", 
     RowBox[{"\[Phi]_", ",", "d\[Phi]_"}], "]"}], ":=", 
    RowBox[{
     RowBox[{"V", "[", "\[Phi]", "]"}], "/", 
     RowBox[{"(", 
      RowBox[{"3.", "-", 
       RowBox[{"KEoverH2", "[", "d\[Phi]", "]"}]}], ")"}]}]}], 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"KE", "[", 
     RowBox[{"\[Phi]_", ",", "d\[Phi]_"}], "]"}], ":=", 
    RowBox[{
     RowBox[{"H2", "[", 
      RowBox[{"\[Phi]", ",", "d\[Phi]"}], "]"}], 
     RowBox[{"KEoverH2", "[", "d\[Phi]", "]"}]}]}], "\[IndentingNewLine]", 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"weqnstate", "[", 
     RowBox[{"\[Phi]_", ",", "d\[Phi]_"}], "]"}], ":=", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{
       RowBox[{"KE", "[", 
        RowBox[{"\[Phi]", ",", "d\[Phi]"}], "]"}], "-", 
       RowBox[{"V", "[", "\[Phi]", "]"}]}], ")"}], "/", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"KE", "[", 
        RowBox[{"\[Phi]", ",", "d\[Phi]"}], "]"}], "+", 
       RowBox[{"V", "[", "\[Phi]", "]"}]}], ")"}]}]}], "\[IndentingNewLine]", 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"\[Epsilon]", "[", "d\[Phi]_", "]"}], ":=", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1.", "/", "2"}], ")"}], 
     RowBox[{"(", 
      RowBox[{"d\[Phi]", ".", "d\[Phi]"}], ")"}]}]}], "\[IndentingNewLine]", 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"a", "[", "\[Alpha]_", "]"}], ":=", 
    RowBox[{"a0", " ", 
     RowBox[{"Exp", "[", "\[Alpha]", "]"}]}]}], 
   "\[IndentingNewLine]"}]}]], "Input",
 CellChangeTimes->{{3.580674681551624*^9, 3.580674760363608*^9}, {
   3.580675130019357*^9, 3.580675133586377*^9}, {3.580675377707961*^9, 
   3.580675399253677*^9}, {3.58067548600621*^9, 3.58067549229268*^9}, {
   3.580676284943541*^9, 3.580676298150949*^9}, {3.580676708968756*^9, 
   3.580676803894329*^9}, {3.580676923743787*^9, 3.580676953972595*^9}, 
   3.580677097783474*^9, {3.580679619145484*^9, 3.580679634540499*^9}, {
   3.580682183859893*^9, 3.580682191118032*^9}, {3.580683168451305*^9, 
   3.580683169847095*^9}, {3.580684776181261*^9, 3.58068477659755*^9}, {
   3.580684827571977*^9, 3.580684828209452*^9}, {3.580684951643981*^9, 
   3.58068495449289*^9}, {3.580758961393017*^9, 3.580759018176546*^9}, {
   3.580764383917068*^9, 3.580764404163507*^9}, 3.580764530079741*^9, {
   3.580771409436032*^9, 3.580771421577495*^9}, {3.580773358504419*^9, 
   3.580773472226832*^9}, {3.580773795241287*^9, 3.580773831822586*^9}, {
   3.580776431957805*^9, 3.580776432398422*^9}}],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{"Background", " ", "equations"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"kgeqn", "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{
       RowBox[{
        RowBox[{
         RowBox[{"d\[Phi]", "'"}], "[", "\[Alpha]", "]"}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{"3", "-", 
           RowBox[{"\[Epsilon]", "[", 
            RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], "]"}]}], ")"}], 
         RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "/", 
           RowBox[{"H2", "[", 
            RowBox[{
             RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], ",", 
             RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}], "]"}]}], ")"}], 
         RowBox[{"dVd\[Phi]", "[", 
          RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], "]"}]}]}], "\[Equal]", 
       "0"}], ",", "\[IndentingNewLine]", " ", 
      RowBox[{
       RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], "==", 
       RowBox[{
        RowBox[{"\[Phi]", "'"}], "[", "\[Alpha]", "]"}]}], ",", 
      "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"\[Phi]", "[", "0", "]"}], "\[Equal]", "\[Phi]0"}], ",", 
      "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"d\[Phi]", "[", "0", "]"}], "\[Equal]", "d\[Phi]0"}]}], 
     "}"}]}], ";"}]}]], "Input",
 CellChangeTimes->{{3.580673182243431*^9, 3.580673191753762*^9}, {
   3.580673990271602*^9, 3.580674009395801*^9}, {3.580674090885081*^9, 
   3.580674091579571*^9}, 3.580675183046353*^9, {3.580675216976145*^9, 
   3.5806752340151*^9}, {3.580675281933912*^9, 3.580675358320064*^9}, {
   3.580675504082036*^9, 3.580675542061685*^9}, {3.580675580760197*^9, 
   3.580675625336132*^9}, {3.580676856923152*^9, 3.580676888361354*^9}, 
   3.580677513208761*^9, {3.580678823068553*^9, 3.580678866143317*^9}, {
   3.580679419654574*^9, 3.580679485282382*^9}, {3.580679568737988*^9, 
   3.580679577978491*^9}, {3.580679638630309*^9, 3.580679639743004*^9}, {
   3.580679712678544*^9, 3.580679728324394*^9}, {3.580680793748465*^9, 
   3.58068081158615*^9}, 3.58068097979736*^9, {3.580681509768722*^9, 
   3.58068151393727*^9}, {3.580681777179855*^9, 3.580681805970356*^9}, {
   3.580682330325386*^9, 3.580682399430808*^9}, 3.580759321245052*^9, 
   3.580760330060726*^9, {3.580762003001998*^9, 3.580762035335165*^9}, {
   3.580762069378211*^9, 3.580762072152883*^9}, 3.580762285598848*^9, {
   3.580764445698*^9, 3.580764449651905*^9}, {3.580765449270306*^9, 
   3.580765484730804*^9}, 3.580771524387938*^9, {3.580773661257434*^9, 
   3.58077367401275*^9}, 3.580773775617428*^9, {3.580773862853437*^9, 
   3.580773864824641*^9}, 3.580776586048234*^9}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"\[IndentingNewLine]", "\[IndentingNewLine]", 
  RowBox[{"(*", 
   RowBox[{
    RowBox[{"--", 
     RowBox[{"--", 
      RowBox[{"--", 
       RowBox[{"--", 
        RowBox[{"--", 
         RowBox[{"--", 
          RowBox[{"--", 
           RowBox[{
            RowBox[{
             RowBox[{
              RowBox[{
               RowBox[{
                RowBox[{
                 RowBox[{
                  RowBox[{
                   RowBox[{
                    RowBox[{
                    RowBox[{
                    RowBox[{"Solver", "--"}], "--"}], "--"}], "--"}], "--"}], 
                  "--"}], "--"}], "--"}], "--"}], "--"}], "--"}], 
            "--"}]}]}]}]}]}]}]}], "-"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{"(*", 
   RowBox[{
   "Evolve", " ", "background", " ", "until", " ", "end", " ", "of", " ", 
    "inflation"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"\[Alpha]end", "=", "10000000"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"background", "=", 
    RowBox[{"NDSolve", "[", 
     RowBox[{"kgeqn", ",", 
      RowBox[{"{", 
       RowBox[{"\[Phi]", ",", "d\[Phi]"}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"\[Alpha]", ",", "0", ",", "\[Alpha]end"}], "}"}], ",", 
      RowBox[{"Method", " ", "\[Rule]", 
       RowBox[{"{", 
        RowBox[{"\"\<EventLocator\>\"", ",", 
         RowBox[{"\"\<Event\>\"", "\[Rule]", 
          RowBox[{"{", 
           RowBox[{
            RowBox[{"weqnstate", "[", 
             RowBox[{
              RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], ",", 
              RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}], "]"}], "+", 
            RowBox[{"1", "/", "3."}]}], "}"}]}]}], "}"}]}], ",", 
      RowBox[{"Method", "\[Rule]", "\"\<BDF\>\""}]}], "]"}]}], 
   "\[IndentingNewLine]", 
   RowBox[{"\[Alpha]end", " ", "=", " ", 
    RowBox[{
     RowBox[{"Flatten", "[", 
      RowBox[{
       RowBox[{"InterpolatingFunctionDomain", "[", 
        RowBox[{"First", "[", 
         RowBox[{"\[Phi]", " ", "/.", "background"}], "]"}], "]"}], ",", 
       "1"}], "]"}], "[", 
     RowBox[{"[", "2", "]"}], "]"}]}]}]}]], "Input",
 CellChangeTimes->{{3.580680277915251*^9, 3.580680280099182*^9}, {
   3.580680400932314*^9, 3.580680418277269*^9}, {3.580680603645959*^9, 
   3.580680626378966*^9}, 3.580680864238819*^9, {3.580680917390745*^9, 
   3.580680918864588*^9}, 3.58068096463529*^9, {3.580681099058112*^9, 
   3.580681103472931*^9}, {3.580681281489835*^9, 3.580681289201413*^9}, {
   3.580682065550435*^9, 3.580682068754941*^9}, {3.580682384534755*^9, 
   3.580682391145324*^9}, {3.580683307031961*^9, 3.58068339785706*^9}, {
   3.580683442052231*^9, 3.580683491670694*^9}, {3.580683637191966*^9, 
   3.580683676842484*^9}, 3.580684437936914*^9, {3.580684527855185*^9, 
   3.580684598291595*^9}, {3.580684920352086*^9, 3.580684935753966*^9}, {
   3.580685158722582*^9, 3.580685159509133*^9}, {3.580762058301507*^9, 
   3.58076206420022*^9}, {3.580771633115419*^9, 3.580771641844033*^9}, {
   3.580771676509135*^9, 3.580771722219262*^9}, {3.580771804139577*^9, 
   3.580771809331786*^9}}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"\[Phi]", "\[Rule]", 
     TagBox[
      RowBox[{"InterpolatingFunction", "[", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"{", 
          RowBox[{"0.`", ",", "69.68134029100207`"}], "}"}], "}"}], 
        ",", "\<\"<>\"\>"}], "]"}],
      False,
      Editable->False]}], ",", 
    RowBox[{"d\[Phi]", "\[Rule]", 
     TagBox[
      RowBox[{"InterpolatingFunction", "[", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"{", 
          RowBox[{"0.`", ",", "69.68134029100207`"}], "}"}], "}"}], 
        ",", "\<\"<>\"\>"}], "]"}],
      False,
      Editable->False]}]}], "}"}], "}"}]], "Output",
 CellChangeTimes->{
  3.580680626880686*^9, 3.580680741921419*^9, 3.580680827173504*^9, 
   3.58068086584886*^9, 3.580680919550401*^9, 3.580680964991293*^9, {
   3.580681090579116*^9, 3.580681103797998*^9}, {3.580681278786426*^9, 
   3.580681296375632*^9}, {3.580681404636532*^9, 3.580681434425297*^9}, 
   3.580681518121801*^9, 3.580681763400054*^9, {3.580682062543356*^9, 
   3.580682068973827*^9}, {3.580682375889111*^9, 3.580682385114364*^9}, 
   3.580683182359728*^9, 3.580683242987722*^9, {3.580683399757669*^9, 
   3.580683492594767*^9}, {3.580683651906568*^9, 3.580683677438596*^9}, 
   3.580683739955552*^9, {3.580684419506109*^9, 3.580684438861991*^9}, {
   3.580684534667762*^9, 3.580684541441492*^9}, {3.58068457656364*^9, 
   3.580684598435437*^9}, {3.580684839948938*^9, 3.580684841150137*^9}, 
   3.58068487754685*^9, {3.580684913482597*^9, 3.580684937136863*^9}, {
   3.580685162023303*^9, 3.58068517210039*^9}, 3.580759052232645*^9, 
   3.580760276519058*^9, {3.580761518161217*^9, 3.580761520182945*^9}, 
   3.580761924667168*^9, {3.580762049273114*^9, 3.580762076584883*^9}, 
   3.580762769243443*^9, {3.580764429531456*^9, 3.580764452739483*^9}, 
   3.580766726862776*^9, 3.58077145195601*^9, 3.580771643114147*^9, {
   3.580771711991204*^9, 3.58077172726698*^9}, {3.580771760091451*^9, 
   3.580771809734677*^9}, {3.580773482006782*^9, 3.580773486545247*^9}, 
   3.580773644272273*^9, 3.58077368821571*^9, 3.580773765225314*^9, 
   3.580773822265229*^9, 3.580773868656042*^9, 3.580777106153455*^9, 
   3.580778393102353*^9, 3.580778551072476*^9}],

Cell[BoxData["69.68134029100207`"], "Output",
 CellChangeTimes->{
  3.580680626880686*^9, 3.580680741921419*^9, 3.580680827173504*^9, 
   3.58068086584886*^9, 3.580680919550401*^9, 3.580680964991293*^9, {
   3.580681090579116*^9, 3.580681103797998*^9}, {3.580681278786426*^9, 
   3.580681296375632*^9}, {3.580681404636532*^9, 3.580681434425297*^9}, 
   3.580681518121801*^9, 3.580681763400054*^9, {3.580682062543356*^9, 
   3.580682068973827*^9}, {3.580682375889111*^9, 3.580682385114364*^9}, 
   3.580683182359728*^9, 3.580683242987722*^9, {3.580683399757669*^9, 
   3.580683492594767*^9}, {3.580683651906568*^9, 3.580683677438596*^9}, 
   3.580683739955552*^9, {3.580684419506109*^9, 3.580684438861991*^9}, {
   3.580684534667762*^9, 3.580684541441492*^9}, {3.58068457656364*^9, 
   3.580684598435437*^9}, {3.580684839948938*^9, 3.580684841150137*^9}, 
   3.58068487754685*^9, {3.580684913482597*^9, 3.580684937136863*^9}, {
   3.580685162023303*^9, 3.58068517210039*^9}, 3.580759052232645*^9, 
   3.580760276519058*^9, {3.580761518161217*^9, 3.580761520182945*^9}, 
   3.580761924667168*^9, {3.580762049273114*^9, 3.580762076584883*^9}, 
   3.580762769243443*^9, {3.580764429531456*^9, 3.580764452739483*^9}, 
   3.580766726862776*^9, 3.58077145195601*^9, 3.580771643114147*^9, {
   3.580771711991204*^9, 3.58077172726698*^9}, {3.580771760091451*^9, 
   3.580771809734677*^9}, {3.580773482006782*^9, 3.580773486545247*^9}, 
   3.580773644272273*^9, 3.58077368821571*^9, 3.580773765225314*^9, 
   3.580773822265229*^9, 3.580773868656042*^9, 3.580777106153455*^9, 
   3.580778393102353*^9, 3.580778551074795*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"Show", "[", 
  RowBox[{
   RowBox[{"ListPlot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{"\[Phi]", "[", "0", "]"}], "}"}], "/.", "background"}], ",", 
     RowBox[{"PlotMarkers", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{"Automatic", ",", "Medium"}], "}"}]}], ",", 
     RowBox[{"PlotStyle", "\[Rule]", "Black"}], ",", 
     RowBox[{"Frame", "\[Rule]", "True"}], ",", 
     RowBox[{"FrameLabel", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"Style", "[", 
         RowBox[{
          SubscriptBox["\[Phi]", "1"], ",", "20"}], "]"}], ",", 
        RowBox[{"Style", "[", 
         RowBox[{
          SubscriptBox["\[Phi]", "2"], ",", "20"}], "]"}]}], "}"}]}]}], "]"}],
    ",", 
   RowBox[{"ParametricPlot", "[", 
    RowBox[{
     RowBox[{"Evaluate", "[", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"\[Phi]", "[", "t", "]"}], "}"}], "/.", "background"}], "]"}],
      ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "\[Alpha]end"}], "}"}], ",", 
     RowBox[{"Frame", "\[Rule]", "True"}], ",", 
     RowBox[{"PlotStyle", "\[Rule]", 
      RowBox[{"{", "Thick", "}"}]}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}]}], "]"}]}], 
  "]"}], "\[IndentingNewLine]", 
 RowBox[{"Show", "[", 
  RowBox[{
   RowBox[{"ParametricPlot", "[", 
    RowBox[{
     RowBox[{"Evaluate", "[", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"d\[Phi]", "[", "t", "]"}], "}"}], "/.", "background"}], 
      "]"}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "\[Alpha]end"}], "}"}], ",", 
     RowBox[{"PlotStyle", "\[Rule]", 
      RowBox[{"{", "Thick", "}"}]}], ",", 
     RowBox[{"PlotRange", "\[Rule]", "All"}], ",", 
     RowBox[{"Frame", "\[Rule]", "True"}], ",", 
     RowBox[{"FrameLabel", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"Style", "[", 
         RowBox[{
          SubscriptBox["d\[Phi]", "1"], ",", "20"}], "]"}], ",", 
        RowBox[{"Style", "[", 
         RowBox[{
          SubscriptBox["d\[Phi]", "2"], ",", "20"}], "]"}]}], "}"}]}]}], 
    "]"}], ",", 
   RowBox[{"ListPlot", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{"d\[Phi]", "[", "0", "]"}], "}"}], "/.", "background"}], ",", 
     RowBox[{"PlotMarkers", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{"Automatic", ",", "Medium"}], "}"}]}], ",", 
     RowBox[{"PlotStyle", "\[Rule]", "Black"}], ",", 
     RowBox[{"Frame", "\[Rule]", "True"}]}], "]"}]}], "]"}]}], "Input",
 CellChangeTimes->{{3.580680893952681*^9, 3.580680898960688*^9}, {
  3.580680933437632*^9, 3.580680959485598*^9}, {3.580681113476451*^9, 
  3.580681166497718*^9}, {3.580681200283665*^9, 3.580681201867835*^9}, {
  3.580681239131625*^9, 3.580681292647827*^9}, {3.580684616378604*^9, 
  3.580684630443976*^9}, {3.580774066100922*^9, 3.580774121185184*^9}, {
  3.580776698024071*^9, 3.580776790641773*^9}, {3.58077682353939*^9, 
  3.580776850916659*^9}}],

Cell[BoxData[
 GraphicsBox[{
   GraphicsComplexBox[{{10.31001, 12.93651}, {10.31001, 12.93651}, {10.31001, 
    12.93651}, {10.31001, 12.93651}}, {
     {GrayLevel[0], InsetBox[
       StyleBox["\<\"\[FilledCircle]\"\>",
        StripOnInput->False,
        FontSize->Medium], 3], InsetBox[
       StyleBox["\<\"\[FilledCircle]\"\>",
        StripOnInput->False,
        FontSize->Medium], 4]}, {}}], {{}, {}, 
    {Hue[0.67, 0.6, 0.6], Thickness[Large], LineBox[CompressedData["
1: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
JS2WRCbwPwL1DFtMsRK+qAte/Q==
      "]]}}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->True,
  AxesOrigin->{0, 0},
  Frame->True,
  FrameLabel->{
    FormBox[
     StyleBox[
      SubscriptBox["\[Phi]", "1"], 20, StripOnInput -> False], 
     TraditionalForm], 
    FormBox[
     StyleBox[
      SubscriptBox["\[Phi]", "2"], 20, StripOnInput -> False], 
     TraditionalForm]},
  ImageSize->{571.908102865759, Automatic},
  PlotRange->{{0., 20.62002}, {0., 25.87302}},
  PlotRangeClipping->True,
  PlotRangePadding->{
    Scaled[0.02], 
    Scaled[0.02]}]], "Output",
 CellChangeTimes->{{3.580680946028036*^9, 3.580680967025308*^9}, {
   3.580681093431184*^9, 3.580681130616556*^9}, 3.580681167196599*^9, 
   3.580681202490477*^9, {3.58068123972228*^9, 3.58068129787571*^9}, 
   3.58068176605229*^9, {3.580684603617066*^9, 3.580684630838601*^9}, 
   3.580759056656705*^9, 3.580766735273832*^9, 3.580771454987028*^9, 
   3.580771646110812*^9, 3.580771732070261*^9, {3.580771791884657*^9, 
   3.580771811325401*^9}, {3.580773484164799*^9, 3.580773487661499*^9}, 
   3.580773870251924*^9, {3.58077407705163*^9, 3.580774121613592*^9}, 
   3.58077432033151*^9, 3.580776629433876*^9, {3.580776716981736*^9, 
   3.580776791098535*^9}, {3.580776832299159*^9, 3.580776851600577*^9}, 
   3.580777107297225*^9, 3.580778552387501*^9}],

Cell[BoxData[
 GraphicsBox[{{{}, {}, 
    {Hue[0.67, 0.6, 0.6], Thickness[Large], LineBox[CompressedData["
1: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
v9xHTFOqnDQ+d690Ig==
      "]]}}, 
   GraphicsComplexBox[{{-0.0015093069460898748`, -0.15339833002044168`}, \
{-0.0015093069460898748`, -0.15339833002044168`}, {-0.0015093069460898748`, \
-0.15339833002044168`}, {-0.0015093069460898748`, -0.15339833002044168`}}, {
     {GrayLevel[0], InsetBox[
       StyleBox["\<\"\[FilledCircle]\"\>",
        StripOnInput->False,
        FontSize->Medium], 3], InsetBox[
       StyleBox["\<\"\[FilledCircle]\"\>",
        StripOnInput->False,
        FontSize->Medium], 4]}, {}}]},
  Axes->True,
  AxesOrigin->{0, 0},
  Frame->True,
  FrameLabel->{
    FormBox[
     StyleBox[
      SubscriptBox["d\[Phi]", "1"], 20, StripOnInput -> False], 
     TraditionalForm], 
    FormBox[
     StyleBox[
      SubscriptBox["d\[Phi]", "2"], 20, StripOnInput -> False], 
     TraditionalForm]},
  ImageSize->{571.908102865759, Automatic},
  PlotRange->{All, All},
  PlotRangeClipping->True,
  PlotRangePadding->{Automatic, Automatic}]], "Output",
 CellChangeTimes->{{3.580680946028036*^9, 3.580680967025308*^9}, {
   3.580681093431184*^9, 3.580681130616556*^9}, 3.580681167196599*^9, 
   3.580681202490477*^9, {3.58068123972228*^9, 3.58068129787571*^9}, 
   3.58068176605229*^9, {3.580684603617066*^9, 3.580684630838601*^9}, 
   3.580759056656705*^9, 3.580766735273832*^9, 3.580771454987028*^9, 
   3.580771646110812*^9, 3.580771732070261*^9, {3.580771791884657*^9, 
   3.580771811325401*^9}, {3.580773484164799*^9, 3.580773487661499*^9}, 
   3.580773870251924*^9, {3.58077407705163*^9, 3.580774121613592*^9}, 
   3.58077432033151*^9, 3.580776629433876*^9, {3.580776716981736*^9, 
   3.580776791098535*^9}, {3.580776832299159*^9, 3.580776851600577*^9}, 
   3.580777107297225*^9, 3.580778552450517*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Show", "[", 
  RowBox[{"Plot", "[", 
   RowBox[{
    RowBox[{
     RowBox[{"\[Epsilon]", "[", 
      RowBox[{"d\[Phi]", "[", "t", "]"}], "]"}], "/.", "background"}], ",", 
    RowBox[{"{", 
     RowBox[{"t", ",", "0", ",", "\[Alpha]end"}], "}"}], ",", 
    RowBox[{"PlotRange", "\[Rule]", "All"}], ",", 
    RowBox[{"Frame", "\[Rule]", "True"}], ",", 
    RowBox[{"FrameLabel", "\[Rule]", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"Style", "[", 
        RowBox[{"\[Alpha]", ",", "20"}], "]"}], ",", 
       RowBox[{"Style", "[", 
        RowBox[{"\[Epsilon]", ",", "20"}], "]"}]}], "}"}]}]}], "]"}], 
  "]"}]], "Input",
 CellChangeTimes->{{3.580773493062811*^9, 3.580773616256361*^9}, {
  3.580773750650866*^9, 3.580773761294973*^9}, {3.580773876699273*^9, 
  3.58077392254605*^9}, {3.580773973636447*^9, 3.580774042792511*^9}, {
  3.58077430699439*^9, 3.580774309367785*^9}, {3.58077435038489*^9, 
  3.580774409479475*^9}, {3.580776865960073*^9, 3.580776884010333*^9}}],

Cell[BoxData[
 GraphicsBox[{{}, {}, 
   {Hue[0.67, 0.6, 0.6], LineBox[CompressedData["
1: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
BTQ=
     "]]}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->True,
  AxesOrigin->{0, 0},
  Frame->True,
  FrameLabel->{
    FormBox[
     StyleBox["\[Alpha]", 20, StripOnInput -> False], TraditionalForm], 
    FormBox[
     StyleBox["\[Epsilon]", 20, StripOnInput -> False], TraditionalForm]},
  ImageSize->{676., Automatic},
  PlotRange->{All, All},
  PlotRangeClipping->True,
  PlotRangePadding->{Automatic, Automatic}]], "Output",
 CellChangeTimes->{{3.580776859812463*^9, 3.580776884717151*^9}, 
   3.580777119533017*^9, 3.580778555370605*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"x", "[", "t_", "]"}], ":=", 
  RowBox[{"Flatten", "[", 
   RowBox[{
    RowBox[{"\[Phi]", "[", "t", "]"}], "/.", "background"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"x", "[", "1.2309005", "]"}]}], "Input",
 CellChangeTimes->{{3.580777179938722*^9, 3.580777233720886*^9}}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"10.30813590891968`", ",", "12.746869126640565`"}], "}"}]], "Output",\

 CellChangeTimes->{{3.580777216241269*^9, 3.580777234083849*^9}}]
}, Open  ]],

Cell[BoxData[
 RowBox[{"\[IndentingNewLine]", "\[IndentingNewLine]", "\[IndentingNewLine]", 
  "\[IndentingNewLine]", "\[IndentingNewLine]", "\[IndentingNewLine]", 
  "\[IndentingNewLine]", "\[IndentingNewLine]", "\[IndentingNewLine]", 
  RowBox[{"(*", 
   RowBox[{
    RowBox[{"Ptb", " ", "ICs", " ", "for", " ", "mode", " ", "matrix", " ", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{
        RowBox[{"See", " ", "Salopek"}], ",", " ", "Bond", ",", " ", 
        "Bardeen"}], ")"}], " ", "--"}]}], "-", " ", "Bunch", "-", "Davies"}],
    "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Clear", "[", 
    RowBox[{"\[Psi]0", ",", "d\[Psi]0"}], "]"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"\[Psi]0", "[", "k_", "]"}], ":=", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1.", "/", 
       RowBox[{"Sqrt", "[", 
        RowBox[{"2", "k"}], "]"}]}], ")"}], 
     RowBox[{"IdentityMatrix", "[", 
      RowBox[{"Length", "[", "\[Phi]0", "]"}], "]"}]}]}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"d\[Psi]0", "[", 
     RowBox[{"k_", ",", "\[Phi]0_", ",", "d\[Phi]0_"}], "]"}], ":=", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{
       RowBox[{"I", "/", "a0"}], "/", " ", 
       RowBox[{"H", "[", 
        RowBox[{"\[Phi]0", ",", "d\[Phi]0"}], "]"}]}], ")"}], 
     RowBox[{"Sqrt", "[", 
      RowBox[{"k", "/", "2."}], "]"}], 
     RowBox[{"IdentityMatrix", "[", 
      RowBox[{"Length", "[", "\[Phi]0", "]"}], "]"}]}]}], 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"kpivot", "=", "1."}], ";"}], "\[IndentingNewLine]"}]}]], "Input",\

 CellChangeTimes->{{3.580677514716848*^9, 3.580677547165261*^9}, {
   3.58067763171722*^9, 3.580677642842213*^9}, {3.580677691201236*^9, 
   3.580677728165931*^9}, {3.580679843069592*^9, 3.580679848817861*^9}, 
   3.580679892823028*^9, {3.580680668070113*^9, 3.580680686535397*^9}, {
   3.580752809607461*^9, 3.580752812456493*^9}, {3.580773489714766*^9, 
   3.580773491346333*^9}, {3.580773928098813*^9, 3.580773939284613*^9}, {
   3.580776910123605*^9, 3.580776920931807*^9}}],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{"Functions", " ", "for", " ", "ptbs"}], "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"MatMul", "[", 
     RowBox[{"vecx_", ",", "vecy_"}], "]"}], ":=", 
    RowBox[{"Outer", "[", 
     RowBox[{"Times", ",", "vecx", ",", "vecy"}], "]"}]}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"ptbmassmatrix", "[", 
     RowBox[{"\[Phi]_", ",", "d\[Phi]_"}], "]"}], ":=", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1.", "/", 
        RowBox[{"H2", "[", 
         RowBox[{"\[Phi]", ",", "d\[Phi]"}], "]"}]}], ")"}], 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"m2matrix", "[", "\[Phi]", "]"}], "+", 
        RowBox[{"MatMul", "[", 
         RowBox[{"d\[Phi]", ",", 
          RowBox[{"dVd\[Phi]", "[", "\[Phi]", "]"}]}], "]"}], "+", 
        RowBox[{"MatMul", "[", 
         RowBox[{
          RowBox[{"dVd\[Phi]", "[", "\[Phi]", "]"}], ",", "d\[Phi]"}], 
         "]"}]}], ")"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"3.", "-", 
        RowBox[{"\[Epsilon]", "[", "d\[Phi]", "]"}]}], ")"}], 
      RowBox[{"MatMul", "[", 
       RowBox[{"d\[Phi]", ",", "d\[Phi]"}], "]"}]}]}]}], 
   "\[IndentingNewLine]", "\[IndentingNewLine]"}]}]], "Input",
 CellChangeTimes->{{3.580677017382814*^9, 3.580677049393753*^9}, {
   3.580677145465765*^9, 3.580677215646185*^9}, {3.580677320178735*^9, 
   3.580677371538552*^9}, {3.580678095791572*^9, 3.580678102750582*^9}, {
   3.580678136011087*^9, 3.580678140068411*^9}, {3.580678318019066*^9, 
   3.580678351319561*^9}, {3.580678383730328*^9, 3.580678392645887*^9}, {
   3.580678632775383*^9, 3.580678764824156*^9}, {3.580679911800038*^9, 
   3.580679912441682*^9}, 3.58068163227066*^9, 3.580681752321387*^9, {
   3.580766527983409*^9, 3.580766528962391*^9}, {3.58077693819866*^9, 
   3.580776968481544*^9}, 3.580776999920064*^9}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{
    RowBox[{"Ptb", " ", 
     RowBox[{"equations", " ", "--"}]}], "-", " ", 
    RowBox[{"Fourier", " ", "transformed"}]}], "*)"}], "\[IndentingNewLine]", 
  
  RowBox[{"ptbeqn", "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{
        RowBox[{"d\[Psi]", "'"}], "[", "\[Alpha]", "]"}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "-", 
          RowBox[{"\[Epsilon]", "[", 
           RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], "]"}]}], ")"}], 
        RowBox[{"d\[Psi]", "[", "\[Alpha]", "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          SuperscriptBox[
           RowBox[{"(", 
            RowBox[{
             RowBox[{"kpivot", "/", 
              RowBox[{"a", "[", "\[Alpha]", "]"}]}], " ", 
             RowBox[{"H", "[", 
              RowBox[{
               RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], ",", 
               RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}], "]"}]}], ")"}], 
           "2"], "-", "2", "+", 
          RowBox[{"\[Epsilon]", "[", 
           RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], "]"}]}], ")"}], 
        RowBox[{"\[Psi]", "[", "\[Alpha]", "]"}]}], "+", 
       RowBox[{
        RowBox[{"ptbmassmatrix", "[", 
         RowBox[{
          RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], ",", 
          RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}], "]"}], ".", 
        RowBox[{"\[Psi]", "[", "\[Alpha]", "]"}]}]}], "\[Equal]", "0"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{
       RowBox[{"\[Psi]", "'"}], "[", "\[Alpha]", "]"}], "\[Equal]", 
      RowBox[{"d\[Psi]", "[", "\[Alpha]", "]"}]}], ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"\[Psi]", "[", "0", "]"}], "\[Equal]", 
      RowBox[{"\[Psi]0", "[", "kpivot", "]"}]}], ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"d\[Psi]", "[", "0", "]"}], "\[Equal]", 
      RowBox[{"d\[Psi]0", "[", 
       RowBox[{"kpivot", ",", "\[Phi]0ptb", ",", "d\[Phi]0ptb"}], "]"}]}]}], 
    "}"}]}]}]], "Input",
 CellChangeTimes->{{3.580678117953779*^9, 3.580678126050189*^9}, {
   3.580678787067428*^9, 3.580678817313135*^9}, {3.580678852283752*^9, 
   3.580678900225189*^9}, {3.580679120123327*^9, 3.580679313894092*^9}, {
   3.580679344230327*^9, 3.580679408985338*^9}, {3.580679696766883*^9, 
   3.580679702875322*^9}, {3.580679754462746*^9, 3.580679834215451*^9}, {
   3.580679885269579*^9, 3.580679886341009*^9}, {3.580680698987757*^9, 
   3.580680720587971*^9}, 3.580680777497466*^9, {3.580776948812899*^9, 
   3.580777017490957*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Outer", "::", "heads"}], "MessageName"], 
  RowBox[{
  ":", " "}], "\<\"Heads \[NoBreak]\\!\\(dVd\[Phi]\\)\[NoBreak] and \
\[NoBreak]\\!\\(d\[Phi]\\)\[NoBreak] at positions \[NoBreak]\\!\\(3\\)\
\[NoBreak] and \[NoBreak]\\!\\(2\\)\[NoBreak] are expected to be the same. \
\\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", ButtonStyle->\\\"Link\\\", \
ButtonFrame->None, ButtonData:>\\\"paclet:ref/message/General/heads\\\", \
ButtonNote -> \\\"Outer::heads\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{
  3.580680778166793*^9, {3.580776983908179*^9, 3.580777008564784*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Outer", "::", "heads"}], "MessageName"], 
  RowBox[{
  ":", " "}], "\<\"Heads \[NoBreak]\\!\\(d\[Phi]\\)\[NoBreak] and \
\[NoBreak]\\!\\(dVd\[Phi]\\)\[NoBreak] at positions \[NoBreak]\\!\\(3\\)\
\[NoBreak] and \[NoBreak]\\!\\(2\\)\[NoBreak] are expected to be the same. \
\\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", ButtonStyle->\\\"Link\\\", \
ButtonFrame->None, ButtonData:>\\\"paclet:ref/message/General/heads\\\", \
ButtonNote -> \\\"Outer::heads\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{
  3.580680778166793*^9, {3.580776983908179*^9, 3.580777008567462*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"IdentityMatrix", "::", "dims"}], "MessageName"], 
  RowBox[{
  ":", " "}], "\<\"Dimension specification \[NoBreak]\\!\\(0\\)\[NoBreak] \
should be a positive machine integer or a pair of positive machine integers. \
\\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", ButtonStyle->\\\"Link\\\", \
ButtonFrame->None, ButtonData:>\\\"paclet:ref/IdentityMatrix\\\", ButtonNote \
-> \\\"IdentityMatrix::dims\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{
  3.580680778166793*^9, {3.580776983908179*^9, 3.580777008569759*^9}}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{"3.`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.5`", " ", 
            RowBox[{
             RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], ".", 
             RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}]}]}], ")"}], " ", 
         RowBox[{"d\[Phi]", "[", 
          RowBox[{"d\[Phi]", "[", 
           SuperscriptBox["\[Alpha]", "2"], "]"}], "]"}]}], "+", 
        FractionBox[
         RowBox[{"2.`", " ", 
          RowBox[{"(", 
           RowBox[{"3.`", "\[VeryThinSpace]", "-", 
            RowBox[{"0.5`", " ", 
             RowBox[{
              RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], ".", 
              RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}]}]}], ")"}], " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{"m2matrix", "[", 
             RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], "]"}], "+", 
            RowBox[{"Outer", "[", 
             RowBox[{"Times", ",", 
              RowBox[{"dVd\[Phi]", "[", 
               RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], "]"}], ",", 
              RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}], "]"}], "+", 
            RowBox[{"Outer", "[", 
             RowBox[{"Times", ",", 
              RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], ",", 
              RowBox[{"dVd\[Phi]", "[", 
               RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], "]"}]}], "]"}]}], 
           ")"}]}], 
         RowBox[{
          RowBox[{"{", 
           RowBox[{"3.776634469324246`*^-11", ",", "3.059073920152644`*^-9"}],
            "}"}], ".", 
          SuperscriptBox[
           RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], "2"]}]]}], ")"}], ".", 
      RowBox[{"\[Psi]", "[", "\[Alpha]", "]"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "-", 
        RowBox[{"0.5`", " ", 
         RowBox[{
          RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], ".", 
          RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}]}]}], ")"}], " ", 
      RowBox[{"d\[Psi]", "[", "\[Alpha]", "]"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "2"}], "+", 
        RowBox[{"0.5`", " ", 
         RowBox[{
          RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], ".", 
          RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}]}], "+", 
        FractionBox[
         RowBox[{"0.5000000000000001`", " ", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{
            RowBox[{"-", "2"}], " ", "\[Alpha]"}]], " ", 
          RowBox[{
           RowBox[{"{", 
            RowBox[{
            "3.776634469324246`*^-11", ",", "3.059073920152644`*^-9"}], "}"}],
            ".", 
           SuperscriptBox[
            RowBox[{"\[Phi]", "[", "\[Alpha]", "]"}], "2"]}]}], 
         RowBox[{"3.`", "\[VeryThinSpace]", "-", 
          RowBox[{"0.5`", " ", 
           RowBox[{
            RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}], ".", 
            RowBox[{"d\[Phi]", "[", "\[Alpha]", "]"}]}]}]}]]}], ")"}], " ", 
      RowBox[{"\[Psi]", "[", "\[Alpha]", "]"}]}], "+", 
     RowBox[{
      SuperscriptBox["d\[Psi]", "\[Prime]",
       MultilineFunction->None], "[", "\[Alpha]", "]"}]}], "\[Equal]", "0"}], 
   ",", 
   RowBox[{
    RowBox[{
     SuperscriptBox["\[Psi]", "\[Prime]",
      MultilineFunction->None], "[", "\[Alpha]", "]"}], "\[Equal]", 
    RowBox[{"d\[Psi]", "[", "\[Alpha]", "]"}]}], ",", 
   RowBox[{
    RowBox[{"\[Psi]", "[", "0", "]"}], "\[Equal]", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"0.7071067811865475`", ",", "0.`"}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"0.`", ",", "0.7071067811865475`"}], "}"}]}], "}"}]}], ",", 
   RowBox[{
    RowBox[{"d\[Psi]", "[", "0", "]"}], "\[Equal]", 
    FractionBox[
     RowBox[{
      RowBox[{"(", 
       RowBox[{"0.`", "\[VeryThinSpace]", "+", 
        RowBox[{"1.`", " ", "\[ImaginaryI]"}]}], ")"}], " ", 
      RowBox[{"IdentityMatrix", "[", "0", "]"}]}], 
     SqrtBox[
      FractionBox[
       RowBox[{
        RowBox[{"{", 
         RowBox[{"3.776634469324246`*^-11", ",", "3.059073920152644`*^-9"}], 
         "}"}], ".", 
        SuperscriptBox["\[Phi]0ptb", "2"]}], 
       RowBox[{"3.`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.5`", " ", 
         RowBox[{"d\[Phi]0ptb", ".", "d\[Phi]0ptb"}]}]}]]]]}]}], 
  "}"}]], "Output",
 CellChangeTimes->{{3.580776986447682*^9, 3.580777008573553*^9}}]
}, Open  ]],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{
   RowBox[{"End", "[", "]"}], "\[IndentingNewLine]", 
   RowBox[{"EndPackage", "[", "]"}]}], "*)"}]], "Input",
 CellChangeTimes->{{3.580778339294511*^9, 3.580778351194725*^9}, {
  3.580778513499149*^9, 3.580778515149114*^9}}]
},
WindowSize->{1621, 1026},
WindowMargins->{{0, Automatic}, {Automatic, 0}},
FrontEndVersion->"8.0 for Linux x86 (32-bit) (October 10, 2011)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[557, 20, 258, 5, 30, "Input"],
Cell[818, 27, 916, 23, 107, "Input"],
Cell[1737, 52, 214, 4, 30, "Input"],
Cell[1954, 58, 654, 20, 107, "Input"],
Cell[2611, 80, 2170, 65, 126, "Input"],
Cell[4784, 147, 1650, 35, 126, "Input"],
Cell[6437, 184, 3544, 91, 278, "Input"],
Cell[9984, 277, 2716, 57, 107, "Input"],
Cell[CellGroupData[{
Cell[12725, 338, 3106, 74, 145, "Input"],
Cell[15834, 414, 2252, 46, 30, "Output"],
Cell[18089, 462, 1613, 23, 30, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19739, 490, 2949, 78, 88, "Input"],
Cell[22691, 570, 24155, 409, 385, "Output"],
Cell[46849, 981, 24771, 417, 401, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[71657, 1403, 1002, 23, 30, "Input"],
Cell[72662, 1428, 12190, 208, 447, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[84889, 1641, 317, 8, 50, "Input"],
Cell[85209, 1651, 176, 4, 30, "Output"]
}, Open  ]],
Cell[85400, 1658, 2079, 49, 316, "Input"],
Cell[87482, 1709, 1889, 46, 107, "Input"],
Cell[CellGroupData[{
Cell[89396, 1759, 2594, 63, 116, "Input"],
Cell[91993, 1824, 625, 12, 24, "Message"],
Cell[92621, 1838, 625, 12, 24, "Message"],
Cell[93249, 1852, 570, 11, 24, "Message"],
Cell[93822, 1865, 4463, 118, 134, "Output"]
}, Open  ]],
Cell[98300, 1986, 266, 6, 50, "Input"]
}
]
*)

(* End of internal cache information *)
