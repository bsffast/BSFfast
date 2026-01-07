(* ::Package:: *)

(*** BSFfast ***)
(*version 4;
  author: Stefan Lederer [TUM, NCBJ];
  created: 2026-01-07*)
BeginPackage["BSFfast`"];

BSFfastSigmavEffBSF::usage = "BSFfast\[Sigma]vtopPlat[ModelString(,alpha)][ mass , x ]\nEffective BSF annihilation cross-section for the specified models:\n     \"QCD-_U\", \"QCD-_D\", \"QCD-_\", \"dQCD-_\", \"dQED-_\", \"dQED-_nT\", \"dQED-F\", \"dQED-FnoTop\".  The \"_\" must be filled by \"S\" (scalar) or \"F\" (fermion). \nQCD refers to SM-QCD and involves running couplings.\n m = constituent mass is here expected in GeV. x = m/T is inverse temperature. alpha defines the gauge coupling strength in \"dQCD\" and \"dQED\" models, and the prescription of the running coupling \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]<1GeV)=const or =0 for \"QCD\" models. (not required for \"QED-_\".)\n Call Print[Information[#,\"Usage\"]]&/@BSFfastFunctionNames for more details.";
BSFfastFunctionNames::usage = "List of BSFfast`Private` interpolation functions of all different models.";


Begin["Private`"];
BSFfastDataDir=DirectoryName[$InputFileName]<>"/BSFfast_Data/";


BSFfastFunctionNames:={BSFfast\[Sigma]vtopPlat ,BSFfast\[Sigma]vtopCut ,BSFfast\[Sigma]vbotPlat ,BSFfast\[Sigma]vbotCut ,BSFfast\[Sigma]vnoTPlat ,BSFfast\[Sigma]vnoTCut ,BSFfast\[Sigma]vdqedConst ,BSFfast\[Sigma]vdqedNoTConst ,BSFfast\[Sigma]vqcdNoTConst ,
					   BSFfast\[Sigma]vtopPlatF,BSFfast\[Sigma]vtopCutF,BSFfast\[Sigma]vbotPlatF,BSFfast\[Sigma]vbotCutF,BSFfast\[Sigma]vnoTPlatF,BSFfast\[Sigma]vnoTCutF,BSFfast\[Sigma]vdqedConstF,BSFfast\[Sigma]vdqedNoTConstF,BSFfast\[Sigma]vqcdNoTConstF
					   BSFfast\[Sigma]vqedInclTopF,BSFfast\[Sigma]vqedExclTopF
					   };


BSFfast\[Sigma]vtopPlat::usage = "BSFfast\[Sigma]vtopPlat[ mass , x ]\nEffective BSF annihilation cross-section for a stop-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(+2\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong running is frozen once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu])==1. (\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) stays 1 at all lower scales.)\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vtopCut::usage  = "BSFfast\[Sigma]vtopCut[ mass , x ]\nEffective BSF annihilation cross-section for a stop-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(+2\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong coupling is set to 0 once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]) > 1.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vbotPlat::usage = "BSFfast\[Sigma]vbotPlat[ mass , x ]\nEffective BSF annihilation cross-section for a sbottom-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(-1\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong running is frozen once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu])==1. (\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) stays 1 at all lower scales.)\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vbotCut::usage  = "BSFfast\[Sigma]vbotCut[ mass , x ]\nEffective BSF annihilation cross-section for a sbottom-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(-1\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong coupling is set to 0 once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]) > 1.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vnoTPlat::usage = "BSFfast\[Sigma]vnoTPlat[ mass , x ]\nEffective BSF annihilation cross-section for a colored scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(0\)]\) multiplet) particle anti-particle pair. (There are no allowed bound-to-bound transitions)\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong running is frozen once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu])==1. (\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) stays 1 at all lower scales.)\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vnoTCut::usage  = "BSFfast\[Sigma]vnoTCut[ mass , x ]\nEffective BSF annihilation cross-section for a colored scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(0\)]\) multiplet) particle anti-particle pair. (There are no allowed bound-to-bound transitions)\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n  Strong coupling is set to 0 once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]) > 1.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vdqedConst::usage = "BSFfast\[Sigma]vdqedConst[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a U(1)-interacting particle (charge 1) including bound-to-bound transitions.\n alpha = the \!\(\*
StyleBox[\"constant\",\nFontVariations->{\"Underline\"->True}]\) gauge coupling-strength g^2/(4Pi) of the SU(3) interaction.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vdqedNoTConst::usage = "BSFfast\[Sigma]vdqedNoTConst[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a U(1)-interacting particle (charge 1) NEGLECTING bound-to-bound transitions.\n alpha = the \!\(\*
StyleBox[\"constant\",\nFontVariations->{\"Underline\"->True}]\) gauge coupling-strength g^2/(4Pi) of the SU(3) interaction.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vqcdNoTConst::usage = "BSFfast\[Sigma]vqcdNoTConst[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a fundamental SU(3)-interacting particle anti-particle pair. (There are no allowed bound-to-bound transitions)\n alpha = the \!\(\*
StyleBox[\"constant\",\nFontVariations->{\"Underline\"->True}]\) gauge coupling-strength g^2/(4Pi) of the SU(3) interaction.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";

(*fermions*)
BSFfast\[Sigma]vtopPlatF::usage = "same as BSFfast\[Sigma]vtopPlat[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vtopCutF::usage  = "same as BSFfast\[Sigma]vtopCut[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vbotPlatF::usage = "same as BSFfast\[Sigma]vbotPlat[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vbotCutF::usage  = "same as BSFfast\[Sigma]vbotCut[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vnoTPlatF::usage = "same as BSFfast\[Sigma]vnoTPlat[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vnoTCutF::usage  = "same as BSFfast\[Sigma]vnoTCut[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vdqedConstF::usage = "same as BSFfast\[Sigma]vdqedConst[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vdqedNoTConstF::usage = "same as BSFfast\[Sigma]vdqedNoTConst[ mass , x ] but for spin-1/2 constituents rathern than scalars.";
BSFfast\[Sigma]vqcdNoTConstF::usage = "same as BSFfast\[Sigma]vqcdNoTConst[ mass , x ] but for spin-1/2 constituents rathern than scalars.";

BSFfast\[Sigma]vqedInclTopF::usage = "BSFfast\[Sigma]vqedInclTop[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a U(1)-interacting fermion (charge 1, spin 1/2) including bound-to-bound transitions and decay into all SM fermions via s-channel photons.\n alpha = the fine structure constant, by default 1/128.9.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vqedExclTopF::usage = "BSFfast\[Sigma]vqedInclTop[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a U(1)-interacting fermion (charge 1, spin 1/2) including bound-to-bound transitions and decay into SM fermions (excluding t\!\(\*OverscriptBox[\(t\), \(_\)]\)) via s-channel photons.\n alpha = the fine structure constant, by default 1/128.9.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";


(*import of interpolation tables*)
SetDirectory[BSFfastDataDir];

	(*SCALAR*)
(*sbottom & stop*)
topPlatS=Import["QCD-SU_plateau.csv","CSV"];
topCutS=Import["QCD-SU_cutoff.csv","CSV"];
botPlatS=Import["QCD-SD_plateau.csv","CSV"];
botCutS=Import["QCD-SD_cutoff.csv","CSV"];
(*no transition*)
noTCutS=Import["QCD-S_cutoff.csv","CSV"];
noTPlatS=Import["QCD-S_plateau.csv","CSV"];
(*constant coupling*)
xQEDS=Import["dQED-S.csv","CSV"];
xNoTQEDS=Import["dQED-SnoTr.csv","CSV"];
xNoTQCDS=Import["dQCD-S.csv","CSV"];

	(*FERMION*)
(*sbottom & stop*)
topPlatF=Import["QCD-FU_plateau.csv","CSV"];
topCutF=Import["QCD-FU_cutoff.csv","CSV"];
botPlatF=Import["QCD-FD_plateau.csv","CSV"];
botCutF=Import["QCD-FD_cutoff.csv","CSV"];
(*no transition*)
noTCutF=Import["QCD-F_cutoff.csv","CSV"];
noTPlatF=Import["QCD-F_plateau.csv","CSV"];
(*constant coupling*)
xdQEDF=Import["dQED-F.csv","CSV"];
xNoTdQEDF=Import["dQED-FnoTr.csv","CSV"];
xNoTQCDF=Import["dQCD-F.csv","CSV"];
(*QED including S1 decays into fermions*)
xQEDFinclTop=Import["QED-F_inclTop.csv","CSV"];
xQEDFexclTop=Import["QED-F_exclTop.csv","CSV"];

ResetDirectory[];

(* define interpolation functions *)
	(*SCALAR*)
	(*QCD-running couplings*)
BSFfast\[Sigma]vtopPlatS[m_,x_]=(Exp@*Evaluate[Interpolation[Log@topPlatS,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vtopCutS[m_,x_]=(Exp@*Evaluate[Interpolation[Log@topCutS,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vbotPlatS[m_,x_]=(Exp@*Evaluate[Interpolation[Log@botPlatS,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vbotCutS[m_,x_]=(Exp@*Evaluate[Interpolation[Log@botCutS,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vnoTCutS[m_,x_]=(Exp@*Evaluate[Interpolation[Log@noTCutS,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vnoTPlatS[m_,x_]=(Exp@*Evaluate[Interpolation[Log@noTPlatS,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
	(*include rescaling to any coupling for asConst tables*)
BSFfast\[Sigma]vdqedConstS[\[Alpha]_][m_,x_]=Block[{r=\[Alpha]/0.1,a0=0.1,m0=xQEDS[[1,1]]},r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xQEDS[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vdqedNoTConstS[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQEDS[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTQEDS[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqcdNoTConstS[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQCDS[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTQCDS[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
	
	(*FERMION*)
	(*QCD-running couplings*)
BSFfast\[Sigma]vtopPlatF[m_,x_]=(Exp@*Evaluate[Interpolation[Log@topPlatF,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vtopCutF[m_,x_]=(Exp@*Evaluate[Interpolation[Log@topCutF,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vbotPlatF[m_,x_]=(Exp@*Evaluate[Interpolation[Log@botPlatF,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vbotCutF[m_,x_]=(Exp@*Evaluate[Interpolation[Log@botCutF,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vnoTCutF[m_,x_]=(Exp@*Evaluate[Interpolation[Log@noTCutF,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vnoTPlatF[m_,x_]=(Exp@*Evaluate[Interpolation[Log@noTPlatF,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
	(*include rescaling to any coupling for asConst tables*)
BSFfast\[Sigma]vdqedConstF[\[Alpha]_][m_,x_]=Block[{r=\[Alpha]/0.1,a0=0.1,m0=xdQEDF[[1,1]]},r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xdQEDF[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vdqedNoTConstF[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTdQEDF[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTdQEDF[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqcdNoTConstF[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQCDF[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTQCDF[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqedInclTopF[\[Alpha]_][m_,x_]=Block[{r=\[Alpha]/0.1,a0=0.1,m0=xQEDFinclTop[[1,1]]},r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xQEDFinclTop[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqedExclTopF[\[Alpha]_][m_,x_]=Block[{r=\[Alpha]/0.1,a0=0.1,m0=xQEDFexclTop[[1,1]]},r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xQEDFexclTop[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];

												
(*map to interface function*)
BSFfastSigmavEffBSF["QCD-SU","Plat"]=BSFfast\[Sigma]vtopPlatS;
BSFfastSigmavEffBSF["QCD-SU","Cut"]=BSFfast\[Sigma]vtopCutS;
BSFfastSigmavEffBSF["QCD-SD","Plat"]=BSFfast\[Sigma]vbotPlatS;
BSFfastSigmavEffBSF["QCD-SD","Cut"]=BSFfast\[Sigma]vbotCutS;
BSFfastSigmavEffBSF["QCD-S","Plat"]=BSFfast\[Sigma]vnoTPlatS;
BSFfastSigmavEffBSF["QCD-S","Cut"]=BSFfast\[Sigma]vnoTCutS;
BSFfastSigmavEffBSF["dQED-S",alpha_?NumericQ]=BSFfast\[Sigma]vdqedConstS[alpha];
BSFfastSigmavEffBSF["dQED-SnoT",alpha_?NumericQ]=BSFfast\[Sigma]vdqedNoTConstS[alpha];
BSFfastSigmavEffBSF["dQCD-S",alpha_?NumericQ]=BSFfast\[Sigma]vqcdNoTConstS[alpha];

BSFfastSigmavEffBSF["QCD-FU","Plat"]=BSFfast\[Sigma]vtopPlatF;
BSFfastSigmavEffBSF["QCD-FU","Cut"]=BSFfast\[Sigma]vtopCutF;
BSFfastSigmavEffBSF["QCD-FD","Plat"]=BSFfast\[Sigma]vbotPlatF;
BSFfastSigmavEffBSF["QCD-FD","Cut"]=BSFfast\[Sigma]vbotCutF;
BSFfastSigmavEffBSF["QCD-F","Plat"]=BSFfast\[Sigma]vnoTPlatF;
BSFfastSigmavEffBSF["QCD-F","Cut"]=BSFfast\[Sigma]vnoTCutF;
BSFfastSigmavEffBSF["dQED-F",alpha_?NumericQ]=BSFfast\[Sigma]vdqedConstF[alpha];
BSFfastSigmavEffBSF["dQED-FnoT",alpha_?NumericQ]=BSFfast\[Sigma]vdqedNoTConstF[alpha];
BSFfastSigmavEffBSF["dQCD-F",alpha_?NumericQ]=BSFfast\[Sigma]vqcdNoTConstF[alpha];

BSFfastSigmavEffBSF["QED-S"]=BSFfast\[Sigma]vdqedConstS[1/128.9];
BSFfastSigmavEffBSF["QED-F"]=BSFfast\[Sigma]vqedInclTopF[1/128.9];
BSFfastSigmavEffBSF["QED-F",alpha_?NumericQ]=BSFfast\[Sigma]vqedInclTopF[alpha];
BSFfastSigmavEffBSF["QED-F_noTop"]=BSFfast\[Sigma]vqedExclTopF[1/128.9];
BSFfastSigmavEffBSF["QED-F_noTop",alpha_?NumericQ]=BSFfast\[Sigma]vqedExclTopF[alpha];



End[];
EndPackage[];
