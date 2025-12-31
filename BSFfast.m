(* ::Package:: *)

(*** BSFfast ***)
(*version 3;
  author: Stefan Lederer [TUM, NCBJ];
  created: 2025-12-28*)
BeginPackage["BSFfast`"];

BSFfastSigmavEffBSF::usage = "BSFfast\[Sigma]vtopPlat[ModelString(,alpha)][ mass , x ]\nEffective BSF annihilation cross-section for the specified models:\n     \"QCD-_U\", \"QCD-_D\", \"QCD-_\", \"dQCD-_\", \"QED-_\", \"dQED-_\", \"dQED-_nT\".  The \"_\" must be filled by \"S\" (scalar) or \"F\" (fermion). \nQCD refers to SM-QCD and involves running couplings.\n m = constituent mass is here expected in GeV. x = m/T is inverse temperature. alpha defines the gauge coupling strength in \"dQCD\" and \"dQED\" models, and the prescription of the running coupling \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]<1GeV)=const or =0 for \"QCD\" models. (not required for \"QED-_\".)\n Call Print[Information[#,\"Usage\"]]&/@BSFfastFunctionNames for more details.";
BSFfastFunctionNames::usage = "List of BSFfast`Private` interpolation functions of all different models.";


Begin["Private`"];
BSFfastDataDir=DirectoryName[$InputFileName]<>"/BSFfast_DataM/";


BSFfastFunctionNames:={BSFfast\[Sigma]vtopPlat,BSFfast\[Sigma]vtopCut,BSFfast\[Sigma]vbotPlat,BSFfast\[Sigma]vbotCut,BSFfast\[Sigma]vnoTPlat,BSFfast\[Sigma]vnoTCut,BSFfast\[Sigma]vqedConst,BSFfast\[Sigma]vqedNoTConst,BSFfast\[Sigma]vqcdNoTConst};


BSFfast\[Sigma]vtopPlat::usage = "BSFfast\[Sigma]vtopPlat[ mass , x ]\nEffective BSF annihilation cross-section for a stop-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(+2\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong running is frozen once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu])==1. (\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) stays 1 at all lower scales.)\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vtopCut::usage  = "BSFfast\[Sigma]vtopCut[ mass , x ]\nEffective BSF annihilation cross-section for a stop-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(+2\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong coupling is set to 0 once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]) > 1.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vbotPlat::usage = "BSFfast\[Sigma]vbotPlat[ mass , x ]\nEffective BSF annihilation cross-section for a sbottom-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(-1\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong running is frozen once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu])==1. (\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) stays 1 at all lower scales.)\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vbotCut::usage  = "BSFfast\[Sigma]vbotCut[ mass , x ]\nEffective BSF annihilation cross-section for a sbottom-like scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(\(-1\)/3\)]\) multiplet) particle anti-particle pair. Only QCD contributes to long-range potentials, BS-capture and -decay. Only U(1\!\(\*SubscriptBox[\()\), \(Y\)]\) can source bound-to-bound transitions.\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong coupling is set to 0 once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]) > 1.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vnoTPlat::usage = "BSFfast\[Sigma]vnoTPlat[ mass , x ]\nEffective BSF annihilation cross-section for a colored scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(0\)]\) multiplet) particle anti-particle pair. (There are no allowed bound-to-bound transitions)\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n Strong running is frozen once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu])==1. (\!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\) stays 1 at all lower scales.)\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vnoTCut::usage  = "BSFfast\[Sigma]vnoTCut[ mass , x ]\nEffective BSF annihilation cross-section for a colored scalar (being a (3,0\!\(\*SubscriptBox[\()\), \(0\)]\) multiplet) particle anti-particle pair. (There are no allowed bound-to-bound transitions)\n \!\(\*SubscriptBox[\(\[Alpha]\), \(1\)]\)=const.\n  Strong coupling is set to 0 once \!\(\*SubscriptBox[\(\[Alpha]\), \(s\)]\)(\[Mu]) > 1.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vqedConst::usage = "BSFfast\[Sigma]vqedConst[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a U(1)-interacting particle (charge 1) including bound-to-bound transitions.\n alpha = the \!\(\*
StyleBox[\"constant\",\nFontVariations->{\"Underline\"->True}]\) gauge coupling-strength g^2/(4Pi) of the SU(3) interaction.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vqedNoTConst::usage = "BSFfast\[Sigma]vqedNoTConst[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a U(1)-interacting particle (charge 1) NEGLECTING bound-to-bound transitions.\n alpha = the \!\(\*
StyleBox[\"constant\",\nFontVariations->{\"Underline\"->True}]\) gauge coupling-strength g^2/(4Pi) of the SU(3) interaction.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";
BSFfast\[Sigma]vqcdNoTConst::usage = "BSFfast\[Sigma]vqcdNoTConst[ alpha ][ mass , x ]\nEffective BSF annihilation cross-section for a fundamental SU(3)-interacting particle anti-particle pair. (There are no allowed bound-to-bound transitions)\n alpha = the \!\(\*
StyleBox[\"constant\",\nFontVariations->{\"Underline\"->True}]\) gauge coupling-strength g^2/(4Pi) of the SU(3) interaction.\n mass = constituent mass [GeV]. \n x = mass / T = inverse temperature parameter.";


(*import of interpolation tables*)
SetDirectory[BSFfastDataDir];
	(*SCALAR*)
(*sbottom & stop*)
topPlatS=Get["QCD-SU_plateau.m"];
topCutS=Get["QCD-SU_cutoff.m"];
botPlatS=Get["QCD-SD_plateau.m"];
botCutS=Get["QCD-SD_cutoff.m"];
(*no transition*)
noTCutS=Get["QCD-S_cutoff.m"];
noTPlatS=Get["QCD-S_plateau.m"];
(*constant coupling*)
xQEDS=Get["dQED-S.m"];
xNoTQEDS=Get["dQED-SnoTr.m"];
xNoTQCDS=Get["dQCD-S.m"];

	(*FERMION*)
(*sbottom & stop*)
topPlatF=Get["QCD-FU_plateau.m"];
topCutF=Get["QCD-FU_cutoff.m"];
botPlatF=Get["QCD-FD_plateau.m"];
botCutF=Get["QCD-FD_cutoff.m"];
(*no transition*)
noTCutF=Get["QCD-F_cutoff.m"];
noTPlatF=Get["QCD-F_plateau.m"];
(*constant coupling*)
xQEDF=Get["dQED-F.m"];
xNoTQEDF=Get["dQED-FnoTr.m"];
xNoTQCDF=Get["dQCD-F.m"];
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
BSFfast\[Sigma]vqedConstS[\[Alpha]_][m_,x_]=Block[{r=\[Alpha]/0.1,a0=0.1,m0=xQEDS[[1,1]]},r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xQEDS[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqedNoTConstS[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQEDS[[1,1]]},r=\[Alpha]/a0;
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
BSFfast\[Sigma]vqedConstF[\[Alpha]_][m_,x_]=Block[{r=\[Alpha]/0.1,a0=0.1,m0=xQEDF[[1,1]]},r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xQEDF[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqedNoTConstF[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQEDF[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTQEDF[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqcdNoTConstF[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQCDF[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTQCDF[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
						
(*map to interface function*)
BSFfastSigmavEffBSF["QCD-SU","Plat"]=BSFfast\[Sigma]vtopPlatS;
BSFfastSigmavEffBSF["QCD-SU","Cut"]=BSFfast\[Sigma]vtopCutS;
BSFfastSigmavEffBSF["QCD-SD","Plat"]=BSFfast\[Sigma]vbotPlatS;
BSFfastSigmavEffBSF["QCD-SD","Cut"]=BSFfast\[Sigma]vbotCutS;
BSFfastSigmavEffBSF["QCD-S","Plat"]=BSFfast\[Sigma]vnoTPlatS;
BSFfastSigmavEffBSF["QCD-S","Cut"]=BSFfast\[Sigma]vnoTCutS;
BSFfastSigmavEffBSF["dQED-S",alpha_?NumericQ]=BSFfast\[Sigma]vqedConstS[alpha];
BSFfastSigmavEffBSF["QED-S"]=BSFfast\[Sigma]vqedConstS[1/128.9];
BSFfastSigmavEffBSF["dQED-SnoT",alpha_?NumericQ]=BSFfast\[Sigma]vqedNoTConstS[alpha];
BSFfastSigmavEffBSF["dQCD-S",alpha_?NumericQ]=BSFfast\[Sigma]vqcdNoTConstS[alpha];

BSFfastSigmavEffBSF["QCD-FU","Plat"]=BSFfast\[Sigma]vtopPlatF;
BSFfastSigmavEffBSF["QCD-FU","Cut"]=BSFfast\[Sigma]vtopCutF;
BSFfastSigmavEffBSF["QCD-FD","Plat"]=BSFfast\[Sigma]vbotPlatF;
BSFfastSigmavEffBSF["QCD-FD","Cut"]=BSFfast\[Sigma]vbotCutF;
BSFfastSigmavEffBSF["QCD-F","Plat"]=BSFfast\[Sigma]vnoTPlatF;
BSFfastSigmavEffBSF["QCD-F","Cut"]=BSFfast\[Sigma]vnoTCutF;
BSFfastSigmavEffBSF["dQED-F",alpha_?NumericQ]=BSFfast\[Sigma]vqedConstF[alpha];
BSFfastSigmavEffBSF["QED-F"]=BSFfast\[Sigma]vqedConstF[1/128.9];
BSFfastSigmavEffBSF["dQED-FnoT",alpha_?NumericQ]=BSFfast\[Sigma]vqedNoTConstF[alpha];
BSFfastSigmavEffBSF["dQCD-F",alpha_?NumericQ]=BSFfast\[Sigma]vqcdNoTConstF[alpha];



End[];
EndPackage[];
