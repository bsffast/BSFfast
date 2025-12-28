(* ::Package:: *)

(*** BSFfast ***)
(*version 2;
  author: Stefan Lederer [TUM, NCBJ];
  created: 2025-12-19*)
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
(*sbottom & stop*)
topPlat=Get["sigBSFeff_stop_asPlateau1.RD.m"];
topCut=Get["sigBSFeff_stop_asCutoff1.RD.m"];
botPlat=Get["sigBSFeff_sbottom_asPlateau1.RD.m"];
botCut=Get["sigBSFeff_sbottom_asCutoff1.RD.m"];
(*no transition*)
noTCut=Get["sigBSFeff-NoTrans_sbottom_asCutoff1.RD.m"];
noTPlat=Get["sigBSFeff-NoTrans_sbottom_asPlateau1.RD.m"];
(*constant coupling*)
xQED=Get["xgridSigBSFeff_QED_asConst0.1_m1GeV.m"];
xNoTQED=Get["xgridSigBSFeff_QEDnoTrans_asConst0.1_m1GeV.m"];
xNoTQCD=Get["xgridSigBSFeff_QCDnoTrans_asConst0.1_m1GeV.m"];
ResetDirectory[];

(* define interpolation functions *)
	(*QCD-running couplings*)
BSFfast\[Sigma]vtopPlat[m_,x_]=(Exp@*Evaluate[Interpolation[Log@topPlat,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vtopCut[m_,x_]=(Exp@*Evaluate[Interpolation[Log@topCut,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vbotPlat[m_,x_]=(Exp@*Evaluate[Interpolation[Log@botPlat,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vbotCut[m_,x_]=(Exp@*Evaluate[Interpolation[Log@botCut,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vnoTCut[m_,x_]=(Exp@*Evaluate[Interpolation[Log@noTCut,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
BSFfast\[Sigma]vnoTPlat[m_,x_]=(Exp@*Evaluate[Interpolation[Log@noTPlat,InterpolationOrder->1]]@@#&@*Log)[{m,x}];
	(*include rescaling to any coupling for asConst tables*)
BSFfast\[Sigma]vqedConst[\[Alpha]_][m_,x_]=Block[{r=\[Alpha]/0.1,a0=0.1,m0=xQED[[1,1]]},r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xQED[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqedNoTConst[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQED[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTQED[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
BSFfast\[Sigma]vqcdNoTConst[\[Alpha]_][m_,x_]=Block[{r,a0=0.1,m0=xNoTQCD[[1,1]]},r=\[Alpha]/a0;
						r^2*(m0/m)^2*(Exp@*Evaluate[Interpolation[Log@xNoTQCD[[;;,{2,3}]],InterpolationOrder->1]]@*Log)[x*r^2]];
						
(*map to interface function*)
BSFfastSigmavEffBSF["QCD-SU","Plat"]=BSFfast\[Sigma]vtopPlat;
BSFfastSigmavEffBSF["QCD-SU","Cut"]=BSFfast\[Sigma]vtopCut;
BSFfastSigmavEffBSF["QCD-SD","Plat"]=BSFfast\[Sigma]vbotPlat;
BSFfastSigmavEffBSF["QCD-SD","Cut"]=BSFfast\[Sigma]vbotCut;
BSFfastSigmavEffBSF["QCD-S","Plat"]=BSFfast\[Sigma]vnoTPlat;
BSFfastSigmavEffBSF["QCD-S","Cut"]=BSFfast\[Sigma]vnoTCut;
BSFfastSigmavEffBSF["dQED-S",alpha_?NumericQ]=BSFfast\[Sigma]vqedConst[alpha];
BSFfastSigmavEffBSF["QED-S"]=BSFfast\[Sigma]vqedConst[1/128.9];
BSFfastSigmavEffBSF["dQED-SnoT",alpha_?NumericQ]=BSFfast\[Sigma]vqedNoTConst[alpha];
BSFfastSigmavEffBSF["dQCD-S",alpha_?NumericQ]=BSFfast\[Sigma]vqcdNoTConst[alpha];



End[];
EndPackage[];
