#ifndef INPUTNAMEDEF_H
#define INPUTNAMEDEF_H


#include "InputNumberDef.h"



/*
TString Rn222BkgNames[nRn222Bkgs] =
{
"^{214}Bi Air",
"^{214}Pb Air",

"^{214}Bi Foil Surface",
"^{214}Pb Foil Surface",

"^{210}Bi Foil Surface",

"^{214}Bi Wire Surface",
"^{214}Pb Wire Surface",

"^{210}Bi Wire Surface",

//"^{214}Bi Scintillator Surface",
//"^{214}Pb Scintillator Surface",
"^{210}Bi Scintillator Surface",

"^{214}Bi Mylar",
"^{214}Pb Mylar"//,

//"^{234m}Pa Scintillator Surface"

};
*/

TString Rn222BkgNames[nRn222Bkgs] =
{
"^{214}Bi Air",
"^{214}Pb Air",

"^{214}Bi Foil Surface IN",
"^{214}Bi Foil Surface OUT",
"^{214}Pb Foil Surface IN",
"^{214}Pb Foil Surface OUT",

"^{210}Bi Foil Surface",

"^{214}Bi Wire Surface L0 IN",
"^{214}Bi Wire Surface L0 OUT",
"^{214}Bi Wire Surface L>0",
"^{214}Pb Wire Surface L0 IN",
"^{214}Pb Wire Surface L0 OUT",
"^{214}Pb Wire Surface L>0",

"^{210}Bi Wire Surface",

//"^{214}Bi Scintillator Surface",
//"^{214}Pb Scintillator Surface",
"^{210}Bi Scintillator Surface",

"^{214}Bi Mylar",
"^{214}Pb Mylar"//,

//"^{234m}Pa Scintillator Surface"

};



TString Rn220BkgNames[nRn220Bkgs] = 
{
"^{208}Tl Air"
};


TString ExternalBkgNames[nExternalBkgs] =
{
"^{214}Bi Fe Shield",

"^{228}Ac Fe shield",
"^{208}Tl Fe Shield",

"^{214}Bi PMT",

"^{228}Ac PMT",
"^{208}Tl PMT",

"^{40}K PMT",

"^{40}K Scintillator IN",
"^{40}K Scintillator OUT",
"^{40}K Scintillator PET",

//"^{234m}Pa sscin",
"^{234m}Pa Scintillator Surface",
"^{60}Co Cu Tower",

//"Neutrons"
};



TString InternalBkgNames[nInternalBkgs] =
{
"^{214}Bi Internal",
"^{214}Pb Internal",
"^{228}Ac Internal",
"^{208}Tl Internal",

"^{212}Bi Internal",
"^{207}Bi Internal",
"^{152}Eu Internal",
"^{154}Eu Internal",

"^{40}K Internal",
"^{234m}Pa Internal"
};



TString Nd150Names[nNd150Samples] =
{
"^{150}Nd 2#nu#beta#beta"
};


TString Nd150SystematicNames[nNd150Systematics] =
{
"^{150}Nd 2#nu#beta#beta nominal thick",
"^{150}Nd 2#nu#beta#beta double thick",
"^{150}Nd 2#nu#beta#beta half thick",
"^{150}Nd 2#nu#beta#beta Hbrem_2vBB",
"^{150}Nd 2#nu#beta#beta Lbrem_2vBB",
"^{150}Nd 2#nu#beta#beta LdEdX_2vBB",
"^{150}Nd 2#nu#beta#beta HdEdX_2vBB"
};



TString NeighbourNames[nNeighbours] =
{
"^{100}Mo 2#nu#beta#beta",
"^{100}Mo int ^{214}Bi",
"^{100}Mo int ^{234m}Pa",
"^{100}Mo int ^{40}K",

"^{48}Ca 2#nu#beta#beta",
"^{48}Ca int ^{214}Bi",
//"^{48}Ca int ^{40}K",
"^{48}Ca int ^{234m}Pa",
"^{48}Ca int ^{90}Y",

"^{96}Zr 2#nu#beta#beta",
"^{96}Zr int ^{214}Bi",
"^{96}Zr int ^{234m}Pa",
"^{96}Zr int ^{40}K"
};


TString DataName = "data";

#endif
