#ifndef INPUTFILES_H
#define INPUTFILES_H

static const int nRn222Bkgs = 8; //7;
//static const int nRn222Bkgs = 8;


// TODO: break out file/name/color definitions to a new header file which
// can be included by both programs
TString Rn222BkgFiles[nRn222Bkgs] =
{
"bi214_sfoil_rot",
"bi214_sfoil",
"pb214_sfoil", // Note: no rot version
"bi210_sfoil",
"bi214_swire",
"pb214_swire",
"bi210_swire",
//,
// NOTE: re-enabled?
//"bi214_sscin", // no events, disabled
//"pb214_sscin" // no events, disabled
"bi210_sscin"
};

TString Rn222BkgNames[nRn222Bkgs] =
{
//"^{214}Bi air",
//"^{214}Pb air",
"^{214}Bi Foil Surface (Rot)",
"^{214}Bi Foil Surface",
"^{214}Pb Foil Surface",
"^{210}Bi Foil Surface",

"^{214}Bi Wire Surface",
"^{214}Pb Wire Surface",
"^{210}Bi Wire Surface",
//,
// NOTE: re-enabled?
//"^{214}Bi Scintillator Surface",
//"^{214}Pb Scintillator Surface"//,
"^{210}Bi Scintillator Surface"
//""
};

static const int nRn220Bkgs = 2;
TString Rn220BkgFiles[nRn220Bkgs] =
{
"tl208_swire",
"tl208_air"
};
TString Rn220BkgNames[nRn220Bkgs] =
{
"^{208}Tl SWire",
"^{208}Tl Air"
};

static const int nExternalBkgs = 12;
TString ExternalBkgFiles[nExternalBkgs] =
{
"bi214_feShield",
"ac228_feShield",
"tl208_feShield",

"bi214_pmt",
"ac228_pmt",
"tl208_pmt",

"k40_pmt",
"k40_scintIN",
"k40_scintOUT",
"k40_scintPET",

"pa234m_sscin",
"co60_cuTower"
};
// Note: mylar moved to interal
TString ExternalBkgNames[nExternalBkgs] =
{
"^{214}Bi Fe shield",
"^{228}Ac Fe shield",
"^{208}Tl Fe shield",

"^{214}Bi PMT",
"^{228}Ac PMT",
"^{208}Tl PMT",

"^{40}K PMT",
"^{40}K Scintillator Surface (IN)",
"^{40}K Scintillator Surface (OUT)",
"^{40}K Scintillator Surface (PET)",

"^{234m}Pa sscin",
"^{60}Co Cu Tower"
};


static const int nInternalBkgs = 12;
TString InternalBkgFiles[nInternalBkgs] =
{
"bi214_int_rot",
"pb214_int_rot",
"ac228_int_rot",
"tl208_int_rot",
"bi212_int_rot",
"bi207_int_rot",
"eu152_int_rot",
"eu154_int_rot",
"k40_int_rot",
"pa234m_int_rot",
"bi214_mylar",
"pb214_mylar"
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
"^{234m}Pa Internal",
"^{214}Bi Mylar",
"^{214}Pb Mylar"
};

static const int nNd150Samples = 1;
TString Nd150Files[nNd150Samples] =
{
//"nd150_rot_2b2n_m4"
"nd150_rot_2n2b_m4"
// changed from nd150_2n2b_rot_m4
};
TString Nd150Names[nNd150Samples] =
{
"^{150}Nd 2#nu#beta#beta"
};



static const int nNeighbours = 13;
TString NeighbourFiles[nNeighbours] =
{
"mo100_99_rot_2n2b_m14",
"mo100_99_rot_bi214",
"mo100_99_rot_pa234m",
"mo100_99_rot_k40",
"ca48_63_rot_2n2b_m4",
"ca48_63_rot_bi214",
"ca48_63_rot_pa234m",
"ca48_63_rot_y90",
"ca48_63_rot_k40",
"zr96_rot_2n2b_m4",
"zr96_rot_bi214",
"zr96_rot_pa234m",
"zr96_rot_k40"
};
TString NeighbourNames[nNeighbours] =
{
"^{100}Mo 2#nu#beta#beta",
"^{100}Mo int ^{214}Bi",
"^{100}Mo int ^{234m}Pa",
"^{100}Mo int ^{40}K",
"^{48}Ca 2#nu#beta#beta",
"^{48}Ca int ^{214}Bi",
"^{48}Ca int ^{234m}Pa",
"^{48}Ca int ^{90}Y",
"^{48}Ca int ^{40}K",
"^{96}Zr 2#nu#beta#beta",
"^{96}Zr int ^{214}Bi",
"^{96}Zr int ^{234m}Pa",
"^{96}Zr int ^{40}K"
};

#endif

