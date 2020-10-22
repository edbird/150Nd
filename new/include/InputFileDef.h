#ifndef INPUTFILEDEF_H
#define INPUTFILEDEF_H


#include "InputNumberDef.h"

// TODO: some of these use links from the 2eNg... directory to the 2e
// directory
// check if the 2e directory actually has any significant number of events
// passing cuts, if not then remove them and remove the links

// this is radon, internal and external
// 214Bi, 214Pb and 210Bi events
TString Rn222BkgFilesInput[nRn222BkgsInput] =
{
"bi214_air",
"pb214_air",

// bad - needs to be split into In/Out
"bi214_sfoil", // not sure what the difference might be between rot/non rot
                // in 2eNg case
"pb214_sfoil", // bad - needs to be split into In/Out

//0 events disabled "bi210_sfoil",
"bi210_sfoil",

"bi214_swire", // bad - needs to be split into In L0/Out L0/L>0
"pb214_swire", //, // bad - needs to be split into In L0/Out L0/L>0

//5 events disabled "bi210_swire",
"bi210_swire",

//0 events disabled "bi214_sscin",
//0 events disabled "pb214_sscin",
//0 events disabled "bi210_sscin",
"bi210_sscin",

"bi214_mylar",
"pb214_mylar"//,

//"pa234m_sscin"

//4 events disabled "bi214_air",
//0 events disabled "pb214_air"
};
// TODO: 214Bi Foil Surface missing

/*
static const int nRn222Bkgs = 9;
TString Rn222BkgFiles[nRn222Bkgs] = {"pb214_air","bi214_sfoil","pb214_sfoil","bi210_sfoil","bi210_sscin","bi214_swire","pb214_swire","bi214_mylar","pb214_mylar"};
TString Rn222BkgNames[nRn222Bkgs] = {"^{214}Pb air","^{214}Bi sfoil","^{214}Pb sfoil","^{210}Bi sfoil","^{210}Bi sscin","^{214}Bi swire","^{214}Pb swire","bi214_mylar","pb214_mylar"};
*/


TString Rn222BkgFilesNew[nRn222Bkgs] =
{
"bi214_air",
"pb214_air",

// bad - needs to be split into In/Out
"bi214_sfoil_IN",
"bi214_sfoil_OUT",
"pb214_sfoil_IN",
"pb214_sfoil_OUT",

//0 events disabled "bi210_sfoil",
"bi210_sfoil",

"bi214_swire_INL0",
"bi214_swire_OUTL0",
"bi214_swire_Lg0",
"pb214_swire_INL0",
"pb214_swire_OUTL0",
"pb214_swire_Lg0",

//5 events disabled "bi210_swire",
"bi210_swire",

//0 events disabled "bi214_sscin",
//0 events disabled "pb214_sscin",
//0 events disabled "bi210_sscin",
"bi210_sscin",

"bi214_mylar",
"pb214_mylar"//,

//"pa234m_sscin"

//4 events disabled "bi214_air",
//0 events disabled "pb214_air"
};



// radon, internal and external, Tl208 related (thoron decays)
TString Rn220BkgFiles[nRn220Bkgs] =
{
// 0 events disabled "tl208_swire", // this is removed because fit result is anticipated to be neglegiable
// but we can test this by adding it back in
// 0 events disabled "tl208_air"
"tl208_air"
};

/*static const int nExternalBkgs = 19;
TString ExternalBkgFiles[nExternalBkgs] = {"ac228_feShield","bi214_feShield","tl208_feShield","co60_cuTower","ac228_steel","bi214_steel","k40_steel","tl208_steel","co60_muMetal","bi214_pmt","k40_pmt","tl208_pmt","ac228_pmt","co60_cuPetals","bi212_steel","pa234m_steel","k40_scintIN","k40_scintOUT","k40_scintPET"};
TString ExternalBkgNames[nExternalBkgs] = {"^{228}Ac Fe shield","^{214}Bi Fe shield","^{208}Tl Fe shield","^{60}Co Cu tower","^{228}Ac steel","^{214}Bi steel","^{40}K steel","^{208}Tl steel","^{60}Co #mu metal","^{214}Bi pmt","^{40}K pmt","^{208}Tl pmt","^{228}Ac pmt","^{60}Co Cu petals","^{212}Bi steel","^{234m}Pa steel","^{40}K sscintIN","^{40}K sscintOUT","^{40}K sscintPET"};
*/

TString ExternalBkgFiles[nExternalBkgs] = 
{
"bi214_feShield",
//1 events disabled "ac228_feShield", // ok
//2 events disabled "tl208_feShield", // ok
"ac228_feShield",

"tl208_feShield",

"bi214_pmt",
//0 events disabled "ac228_pmt", // ok

"ac228_pmt",
"tl208_pmt",

//5 events disabled "k40_pmt", // ok
"k40_pmt",
//0 events disabled "k40_scintIN", // ok
//1 events disabled "k40_scintOUT",
//1 events disabled "k40_scintPET",
"k40_scintIN",
"k40_scintOUT",
"k40_scintPET",

//0 events disabled "pa234m_sscin", // ok
//1 events disabled "co60_cuTower" // ok
"pa234m_sscin",

"co60_cuTower"//,
//"neutrons"
};
// 214Bi Air is missing (no it isn't just no events expect in P2, so activity is zero)
// NOTE: mylar moved to internal

/*
static const int nExternalBkgs = 13;
TString ExternalBkgFiles[nExternalBkgs] = {"ac228_feShield","tl208_feShield","co60_cuTower","bi214_steel","tl208_steel","co60_muMetal","tl208_pmt","co60_cuPetals","bi212_steel","pa234m_steel","k40_scintIN","k40_scintOUT","k40_scintPET"};
TString ExternalBkgNames[nExternalBkgs] = {"^{228}Ac Fe shield","^{208}Tl Fe shield","^{60}Co Cu tower","^{214}Bi steel","^{208}Tl steel","^{60}Co #mu metal","^{208}Tl pmt","^{60}Co Cu petals","^{212}Bi steel","^{234m}Pa steel","^{40}K sscintIN","^{40}K sscintOUT","^{40}K sscintPET"};*/


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
"pa234m_int_rot"//,
};


// note: use symlink to fix issue
//#if MODE == 1
//TString Nd150Files[nNd150Samples] = 
//{
////"nd150_rot_2n2b_m4"
//"nd150_rot_2b2n_m4"
//};
//#elif MODE == 2
TString Nd150Files[nNd150Samples] =
{
//"nd150_rot_2b2n_m4"
"nd150_rot_2n2b_m4"
};
////TString Nd150Files[nNd150Samples] = {"nd150_rot_2n2b_m4"};
//#endif

TString Nd150Systematics[nNd150Systematics] =
{
"nd150_nominalThick",
"nd150_doubleThick",
"nd150_halfThick",
"Hbrem_2vBB",
"Lbrem_2vBB",
"LdEdX_2vBB",
"HdEdX_2vBB"
};

TString NeighbourFiles[nNeighbours] =
{
"mo100_99_rot_2n2b_m14",
"mo100_99_rot_bi214",       // 5 events disabled
"mo100_99_rot_pa234m",      // 4 evnets disabled
"mo100_99_rot_k40",         // 0 events disabled

"ca48_63_rot_2n2b_m4",
"ca48_63_rot_bi214",       // 1 events disabled
//"ca48_63_rot_k40",         // 1 events disabled
"ca48_63_rot_pa234m",      // 5 events disabled
"ca48_63_rot_y90",

"zr96_rot_2n2b_m4",
"zr96_rot_bi214",           // 14 events disabled
"zr96_rot_pa234m",          // 24 events disabled
"zr96_rot_k40"              // 12 events disabled
};
// 100Mo 40K missing


TString DataFile = "data_2e";
//TString FakeDataFile = "data_2e_fake";

#endif
