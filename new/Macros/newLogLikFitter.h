#ifndef NEWLOGLIKFITTER_H
#define NEWLOGLIKFITTER_H

#include "TH1.h"
#include "TFractionFitter.h"
#include "TArray.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TPaveText.h"

#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <string>
#include "TLegend.h"
#include "TColor.h"
#include "TVector.h"
#include <vector>
#include <map>

//#include "InputFiles.h"

#include "../include/InputNumberDef.h"
#include "../include/InputFileDef.h"
#include "../include/InputNameDef.h"
#include "../include/InputColorDef.h"

const Int_t MODE_PARALLEL = 0;

// the phase, as string either "1" for Phase 1 or "2" for Phase 2
//Int_t thePhase = -1;
//TString Phase = "BOTH"; // remove these now, fit both phases simultaniously
//const Int_t gMODE_thePhase = -1; // -1=both, 0=P1, 1=P2
//const TString gMODE_Phase = "BOTH"; // BOTH/P1/P2
bool gEnablePhase1;
bool gEnablePhase2;
//std::map<int, int> map_1d_channel_to_phase;
//std::map<int, int> map_2d_channel_to_phase;

parameter_group g_pg;

//bool g_mode_fake_data = false;
//bool g_mode_fake_data = true;

// true = fit to fake data
// false = fit to real data
// fit is done using MinimizeFCNAxialVector
// the exact fit mode / choice of algorithm is set elsewhere
//bool g_mode_fit_fake_data = false;
bool g_mode_fake_data = false;
// this might change during program execution
// for example, it is necessary to set this to true to fit the fake data
// to measure systematic effects, before switching back to false to fit
// real data for the rest of the algorithm

//std::string g_datetimestamp_string;

// globals required in logLikelihood function but cannot be passed as
// parameters easily
double xi_31_baseline; // remove this? TODO
TH2D *h_nEqNull;
TH2D *h_nEqTwo;
double psiN0;
double psiN2;
double bb_Q;

// This does not change between P1 and P2, so I did
// not add a variable for each phase
//double xi_31_init_value;
//double xi_31_init_error;
//int xi_31_int_param_index;
//int xi_31_ext_param_index;

// store last xi_31 value used in logliklihood function
//double last_xi_31_parameter_value;
// store last parameter values used in logliklihood function
//int numberParams; // 45;

// Phase 1
//double paramLastValueP1Map[numberParams];
//double paramInitValueP1Map[numberParams];

// Phase 2
//double paramLastValueP2Map[numberParams];
//double paramInitValueP2Map[numberParams];
//std::vector<double> paramLastValueMap;
//std::vector<double> paramInitValueMap; // TODO: can I remove these?

//std::vector<bool> paramP1EnabledMap;
//std::vector<bool> paramP2EnabledMap;

//std::vector<file_parameter> file_params;

// Note: commented out code in fitBackgrounds.h
// where these were inialized using new and num_params
// which is the number of free params
// however here we just use numberParams
// this will work because num_free_params <= numberParams

// Phase 1
//double minuitParamCurrentP1[numberParams];
//double minuitParamInitP1[numberParams];
//double minuitParamLastP1[numberParams];

// Phase 2
//double minuitParamCurrentP2[numberParams];
//double minuitParamInitP2[numberParams];
//double minuitParamLastP2[numberParams];
//double minuitParamCurrent[numberParams];
//double minuitParamInit[numberParams];
//double minuitParamLast[numberParams]; // TODO: try and get rid of these as well


// TODO: set numberParams dynamically
// using input file
//static const int numberParams = 40; // NOTE: moved above
//int numberEnabledParams = 0;
//int minuitParamNumberCounter;
// 2 options here - either add 1 to the number of params
// and ensure final parameter is a copy of parameter 0
// adjusted for spectral shape?
// this isn't right, think about it...
//static const int numberParams_spectral = 1;

// map parameter number to Minuit parameter number (which has less parameters
// in the case of disabled parameters)
//
//std::map<int, int> paramNumberToMinuitParamNumberMap;
// skip over parameters which are not enabled

// TODO: this is now wrong because we have 2 phases
// but there is a method (a hack) to fix it
// if you want the P2 parameter instead of the P1 parameter
// then add numberEnabledParams to obtain the minuit internal index
// remember that the external index is the same for both P1 and P2
// whereas the internal index is not

// and reverse mapping
//std::map<int, int> minuitParamNumberToParamNumberMap;
// TODO: similarly this map is now wrong in the case of P2 paramter,
// but in this case you subtract numberEnabledParams from the minuit
// internal index before converting to external index

// TODO: remove vectorness of all these objects, no longer required

// appears to be a list of strings corresponding to parameter numbers
// eg, parameter 0 maps to "0" (integer to string)
// doesn't appear to have any actual use in code
// Map: parameter number/index (integer) to parameter number string (string)
// eg: 0 -> "0"
//TString sample_names[numberParams];
// TODO: remove

// maps parameter index (integer) to list of names
// fairly sure this should be changed to std::string
// NOTE: done
// Map: parameter number/index (integer) to parameter name string (string)
// eg: 1 -> "ac228_int_rot,bi212_int_rot"
// parameter name string is a unique string corresponding to each parameter
// it is constructed from joining the names of each MC sample with comma ","
//std::string paramNameMap[numberParams];


// list of each individual MC sample
// not sure if should be TString or std::string
//std::vector<TString> paramMCNameMap[numberParams];
// map a single MC file name to the corresponding parameter name,
// which consists of all included single MC names separated by comma ','
// Map: single MC sample name / file name (string) to parameter name string
// (string)
// eg: "ac228_int_rot" -> "ac228_int_rot_bi212_int_rot"
//std::map<std::string, std::string> MCNameToParamNameMap; // TODO: try to remove
// same as above but map to parameter index/number
//std::map<std::string, int> MCNameToParamNumberMap; // TODO: try to remove

// TODO: may want to remove MCNameToParamNameMap
// TODO: definitly want to clean up code in .C file which searches for
// correct MC sample in TCanvas loop
// TODO: introduce a paramEnabledMap[] (perhaps for P1 and P2)
// TODO: introduce a paramFixedP1Map and paramFixedP2Map [DONE]
// and remove the fixedParams list

// TODO: from now on, use name as described below
// maps parameter index (integer) to list of activity/error
//std::vector<double> paramActMap[numberParams];
//std::vector<double> paramActErrMap[numberParams];
//std::vector<double> paramActP2Map[numberParams];
//std::vector<double> paramActErrP2Map[numberParams];

// old names above are confusing
//
// new names below, use instead

// initial values, read from file, and initial errors, Phase 1
//double paramInitValueP1Map[numberParams];
//double paramInitErrorP1Map[numberParams];

// initial values, read from file, and initial errors, Phase 2
//double paramInitValueP2Map[numberParams];
//double paramInitErrorP2Map[numberParams];

// constraint values, and errors, phase 1
//double paramConstraintValueP1Map[numberParams];
//double paramConstraintErrorP1Map[numberParams];

// constraint values, and errors, phase 2
//double paramConstraintValueP2Map[numberParams];
//double paramConstraintErrorP2Map[numberParams];

// used to convert param name to index
// TODO: from now on, param names are given by
// single parameter: "nd150_rot_2n2b_m4"
// multiple parameter: as list "bi214_mylar,pb214_mylar"
//std::map<std::string, int> paramNameToNumberMap;

// list of parameter index for fixed parameters
//std::vector<int> fixed_params;
// now, parameter "fixed" is dependent on the phase
//std::vector<int> fixed_params_P1;
//std::vector<int> fixed_params_P2;
// and names
////std::vector<TString> fixed_param_names; // TODO unused
// list of parameter index for free parameters
//std::vector<int> free_params;
// now, parameter "free" is dependent on the phase
//std::vector<int> free_params_P1;
//std::vector<int> free_params_P2;
// and names
//std::vector<TString> free_param_names; // TODO: is this used - answer yes
// MARKER
//std::map<TString, TString> paramNameToHumanParamNameMap; // TODO: is this used?
// NOTE: required for filling histogram with names as labels
// NOTE: changed to map, now includes all parameters
// list of strings where each free parameter index is converted to string
//std::vector<std::string> index_free_params;
////std::vector<TString> index_free_params; // TODO: is this used?
// NOTE: removed and replaced with 2x Maps below

//std::vector<std::string> paramMCList[numberParams];
//std::map<std::string, std::string> paramMCNameToHumanMCNameMap;
// I think I need two additional maps
// one to convert between MC Sample Name and human readable MC Sample Name
// another to convert between Parameter Name and human readable Parameter Name
//std::map<TString, TString> MCSampleNameToHumanReadableMCSampleNameMap;
//std::map<TString, TString> paramNameToHumanReadableParamNameMap;

// new map to avoid conflicts, consider removing above 2 maps
//std::map<int, TString> paramNumberToHumanReadableParamNameMap;


//int paramConstrainModeP1Map[numberParams];
//int paramConstrainModeP2Map[numberParams];

/*
const int MODE_PARAM_LOCK_UNDEFINED = -1;
const int MODE_PARAM_LOCKED = 1; // use only P1 values for everything
const int MODE_PARAM_UNLOCKED = 0; // user P1 and P2 values separatly

int paramLockedModeMap[numberParams];
*/

// list of parameter index for enabled parameters
// enabled parameters are drawn in the final output
//std::vector<int> enabled_params;
// and names
////std::vector<TString> enabled_param_names; // TODO: is this used?
// list of parameter index for disabled parameters
//std::vector<int> disabled_params;
// and names
////std::vector<TString> disabled_param_names; // TODO unused
// list of strings where each enabled parameter index is converted to string
////std::vector<TString> index_enabled_params; // TODO: is this used?

//const int numberPhases = 2;
//const int number1DHists_perphase = 6;
//const int number2DHists_perphase = 1;
const int number1DHists = 6; //number1DHists_perphase * numberPhases;
const int number2DHists = 1; //number2DHists_perphase * numberPhases;


///////////////////////////////////////////////////////////////////////////////
// channel enable flags

const double channel_enable_1D[number1DHists] =
{
0, // ch  0 = hTotalE        (P1&2)
1, // ch  1 = hSingleEnergy  (P1&2)
0, // ch  2 = hHighEnergy    (P1&2)
0, // ch  3 = hLowEnergy     (P1&2)
0, // ch  4 = hEnergySum     (P1&2)
0  // ch  5 = hEnergyDiff    (P1&2)
};

const double channel_enable_2D[number2DHists] =
{
0  // ch  0 = hHighLowEnergy    (P1&2)
};


///////////////////////////////////////////////////////////////////////////////
// histogram names

const std::string channel_histname_1D[number1DHists] =
{
    "hTotalE_",
    "hSingleEnergy_",
    "hHighEnergy_",
    "hLowEnergy_",
    "hEnergySum_",
    "hEnergyDiff_"
};

const std::string channel_histname_2D[number2DHists] =
{
    "hHighLowEnergy_"
};


///////////////////////////////////////////////////////////////////////////////
// draw flags

const bool channel_enable_draw_1D[number1DHists] =
{
    true,
    true,
    false,
    false,
    false,
    false
};

const bool channel_enable_draw_2D[number2DHists] =
{
    false
};

// global variable to hold chisquare result
// there is probably a better way to do this using Eval
// but I am too lazy to implement it properly
//Double_t global_chisquare; // TODO: removed, no longer used? 

TObjArray *allMCSamples1D[number1DHists];
TObjArray *allMCSamples2D[number2DHists];
TObjArray *allDataSamples1D;
TObjArray *allDataSamples2D;
TObjArray *allFakeDataSamples1D;
TObjArray *allFakeDataSamples2D;

// stop printing unneccessary debug statements
bool rebuild_fake_data_first_run = true;

// select between data / fake data in transparent way
//TObjArray *allDataOrFakeDataSamples1D;
//TObjArray *allDataOrFakeDataSamples2D;


std::vector<std::pair<double,double>> ll_walk;
std::vector<std::pair<double,double>> ll_walk_save;


///////////////////////////////////////////////////////////////////////////////
// minimization points
///////////////////////////////////////////////////////////////////////////////

// systematic: SYS1
// h = high
// l = low
double min_point[2] = {0.0, 0.0}; // minimum point found, all parameter fit
double min_point_sys1_h[2] = {0.0, 0.0};
double min_point_sys1_l[2] = {0.0, 0.0};

///////////////////////////////////////////////////////////////////////////////
// systematic objects - Phase 1
///////////////////////////////////////////////////////////////////////////////

std::vector<double> *systematic_offset_nominal_1D_P1[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_offset_low_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_offset_high_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_offset_V_MATRIX_coeff_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};


///////////////////////////////////////////////////////////////////////////////
// systematic objects - Phase 2
///////////////////////////////////////////////////////////////////////////////

std::vector<double> *systematic_offset_nominal_1D_P2[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_offset_low_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_offset_high_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_offset_V_MATRIX_coeff_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};


///////////////////////////////////////////////////////////////////////////////
// chi2 objects for V MATRIX method, Phase 1
///////////////////////////////////////////////////////////////////////////////


const bool V_ENABLE_SYS1 = false;
const bool V_ENABLE_STAT = true; // leave on

TMatrixD *V_PHYS_1D_P1_MATHMORE[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<bool> *V_ENABLE_BIN_1D_P1[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*
TH2D *V_PHYS_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *V_PHYS_STAT_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_STAT_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *V_PHYS_SYS1_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_SYS1_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *V_PHYS_SYS2_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_SYS2_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};


/*TH2D *D_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *D_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *D_minus_M_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *D_minus_M_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *M_1D_P1[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *M_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

///////////////////////////////////////////////////////////////////////////////
// chi2 objects for V MATRIX method, Phase 2
///////////////////////////////////////////////////////////////////////////////

TMatrixD *V_PHYS_1D_P2_MATHMORE[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

bool recalculate_V_PHYS_xD_Px_MATHMORE = true;

std::vector<bool> *V_ENABLE_BIN_1D_P2[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *V_PHYS_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *V_PHYS_STAT_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_STAT_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *V_PHYS_SYS1_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_SYS1_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *V_PHYS_SYS2_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *V_PHYS_SYS2_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};


/*TH2D *D_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *D_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *D_minus_M_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *D_minus_M_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*TH2D *M_1D_P2[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};*/
std::vector<double> *M_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

///////////////////////////////////////////////////////////////////////////////
// chi2 objects for V MATRIX method, Phase 1
///////////////////////////////////////////////////////////////////////////////

/*
TH2D *V_PHYS_2D_P1[number2DHists] = 
{
    nullptr
};

TH2D *V_PHYS_STAT_2D_P1[number2DHists] = 
{
    nullptr
};

TH2D *V_PHYS_SYS1_2D_P1[number2DHists] = 
{
    nullptr
};

TH2D *V_PHYS_SYS2_2D_P1[number2DHists] = 
{
    nullptr
};


TH2D *D_2D_P1[number2DHists] = 
{
    nullptr
};

TH2D *D_minus_M_2D_P1[number2DHists] = 
{
    nullptr
};

TH2D *M_2D_P1[number2DHists] = 
{
    nullptr
};
*/

///////////////////////////////////////////////////////////////////////////////
// chi2 objects for V MATRIX method, Phase 2
///////////////////////////////////////////////////////////////////////////////

/*
TH2D *V_PHYS_2D_P2[number2DHists] = 
{
    nullptr
};

TH2D *V_PHYS_STAT_2D_P2[number2DHists] = 
{
    nullptr
};

TH2D *V_PHYS_SYS1_2D_P2[number2DHists] = 
{
    nullptr
};

TH2D *V_PHYS_SYS2_2D_P2[number2DHists] = 
{
    nullptr
};


TH2D *D_2D_P2[number2DHists] = 
{
    nullptr
};

TH2D *D_minus_M_2D_P2[number2DHists] = 
{
    nullptr
};

TH2D *M_2D_P2[number2DHists] = 
{
    nullptr
};
*/

// chisquare objects
/*
TH2D *V_CHEN = nullptr; // V channel enable (bin enable)
TH2D *V_PHYS = nullptr; // V physics
TH2D *V_PHYS_STAT = nullptr; // V physics, statistical
TH2D *V_PHYS_SYS = nullptr; // V physics, systematic
TH2D *D = nullptr; // data
TH2D *M = nullptr; // MC (model)
TH2D *D_minus_M = nullptr; // D - M
TH2D *V_SUPER = nullptr; // containing all information from V matricies
*/


#if 0

static const int nRn222Bkgs = 6;
//static const int nRn222Bkgs = 8;


// TODO: break out file/name/color definitions to a new header file which
// can be included by both programs
TString Rn222BkgFiles[nRn222Bkgs] =
{
"bi214_sfoil_rot",
"pb214_sfoil", // Note: no rot version
"bi210_sfoil",
"bi214_swire",
"pb214_swire",
"bi210_swire"
//,
// NOTE: re-enabled?
//"bi214_sscin", // no events, disabled
//"pb214_sscin" // no events, disabled
};

TString Rn222BkgNames[nRn222Bkgs] =
{
//"^{214}Bi air",
//"^{214}Pb air",
"^{214}Bi Foil Surface",
"^{214}Pb Foil Surface",
"^{210}Bi Foil Surface",

"^{214}Bi Wire Surface",
"^{214}Pb Wire Surface",
"^{210}Bi Wire Surface"//,
//,
// NOTE: re-enabled?
//"^{214}Bi Scintillator Surface",
//"^{214}Pb Scintillator Surface"//,
//"^{210}Bi Scintillator Surface"
//""
};

static const int nRn220Bkgs = 1;
TString Rn220BkgFiles[nRn220Bkgs] =
{
"tl208_swire"
};
TString Rn220BkgNames[nRn220Bkgs] =
{
"^{208}Tl SWire"
};

static const int nExternalBkgs = 11;
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
"pa234m_sscin",
//"bi210_sfoil", // Note: no rot version // moved
//"bi210_swire", // moved
//"bi210_sscin", // no events, disabled
"bi214_air",
//"pb214_air", // disabled, no events in P2
//"tl208_air",
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
"^{40}K sscintIN",
"^{234m}Pa sscin",
//"^{210}Bi Foil Surface",
//"^{210}Bi Wire Surface",
//"^{210}Bi Scintillator Surface",
"^{214}Bi Air",
//"^{214}Pb Air",
//"^{208}Tl Air",
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



static const int nNeighbours = 12;
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
"^{96}Zr 2#nu#beta#beta",
"^{96}Zr int ^{214}Bi",
"^{96}Zr int ^{234m}Pa",
"^{96}Zr int ^{40}K"
};

#endif

#endif
