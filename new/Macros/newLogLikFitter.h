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

const Int_t MODE_PARALLEL = 1;

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

int channel_enable_1D[number1DHists] =
{
0, // ch  0 = hTotalE        (P1&2)
1, // ch  1 = hSingleEnergy  (P1&2)
0, // ch  2 = hHighEnergy    (P1&2)
0, // ch  3 = hLowEnergy     (P1&2)
0, // ch  4 = hEnergySum     (P1&2)
0  // ch  5 = hEnergyDiff    (P1&2)
};

int channel_enable_2D[number2DHists] =
{
0  // ch  0 = hHighLowEnergy    (P1&2)
};


///////////////////////////////////////////////////////////////////////////////
// number of free params
// this is used in the draw code
// must be updated manually
// it is a combination of two things
// in the gParameterGroup, params can be HARD/SOFT constrained or FREE
// they are only free parameters if this mode is FREE
// then some code blocks call Fix() on MnUserParameterState which effectively
// sets the mode to HARD
// i could find a way to do this automatically but it isn't worth the time
// at this late stage
// 2020-11-03
int gNumberFreeParams = 0;



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
// Fix all background channels (forced all background channels to be fit
// as HARD parameters)
///////////////////////////////////////////////////////////////////////////////

bool FORCE_BKG_HARD = false;



///////////////////////////////////////////////////////////////////////////////
// state machine
///////////////////////////////////////////////////////////////////////////////

bool V_ENABLE_STAT = true; // leave on
// TODO: might need to update code which enables/disables bins
// now that SYS are included

// enable/disable systematics
bool V_ENABLE_SYSALL = true;

const int N_SYSTEMATICS = 11;
bool V_ENABLE_SYSn[N_SYSTEMATICS] =
{
    false, // constant 1.0 MeV shift
    true, // scale factor: m = 1 % + 0.2 %
    true, // +- 5.55 % efficiency
    true, // +- 0.50 % enrichment
    true, // +- 3 keV
    true, // foil thickness
    true, // dE/dX
    false, // brem
    false, // foil thickness (nominal)
    false, // dE/dX (nominal)
    true  // brem (nominal)
};

//bool V_ENABLE_SYS1 = false; // constant 1.0 MeV shift
//bool V_ENABLE_SYS2 = true; // scale factor: m = 1 % + 0.2 %
//bool V_ENABLE_SYS3 = true; // +- 5.55 % efficiency
//bool V_ENABLE_SYS4 = true; // +- 0.50 % enrichment
//bool V_ENABLE_SYS5 = true; // +- 3 keV


// set to true to recalculate V_PHYS_MATHMORE
// ... and maybe some other stuff I got confused about it
bool recalculate_V_PHYS_MATHMORE = true;

// set to true after changing which systematics are enabled
// related to functions inside MinimizeFCNAxialVector::operator()
bool recalculate_V_PHYS_SYS = true;


const bool ENABLE_MIN_POINT_SYSn[N_SYSTEMATICS] =
{
    false, // +- 0.1 MeV
    true, // +- 1.2 % scale
    true, // +- 5.55 % efficiency
    false, // +- 0.50 % enrichment (overlays 3 keV)
    true, // +- 3 keV
    true, // foil thickness
    true, // energy loss (overlays bremsstrahlung)
    false, // bremsstrahlung
    false, // foil thickness (nominal)
    false, // dE/dX (nominal)
    true  // brem (nominal)
};
// enable/disable drawing of minimum points for each systematic fit
//const bool ENABLE_MIN_POINT_SYS1 = true; // +- 0.1 MeV
//const bool ENABLE_MIN_POINT_SYS2 = true; // +- 1.2 % scale
//const bool ENABLE_MIN_POINT_SYS3 = true; // +- 5.55 % efficiency
//const bool ENABLE_MIN_POINT_SYS4 = true; // +- 0.50 % enrichment
//const bool ENABLE_MIN_POINT_SYS5 = true; // +- 3 keV





// TODO: get rid of the stack
std::vector<bool> V_ENABLE_SYSALL_stack;
std::vector<bool> V_ENABLE_SYSn_stack[N_SYSTEMATICS]; //TODO anything here to init?
//std::vector<bool> V_ENABLE_SYS1_stack;
//std::vector<bool> V_ENABLE_SYS2_stack;
//std::vector<bool> V_ENABLE_SYS3_stack;
//std::vector<bool> V_ENABLE_SYS4_stack;
//std::vector<bool> V_ENABLE_SYS5_stack;

// 2020-10-04
// wanted to implement this as a debugging step but turned out to be not
// trivial so did not bother
bool DRAW_V_PHYS_MATRIX = false;

///////////////////////////////////////////////////////////////////////////////
// state machine support functions

void check_V_ENABLE_SYS_stack()
{

    std::size_t V_ENABLE_SYSALL_stack_size = V_ENABLE_SYSALL_stack.size();

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        std::size_t stack_size = V_ENABLE_SYSn_stack[i].size();
        if(V_ENABLE_SYSALL_stack_size != stack_size)
        {
            std::cout << "ERROR: " << __func__ << " i=" << i << std::endl;
            std::cout << "V_ENABLE_SYSALL_stack_size=" << V_ENABLE_SYSALL_stack_size << std::endl;
            std::cout << "stack_size=" << stack_size << std::endl;
            throw "V_ENABLE_SYSx_stack size error";
        }
    }

    /*
    if(V_ENABLE_SYSALL_stack.size() != V_ENABLE_SYS1_stack.size())
    {
        throw "V_ENABLE_SYSx_stack size error";
    }
    if(V_ENABLE_SYSALL_stack.size() != V_ENABLE_SYS2_stack.size())
    {
        throw "V_ENABLE_SYSx_stack size error";
    }
    if(V_ENABLE_SYSALL_stack.size() != V_ENABLE_SYS3_stack.size())
    {
        throw "V_ENABLE_SYSx_stack size error";
    }
    if(V_ENABLE_SYSALL_stack.size() != V_ENABLE_SYS4_stack.size())
    {
        throw "V_ENABLE_SYSx_stack size error";
    }
    if(V_ENABLE_SYSALL_stack.size() != V_ENABLE_SYS5_stack.size())
    {
        throw "V_ENABLE_SYSx_stack size error";
    }
    */

    return;
}

void V_ENABLE_SYS_stack_push()
{
    check_V_ENABLE_SYS_stack();

    V_ENABLE_SYSALL_stack.push_back(V_ENABLE_SYSALL);

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        V_ENABLE_SYSn_stack[i].push_back(V_ENABLE_SYSn[i]);
    }

    check_V_ENABLE_SYS_stack();
}

void V_ENABLE_SYS_stack_pop()
{
    check_V_ENABLE_SYS_stack();

    if(V_ENABLE_SYSALL_stack.size() > 0)
    {
        V_ENABLE_SYSALL = V_ENABLE_SYSALL_stack.back();
        V_ENABLE_SYSALL_stack.pop_back();

        for(int i = 0; i < N_SYSTEMATICS; ++ i)
        {
            V_ENABLE_SYSn[i] = V_ENABLE_SYSn_stack[i].back();
            V_ENABLE_SYSn_stack[i].pop_back();
        }
    }
    else
    {
        std::cout <<  "V_ENABLE_SYS_stack_pop ERROR" << std::endl;
        throw "V_ENABLE_SYS_stack_pop ERROR";
    }

    check_V_ENABLE_SYS_stack();
}




///////////////////////////////////////////////////////////////////////////////
// minimization points
///////////////////////////////////////////////////////////////////////////////

// systematic: SYS1
// h = high
// l = low
double min_point[2] = {0.0, 0.0}; // minimum point found, all parameter fit
double min_point_fval = 0.0;
double min_point_fake_data[2] = {0.0, 0.0};
double min_point_fake_data_fval = 0.0;

// +- 0.1 MeV
// +- 1.2 % scale
// +- 5.55 % efficiency
// +- 0.50 % enrichment
// +- 3 keV

double min_point_sysn_h[N_SYSTEMATICS][2] =
{
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0}
};

double min_point_sysn_h_fval[N_SYSTEMATICS] =
{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

double min_point_sysn_l[N_SYSTEMATICS][2] =
{
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0},
    {0.0, 0.0}
};

double min_point_sysn_l_fval[N_SYSTEMATICS] =
{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};


///////////////////////////////////////////////////////////////////////////////
// systematic objects - Phase 1
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// nominal

std::vector<double> *systematic_nominal_1D_P1[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_n_low_1D_P1[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_high_1D_P1[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_V_MATRIX_coeff_1D_P1[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};


///////////////////////////////////////////////////////////////////////////////
// systematic objects - Phase 2
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// nominal

std::vector<double> *systematic_nominal_1D_P2[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *systematic_n_low_1D_P2[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_high_1D_P2[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *systematic_n_V_MATRIX_coeff_1D_P2[N_SYSTEMATICS][number1DHists] = 
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};



///////////////////////////////////////////////////////////////////////////////
// chi2 objects for V MATRIX method, Phase 1
///////////////////////////////////////////////////////////////////////////////

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

// not used
/*
std::vector<double> *V_PHYS_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};
*/

// changes depending on M
std::vector<double> *V_PHYS_STAT_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *V_PHYS_SYSn_1D_P1_data[N_SYSTEMATICS][number1DHists] =
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};

std::vector<double> *V_PHYS_SYSALL_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *D_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *D_minus_M_1D_P1_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

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

std::vector<bool> *V_ENABLE_BIN_1D_P2[number1DHists] =
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

/*
std::vector<double> *V_PHYS_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};
*/

std::vector<double> *V_PHYS_STAT_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *V_PHYS_SYSn_1D_P2_data[N_SYSTEMATICS][number1DHists] =
{
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
    {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}
};


std::vector<double> *V_PHYS_SYSALL_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};


std::vector<double> *D_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *D_minus_M_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

std::vector<double> *M_1D_P2_data[number1DHists] = 
{
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

#endif
