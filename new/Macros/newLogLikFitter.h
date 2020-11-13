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


///////////////////////////////////////////////////////////////////////////////
// add random variance to 150Nd MC
//
// this is used to test the SYSALL contours for fake data
// they are much larger in real data - this may be due to variance
// NOTE: removed. issue was due to fitting data with systematics enabled
// but not using different fval in drawing code.


// used to simulate fake_data for min point of (0.7, 1.15) which is close
// to the fit point for real data
double ND150_FAKEDATA_SCALE_FACTOR = 1.0;

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
bool V_ENABLE_SYSALL = false;

const int N_SYSTEMATICS = 11;
bool V_ENABLE_SYSn[N_SYSTEMATICS] =
{
    false, // constant 1.0 MeV shift
    true,//true, // scale factor: m = 1 % + 0.2 %
    true,//true, // +- 5.55 % efficiency
    true,//true, // +- 0.50 % enrichment
    true,//true, // +- 3 keV
    true,//true, // foil thickness
    true,//true, // dE/dX
    false, // brem
    false, // foil thickness (nominal)
    false, // dE/dX (nominal)
    true//true  // brem (nominal)
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

// min point data, xi free
double min_point_data[2] = {0.0, 0.0};
double min_point_data_fval = 0.0;
double min_point_data_err[2] = {0.0, 0.0};

// min point fake xi free
double min_point_fake[2] = {0.0, 0.0};
double min_point_fake_fval = 0.0;
double min_point_fake_err[2] = {0.0, 0.0};

// min point data with systematics, xi free
double min_point_data_SYSALL[2] = {0.0, 0.0};
double min_point_data_SYSALL_fval = 0.0;
double min_point_data_SYSALL_err[2] = {0.0, 0.0};

// min point fake with systematics, xi free
double min_point_fake_SYSALL[2] = {0.0, 0.0};
double min_point_fake_SYSALL_fval = 0.0;
double min_point_fake_SYSALL_err[2] = {0.0, 0.0};

// min point data, xi=HSD, total electron energy
double min_point_data_HSD_CH0[2] = {0.0, 0.0};
double min_point_data_HSD_CH0_fval = 0.0;
double min_point_data_HSD_CH0_err[2] = {0.0, 0.0};

// min point data, xi=HSD
double min_point_data_HSD[2] = {0.0, 0.0};
double min_point_data_HSD_fval = 0.0;
double min_point_data_HSD_err[2] = {0.0, 0.0};

// min point data, xi=SSD
double min_point_data_SSD[2] = {0.0, 0.0};
double min_point_data_SSD_fval = 0.0;
double min_point_data_SSD_err[2] = {0.0, 0.0};

// +- 0.1 MeV
// +- 1.2 % scale
// +- 5.55 % efficiency
// +- 0.50 % enrichment
// +- 3 keV

double min_point_fake_sysn_h[N_SYSTEMATICS][2] =
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
double min_point_fake_sysn_h_err[N_SYSTEMATICS][2] =
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

double min_point_fake_sysn_h_fval[N_SYSTEMATICS] =
{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

double min_point_fake_sysn_l[N_SYSTEMATICS][2] =
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
double min_point_fake_sysn_l_err[N_SYSTEMATICS][2] =
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

double min_point_fake_sysn_l_fval[N_SYSTEMATICS] =
{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};


void min_point_save(
    const std::string &min_point_fname,
    double *const min_point,
    double *const min_point_err,
    const double &min_point_fval)
{
    std::ofstream ofs(min_point_fname);
    ofs << min_point[0] << " " << min_point[1] << std::endl;
    ofs << min_point_fval << std::endl;
    ofs << min_point_err[0] << " " << min_point_err[1] << std::endl;
    ofs.close();
}


int min_point_load(
    const std::string &min_point_fname,
    double *const min_point,
    double *const min_point_err,
    double &min_point_fval)
{
    int ret = 0;
    std::ifstream ifs(min_point_fname);
    if(ifs.is_open())
    {
        ifs >> min_point[0] >> min_point[1];
        ifs >> min_point_fval;
        ifs >> min_point_err[0] >> min_point_err[1];
        ifs.close();
    }
    else
    {
        ret = -1;
    }
    return ret;
}


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



void
draw_V_PHYS_init_helper(TH2D* histo, const TString &latexstring, const TString& xtitle, const TString& ytitle)
{
    //histo->SetTitle(0);
    histo->SetTitle(latexstring);
    histo->SetStats(0);
    histo->SetContour(1000);

    histo->GetXaxis()->SetTitle(xtitle);
    histo->GetYaxis()->SetTitle(ytitle);

    histo->GetZaxis()->SetLabelOffset(0.005);
    histo->GetXaxis()->SetLabelSize(17.0);
    histo->GetXaxis()->SetLabelFont(43);
    histo->GetYaxis()->SetLabelSize(17.0);
    histo->GetYaxis()->SetLabelFont(43);
    histo->GetZaxis()->SetLabelSize(17.0);
    histo->GetZaxis()->SetLabelFont(43);
    histo->GetXaxis()->SetTitleSize(18.0);
    histo->GetXaxis()->SetTitleFont(43);
    histo->GetYaxis()->SetTitleSize(18.0);
    histo->GetYaxis()->SetTitleFont(43);
    histo->GetXaxis()->SetTitleOffset(1.5);
    histo->GetYaxis()->SetTitleOffset(1.2);
    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetLabelOffset(0.01);

}


void 
draw_V_PHYS_canvas_helper(TCanvas *canvas, TH2D* histo, Color_t textcolor)
{
    canvas->cd();
    canvas->SetTicks(2, 2);
    canvas->SetRightMargin(0.15);
    canvas->SetBottomMargin(0.15);

    // TODO: try something like
    // https://root-forum.cern.ch/t/drawing-histograms-with-colz-how-to-move-the-z-axis-color-legend/41915

    canvas->Update();
    TPaletteAxis *palette = (TPaletteAxis*)histo->GetListOfFunctions()->FindObject("palette");
    if(palette != nullptr)
    {
        palette->SetX1NDC(0.88 + 0.03);
        palette->SetX2NDC(0.92 + 0.03);
        palette->SetY1NDC(0.15);
        palette->SetY2NDC(0.9);
        palette->Draw();
    }
    gPad->Modified();
    gPad->Update();
    canvas->Modified();

    TString latexstring = histo->GetTitle();
    histo->SetTitle(0);

    // latex top right label
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(43);
    latex.SetTextSize(20);
    latex.SetTextColor(textcolor);
    latex.DrawLatex(0.2, 0.8, latexstring);
}



bool DRAW_V_PHYS_ = false;

void draw_V_PHYS_(
    const TString &dir
    )
{

    const int channel = 1;

    ///////////////////////////////////////////////////////////////////////////
    // Phase 1 - initialize pointers: set to null

    // stat
    TH2D *draw_V_PHYS_STAT_P1 = nullptr;

    // sys N
    TH2D *draw_V_PHYS_SYSn_P1[N_SYSTEMATICS];
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        draw_V_PHYS_SYSn_P1[i] = nullptr;
    }

    // sys all
    TH2D *draw_V_PHYS_SYSALL_P1 = nullptr;
    
    // PHYS
    TH2D *draw_V_PHYS_P1 = nullptr;
    TH2D *draw_V_PHYSINV_P1 = nullptr;
    TH2D *draw_V_DmMPHYSINVDmM_P1 = nullptr;
    
    ///////////////////////////////////////////////////////////////////////////
    // Phase 2 - initialize pointers: set to null

    // stat
    TH2D *draw_V_PHYS_STAT_P2 = nullptr;

    // sys N
    TH2D *draw_V_PHYS_SYSn_P2[N_SYSTEMATICS]; 
    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        draw_V_PHYS_SYSn_P2[i] = nullptr;
    }

    // sys all
    TH2D *draw_V_PHYS_SYSALL_P2 = nullptr;

    // PHYS
    TH2D *draw_V_PHYS_P2 = nullptr;
    TH2D *draw_V_PHYSINV_P2 = nullptr;
    TH2D *draw_V_DmMPHYSINVDmM_P2 = nullptr;


    ///////////////////////////////////////////////////////////////////////////
    // initialize pointers: allocate memory

    ///////////////////////////////////////////////////////
    // PHASE 1
    ///////////////////////////////////////////////////////

    draw_V_PHYS_STAT_P1 =
        new TH2D("draw_V_PHYS_STAT_P1", "draw_V_PHYS_STAT_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_STAT_P1,
        "V^{STAT} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("draw_V_PHYS_SYS%d_P1", i);
        draw_V_PHYS_SYSn_P1[i] = new TH2D(name, name, 50, 0.0, 5.0, 50, 0.0, 5.0);
        TString title;
        title.Form("V^{SYS%d} (Phase 1)", i);
        draw_V_PHYS_init_helper(
            draw_V_PHYS_SYSn_P1[i], title,
            "V_{ij} bin i Energy (MeV)",
            "V_{ij} bin j Energy (MeV)");
    }

    draw_V_PHYS_SYSALL_P1
        = new TH2D("draw_V_PHYS_SYSALL_P1", "draw_V_PHYS_SYSALL_P1",
                   50, 0.0, 5.0,
                   50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_SYSALL_P1,
        "V^{SYS} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYS_P1 =
        new TH2D("draw_V_PHYS_P1", "draw_V_PHYS_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_P1,
        "V (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYSINV_P1 = 
        new TH2D("draw_V_PHYSINV_P1", "draw_V_PHYSINV_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYSINV_P1,
        "V^{-1} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");
    
    draw_V_DmMPHYSINVDmM_P1 = 
        new TH2D("draw_V_DmMPHYSINVDmM_P1", "draw_V_DmMPHYSINVDmM_P1",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_DmMPHYSINVDmM_P1,
        "V^{-1}_{RES^{2}} (Phase 1)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");
    

    ///////////////////////////////////////////////////////
    // PHASE 2
    ///////////////////////////////////////////////////////
    
    draw_V_PHYS_STAT_P2 = new TH2D(
        "draw_V_PHYS_STAT_P2", "draw_V_PHYS_STAT_P2",
        50, 0.0, 5.0,
        50, 0.0, 5.0);
    draw_V_PHYS_init_helper(draw_V_PHYS_STAT_P2,
        "V^{STAT} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("draw_V_PHYS_SYS%d_P2", i);
        draw_V_PHYS_SYSn_P2[i] =
            new TH2D(name, name, 50, 0.0, 5.0, 50, 0.0, 5.0);

        TString title;
        title.Form("V^{SYS%d} (Phase 2)", i);
        draw_V_PHYS_init_helper(
            draw_V_PHYS_SYSn_P2[i],
            title,
            "V_{ij} bin i Energy (MeV)",
            "V_{ij} bin j Energy (MeV)");
    }

    draw_V_PHYS_SYSALL_P2 =
        new TH2D("draw_V_PHYS_SYSALL_P2", "draw_V_PHYS_SYSALL_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_SYSALL_P2,
        "V^{SYS}_{ij} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYS_P2 =
        new TH2D("draw_V_PHYS_P2", "draw_V_PHYS_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYS_P2,
        "V (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_PHYSINV_P2 =
        new TH2D("draw_V_PHYSINV_P2", "draw_V_PHYSINV_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_PHYSINV_P2,
        "V^{-1} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    draw_V_DmMPHYSINVDmM_P2 =
        new TH2D("draw_V_DmMPHYSINVDmM_P2", "draw_V_DmMPHYSINVDmM_P2",
                 50, 0.0, 5.0,
                 50, 0.0, 5.0);
    draw_V_PHYS_init_helper(
        draw_V_DmMPHYSINVDmM_P2,
        "V^{-1}_{RES^{2}} (Phase 2)",
        "V_{ij} bin i Energy (MeV)",
        "V_{ij} bin j Energy (MeV)");

    //std::cout << "ALLOC " << __func__ << std::endl;
                        

    ///////////////////////////////////////////////////////////////////////////
    // set contents of TH2D objects


    double fval_fake = 0.0;

    ///////////////////////////////////////////////////////////
    // Phase 1

    int j_counter = 0;
    for(Int_t j = 0; j < 50; ++ j)
    {
        if(V_ENABLE_BIN_1D_P1[channel]->at(j) == true)
        {
            // do nothing
        }
        else
        {
            continue;
        }

        int i_counter = 0;
        for(Int_t i = 0; i < 50; ++ i)
        {
            if(V_ENABLE_BIN_1D_P1[channel]->at(i) == true)
            {
                // do nothing
            }
            else
            {
                continue;
            }

            double cstat = 0.0;
            double csys = 0.0;
            double sysn[N_SYSTEMATICS];
            double content_MATHMORE = 0.0;
            double D_minus_M_content_1 = 0.0;
            double D_minus_M_content_2 = 0.0;

            if(V_ENABLE_STAT == true)
            {
                if(i == j)
                {
                    cstat = V_PHYS_STAT_1D_P1_data[channel]->at(i + j * 50);
                }
            }
            if(V_ENABLE_SYSALL == true)
            {
                csys = V_PHYS_SYSALL_1D_P1_data[channel]->at(i + j * 50);
            }

            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                sysn[n] = V_PHYS_SYSn_1D_P1_data[n][channel]->at(i + j * 50);
            }
            
            // assumes inversion done - which it should be if det != 0.0
            D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->at(i);
            content_MATHMORE = V_PHYS_1D_P1_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);
            D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(j);

            // set

            draw_V_PHYS_STAT_P1->SetBinContent(i + 1, j + 1, cstat);

            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                draw_V_PHYS_SYSn_P1[n]->SetBinContent(i + 1, j + 1, sysn[n]);
            }

            draw_V_PHYS_SYSALL_P1->SetBinContent(i + 1, j + 1, csys);
            draw_V_PHYS_P1->SetBinContent(i + 1, j + 1, cstat + csys);
            draw_V_PHYSINV_P1->SetBinContent(i + 1, j + 1, content_MATHMORE);
            draw_V_DmMPHYSINVDmM_P1->SetBinContent(i + 1, j + 1, D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2);

            fval_fake += D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2;

            ++ i_counter;
        }

        ++ j_counter;
    }

    ///////////////////////////////////////////////////////////
    // Phase 2

    j_counter = 0;
    for(Int_t j = 0; j < 50; ++ j)
    {
        if(V_ENABLE_BIN_1D_P2[channel]->at(j) == true)
        {
            // do nothing
        }
        else
        {
            continue;
        }

        int i_counter = 0;
        for(Int_t i = 0; i < 50; ++ i)
        {
            if(V_ENABLE_BIN_1D_P2[channel]->at(i) == true)
            {
                // do nothing
            }
            else
            {
                continue;
            }

            double cstat = 0.0;
            double csys = 0.0;
            double sysn[N_SYSTEMATICS];
            double content_MATHMORE = 0.0;
            double D_minus_M_content_1 = 0.0;
            double D_minus_M_content_2 = 0.0;

            if(V_ENABLE_STAT == true)
            {
                if(i == j)
                {
                    cstat = V_PHYS_STAT_1D_P2_data[channel]->at(i + j * 50);
                }
            }
            if(V_ENABLE_SYSALL == true)
            {
                csys = V_PHYS_SYSALL_1D_P2_data[channel]->at(i + j * 50);
            }

            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                sysn[n] = V_PHYS_SYSn_1D_P2_data[n][channel]->at(i + j * 50);
                csys += sysn[n];
            }
            
            // assumes inversion done - which it should be if det != 0.0
            D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->at(i);
            content_MATHMORE = V_PHYS_1D_P2_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);
            D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->at(j);
 
            // set

            draw_V_PHYS_STAT_P2->SetBinContent(i + 1, j + 1, cstat);
            for(int n = 0; n < N_SYSTEMATICS; ++ n)
            {
                draw_V_PHYS_SYSn_P2[n]->SetBinContent(i + 1, j + 1, sysn[n]);
            }

            draw_V_PHYS_SYSALL_P2->SetBinContent(i + 1, j + 1, csys);
            draw_V_PHYS_P2->SetBinContent(i + 1, j + 1, cstat + csys);
            draw_V_PHYSINV_P2->SetBinContent(i + 1, j + 1, content_MATHMORE);
            draw_V_DmMPHYSINVDmM_P2->SetBinContent(i + 1, j + 1, D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2);

            fval_fake += D_minus_M_content_1 * content_MATHMORE * D_minus_M_content_2;

            ++ i_counter;
        }

        ++ j_counter;
    }





    ///////////////////////////////////////////////////////////////////////////
    // draw canvas

    ///////////////////////////////////////////////////////
    // PHASE 1
    ///////////////////////////////////////////////////////

    TCanvas *canvas_V_PHYS_STAT_P1 = new TCanvas("canvas_V_PHYS_STAT_P1", "canvas_V_PHYS_STAT_P1");
    draw_V_PHYS_STAT_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_STAT_P1, draw_V_PHYS_STAT_P1, kBlack);
    canvas_V_PHYS_STAT_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_STAT_P1.png");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("canvas_V_PHYS_SYS%d_P1", i);
        TString oname;
        oname = name + ".png";
        TCanvas *canvas_V_PHYS_SYSn_P1 = new TCanvas(name, name);
        draw_V_PHYS_SYSn_P1[i]->Draw("colz");
        if(i == 0 || i == 4)
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P1, draw_V_PHYS_SYSn_P1[i], kWhite);
        }
        else
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P1, draw_V_PHYS_SYSn_P1[i], kBlack);
        }
        canvas_V_PHYS_SYSn_P1->SaveAs(TString("./") + dir + "/" + oname);
    }

    TCanvas *canvas_V_PHYS_SYSALL_P1 = new TCanvas("canvas_V_PHYS_SYSALL_P1", "canvas_V_PHYS_SYSALL_P1");
    draw_V_PHYS_SYSALL_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSALL_P1, draw_V_PHYS_SYSALL_P1, kWhite);
    canvas_V_PHYS_SYSALL_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_SYSALL_P1.png");

    TCanvas *canvas_V_PHYS_P1 = new TCanvas("canvas_V_PHYS_P1", "canvas_V_PHYS_P1");
    draw_V_PHYS_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_P1, draw_V_PHYS_P1, kWhite);
    canvas_V_PHYS_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_P1.png");

    TCanvas *canvas_V_PHYSINV_P1 = new TCanvas("canvas_V_PHYSINV_P1", "canvas_V_PHYSINV_P1");
    draw_V_PHYSINV_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYSINV_P1, draw_V_PHYSINV_P1, kWhite);
    canvas_V_PHYSINV_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYSINV_P1.png");

    TCanvas *canvas_V_DmMPHYSINVDmM_P1 = new TCanvas("canvas_V_DmMPHYSINVDmM_P1", "canvas_V_DmMPHYSINVDmM_P1");
    draw_V_DmMPHYSINVDmM_P1->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_DmMPHYSINVDmM_P1, draw_V_DmMPHYSINVDmM_P1, kWhite);
    canvas_V_DmMPHYSINVDmM_P1->SaveAs(TString("./") + dir + "/" + "canvas_V_DmMPHYSINVDmM_P1.png");

    ///////////////////////////////////////////////////////
    // PHASE 2
    ///////////////////////////////////////////////////////

    TCanvas *canvas_V_PHYS_STAT_P2 = new TCanvas("canvas_V_PHYS_STAT_P2", "canvas_V_PHYS_STAT_P2");
    draw_V_PHYS_STAT_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_STAT_P2, draw_V_PHYS_STAT_P2, kBlack);
    canvas_V_PHYS_STAT_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_STAT_P2.png");

    for(int i = 0; i < N_SYSTEMATICS; ++ i)
    {
        TString name;
        name.Form("canvas_V_PHYS_SYS%d_P2", i);
        TString oname;
        oname = name + ".png";
        TCanvas *canvas_V_PHYS_SYSn_P2 = new TCanvas(name, name);
        draw_V_PHYS_SYSn_P2[i]->Draw("colz");
        if(i == 0 || i == 4)
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P2, draw_V_PHYS_SYSn_P2[i], kWhite);
        }
        else
        {
            draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSn_P2, draw_V_PHYS_SYSn_P2[i], kBlack);
        }
        canvas_V_PHYS_SYSn_P2->SaveAs(TString("./") + dir + "/" + oname);
    }

    TCanvas *canvas_V_PHYS_SYSALL_P2 = new TCanvas("canvas_V_PHYS_SYSALL_P2", "canvas_V_PHYS_SYSALL_P2");
    draw_V_PHYS_SYSALL_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_SYSALL_P2, draw_V_PHYS_SYSALL_P2, kWhite);
    canvas_V_PHYS_SYSALL_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_SYSALL_P2.png");

    TCanvas *canvas_V_PHYS_P2 = new TCanvas("canvas_V_PHYS_P2", "canvas_V_PHYS_P2");
    draw_V_PHYS_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYS_P2, draw_V_PHYS_P2, kWhite);
    canvas_V_PHYS_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYS_P2.png");

    TCanvas *canvas_V_PHYSINV_P2 = new TCanvas("canvas_V_PHYSINV_P2", "canvas_V_PHYSINV_P2");
    draw_V_PHYSINV_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_PHYSINV_P2, draw_V_PHYSINV_P2, kWhite);
    canvas_V_PHYSINV_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_PHYSINV_P2.png");

    TCanvas *canvas_V_DmMPHYSINVDmM_P2 = new TCanvas("canvas_V_DmMPHYSINVDmM_P2", "canvas_V_DmMPHYSINVDmM_P2");
    draw_V_DmMPHYSINVDmM_P2->Draw("colz");
    draw_V_PHYS_canvas_helper(canvas_V_DmMPHYSINVDmM_P2, draw_V_DmMPHYSINVDmM_P2, kWhite);
    canvas_V_DmMPHYSINVDmM_P2->SaveAs(TString("./") + dir + "/" + "canvas_V_DmMPHYSINVDmM_P2.png");

    std::cout << "printed all canvas" << std::endl;
}

#endif
