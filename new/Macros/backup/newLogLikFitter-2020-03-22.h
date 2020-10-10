#include "TH1F.h"
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
#include "TH1F.h"
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

static const int numberParams = 37;
static const int number1DHists = 1;
static const int number2DHists = 0;
TString sample_names[numberParams];

std::vector<int> fixed_params;
std::vector<int> free_params;
std::vector<TString> index_free_params;

std::vector<TString> paramNameMap[numberParams];
std::vector<double> paramActMap[numberParams];
std::vector<double> paramActErrMap[numberParams];

TObjArray *allMCSamples1D[number1DHists];
TObjArray *allMCSamples2D[number2DHists];
TObjArray *allDataSamples1D;
TObjArray *allDataSamples2D;

static const int nRn222Bkgs = 6;
TString Rn222BkgFiles[nRn222Bkgs] =
{
"bi214_sfoil",
"pb214_sfoil",
"bi214_swire",
"pb214_swire",
"bi214_sscin",
"pb214_sscin"
};

static const int nRn220Bkgs = 1;
TString Rn220BkgFiles[nRn220Bkgs] =
{
"tl208_swire"
};

static const int nExternalBkgs_oneParam = 18;
TString ExternalBkgOneParamFiles[nExternalBkgs_oneParam] =
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
"bi210_sfoil",
"bi210_swire",
"bi210_sscin",
"bi214_air",
"pb214_air",
"tl208_air",
"bi214_mylar",
"pb214_mylar",
"co60_cuTower"
};

static const int nExternalBkgs_twoParam = 0;
//TString ExternalBkgTwoParamFiles[nExternalBkgs_twoParam] = {"co60_cuTower","co60_steel","co60_muMetal","co60_cuPetals"};
//TString ExternalBkgTwoParamFiles[nExternalBkgs_twoParam] = {""};

static const int nInternalBkgs = 10;
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
"pa234m_int_rot"
};

static const int nNd150Samples = 1;
TString Nd150Files[nNd150Samples] =
{
//"nd150_2b2n_rot_m4"
"nd150_2n2b_rot_m4"
// changed from nd150_2n2b_rot_m4
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
