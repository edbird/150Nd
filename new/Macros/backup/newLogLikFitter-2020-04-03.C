#include "TH1F.h"
#include "TArray.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"

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
#include <vector>
#include <TMinuit.h>

#include "newLogLikFitter.h"

//-----------------------------------------
//    Functions
//----------------------------------------
void loadFiles();

void  book1DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);
void  book2DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);
void logLikelihood(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */); 
Double_t getNumberMC1D(Int_t channel, Int_t binx, Double_t *p);
Double_t getNumberMC2D(Int_t channel, Int_t binx, Int_t biny, Double_t *p);

void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, Int_t);


// helper functions to set values read from file
//


// this anonymous function
// (I couldn't think of a name)
// copies the value from input to output if
// string_input == string_compare
// if string_input != string_compare,
// then if string_input == "none", output is
// set to zero
// else
// an attempt to set output using std::stod
// is made, if this fails, an exception is
// caught and an error is printed
void function_A(const std::string& string_compare,
                const double input,
                const std::string& string_input,
                double& output)
{
    if(string_input == string_compare)
    {
        output = input;
    }
    else if(string_input == "none")
    {
        output = 0.;
    }
    else
    {
        try
        {
            output = std::stod(string_input);
        }
        catch(std::exception& e)
        {
            std::cout << e.what() << std::endl;
            std::cout << "string_input=" << string_input << std::endl;
        }
    }
}

// function for setting either the Value or Error
// of the constraint parameter for Phase 1, using
// initial value data from Phase 1
//
// only for Phase 1, because "useinit" is used to
// obtain value copied from Phase 1 initial value
//void function_paramConstrainValueErrorP1(const double paramInitValueP1,
//                                         const std::string& paramConstrainValueP1_str,
//                                         double& paramConstrainValueP1);
//void function_paramConstrainErrorP1();

// function for setting either the Value or Error
// of the constraint parameter for Phase 2, using
// constraint data from Phase 1
//
// only for Phase 2, because "same" is used to
// obtain value copied from Phase 1
//void function_paramConstrainValueErrorP2(paramConstrainValueErrorP1,
//                                         paramConstrainValueErrorP2_str,
//                                         paramConstrainValueErrorP2);

// TODO: change to generic name, both functions the same
//void function_paramConstrainErrorP2(paramConstrainErrorP1,
//                                    paramConstrainErrorP2_str,
//                                    paramConstrainErrorP2);


void function_paramConstrainMode(const std::string& paramConstrainModePX_str, int &paramConstrainModePX)
{

    // Phase 1 / Phase 2 Constran Mode
    if(paramConstrainModePX_str == "free")
    {
        paramConstrainModePX = MODE_PARAM_FREE;
    }
    else if(paramConstrainModePX_str == "soft")
    {
        paramConstrainModePX = MODE_PARAM_SOFT;
    }
    else if(paramConstrainModePX_str == "hard")
    {
        paramConstrainModePX = MODE_PARAM_HARD;
    }
    else
    {
        std::cout << "ERROR: Unrecognized value: paramConstrainModePX_str=" << paramConstrainModePX_str << std::endl;
        paramConstrainModePX = MODE_PARAM_UNDEFINED;
    }
}


void loadFiles() {

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    allDataSamples1D = new TObjArray();
    allDataSamples2D = new TObjArray();



    // read parameter_name.lst file
    for(int i = 0; i < numberParams; i++)
    {
        // paramNameMap[i] = std::vector<TString>();
        //paramActMap[i] = std::vector<double>();
        //paramActErrMap[i] = std::vector<double>();

        /*paramNameMap[i].clear();*/
        //paramActMap[i].clear();
        //paramActErrMap[i].clear();
        //paramActP2Map[i].clear();
        //paramActErrP2Map[i].clear();
        // old names above, new below
        /*paramInitValueP1Map[i].clear();
        paramInitErrorP1Map[i].clear();
        paramInitValueP2Map[i].clear();
        paramInitErrorP2Map[i].clear();
        paramConstraintValueP1Map[i].clear();
        paramConstraintErrorP1Map[i].clear();
        paramConstraintValueP2Map[i].clear();
        paramConstraintErrorP2Map[i].clear();*/

        MCNameToParamNameMap.clear();
        MCNameToParamNumberMap.clear();

        paramNameToNumberMap.clear();

        fixed_params.clear();
        free_params.clear();
        index_free_params.clear();

        enabled_params.clear();
        disabled_params.clear();
        index_enabled_params.clear();

        // TODO: this is done by calling "new"
        /*
        for(int i = 0; i < number1DHists; ++ i)
        {
            allMCSamples1D[i]->Clear();
        }
        for(int i = 0; i < number2DHists; ++ i)
        {
            allMCSamples2D[i]->Clear();
        }
        allDataSamples1D->Clear();
        allDataSamples2D->Clear();
        */
    }

    std::ifstream paramFile;
    paramFile.open("parameter_names.lst");
    // TODO: note that params file contains "same" in many backgrounds
    // which clearly should not be the same in P1 and P2
    // is tl208 sfoil and sscin missing?
    // NOTE: this data is not used in the fit unless constraints are
    // introduced
    // TODO: understand fully how this works, the parameter numbers
    // may not be arbitrary and the _P1 / _P2 strings may be necessary
    // whereas I removed them
    while(!paramFile.eof())
    {
        
        //std::string paramName_str;
        //int paramName;

        std::string paramName;

        std::string paramNumber_str;
        int paramNumber;
        
        std::string paramEnabled_str;
        
        //double activity;
        //double activity_error;
        //double activity_P2;
        //double activity_error_P2;
        // old naming convention above, new below
        //double activity_init_value_P1

        double paramInitValueP1 = 0.;
        double paramInitErrorP1 = 0.;
        double paramInitValueP2 = 0.;
        double paramInitErrorP2 = 0.;

        std::string paramInitValueP1_str;
        std::string paramInitErrorP1_str;
        std::string paramInitValueP2_str;
        std::string paramInitErrorP2_str;

        double paramConstraintValueP1 = 0.;
        double paramConstraintErrorP1 = 0.;
        double paramConstraintValueP2 = 0.;
        double paramConstraintErrorP2 = 0.;

        std::string paramConstraintValueP1_str;
        std::string paramConstraintErrorP1_str;
        std::string paramConstraintValueP2_str;
        std::string paramConstraintErrorP2_str;

        int paramConstrainModeP1 = MODE_PARAM_UNDEFINED;
        int paramConstrainModeP2 = MODE_PARAM_UNDEFINED;

        std::string paramConstrainModeP1_str;
        std::string paramConstrainModeP2_str;

        //std::string activity_P2_str;
        //std::string activity_error_P2_str;
        //TString fixed_str;
        //TString enabled_str;

        std::stringstream ss;
        std::string s;
        std::getline(paramFile, s);

        /*
        if(s.size() > 0)
        {
            if(s[0] == '#')
            {
                continue;
            }
        }
        else
        {
            continue;
        }
        */
        
        // ignore blank line
        if((!(s.size() > 0)) || s[0] == '#')
        {
            continue;
        }
        //if(s.size() > 0)
        else
        {

            ss << s;

            //ss >> paramNumber >> name
            //   >> activity >> activity_error
            //   >> activity_P2_str >> activity_error_P2_str
            //   >> fixed_str >> enabled_str;
            /*
            ss >> paramNumber >> paramName
               >> paramInitValueP1_str >> paramInitErrorP1_str
               >> paramInitValueP2_str >> paramInitErrorP2_str
               >> paramConstraintValueP1_str >> paramConstraintErrorP1_str
               >> paramConstraintValueP2_str >> paramConstraintErrorP2_str
               >> fixed_str >> enabled_str;
            */
            // read in fixed width entries
            ss >> paramNumber_str
               >> paramEnabled_str
               >> paramInitValueP1_str >> paramInitErrorP1_str
               >> paramInitValueP2_str >> paramInitErrorP2_str
               >> paramConstraintValueP1_str >> paramConstraintErrorP1_str
               >> paramConstraintValueP2_str >> paramConstraintErrorP2_str
               >> paramConstrainModeP1_str >> paramConstrainModeP2_str;

            //std::cout << "read paramNumber_str=" << paramNumber_str << std::endl;

            // read in non fixed width entries
            // marked by END
            // construct the parameter name
            std::vector<std::string> unique_mc_names;
            for(;;)
            {
                std::string nextstring;
                ss >> nextstring;
                //std::cout << "nextstring=" << nextstring << std::endl;
                if(nextstring == "END")
                {
                    break;
                }
                else
                {
                    // construct name
                    if(paramName.size() > 0)
                    {
                        paramName += ",";
                    }
                    paramName += nextstring;
                    //std::cout << "paramName=" << paramName << std::endl;

                    // add MC to list of MC
                    //MCNameToParamNameMap[paramNumber].push_back(nextstring);
                    unique_mc_names.push_back(nextstring);
                }
            }

            // setup the MCNameToParamNameMap
            for(int i = 0; i < unique_mc_names.size(); ++ i)
            {
                MCNameToParamNameMap[unique_mc_names.at(i)] = paramName;
                MCNameToParamNumberMap[unique_mc_names.at(i)] = paramNumber;
            }


            //std::cout << paramNumber << "\t" << name
            //          << "\t" << activity << "\t" << activity_error
            //          << "\t" << activity_P2_str << "\t" << activity_error_P2_str
            //          << " fixed_str=" << fixed_str << " enabled_str=" << enabled_str << std::endl;
            //std::cout << paramNumber << paramName
            //          << paramInitValueP1_str << paramInitErrorP1_str
            //          << paramInitValueP2_str << paramInitErrorP2_str
            //          << paramConstraintValueP1_str << paramConstraintErrorP1_str
            //          << paramConstraintValueP2_str << paramConstraintErrorP2_str
            //          << fixed_str << enabled_str;
            // moved below

            // phase 1 initial values
            // read value from file, no other option
            //std::cout << "paramInitValueP1_str=" << paramInitValueP1_str << "END" << std::endl;
            //paramInitValueP1 = std::stod("3.535e-04");
            //std::cout << "worked!" << std::endl;
            //std::cout << paramInitValueP1_str << std::endl;
            //std::cout << paramInitValueP2_str << std::endl;
            paramInitValueP1 = std::stod(paramInitValueP1_str);
            paramInitErrorP1 = std::stod(paramInitErrorP1_str);
            //paramInitValueP1 = atof(paramInitValueP1_str.c_str());
            //paramInitValueP2 = atof(paramInitValueP2_str.c_str());

            // phase 2 intial values
            // check if string is "same", if so use values from phase 1
            // else read value from file
            function_A("same", paramInitValueP1, paramInitValueP2_str, paramInitValueP2);
            function_A("same", paramInitErrorP1, paramInitErrorP2_str, paramInitErrorP2);
            /*
            if(paramInitValueP2_str == "same")
            {
                paramInitValueP2 = paramInitValueP1;
            }
            else
            {
                //std::cout << paramInitValueP2_str << std::endl;
                paramInitValueP2 = std::stod(paramInitValueP2_str);
                //paramInitValueP2 = atof(paramInitValueP2_str.c_str());
            }
            if(paramInitErrorP2_str == "same")
            {
                paramInitErrorP2 = paramInitErrorP1;
            }
            else
            {
                //std::cout << paramInitErrorP2_str << std::endl;
                paramInitErrorP2 = std::stod(paramInitErrorP2_str);
                //paramInitErrorP2 = atof(paramInitErrorP2_str.c_str());
            }
            */

            // phase 1 constraints
            // check if string is "useinit", if so use values from phase 1
            // as constraints for phase 1
            // else read value from file
            function_A("useinit", paramInitValueP1, paramConstraintValueP1_str, paramConstraintValueP1);
            //function_paramConstrainValueErrorPX(paramInitValueP1, paramConstrainValueP1_str, paramConstrainValueP1);
            /*
            if(paramConstraintValueP1_str == "useinit")
            {
                paramConstraintValueP1 = paramInitValueP1;
            }
            else if(paramConstraintValueP1_str == "none")
            {
                paramConstraintValueP1 = 0.;
            }
            else
            {
                //std::cout << paramConstraintValueP1_str << std::endl;
                paramConstraintValueP1 = std::stod(paramConstraintValueP1_str);
                //paramConstraintValueP1 = atof(paramConstraintValueP1_str.c_str());
            }
            */

            function_A("useinit", paramInitErrorP1, paramConstraintErrorP1_str, paramConstraintErrorP1);
            //function_paramConstrainValueErrorPX(paramInitValueP2, paramConstrainValueP2_str, paramConstrainValueP2);
            /*
            if(paramConstraintErrorP1_str == "useinit")
            {
                paramConstraintErrorP1 = paramInitErrorP1;
            }
            else if(paramConstraintErrorP1_str == "none")
            {
                paramConstraintValueP1 = 0.;
            }
            else
            {
                //std::cout << paramConstraintErrorP1_str << std::endl;
                paramConstraintErrorP1 = std::stod(paramConstraintErrorP1_str);
                //paramConstraintErrorP1 = atof(paramConstraintErrorP1_str.c_str());
            }
            */


            // phase 2 initial constraints
            // check if string is "same", if so use constraints from phase 1
            // else read value from file
            //if(activity_P2_str == "same")
            //useinit should be same
            function_A("same", paramInitValueP2, paramConstraintValueP2_str, paramConstraintValueP2);
            function_A("same", paramInitErrorP2, paramConstraintErrorP2_str, paramConstraintErrorP2);
            //function_paramConstrainValueErrorP2(paramConstrainValueP1,
            //                                    paramConstrainValueP2_str,
            //                                    paramConstrainValueP2);

            //function_paramConstrainValueErrorP2(paramContrainErrorP1,
            //                                    paramConstrainErrorP2_str,
            //                                    paramConstrainErrorP2);
            /*
            if(paramConstraintValueP2_str == "same")
            {
                //activity_P2 = activity;
                paramConstraintValueP2 = paramConstraintValueP1;
            }
            else if(paramConstraintValueP2_str == "none")
            {
                paramConstraintValueP2 = 0.;
            }
            else
            {
                //std::cout << paramConstraintValueP2_str << std::endl;
                paramConstraintValueP2 = std::stod(paramConstraintValueP2_str);
                //paramConstraintValueP2 = atof(paramConstraintValueP2_str.c_str());
            }
            */
            //if(activity_error_P2_str == "same")
            /*
            if(paramConstraintErrorP2_str == "same")
            {
                //activity_error_P2 = activity_error;
                paramConstraintErrorP2 = paramConstraintErrorP1;
            }
            else if(paramConstraintErrorP2_str == "none")
            {
                paramConstraintErrorP2 = 0.;
            }
            else
            {
                //std::cout << paramConstraintErrorP2_str << std::endl;
                paramConstraintErrorP2 = std::stod(paramConstraintErrorP2_str);
                //paramConstraintErrorP2 = atof(paramConstraintErrorP2_str.c_str());
            }
            */


            // now check if initial value strs are "useconstraint"
            // and copy values from constraints
            // also need to check for "same" again, for P2 init value and error
        
            paramNumber = std::stoi(paramNumber_str);

            // constrain mode
            // Phase 1
            /*
            if(fixed_str.CompareTo("fixed") == 0)
            {
                //std::cout << "detected fixed parameter" << std::endl;
                fixed_params.push_back(paramNumber);
            }
            else if(fixed_str.CompareTo("notfixed") == 0 || fixed_str.CompareTo("free") == 0)
            {
                //std::cout << "detected free parameter" << std::endl;
                // TODO: move code which picks up free params to this block
            }*/
            function_paramConstrainMode(paramConstrainModeP1_str, paramConstrainModeP1);
            /*
            if(paramConstrainModeP1_str == "free")
            {
                paramConstrainModeP1 = 0;
            }
            else if(paramConstrainModeP1_str == "soft")
            {
                paramConstrainModeP1 = 1;
            }
            else if(paramConstrainModeP1_str == "hard")
            {
                paramConstrainModeP1 = 2;
            }
            else
            {
                std::cout << "ERROR: Unrecognized value: paramConstrainModeP1_str=" << paramConstrainModeP1_str << std::endl;
                //std::cout << "unknown fixed/free parameter specification" << std::endl;
                //std::cout << "fixed_str=" << fixed_str << std::endl;
            }
            */

            // Phase 2
            function_paramConstrainMode(paramConstrainModeP2_str, paramConstrainModeP2);
            /*
            if(paramConstrainModeP2_str == "free")
            {
                paramConstrainModeP2 = MODE_PARAM_FREE;
            }
            else if(paramConstrainModeP2_str == "soft")
            {
                paramConstrainModeP2 = MODE_PARAM_SOFT;
            }
            else if(paramConstrainModeP2_str == "hard")
            {
                paramConstrainModeP2 = MODE_PARAM_HARD;
            }
            else
            {
                std::cout << "ERROR: Unrecognized value: paramConstrainModeP2_str=" << paramConstrainModeP2_str << std::endl;
                paramConstrainModeP2 = MODE_PARAM_UNDEFINED;
            }
            */
            // TODO
            
            // print out using actual values, not strings
            //std::cout << paramNumber << ", " << paramName << ", "
            //          << paramEnabled_str << ", "
            //          << paramInitValueP1 << ", " << paramInitErrorP1 << ", "
            //          << paramInitValueP2 << ", " << paramInitErrorP2 << ", "
            //          << paramConstraintValueP1 << ", " << paramConstraintErrorP1 << ", "
            //          << paramConstraintValueP2 << ", " << paramConstraintErrorP2 << ", "
            //          << paramConstrainModeP1 << ", " << paramConstrainModeP2 << ", ";
            // print out the MCNameToParamNameMap
            //for(auto it = MCNameToParamNameMap.cbegin(); it != MCNameToParamNameMap.cend(); ++ it)
            //{
            //    std::cout << it->first << " -> " << it->second << std::endl;
            //}

            /*
            if(enabled_str.CompareTo("enabled") == 0)
            {
                enabled_params.push_back(paramNumber);
            }
            else if(enabled_str.CompareTo("disabled") == 0)
            {
                // do nothing
            }
            */
            if(paramEnabled_str == "enabled")
            {
                auto count = std::count(enabled_params.begin(), enabled_params.end(), paramNumber);
                if(count > 0)
                {
                    std::cout << "ERROR: enabled_params already contains paramNumber=" << paramNumber << std::endl;
                    std::cout << "HALT" << std::endl;
                    std::cin.get();
                }
                else
                {
                    enabled_params.push_back(paramNumber);
                }
            }
            else if(paramEnabled_str == "disabled")
            {
                // do nothing
            }
            else
            {
                std::cout << "ERROR: Unrecognized value: paramEnabled_str=" << paramEnabled_str << std::endl;
                //std::cout << "unknown enabled/disabled parameter specification" << std::endl;
                //std::cout << "enabled_str=" << enabled_str << std::endl;
                continue;
            }

            //if(paramEnabled_str == "enabled")
            //{

                // TODO: consider removing later, this was added after
                // numberParams was changed to 1
                if(paramNumber < numberParams)
                {
                    paramNameMap[paramNumber] = paramName;

                    //paramActMap[paramNumber].push_back(activity);
                    //paramActErrMap[paramNumber].push_back(activity_error);
                    //paramActP2Map[paramNumber].push_back(activity_P2);
                    //paramActErrP2Map[paramNumber].push_back(activity_error_P2);
                    // old names above, new below
                    
                    paramInitValueP1Map[paramNumber] = paramInitValueP1;
                    paramInitErrorP1Map[paramNumber] = paramInitErrorP1;
                    paramInitValueP2Map[paramNumber] = paramInitValueP2;
                    paramInitErrorP2Map[paramNumber] = paramInitErrorP2;

                    paramConstraintValueP1Map[paramNumber] = paramConstraintValueP1;
                    paramConstraintErrorP1Map[paramNumber] = paramConstraintErrorP1;
                    paramConstraintValueP2Map[paramNumber] = paramConstraintValueP2;
                    paramConstraintErrorP2Map[paramNumber] = paramConstraintErrorP2;

                    paramConstrainModeP1Map[paramNumber] = paramConstrainModeP1;
                    paramConstrainModeP2Map[paramNumber] = paramConstrainModeP2;

                    if(thePhase == 0)
                    {
                        if(paramConstrainModeP1 == MODE_PARAM_HARD)
                        {
                            fixed_params.push_back(paramNumber);
                        }
                    }
                    else if(thePhase == 1)
                    {
                        if(paramConstrainModeP2 == MODE_PARAM_HARD)
                        {
                            fixed_params.push_back(paramNumber);
                        }
                    }
                    else
                    {
                        std::cout << "ERROR: thePhase=" << thePhase << " invalid value" << std::endl;
                    }

                    paramNameToNumberMap[paramName] = paramNumber;
            
                    paramMCList[paramNumber] = unique_mc_names;

                }

            //}

            // write back out
            std::cout << std::endl;
            std::cout << "parameter number: " << paramNumber << " " << paramEnabled_str << std::endl;
            std::cout << "Parameter name: " << paramName << std::endl;
            std::cout << "List of MC datafiles: ";
            for(auto it = paramMCList[paramNumber].cbegin(); it != paramMCList[paramNumber].cend(); )
            {
                std::cout << *it;
                if(++ it != paramMCList[paramNumber].cend())
                {
                    std::cout << ", ";
                }
                else
                {
                    break;
                }
            }
            std::cout << std::endl;
            std::cout << "Phase 1: Initial Value: " << paramInitValueP1Map[paramNumber] << " +- " << paramInitErrorP1Map[paramNumber] << std::endl;
            std::cout << "         Constraint: " << paramConstraintValueP1Map[paramNumber] << " +- " << paramConstraintErrorP1Map[paramNumber] << std::endl;
            std::cout << "         Mode: " << paramConstrainModeP1Map[paramNumber];
            if(paramConstrainModeP1Map[paramNumber] == MODE_PARAM_FREE)
            {
                std::cout << " (free)";
            }
            else if(paramConstrainModeP1Map[paramNumber] == MODE_PARAM_SOFT)
            {
                std::cout << " (soft)";
            }
            else if(paramConstrainModeP1Map[paramNumber] == MODE_PARAM_HARD)
            {
                std::cout << " (hard)";
            }
            else
            {
                std::cout << " (unknown/ERROR)";
            }
            std::cout << std::endl;
            std::cout << "Phase 2: Initial Value: " << paramInitValueP2Map[paramNumber] << " +- " << paramInitErrorP2Map[paramNumber] << std::endl;
            std::cout << "         Constraint: " << paramConstraintValueP2Map[paramNumber] << " +- " << paramConstraintErrorP2Map[paramNumber] << std::endl;
            std::cout << "         Mode: " << paramConstrainModeP2Map[paramNumber];
            if(paramConstrainModeP2Map[paramNumber] == MODE_PARAM_FREE)
            {
                std::cout << " (free)";
            }
            else if(paramConstrainModeP2Map[paramNumber] == MODE_PARAM_SOFT)
            {
                std::cout << " (soft)";
            }
            else if(paramConstrainModeP2Map[paramNumber] == MODE_PARAM_HARD)
            {
                std::cout << " (hard)";
            }
            else
            {
                std::cout << " (unknown/ERROR)";
            }
            std::cout << std::endl;

        }

    }

    std::cout << "read: parameter_names.lst -> done" << std::endl;



    // First, create a root file to hold all of the histograms
    TFile *myFile = TFile::Open("Nd150_loglikResults.root", "RECREATE");
    //myFile->Close();
    // 2020-04-02: removed call to close, no idea why this was here
    // or if it is supposed to be here for some reason?

    // Read in all the histograms we want
    //book2DHistograms(0,"1e1g_","P1","hEeVEg_");
    //book2DHistograms(0,"1e1g_","P2","hEeVEg_");
    //book1DHistograms(0,"1e1g_","P1","hClusterEnergy_");
    //book1DHistograms(1,"1e1g_","P2","hClusterEnergy_");
    
    //book1DHistograms(0,"1e_","P1","hEe_");
    //book1DHistograms(0,"1e_","P2","hEe_");
    //book1DHistograms(2,"1e1g_","P1","hTotalE_");
    //book1DHistograms(3,"1e1g_","P2","hTotalE_");
    //book1DHistograms(4,"1e2g_","P1","hClusterEnergyMax_");
    //book1DHistograms(5,"1e2g_","P2","hClusterEnergyMax_");
    //book1DHistograms(2,"1e1gExt_","P1","hClusterEnergy_");
    //book1DHistograms(3,"1e1gExt_","P2","hClusterEnergy_");
    //book1DHistograms(6,"2e_","P1","hTotalE_");
    //book1DHistograms(7,"2e_","P2","hTotalE_");
    //book1DHistograms(4,"OCE_","P1","hTotalE_");
    //book1DHistograms(5,"OCE_","P2","hTotalE_");
    //book1DHistograms(0, "2e_", "P1", "hTotalE_");
    //book1DHistograms(0, "2e_", "P2", "hTotalE_");
    book1DHistograms(0, "2e_", "P" + Phase, "hTotalE_");
  
    // Array to hold activity adjustment parameters
    Double_t AdjustActs[numberParams];
    Double_t AdjustActs_Err[numberParams];
    for(Int_t ix{0}; ix < numberParams; ++ ix)
    {
        AdjustActs[ix] = 1.0;
        AdjustActs_Err[ix] = 1.0;
    }
  
    fitBackgrounds(AdjustActs, AdjustActs_Err, thePhase);

    std::cout << "The following adjustments should be made:" << std::endl;
    for(int i = 0; i < numberParams; i++)
    {
        std::cout << i << " :\t" << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;  
        std::ofstream myFile("fit_results.txt", std::ios::out | std::ios::app);
        myFile << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;

        myFile.close();
    }
  
    //Now let's draw the results...

    THStack *stacks1D[number1DHists];
    TCanvas *c;
    TH1F *data1D[number1DHists];

    // TODO: copy over legend code from fit_2e.C

    std::cout << "debug: number of data samples: " << allDataSamples1D->GetEntries() << std::endl;
    std::cout << "debug: number of MC samples: " << allMCSamples1D->GetEntries() << std::endl;


    // each channel 1D hists
    for(int i = 0; i < allDataSamples1D->GetEntries(); i++)
    {

        data1D[i] = (TH1F*)allDataSamples1D->At(i)->Clone();

        TString i_str;
        i_str.Form("%i", i);
    
        c = new TCanvas("c" + i_str);
        c->SetFillColor(kWhite);

        
        stacks1D[i] = new THStack("stacks1D" + i_str, i_str);


        TH1F *tmpHist_draw1D;

        // j list MC samples for this channel i
        std::cout << "debug: number of MC samples (i=" << i << "): " << allMCSamples1D[i]->GetEntries() << std::endl;

        std::cout << "list of all objects in allMCSamples1D[i=" << i << "]" << std::endl;
        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {
            //std::cout << (((TH1F*)allMCSamples1D[i])->At(j))->GetName() << std::endl;
            tmpHist_draw1D = (TH1F*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw1D->GetName();
            std::cout << tmpName << std::endl;
        }


        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {

            std::cout << "j=" << j << std::endl;

            TString j_str;
            j_str.Form("%i", j);

            tmpHist_draw1D = (TH1F*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw1D->GetName();

            std::cout << "(1) tmpName=" << tmpName << std::endl;

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param = -1;
            bool foundParam = false;

            // search through parameters to find right one
            for(int k = 0; (k < numberParams) && !foundParam; k++)
            {

                //std::cout << "k=" << k << std::endl;

	            // match up the isotope with parm
                for(int n = 0; (n < paramNameMap[k].size()) && !foundParam; n++)
                {
                    //std::cout << "n=" << n << std::endl;

                    // summer
	                //std::cout << "is it...  " << paramNameMap[k].at(j) << std::endl;
                    // me
	                //std::cout << "is it...  " << paramNameMap[k].at(n) << std::endl;

	                if(tmpName.Contains(paramNameMap[k].at(n)))
                    {
                        //std::cout << tmpName << " contains " << paramNameMap[k].at(n) << std::endl;

	                    foundParam = true;
                        which_param = k;
	                }
                    else
                    {
                        //std::cout << tmpName << " does not contain " << paramNameMap[k].at(n) << std::endl;
                    }
	            }
            }

            if(foundParam)
            {
                std::cout << "found histogram: tmpName=" << tmpName << std::endl;

                // scale histogram to correct size using output parameter
                // from fit
                tmpHist_draw1D->Scale(AdjustActs[which_param]);
                
                if(tmpHist_draw1D->Integral() > 0)
                {
                    stacks1D[i]->Add(tmpHist_draw1D);
	            }
                else
                {
                    std::cout << "not adding to stack, Integral() <= 0" << std::endl;
                }
            }
            else
            {
                std::cout << "error could not find histograms" << std::endl;
                std::cout << "tmpName=" << tmpName << std::endl;
            } 

        }


        stacks1D[i]->SetMaximum(350.);
        stacks1D[i]->Draw("hist");
        data1D[i]->SetLineWidth(2);
        data1D[i]->SetMarkerStyle(20);
        data1D[i]->SetMarkerSize(0.5);
        data1D[i]->Draw("PEsames");

        std::cout << "saving to file (canvas)" << std::endl;
        c->SaveAs("finalHisto1D_" + i_str + ".C");
        c->SaveAs("finalHisto1D_" + i_str + ".eps");
        c->SaveAs("finalHisto1D_" + i_str + ".png");
    
    }

    // each channel 1D hists
    THStack *stacks2D[number2DHists];
    TH2F *data2D[number2DHists];
    for(int i = 0; i < allDataSamples2D->GetEntries(); i++)
    {
        data2D[i] = (TH2F*)allDataSamples2D->At(i)->Clone();
        TString i_str;
        i_str.Form("%i",i);
        c = new TCanvas("c_"+i_str);
        c->SetFillColor(kWhite);
        stacks2D[i] = new THStack("stacks2D"+i_str,i_str);

        // j list MC samples for this channel i
        TH2F *tmpHist_draw2D;
        for(int j = 0; j < allMCSamples2D[i]->GetEntries(); j++)
        {
            TString j_str;
            j_str.Form("%i",j);
            tmpHist_draw2D = (TH2F*)allMCSamples2D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw2D->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            
            // k searching through array of params for the right one
            int which_param;
            bool foundParam = false;
            for(int k = 0; (k < numberParams) && !foundParam; k++)
            {
            
                // match up the isotope with parm
                for(int n = 0; (n < paramNameMap[k].size()) && !foundParam; n++)
                {
	                //std::cout <<"is it...  " <<paramNameMap[k].at(j) << std::endl;

	                if(tmpName.Contains(paramNameMap[k].at(n)))
                    {
	                    foundParam = true;
                        which_param = k;
	                }
	            }
            }

            if(foundParam)
            {
                tmpHist_draw2D->Scale(AdjustActs[which_param]);
     
                if(tmpHist_draw2D->Integral() > 0)
                {
                    stacks2D[i]->Add(tmpHist_draw2D);
                }
            }
            else
            {
                std::cout << "error could not find histograms" << std::endl;
            }

        }

        stacks2D[i]->Draw("");
        data2D[i]->SetLineWidth(2);
        data2D[i]->SetMarkerStyle(20);
        data2D[i]->SetMarkerSize(0.5);
        data2D[i]->Draw("PEsames");

        c->SaveAs("finalHisto2D_"+i_str+".C");
        c->SaveAs("finalHisto2D_"+i_str+".eps");
        c->SaveAs("finalHisto2D_"+i_str+".png");
        std::cout << "saved 2D histogram to file" << std::endl;

    }

}


///////////////////////////////////////////////////////////////////////////////
// book1DHistograms
///////////////////////////////////////////////////////////////////////////////

void book1DHistograms_helper(Int_t channel_counter, TString theChannel, TString thePhase, TString theHistogram,
    const int nBkgs, TString *BkgFiles, TH1F *tmpHist)
{

    for(int i = 0; i < nBkgs; i++)
    {
        // check if parameter is enabled
        // convert parameter string name to index
        //std::cout << "searching for string " << ExternalBkgOneParamFiles[i] << " in paramNameMap" << std::endl;
        //TString search_object = ExternalBkgOneParamFiles[i];
        std::string mc_name = std::string(BkgFiles[i]);
        std::string search_object = MCNameToParamNameMap.at(mc_name);
        if(paramNameToNumberMap.count(search_object) > 0)
        {
            int param_number = paramNameToNumberMap.at(search_object);
            //std::cout << "parameber number " << param_number << " is in the paramNameToNumberMap" << std::endl;
            if(std::find(enabled_params.begin(), enabled_params.end(), param_number) != enabled_params.end())
            {
                std::cout << "parameter number " << param_number << " is enabled" << std::endl;

                std::string name(theHistogram + BkgFiles[i] + "_fit");
                if(gDirectory->GetListOfKeys()->Contains(theHistogram + BkgFiles[i] + "_fit"))
                {
                    //check if the histograms exists 
                    std::string hist_name(BkgFiles[i] + "_" + theChannel + thePhase);
                    //std::cout << "Get() : " << name << " from file, Clone() : " << hist_name << std::endl;
                    tmpHist = (TH1F*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
                    allMCSamples1D[channel_counter]->Add(tmpHist);
                    // TODO: does this work as expected for secular equlibrium samples?

                    //std::cout << tmpHist->GetName() << std::endl;

                }
                else
                {
                    std::cout << "gDirectory->GetListOfKeys() does not contain " << name << " - disabling parameter number " << param_number << std::endl;
                    // cannot find histogram input data, so disable parameter
                    std::remove(enabled_params.begin(), enabled_params.end(), param_number);
                }
            }
            else
            {
                std::cout << "parameter number " << param_number << " is not enabled (not found in vector)" << std::endl;
            }
        }
        else
        {
            std::cout << "!!!!! ERROR: search_object=" << search_object << " not found in paramNameToNumberMap" << std::endl;
            std::cout << "mc_name=" << mc_name << std::endl;

            std::cout << "contents of map paramNameToNumberMap:" << std::endl;
            for(auto it = paramNameToNumberMap.cbegin(); it != paramNameToNumberMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
            std::cout << "contents of map MCNameToParamNameMap:" << std::endl;
            for(auto it = MCNameToParamNameMap.cbegin(); it != MCNameToParamNameMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
        }
    }
}

// channel_counter = 0
// theChannel = "2e_"
// thePhase = "P1"
// theHistogram = "hTotalE_"
void book1DHistograms(Int_t channel_counter, TString theChannel, TString thePhase, TString theHistogram) {

    std::cout << "booking 1D hists for " << theChannel << " " << thePhase << std::endl;
    allMCSamples1D[channel_counter] = new TObjArray();

    TFile *aFile = TFile::Open("/home/ebirdsall/NEMO3/Nd150_analysis/MeasureStuff/new/Macros/Nd150_" + theChannel + thePhase + ".root");
    //gDirectory->cd("singleHistos");
    //gDirectory->ls();


    TH1F *tmpHist  = new TH1F("tmpHist_" + theChannel + thePhase, "" , 1, 0, 1);
    
    book1DHistograms_helper(channel_counter, theChannel,
                            thePhase, theHistogram,
                            nExternalBkgs_oneParam,
                            ExternalBkgOneParamFiles,
                            tmpHist);

    // TODO: does this work as expected for secular equlibrium samples?

    book1DHistograms_helper(channel_counter, theChannel,
                            thePhase, theHistogram,
                            nInternalBkgs,
                            InternalBkgFiles,
                            tmpHist);

    book1DHistograms_helper(channel_counter, theChannel,
                            thePhase, theHistogram,
                            nRn222Bkgs,
                            Rn222BkgFiles,
                            tmpHist);

    book1DHistograms_helper(channel_counter, theChannel,
                            thePhase, theHistogram,
                            nRn220Bkgs,
                            Rn220BkgFiles,
                            tmpHist);

    book1DHistograms_helper(channel_counter, theChannel,
                            thePhase, theHistogram,
                            nNd150Samples,
                            Nd150Files,
                            tmpHist);

    book1DHistograms_helper(channel_counter, theChannel,
                            thePhase, theHistogram,
                            nNeighbours,
                            NeighbourFiles,
                            tmpHist);

    // TODO here
    // what is name in other section of code
    //std::string name(theHistogram + "data_2e");
    std::string name(theHistogram + "data");
    if(gDirectory->GetListOfKeys()->Contains(name.c_str()))
    {
        std::string hist_name("data_" + theChannel + thePhase);
        //std::cout << "Get() : " << name << " from file, Clone() : " << hist_name << std::endl;
        tmpHist = (TH1F*)gDirectory->Get(name.c_str())->Clone(hist_name.c_str());
        allDataSamples1D->Add((TH1F*)tmpHist);
    }
    else
    {
        std::cout << "gDirectory->GetListOfKeys() does not contain " << name << std::endl;
    }
    /*
    if(gDirectory->GetListOfKeys()->Contains(theHistogram + "Data"))
    {
        std::string name(theHistogram + "Data");
        std::cout << "Get() : " << name << " from file, Clone() : " << "Data_" + theChannel + thePhase << std::endl;
        tmpHist = (TH1F*)gDirectory->Get(name.c_str())->Clone("Data_" + theChannel + thePhase);
        allDataSamples1D->Add((TH1F*)tmpHist);
    }
    else
    {
        std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + "Data" << std::endl;
    }
    */

    // std::cout << tmpHist->GetName() << std::endl;
    // tmpHist->Delete();
    // aFile->Close();
    // aFile->Delete();
}


///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////

void book2DHistograms(Int_t channel_counter, TString theChannel, TString thePhase, TString theHistogram) {

  std::cout << "booking 2D hists for "<< theChannel<<" " <<thePhase<<std::endl;
  allMCSamples2D[channel_counter] = new TObjArray();

  TFile *aFile = TFile::Open("/home/blotsd/NEMO3/Nd150_analysis/MeasureStuff/Macros/Nd150_"+theChannel+thePhase+".root");
  gDirectory->cd("singleHistos");
  // gDirectory->ls();
  TH2F *tmpHist2  = new TH2F("tmpHist2_"+theChannel+thePhase,"",1,0,1,1,0,1);
  for ( int i = 0; i < nExternalBkgs_oneParam; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+ExternalBkgOneParamFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+ExternalBkgOneParamFiles[i]+"_fit")->Clone(ExternalBkgOneParamFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);

      //std::cout << tmpHist->GetName() << std::endl;

    }
  }
  /*
  for ( int i = 0; i < nExternalBkgs_twoParam; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+ExternalBkgTwoParamFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+ExternalBkgTwoParamFiles[i]+"_fit")->Clone(ExternalBkgTwoParamFiles[i]+"_"+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);

      //std::cout << tmpHist->GetName() << std::endl;

    }
  }*/

  for ( int i = 0; i < nInternalBkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+InternalBkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+InternalBkgFiles[i]+"_fit")->Clone(InternalBkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nRn222Bkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Rn222BkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Rn222BkgFiles[i]+"_fit")->Clone(Rn222BkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nRn220Bkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Rn220BkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Rn220BkgFiles[i]+"_fit")->Clone(Rn220BkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nNd150Samples; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Nd150Files[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Nd150Files[i]+"_fit")->Clone(Nd150Files[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nNeighbours; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+NeighbourFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+NeighbourFiles[i]+"_fit")->Clone(NeighbourFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  if ( gDirectory->GetListOfKeys()->Contains(theHistogram+"Data") ) {
    tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+"Data")->Clone("Data_"+theChannel+thePhase);
    allDataSamples2D->Add((TH2F*)tmpHist2);
  }
  // NOTE: Data changed to data but not updated here

  std::cout << tmpHist2->GetName() << std::endl;
  // tmpHist2->Delete();
  //aFile->Close();
  // aFile->Delete();
}


///////////////////////////////////////////////////////////////////////////////
// fitBackgrounds
///////////////////////////////////////////////////////////////////////////////

void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, Int_t thePhase)
{

    std::cout << ">>>>> fitBackgrounds()" << std::endl;

    //TVirtualFitter::SetDefaultFitter("Minuit2");
    TMinuit * minuit = new TMinuit(numberParams);

    std::cout << "Fit created" << std::endl;

    // moved reading of parameter list file to before book histogram function
    // calls
  
    TString sample_names[numberParams];
    // read samples names out again just to make sure its correct (useful when changing the number of params)
    for(int i = 0; i < numberParams; i++)
    {
        // tmpStr doesn't appear to do anything ?
        // TODO
        TString tmpStr = "";

        TString i_str;
        i_str.Form("%i", i);

        sample_names[i] = i_str;
        
        std::cout << "param " << i << ": sample_names[" << i << "]=" << sample_names[i] << " -> ";
        std::cout << paramNameMap[i] << " : ";
        bool first = true;
        for(auto it = MCNameToParamNameMap.cbegin(); it != MCNameToParamNameMap.cend(); ++ it)
        {
            if(it->second == paramNameMap[i])
            {
                // match for MC name
                if(!first) std::cout << ", ";
                std::cout << it->first;
                first = false;
            }
        }
        std::cout << std::endl;

        /*
        for(int j = 0; j < paramNameMap[i].size(); j++)
        { 
            std::cout << paramNameMap[i].at(j);
            tmpStr += paramNameMap[i].at(j);
            
            if(j + 1 < paramNameMap[i].size())
            {
                std::cout << ", ";
                tmpStr += ",";
            }
        }
        */
        
        std::cout << std::endl;
    }

    std::cout << "set errors" << std::endl;

    // Set the parameters for the fit, and give them an arbitrary 10% error to start.
    for(int i = 0; i < numberParams; i++)
    {
        // AdjustActs[i] = 1.;
        //AdjustActs_Err[i] = 1.;
        TString i_str;
        i_str.Form("%i",i);
        std::cout << "DefineParameter: i=" << i << std::endl;
        //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 1000.0);

        // get initial parameter values and error
        Double_t param_init_value = 0.;
        Double_t param_init_error = 0.; 
        if(thePhase == 0)
        {
            /*
            for(int j = 0; j < paramInitValueP1Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP1Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP1Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP1Map[i].at(j);
            }
            */

            param_init_value = paramInitValueP1Map[i];
            param_init_error = paramInitErrorP1Map[i];
        }
        else if(thePhase == 1)
        {
            /*
            for(int j = 0; j < paramInitValueP2Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP2Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP2Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP2Map[i].at(j);
            }
            */

            param_init_value = paramInitValueP2Map[i];
            param_init_error = paramInitErrorP2Map[i];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }


        if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            // define parameter using constrained value if hard constrained

            std::cout << "minuit: fixing parameter " << i << std::endl;
            
            //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 2.0);
            minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 2.0);
            minuit->FixParameter(i);
        }
        // TODO: not sure if this is the correct thing to do
        // I think that disabled params do not necessarily appear in the
        // fixed params list, and therefore I have to fix them (not sure
        // about that) here if they are disabled
        else if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
        {
            // NOT enabled, therefore fix
            
            std::cout << "fixing parameter i=" << i << std::endl;
         
            // leave parameter (amplitude) at default value of 1.0, which
            // is more likely to show up any errors in code
            minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 2.0);
            //minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 2.0);
            minuit->FixParameter(i);
        }
        else
        {
            // define parameter using initial value if free/soft constrained
            
            std::cout << "parameter i=" << i << " is not fixed, leaving free" << std::endl;

            minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 10000.0);
            
        }
        // TODO: what about disabled parameters

        //if ((i != 17 ) && (i < 32) || (i > 38)  ) {
        //if ( (i > 32) ) {
        //  minuit->FixParameter(i);
        //  fixed_params.push_back(i);
        //}

        //if ( i == 24 ) {
        //  minuit->FixParameter(i);
        //  fixed_params.push_back(i);
        //}
    }

    std::cout << "all parameters fixed" << std::endl;

    // Fix the tracker radon activities
    //minuit->FixParameter(34);
    //minuit->FixParameter(35);
    //minuit->FixParameter(23);
    // minuit->FixParameter(24);
    /*
    minuit->FixParameter(1);
    minuit->FixParameter(3);
    minuit->FixParameter(7);
    minuit->FixParameter(11);
    minuit->FixParameter(21);
    minuit->FixParameter(22);
    minuit->FixParameter(29);
    minuit->FixParameter(30);
    //minuit->FixParameter(31);
    minuit->FixParameter(34);
    minuit->FixParameter(35);
    //minuit->FixParameter(48);
    minuit->FixParameter(49);
    minuit->FixParameter(50);
    minuit->FixParameter(51);
    minuit->FixParameter(52);
    minuit->FixParameter(53);
    minuit->FixParameter(54);

    fixed_params.push_back(1);
    fixed_params.push_back(3);
    fixed_params.push_back(7);
    fixed_params.push_back(11);
    fixed_params.push_back(21);
    fixed_params.push_back(22);
    fixed_params.push_back(29);
    fixed_params.push_back(30);
    //fixed_params.push_back(31);
    fixed_params.push_back(34);
    fixed_params.push_back(35);
    fixed_params.push_back(49);
    fixed_params.push_back(50);
    fixed_params.push_back(51);
    fixed_params.push_back(52);
    fixed_params.push_back(53);
    fixed_params.push_back(54);
    */
    //fixed_params.push_back(48);

    // Sort out fixed vs free params (this can be improved, its too hackish now)

    // data stored here used as name for each item in correlation matrix
    // TODO other names should be added not just paramNameMap[i].front()
    // loop over paramNameMap[i], adding each name to the free_param_names
    // vector
    std::vector<TString> free_param_names;
    std::vector<TString> enabled_param_names;

    for(int i = 0; i < numberParams; i++)
    {
        TString i_str;
        i_str.Form("%i", i);

        bool fixed = false;
        if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            fixed = true;
        }
        /*
        for(int j = 0; j < fixed_params.size(); j++)
        {
            if(i == fixed_params.at(j))
            {
                fixed = true;
            }
        }
        */

        if(fixed)
        {
            //continue;
            // do nothing
        }
        else
        {
            free_params.push_back(i);
            index_free_params.push_back(i_str);
            free_param_names.push_back(paramNameMap[i].front());
            //for(int ix = 0; ix < paramNameMap[i].size(); ++ ix)
            //{
            //    tmpStr += paramNameMap[i][ix];
            //    if(ix + 1 < paramNameMap[i].size())
            //    {
            //        tmpStr += ", ";
            //    }
            //}
            //TODO: enable this code
        }

        // enabled/disabled
        bool enabled = false;
        if(std::find(enabled_params.begin(), enabled_params.end(), i) != enabled_params.end())
        {
            enabled = true;
        }
        if(enabled)
        {
            // do nothing
            index_enabled_params.push_back(i_str);
            enabled_param_names.push_back(paramNameMap[i].front());
        }
        else
        {
            disabled_params.push_back(i);
            //disabled_param_names.push_back(paramNameMap[i].front());
        }
    }

    // TODO: numbers appear multiple times
    // print out a list of fixed and free parameters
    std::cout << "List of parameters: free/fixed" << std::endl;
    for(int i = 0; i < free_params.size(); ++ i)
    {
        std::cout << "free parameter: " << free_params[i] << std::endl;
    }
    for(int i = 0; i < fixed_params.size(); ++ i)
    {
        std::cout << "fixed parameter: " << fixed_params[i] << std::endl;
    }
    std::cout << "List of parameters: enabled/disabled" << std::endl;
    for(int i = 0; i < enabled_params.size(); ++ i)
    {
        std::cout << "enabled parameter: " << enabled_params[i] << std::endl;
    }
    for(int i = 0; i < disabled_params.size(); ++ i)
    {
        std::cout << "diabled parameter: " << disabled_params[i] << std::endl;
    }



    // Array to pass arguments to fitter
    //double arglist[100];
    //arglist[0] = 0;

    // Set print level
    //  minuit->ExecuteCommand("SET PRINT",arglist,2);
    std::cout << "lets do the fit" << std::endl;

    // Do minimisation
    // arglist[0] = 50000;  // number of function calls
    // arglist[1] = 0.1;  // tolerance

    minuit->SetErrorDef(0.5);
    minuit->SetMaxIterations(50000);
    //give it the function to use
    minuit->SetFCN(logLikelihood);
    //std::cout << "calling: minuit->mnsimp()" << std::endl;
    //minuit->mnsimp();
    std::cout << "calling: minuit->Migrad()" << std::endl;
    minuit->Migrad();
    //minuit->ExecuteCommand("SIMPLEX",arglist,2);

    // Then get results
    for(int i = 0; i < numberParams; i++)
    {
        minuit->GetParameter(i, AdjustActs[i], AdjustActs_Err[i]);
        //AdjustActs_Err[i] = minuit->GetParError(i);
    }


    Double_t CovMatrix[free_params.size()][free_params.size()];
    minuit->mnemat(&CovMatrix[0][0],free_params.size());

  
    /*
    TH2D* hCorrMatrix = new TH2D("hCorrMatrix","Correlation Matrix", free_params.size(), 0, free_params.size(), free_params.size(), 0,free_params.size());

    for(int i = 0; i < free_params.size(); i++)
    {
        for(int j = 0; j < free_params.size(); j++)
        {
          
            //std::cout << index_free_params.at(i) << "  " << index_free_params.at(j)<< std::endl;
            //CovMatrix[i][j] = minuit->GetCovarianceMatrixElement(i,j);

            // std::cout << CovMatrix[i][j] << std::endl;
            Double_t value = CovMatrix[i][j]/(sqrt(CovMatrix[i][i])*sqrt(CovMatrix[j][j]));
            hCorrMatrix->Fill(free_params_names.at(i), free_params_names.at(j), value);

        }
    }
  
    */
    // these were commented out for some reason, I did not remove them
    /*
    double chi2, edm, errdef; 
    int nvpar, nparx;
    minuit->GetStats(chi2, edm, errdef, nvpar, nparx);
    int ndf = npfits-nvpar;
    */
    
    // create a function to draw the final stacked histograms
    /*
    hCorrMatrix->LabelsDeflate();
    hCorrMatrix->GetZaxis()->SetRangeUser(-1,1);
    hCorrMatrix->SetTitle("Correlation Matrix");
    gStyle->SetGridStyle(0);
    hCorrMatrix->GetXaxis()->LabelsOption("v");

    TCanvas *c = new TCanvas("c");
    c->SetFillColor(kWhite);
    c->SetGrid();
    hCorrMatrix->Draw("colz");
    */

}



///////////////////////////////////////////////////////////////////////////////
// logLikelihood
///////////////////////////////////////////////////////////////////////////////

void logLikelihood(Int_t & /*nPar*/, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{


    double loglik = 0.; 
    //double tmp;


//   std::cout << "getting 1D histograms" << std::endl;

    TH1F *tmpData1D;
    // std::cout << allDataSamples1D->GetEntries()  << std::endl;

    // there are i samples for each channel
    for(int i = 0; i < allDataSamples1D->GetEntries(); ++ i)
    {
        TString i_str;
        i_str.Form("%i", i);
        //  std::cout << i << std::endl;

        tmpData1D = (TH1F*)allDataSamples1D->At(i)->Clone("tmpData1D" + i_str + "_");

        // std::cout << tmpData1D->Integral() << std::endl;

        int nBinsX = tmpData1D->GetNbinsX();
        for(int ix = 1; ix <= nBinsX; ++ ix)
        {
            Int_t nData = (Int_t)tmpData1D->GetBinContent(ix);
            // i is the index of the sample ?
            // ix is the bin index
            // p is a pointer to an array of parameter values
            double nMC = getNumberMC1D(i, ix, p);
            //Int_t new_i = -1;
            //TString name = names.at(i);
            //new_i = paramNameToNumberMap[name];
            //std::cout << "the new i value is new_i=" << new_i << std::endl;
            //double nMC = getNumberMC1D(new_i, ix, p);

            //std::cin.get();

            if(nMC > 0. && TMath::Poisson(nData, nMC) > 0.)
            {
	            loglik += TMath::Log(TMath::Poisson(nData, nMC));
            }
            else
            {
	            loglik -= 10.;
            }
        
        } //~bins
    } //channels
    
    
    //  std::cout << "getting 2D histograms" << std::endl;

    TH2F *tmpData2D;
    // std::cout << allDataSamples->GetEntries()  << std::endl;
    for(int i = 0; i < allDataSamples2D->GetEntries(); ++ i)
    {
        // there are i samples for each channel
        TString i_str;
        i_str.Form("%i", i);

        tmpData2D = (TH2F*)allDataSamples2D->At(i)->Clone("tmpData2D" + i_str + "_");

        int nBinsX = tmpData2D->GetNbinsX();
        int nBinsY = tmpData2D->GetNbinsY();

        for(int ix = 1; ix <= nBinsX; ++ ix)
        {
            for(int iy = 1; iy <= nBinsY; ++ iy)
            {
                Int_t nData = (Int_t)tmpData2D->GetBinContent(ix, iy);
	            double nMC = getNumberMC2D(i, ix, iy, p);

	            if(nMC > 0 && TMath::Poisson(nData, nMC) > 0)
                {
	                loglik += TMath::Log(TMath::Poisson(nData, nMC));
	            }
                else
                {
	                loglik -= 10.;
	            }
            } //~ybins
        } //~xbins
    } //channels

 
    // add constraints to improve likelihood
    //will eventually add gaussian constraint
    // disabled
    /*
    for ( int i = 0; i < 14; i++ ) {
        double constraint = 0.;
        //constraint = pow( (p[i] - 1) / (0.01), 2);
        // std::cout << i << "   " << p[i];
        bool fixed = false;
        for ( int j = 0; j < fixed_params.size(); j++ ) {
            if ( i == fixed_params.at(j) ) fixed = true;
        }

        //if ( (i <2) || (i == 7 ) || (i==10) || (i == 14) || (i==15) || (i==16) || (i == 18) || (i == 19) || (i == 21) || 9i==22) ) {
        // constraint = pow( (p[i] - 1.) / (paramActErrMap[i].front() / paramActMap[i].front()), 2);
        constraint = TMath::Gaus(p[i],1.,(paramActErrMap[i].front() / paramActMap[i].front()),true);
      

        // std::cout << "   " << constraint << std::endl;
        //  loglik += log(constraint);  
	    // }
    }
    */

    for(int i = 0; i < numberParams; ++ i)
    {

        if(thePhase == 0)
        {
            if(paramConstrainModeP1Map[i] == MODE_PARAM_SOFT)
            {
                // do nothing, soft constraint will be applied below
            }
            else if(paramConstrainModeP1Map[i] == MODE_PARAM_HARD)
            {
                // parameter fixed by minuit, continue to next param
                continue;
            }
            else if(paramConstrainModeP1Map[i] == MODE_PARAM_FREE)
            {
                // no constraint to apply, continue to next param
                continue;
            }
            else
            {
                std::cout << "ERROR: Invalid value in paramConstrainModeP1Map: paramConstrainModeP1Map[" << i << "]=" << paramConstrainModeP1Map[i] << std::endl;
            }
        }
        else if(thePhase == 1)
        {
            if(paramConstrainModeP2Map[i] == MODE_PARAM_SOFT)
            {
                // do nothing, soft constraint will be applied below
            }
            else if(paramConstrainModeP2Map[i] == MODE_PARAM_HARD)
            {
                // parameter fixed by minuit, continue to next param
                continue;
            }
            else if(paramConstrainModeP2Map[i] == MODE_PARAM_FREE)
            {
                // no constraint to apply, continue to next param
                continue;
            }
            else
            {
                std::cout << "ERROR: Invalid value in paramConstrainModeP2Map: paramConstrainModeP2Map[" << i << "]=" << paramConstrainModeP2Map[i] << std::endl;
            }
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        
        // soft constraint is applied here
        double constraint = 0.;
        double error = 0.;

        /*
        if(thePhase == 0)
        {
            for(int j = 0; j < paramInitValueP1Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP1Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP1Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP1Map[i].at(j);
            }
        }
        else if(thePhase == 1)
        {
            for(int j = 0; j < paramInitValueP2Map[i].size(); ++ j)
            {
                param_init_value += paramInitValueP2Map[i].at(j);
            }
            for(int j = 0; j < paramInitErrorP2Map[i].size(); ++ j)
            {
                param_init_error += paramInitErrorP2Map[i].at(j);
            }
        }
        */
        if(thePhase == 0)
        {
            //for(int j = 0; j < paramConstraintValueP1Map[i].size(); ++ j)
            //{
                //constraint += paramConstraintValueP1Map[i].at(j);
            //}
            //for(int j = 0; j < paramConstraintErrorP1Map[i].size(); ++ j)
            //{
                //error += paramConstraintErrorP1Map[i].at(j);
            //}
            constraint = paramConstraintValueP1Map[i];
            error = paramConstraintErrorP1Map[i];
            // TODO: this is a bit weird and perhaps not what would
            // be expected looking at the input file
            // perhaps it would be better to use .at(0)
            // and ignore all later values in parameters
            // in particular, errors do not usually add, but usually
            // add in quadrature, which is not what happens here
            // TODO: consider adding errors in quadrature, or
            // changing input file to a list of samples as a single parameter
            // rather than having multiple parameter numbers, each with a
            // single sample
            //
            //constraint = paramConstraintValueP1Map[i];
            //error = paramConstraintErrorP1Map[i];
        }
        else if(thePhase == 1)
        {
            //for(int j = 0; j < paramConstraintValueP2Map[i].size(); ++ j)
            //{
            //    constraint += paramConstraintValueP2Map[i].at(j);
            //}
            //for(int j = 0; j < paramConstraintErrorP2Map[i].size(); ++ j)
            //{
            //    error += paramConstraintErrorP2Map[i].at(j);
            //}
            constraint = paramConstraintValueP2Map[i];
            error = paramConstraintErrorP2Map[i];

            //constraint = paramConstraintValueP2Map[i];
            //error = paramConstraintErrorP2Map[i];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        // TODO: can optimize this code

        double param_value = p[i];
        double penalty = std::pow((param_value - constraint) / error, 2.0);

        fval += penalty;
    }
  
    fval = -2.0 * loglik; 
    //tmpData->Delete();

}


Double_t getNumberMC1D(Int_t channel, Int_t binx, Double_t *p) {


    //std::cout << "getNumberMC1D() called with channel=" << channel << " binx=" << binx << std::endl;
    //std::cout << "printing contents of parameter array" << std::endl;
    //for(int i = 0; i < numberParams; ++ i)
    //{
    //    std::cout << "i=" << i << " p[i]=" << p[i] << std::endl;
    //}


    // std::cout <<"getting number of1D MC... "  <<channel << std::endl;

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1F *tmpHist;
    int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;

    //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

    for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    {
        tmpHist = (TH1F*)allMCSamples1D[channel]->At(k);

        //if ( tmpHist->Integral() == 0 ) continue;
       
        TString tmpName = tmpHist->GetName();

        //std::cout << "looking for " << tmpName << std::endl;
        bool foundParam = false;

        for(int i = 0; (i < numberParams) && !foundParam; i++)
        {
            // std::vector<int>::iterator it = find(fixed_params.begin(), fixed_params.end(), i);
            //if ( it != fixed_params.end() ) continue;      

            for(int j = 0; (j < paramNameMap[i].size()) && !foundParam; j++)
            {
                //std::cout << "is it: " << paramNameMap[i].at(j) << std::endl;

	            if(tmpName.Contains(paramNameMap[i].at(j)))
                {
	                foundParam = true;
                    which_param = i;
                    //std::cout << "match found: " << tmpName << " -> " << paramNameMap[i].at(j) << ", which_param=" << which_param << std::endl;
	            }

            } //~j searching through array of params for the right one
        } //~i list of isotopes with same parameter

        if(foundParam)
        {
            //std::cout << "adding to nMC with index of which_param=" << which_param << std::endl;
            // TODO: think this is collecting the wrong parameter? or is it?
            nMC += p[which_param] * tmpHist->GetBinContent(binx);
        }
        else
        {
            std::cout << "error could not find histogram: " << tmpName << std::endl;
        }	 
    } //each histogram

    return nMC;
}


Double_t getNumberMC2D(Int_t channel, Int_t binx, Int_t biny, Double_t *p) {

  // std::cout <<"getting number of 2D MC... "  <<channel << std::endl;

  double nMC = 0;

  // (1) grab a hist from the sample list of this channel
  // (2) figure out which parameter it corresponds to
  TH2F *tmpHist2;
  int which_param;
  //  std::cout <<"getting number of MC... "  <<channel << std::endl;

  //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

  for ( int k = 0; k < allMCSamples2D[channel]->GetEntries(); k++ ) {
    tmpHist2 = (TH2F*)allMCSamples2D[channel]->At(k);
    TString tmpName = tmpHist2->GetName();

    //  std::cout << "looking for " << tmpName << std::endl;
    bool foundParam = false;

    for ( int i = 0; (i < numberParams) && !foundParam; i++ ) {

      std::vector<int>::iterator it = find(fixed_params.begin(), fixed_params.end(), i);
      if ( it != fixed_params.end() ) continue;    

      for ( int j = 0; (j < paramNameMap[i].size()) && !foundParam; j++ ) {
	//	   std::cout <<"is it...  " <<paramNameMap[i].at(j) << std::endl;

	if ( tmpName.Contains(paramNameMap[i].at(j)) ) {
	  foundParam = true;
          which_param = i;
	}

      }//~j searching through array of params for the right one
    }//~i list of isotopes with same parameter

    if ( foundParam ) {nMC += p[which_param]*tmpHist2->GetBinContent(binx,biny);}
    // else {std::cout << "error could not find histogram: " << tmpName << std::endl;}	 
  }//each histogram

  return nMC;
}
