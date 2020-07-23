#ifndef NEWLOGLIKFITTER_READ_PARAMETERNAMES_LST_H
#define NEWLOGLIKFITTER_READ_PARAMETERNAMES_LST_H

#include "TH1.h"
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
#include <vector>
#include <TMinuit.h>


// 
#include <time.h>


#include "newLogLikFitter.h"





//-----------------------------------------
//    Functions
//----------------------------------------







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








void clear_maps()
{

    MCNameToParamNameMap.clear();
    MCNameToParamNumberMap.clear();

    paramNameToNumberMap.clear();

    paramNumberToMinuitParamNumberMap.clear();
    minuitParamNumberToParamNumberMap.clear();
    minuitParamNumberCounter = 0;

    numberEnabledParams = 0;

    fixed_params.clear();
    free_params.clear();

    //paramNameToHumanReadableParamNameMap.clear();
    MCSampleNameToHumanReadableMCSampleNameMap.clear();



    enabled_params.clear();
    disabled_params.clear();

}


void randomize_parameters()
{

    std::cout << "randomizing initial parameters" << std::endl;
    
    TRandom3 rng(1);

    for(int i = 0; i < numberParams; i++)
    {
        Double_t randomnumber = rng.Gaus();
        Double_t activity = 0.;
        Double_t uncertainty = 0.;
        Double_t randomactivity = 0.;

        if(thePhase == 0)
        {
            activity = paramInitValueP1Map[i];
            uncertainty = paramInitErrorP1Map[i];
            randomactivity = activity + uncertainty * randomnumber;

            paramInitValueP2Map[i] = randomactivity;

            std::cout << "parameter number: " << i
                      << " randomnumber=" << randomnumber
                      << " activity=" << activity
                      << " uncertainty=" << uncertainty
                      << " randomactivity=" << randomactivity << std::endl;
        }
        else if(thePhase == 1)
        {
            activity = paramInitValueP2Map[i];
            uncertainty = paramInitErrorP2Map[i];
            randomactivity = activity + uncertainty * randomnumber;

            paramInitValueP2Map[i] = randomactivity;

            std::cout << "parameter number: " << i
                      << " randomnumber=" << randomnumber
                      << " activity=" << activity
                      << " uncertainty=" << uncertainty
                      << " randomactivity=" << randomactivity << std::endl;
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
    }

    std::cout << "randomization done" << std::endl;

}


void process_line_NAME(std::string& s, std::stringstream &ss)
{

    ss << s;

    std::string NAME_str;
    std::string param_number_str;
    std::string param_name_str;
    int param_number;
    std::string param_name;

    ss >> NAME_str >> param_number_str;
    param_number = std::stoi(param_number_str);
    
    std::string::size_type found_ix = s.find(param_number_str);
    if(found_ix != std::string::npos)
    {
        param_name_str = s.substr(found_ix + param_number_str.size());
        std::string::size_type whitespace_ix = 0;
        for(;;)
        {
            if(whitespace_ix >= param_name_str.size())
            {
                break;
            }

            char whitespace_char = param_name_str.at(whitespace_ix);
            if(std::isspace(whitespace_char))
            {
                ++ whitespace_ix;
            }
            else
            {
                break;
            }
        }
        param_name = param_name_str.substr(whitespace_ix);
    }
    else
    {
        std::cout << "Error: found_ix == std::string::npos, string " << param_number_str << " not found in input" << std::endl;
    }

    //std::cout << "NAME: param_number=" << param_number << " param_name=" << param_name << std::endl;

    // set
    // use a new map to avoid conflicts
    paramNumberToHumanReadableParamNameMap.insert(std::make_pair(param_number, TString(param_name)));

}



// paramName is an output
// this value is constructed from all the mc name strings
// by appending them separated by comma ','
// the global variables
// MCNameToParamNameMap
// MCNameToParamNumberMap
// are set by this function
void process_unique_mc_names(
    std::string& paramName,
    std::vector<std::string> &unique_mc_names_ret,
    std::stringstream &ss,
    const int paramNumber)
{

    // read in non fixed width entries
    // marked by END
    // construct the parameter name
    std::vector<std::string> unique_mc_names;
    for(;;)
    {

        std::string nextstring;
        ss >> nextstring;


        //std::cout << "nextstring=" << nextstring << std::endl;
        //std::cin.get(); 
        
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
                paramName += std::string(",");
            }
            paramName += nextstring;
            //std::cout << "paramName=" << paramName << std::endl;

            // add MC to list of MC
            //MCNameToParamNameMap[paramNumber].push_back(nextstring);
            unique_mc_names.push_back(nextstring);

            // log the unique mc names in the MC name to human MC name
            // map
            //MCNameToHumanMCNameMap.insert(std::make_pair(nextstring, ));
            // NOTE: don't do this here, this was already the map
            // I constructed before using the arrays of MC sample
            // names and MC human readable names
            // TODO: remove this commented out block
        }
    }


    // setup the MCNameToParamNameMap
    for(int i = 0; i < unique_mc_names.size(); ++ i)
    {
        MCNameToParamNameMap[unique_mc_names.at(i)] = paramName;
        MCNameToParamNumberMap[unique_mc_names.at(i)] = paramNumber;
    }


    unique_mc_names_ret = unique_mc_names;
}



void process_line_else(std::string& s, std::stringstream &ss)
{

    std::string paramName = "";
    std::string paramNumber_str = "";
    int paramNumber = -1;
    std::string paramEnabled_str = "";

    // reset
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

    // moved from below, required in below block
    paramNumber = std::stoi(paramNumber_str);

    std::vector<std::string> unique_mc_names;
    process_unique_mc_names(paramName, unique_mc_names, ss, paramNumber);
    //process_unique_mc_names(paramName, /*unique_mc_names,*/ ss, paramNumber);






    // phase 1 initial values
    // read value from file, no other option
    //std::cout << "paramInitValueP1_str=" << paramInitValueP1_str << "END" << std::endl;
    //std::cout << paramInitValueP1_str << std::endl;
    //std::cout << paramInitValueP2_str << std::endl;
    paramInitValueP1 = std::stod(paramInitValueP1_str);
    paramInitErrorP1 = std::stod(paramInitErrorP1_str);

    // phase 2 intial values
    // check if string is "same", if so use values from phase 1
    // else read value from file
    function_A("same", paramInitValueP1, paramInitValueP2_str, paramInitValueP2);
    function_A("same", paramInitErrorP1, paramInitErrorP2_str, paramInitErrorP2);

    // now check if initial value strs are "useconstraint"
    // and copy values from constraints
    // also need to check for "same" again, for P2 init value and error
    // (this was an old comment - update TODO)

    // phase 1 constraints
    // check if string is "useinit", if so use values from phase 1
    // as constraints for phase 1
    // else read value from file
    function_A("useinit", paramInitValueP1, paramConstraintValueP1_str, paramConstraintValueP1);
    function_A("useinit", paramInitErrorP1, paramConstraintErrorP1_str, paramConstraintErrorP1);

    // phase 2 initial constraints
    // check if string is "same", if so use constraints from phase 1
    // else read value from file
    //if(activity_P2_str == "same")
    //useinit should be same
    function_A("same", paramConstraintValueP1, paramConstraintValueP2_str, paramConstraintValueP2);
    function_A("same", paramConstraintErrorP1, paramConstraintErrorP2_str, paramConstraintErrorP2);

    //paramNumber = std::stoi(paramNumber_str);

    // constrain mode
    // Phase 1
    function_paramConstrainMode(paramConstrainModeP1_str, paramConstrainModeP1);
    // Phase 2
    function_paramConstrainMode(paramConstrainModeP2_str, paramConstrainModeP2);



    bool is_enabled = false;
    if(paramEnabled_str == "enabled")
    {
        //auto count = std::count(enabled_params.begin(), enabled_params.end(), paramNumber);
        //if(count > 0)
        if(std::find(enabled_params.begin(), enabled_params.end(), paramNumber) != enabled_params.end())
        //if(enabled_params.find(paramNumber) != enabled_params.end())
        // TODO: change all std::find to enabled_params.find()
        {
            std::cout << "ERROR: enabled_params already contains paramNumber=" << paramNumber << std::endl;
            std::cout << "HALT" << std::endl;
            std::cin.get();
        }
        else
        {
            is_enabled = true; // TODO: can remove this and replace
                               // code which uses it below with the
                               // usual std::find call
            enabled_params.push_back(paramNumber);
            ++ numberEnabledParams; // TODO: can just set this at
                                    // end by calling enabled_params.size()
            // moved into block (below) where other map data is set
            //paramNumberToMinuitParamNumberMap[paramNumber] = minuitParamNumberCounter;
            //++ minuitParamNumberCounter;
            // could also have implemented this using .size() - 1
            // or using variable numberEnabledParams
            // but I prefered this method
            // NOTE: enabled_params contains essentially the same
            // data as paramNumberToMinuitParamNumberMap
            // TODO: remove one
        }
    }
    else if(paramEnabled_str == "disabled")
    {
        // do nothing
        disabled_params.push_back(paramNumber);
    }
    else
    {
        std::cout << "ERROR: Unrecognized value: paramEnabled_str=" << paramEnabled_str << std::endl;
        //std::cout << "unknown enabled/disabled parameter specification" << std::endl;
        //std::cout << "enabled_str=" << enabled_str << std::endl;
        //continue;
        return;
    }

    //if(paramEnabled_str == "enabled")
    //{

    // TODO: consider removing later, this was added after
    // numberParams was changed to 1
    if(paramNumber < numberParams)
    {
        if(is_enabled == true)
        {
            paramNumberToMinuitParamNumberMap[paramNumber] = minuitParamNumberCounter;
            minuitParamNumberToParamNumberMap[minuitParamNumberCounter] = paramNumber;
            ++ minuitParamNumberCounter;
        }

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

        // fixed_params set here
        if(is_enabled)
        {
            if(thePhase == 0)
            {
                if(paramConstrainModeP1 == MODE_PARAM_HARD)
                {
                    fixed_params.push_back(paramNumber);
                }
                else
                {
                    // mode is either SOFT or FREE
                    free_params.push_back(paramNumber);
                }
            }
            else if(thePhase == 1)
            {
                if(paramConstrainModeP2 == MODE_PARAM_HARD)
                {
                    fixed_params.push_back(paramNumber);
                }
                else
                {
                    // mode is either SOFT or FREE
                    free_params.push_back(paramNumber);
                }
            }
            else
            {
                std::cout << "ERROR: thePhase=" << thePhase << " invalid value" << std::endl;
            }
        }

        paramNameToNumberMap[paramName] = paramNumber;

        paramMCList[paramNumber] = unique_mc_names;

    }

    //}


    ///////////////////////////////////////////////////////////////////////////
    // write back out - debugging / useful check for user to see if his values
    // have been read by the program correctly
    ///////////////////////////////////////////////////////////////////////////

    /*
    std::cout << "values read:" << std::endl;
    std::cout << paramInitValueP1 << " +- " << paramInitErrorP1 << std::endl;
    std::cout << paramInitValueP2 << " +- " << paramInitErrorP2 << std::endl;
    std::cout << paramConstraintValueP1 << " +- " << paramConstraintErrorP1 << std::endl;
    std::cout << paramConstraintValueP2 << " +- " << paramConstraintErrorP2 << std::endl;
    std::cout << "mode: " << paramConstrainModeP1 << ", " << paramConstrainModeP2 << std::endl;
    std::cout << "check these numbers, waiting..." << std::endl;
    std::cin.get();
    */

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

void process_MCSampleNameToHumanReadableMCSampleNameMap(
    const int nSampleBkgs,
    const TString *const SampleBkgFiles,
    const TString *const SampleBkgNames)
{

    for(int i = 0; i < nSampleBkgs; ++ i)
    {
        //paramNameToHumanReadableParamNameMap.insert(std::make_pair(ExternalBkgFiles[i], ExternalBkgNames[i]));
        MCSampleNameToHumanReadableMCSampleNameMap.insert(std::make_pair(SampleBkgFiles[i], SampleBkgNames[i]));
    }
}


void read_parameter_list_file()
{


    clear_maps();



    ///////////////////////////////////////////////////////////////////////////
    // parameter_names.lst
    //
    // loading of data from parameters file
    ///////////////////////////////////////////////////////////////////////////

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
    std::size_t line_count = 1;
    while(!paramFile.eof())
    {
        std::cout << line_count << std::endl;
        

        std::stringstream ss;
        std::string s;
        std::getline(paramFile, s);
        ++ line_count;

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
        else if((s.size() >= 5) && (s[0] == 'B') && (s[1] == 'R') && (s[2] == 'E') && (s[3] == 'A') && (s[4] == 'K'))
        {
            // this is to fix a bug where lines containing whitespace which
            // isn't visible in the file
            continue;
        }
        else if((s.size() >= 4) && (s.substr(0, 4) == std::string("NAME")))
        {
            process_line_NAME(s, ss);
        }
        //if(s.size() > 0)
        else
        {
            process_line_else(s, ss);
            // there is a fail event in this function that calls return
            // same as block of code here calling continue
            // because there is no code which follows this function call
            // NOTE: if code is inserted below this function call, may need
            // to add return statement and check value, then call continue
        }

    }

    std::cout << "read: parameter_names.lst -> done" << std::endl;


    ///////////////////////////////////////////////////////////////////////////
    // end of loading data from parameter list file
    ///////////////////////////////////////////////////////////////////////////




    ///////////////////////////////////////////////////////////////////////////
    // initialize human readable name map

    
    // first: initialize map to convert from MC Sample Name to
    // Human Readable MC Sample Name


    //paramNameToHumanReadableParamNameMap.clear();
    // does not work
    // these objects are the individual MC samples and corresponding human
    // readable names
    // NOTE: can be made to work with 2 maps
    

    /*
    for(int i = 0; i < nExternalBkgs; ++ i)
    {
        //paramNameToHumanReadableParamNameMap.insert(std::make_pair(ExternalBkgFiles[i], ExternalBkgNames[i]));
        MCSampleNameToHumanReadableMCSampleNameMap.insert(std::make_pair(ExternalBkgFiles[i], ExternalBkgNames[i]));
    }

    for(int i = 0; i < nInternalBkgs; ++ i)
    {
        //paramNameToHumanReadableParamNameMap.insert(std::make_pair(InternalBkgFiles[i], InternalBkgNames[i]));
        MCSampleNameToHumanReadableMCSampleNameMap.insert(std::make_pair(InternalBkgFiles[i], InternalBkgNames[i]));
    }
    
    for(int i = 0; i < nRn222Bkgs; ++ i)
    {
        //paramNameToHumanReadableParamNameMap.insert(std::make_pair(Rn222BkgFiles[i], Rn222BkgNames[i]));
        //MCSampleNameToHumanReadableMCSampleNameMap.insert(std::make_pair(Rn222BkgFiles[i], Rn222BkgNames[i]));
        MCSampleNameToHumanReadableMCSampleNameMap.insert(std::make_pair(Rn222BkgFilesNew[i], Rn222BkgNames[i]));
    }
    
    for(int i = 0; i < nRn220Bkgs; ++ i)
    {
        //paramNameToHumanReadableParamNameMap.insert(make_pair(Rn220BkgFiles[i], Rn220BkgNames[i]));
        MCSampleNameToHumanReadableMCSampleNameMap.insert(make_pair(Rn220BkgFiles[i], Rn220BkgNames[i]));
    }
    
    for(int i = 0; i < nNeighbours; ++ i)
    {
        //paramNameToHumanReadableParamNameMap.insert(make_pair(NeighbourFiles[i], NeighbourNames[i]));
        MCSampleNameToHumanReadableMCSampleNameMap.insert(make_pair(NeighbourFiles[i], NeighbourNames[i]));
    }
    
    for(int i = 0; i < nNd150Samples; ++ i)
    {
        //paramNameToHumanReadableParamNameMap.insert(make_pair(Nd150Files[i], Nd150Names[i]));
        MCSampleNameToHumanReadableMCSampleNameMap.insert(make_pair(Nd150Files[i], Nd150Names[i]));
    }
    */

    // regular parameters
    process_MCSampleNameToHumanReadableMCSampleNameMap(nExternalBkgs, ExternalBkgFiles, ExternalBkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nInternalBkgs, InternalBkgFiles, InternalBkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nRn222Bkgs, Rn222BkgFilesNew, Rn222BkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nRn220Bkgs, Rn220BkgFiles, Rn220BkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nNeighbours, NeighbourFiles, NeighbourNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nNd150Samples, Nd150Files, Nd150Names);
    // add additional parameters
    MCSampleNameToHumanReadableMCSampleNameMap.insert(make_pair("axial_vector_parameter_0", "^{150}Nd 2#nu#beta#beta g_A #xi_{31}"));

    // done

    // second: convert parameter numbers to parameter names and create map
    // to convert between parameter names and human readable parameter names

    // NOTE: this block of code depends on the map
    // MCSampleNameToHumanReadableMCSampleNameMap
    // containing all the MC Sample Names which appear in the parameter list
    // file.
    // Because some parameters MC samples may have been disabled when pre-
    // processing the data, and are therefore set to "disabled" in the
    // parameter list file, need to check that parameter is enabled before
    // processing any names in this code block

    for(int param_number = 0; param_number < numberParams; ++ param_number)
    {
        
        // check if parameter is enabled
        if(std::find(enabled_params.begin(), enabled_params.end(), param_number) != enabled_params.end())
        {
            // do nothing
        }
        else
        {
            // skip
            continue;
        }

        std::cout << "param_number=" << param_number << std::endl;
        std::string parameter_name = paramNameMap[param_number];
        std::string human_readable_parameter_name;

        std::cout << "parameter_name=" << parameter_name << std::endl;

        std::string::size_type find_ix{0};
        for(;;)
        {
            // split name by comma
            // find corresponding MC sample name and corresponding MC sample
            // human name
            std::string::size_type found_ix = parameter_name.find(",", find_ix);
            //std::cout << "found ',' at found_ix=" << found_ix << std::endl;


            std::string mc_sample_name;
            //if(found_ix == std::string::npos)
            //{
                mc_sample_name = parameter_name.substr(find_ix, found_ix - find_ix);
            //}
            //else
            //{
            //    // TODO: this if is not necessary?
            //    mc_sample_name = parameter_name.substr(find_ix, std::string::npos);
            //}

            std::cout << "mc_sample_name=" << mc_sample_name << std::endl;

            //std::cout << "mc_sample_name=" << mc_sample_name << std::endl;
            std::string human_readable_mc_sample_name = std::string(MCSampleNameToHumanReadableMCSampleNameMap.at(mc_sample_name).Data());
            //std::cout << "human_readable_mc_sample_name=" << human_readable_mc_sample_name << std::endl;
         
            // construct human readable parameter name
            if(human_readable_parameter_name.size() > 0)
            {
                human_readable_parameter_name += std::string(",");
            }
            human_readable_parameter_name += human_readable_mc_sample_name;

            //TString human_name = paramMCNameToHumanMCNameMap.at(name);
            // add human names by comma


            // TODO: in case of no further comma, last string is until end
            if(found_ix != std::string::npos)
            {
                // do nothing, code moved above
                // NOTE: may have solved above TODO, check
                find_ix = found_ix + 1;
            }
            else
            {
                //std::cout << "break" << std::endl;
                break;
            }

        }

        std::cout << "human_readable_parameter_name=" << human_readable_parameter_name << std::endl;

        paramNameToHumanReadableParamNameMap.insert(std::make_pair(parameter_name, human_readable_parameter_name));
    }

    // human readable name map
    ///////////////////////////////////////////////////////////////////////////

    //std::cin.get();

    // Some notes after making the above code work.
    // The above code is: The code that builds the names for the maps used
    // in the correlation matrix labels.
    // Map name is: paramNameToHumanReadableParamNameMap
    //
    // some comments:
    // 
    // there are some blank rows/cols. First one at 12 (?)
    // looks like some disabled/fixed parameter is being included in the
    // corralation matrix where it should not be
    // alternatively these could be exactly zero entries?
    // note: these are (probably) due to zero entries
    //
    // the neighbour foils are set to HARD constrained, and yet they appear
    // in the correlation matrix
    //
    // 22 parameters appear in the correlation matrix, which is wrong
    //
    // check how Summer built her free_param_names vector
    //
    //
    //



    // randomize parameters

    // TODO: re-enable
    if(false)
    {
        randomize_parameters();
    }
    // randomize parameters


    // record enabled/disabled fixed/free parameters
    //
    //
    
    
    
    // Sort out fixed vs free params (this can be improved, its too hackish now)

    // data stored here used as name for each item in correlation matrix
    // TODO other names should be added not just paramNameMap[i].front()
    // loop over paramNameMap[i], adding each name to the free_param_names
    // vector
    
    // moved to header
    //std::vector<TString> free_param_names;
    //std::vector<TString> enabled_param_names;
    #if 0
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
            
            //free_params.push_back(i);
            // moved into loop which reads parameter list file

            // TODO: what to do with these? are they used?
            index_free_params.push_back(i_str);
            // TODO: .front() is bug?
            ////free_param_names.push_back(paramNameMap[i].front());
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
            // TODO: .front() is bug?
            ////enabled_param_names.push_back(paramNameMap[i].front());
        }
        else
        {
            //disabled_params.push_back(i);
            // moved into loop which reads parameter list file

            //disabled_param_names.push_back(paramNameMap[i].front());
        }
        // TODO: this may no longer be used
        // TODO: can probably remove? (check)
    }
    #endif

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
    // TODO: this should probably be moved to just after the section which
    // reads parameter list


    // 2020-06-17
    // set last values map
    for(int i = 0; i < numberParams; ++ i)
    {
        if(thePhase == 0)
        {
            paramLastValueMap[i] = paramInitValueP1Map[i];
            paramInitValueMap[i] = paramInitValueP1Map[i];

        }
        else if(thePhase == 1)
        {
            paramLastValueMap[i] = paramInitValueP2Map[i];
            paramInitValueMap[i] = paramInitValueP2Map[i];
        }
        else
        {
            std::cerr << "Error: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
            throw "problem";
        }

        int axial_vector_parameter_0_param_number = get_axial_vector_parameter_index(); 
        if(i == axial_vector_parameter_0_param_number)
        {
            paramLastValueMap[i] = 0.0; // HSD, forces reweighting
        }
    }
    // set initial values map
    /*
    for(int i = 0; i < numberParams; ++ i)
    {
        if(thePhase == 0)
        {
            paramInitValueMap[i] = paramInitValueP1Map[i];
        }
        else if(thePhase == 1)
        {
            paramInitValueMap[i] = paramInitValueP2Map[i];
        }
    }
    */
    // TODO: implement this doing away with if statements to select phase
    // in other parts of code


}



#endif // NEWLOGLIKFITTER_READ_PARAMETERNAMES_LST_H
