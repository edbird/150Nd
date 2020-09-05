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
int function_A(const std::string& string_compare,
                const double input,
                const std::string& string_input,
                double& output)
{
    int ret = 0;

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
            ret = -1;
        }
    }

    return ret;
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


int function_paramConstraintMode(const std::string& paramConstraintMode_str, int &paramConstraintMode)
{
    int ret = 0;

    // Phase 1 / Phase 2 Constran Mode
    if(paramConstraintMode_str == "free")
    {
        paramConstraintMode = MODE_CONSTRAINT_FREE;
    }
    else if(paramConstraintMode_str == "soft")
    {
        paramConstraintMode = MODE_CONSTRAINT_SOFT;
    }
    else if(paramConstraintMode_str == "hard")
    {
        paramConstraintMode = MODE_CONSTRAINT_HARD;
    }
    else
    {
        std::cout << "ERROR: Unrecognized value: paramConstraintMode_str=" << paramConstraintMode_str << std::endl;
        paramConstraintMode = MODE_CONSTRAINT_UNDEFINED;
        ret = -1;
    }

    return ret;
}

/*
void function_paramLockMode(const std::string &paramLockMode_str, int &paramLockMode)
{
    if(paramLockMode_str == "locked")
    {
        paramLockMode = MODE_PARAM_LOCKED;
    }
    else if(paramLockMode_str == "unlocked")
    {
        paramLockMode = MODE_PARAM_UNLOCKED;
    }
    else
    {
        std::cout << "ERROR: Unrecognized value: paramLockMode_str=" << paramLockMode_str << std::endl;
        paramLockMode = MODE_PARAM_LOCK_UNDEFINED;
    }
}
*/






/*
void clear_maps()
{

    //numberParams = 0;
    numberEnabledParams = 0;
    minuitParamNumberCounter = 0;

    MCNameToParamNameMap.clear();
    MCNameToParamNumberMap.clear();

    paramNameToNumberMap.clear();

    paramNumberToMinuitParamNumberMap.clear();
    minuitParamNumberToParamNumberMap.clear();
//    minuitParamNumberCounter = 0;

//    numberEnabledParams = 0;

    fixed_params.clear();
    free_params.clear();

    //paramNameToHumanReadableParamNameMap.clear();
    MCSampleNameToHumanReadableMCSampleNameMap.clear();



    enabled_params.clear();
    disabled_params.clear();

}
*/


/*
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
*/


int process_line_NAME(std::string& s, std::stringstream &ss)
{
    int ret = 0;

    ss << s;

    std::string NAME_str; // contains "NAME"
    std::string paramNumber_str;
    std::string paramHumanReadableName_str;
    int paramNumber;
    std::string paramHumanReadableName;

    ss >> NAME_str >> paramNumber_str;
    paramNumber = std::stoi(paramNumber_str);
    
    std::string::size_type found_ix = s.find(paramNumber_str);
    if(found_ix != std::string::npos)
    {
        paramHumanReadableName_str = s.substr(found_ix + paramNumber_str.size());
        std::string::size_type whitespace_ix = 0;
        for(;;)
        {
            if(whitespace_ix >= paramHumanReadableName_str.size())
            {
                break;
            }

            char whitespace_char = paramHumanReadableName_str.at(whitespace_ix);
            if(std::isspace(whitespace_char))
            {
                ++ whitespace_ix;
            }
            else
            {
                break;
            }
        }
        paramHumanReadableName = paramHumanReadableName_str.substr(whitespace_ix);
    }
    else
    {
        std::cout << "Error: found_ix == std::string::npos, string " << paramNumber_str << " not found in input" << std::endl;
        ret = -1;
    }

    //std::cout << "NAME: param_number=" << param_number << " param_name=" << param_name << std::endl;

    // set
    // use a new map to avoid conflicts
    //paramNumberToHumanReadableParamNameMap.insert(std::make_pair(param_number, TString(param_name)));

//working here
    std::map<int, file_parameter>::iterator it = g_pg.file_params.find(paramNumber);
    if(it != g_pg.file_params.end())
    {
        it->second.paramHumanReadableName = paramHumanReadableName;
    }
    else
    {
        std::cout << __func__ << " " << " Error: could not find parameter in g_pg with paramNumber=" << paramNumber << std::endl;
        ret = -1;
    }

    return ret;
}



// paramName is an output
// this value is constructed from all the mc name strings
// by appending them separated by comma ','
// the global variables
// MCNameToParamNameMap
// MCNameToParamNumberMap
// are set by this function
int process_unique_mc_names(
    std::string& paramName,
    std::vector<std::string> &unique_mc_names_ret,
    std::stringstream &ss,
    const int paramNumber) // paramNumber is external format
{
    int ret = 0;

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
        else if(nextstring == "NOMC")
        {
            // do nothing
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
        //MCNameToParamNameMap[unique_mc_names.at(i)] = paramName;
        //MCNameToParamNumberMap[unique_mc_names.at(i)] = paramNumber;
        g_pg.MCNameToExtParamNumberMap[unique_mc_names.at(i)] = paramNumber;
    }


    unique_mc_names_ret = unique_mc_names;

    return ret;
}


int stringToBool(std::string &valuestring, bool &output, const std::string& value_if_true, const std::string& value_if_false)
{
    int ret = 0;
    if(valuestring == value_if_true)
    {
        output = true;
    }
    else if(valuestring == value_if_false)
    {
        output = false;
    }
    else
    {
        std::cout << "Error: " << valuestring << " not \"" << value_if_true << "\" or \"" << value_if_false << "\"" << std::endl;
        ret = -1;
    }
    return ret;
}


void stringToBool(std::string &valuestring, bool &output)
{
    int ret = 0;
    if(valuestring == "true")
    {
        output = true;
    }
    else if(valuestring == "false")
    {
        output = false;
    }
    else
    {
        std::cout << "Error: " << valuestring << " not \"true\" or \"false\"" << std::endl;
        ret = -1;
    }
    return ret;
}


int process_line_else(std::string& s, std::stringstream &ss)
{
    int ret = 0;

    std::string paramNumber_str = "";
    std::string paramEnabled_str = "";
    std::string paramEnabledP1_str = "";
    std::string paramEnabledP2_str = "";

    std::string paramInitValue_str;
    std::string paramInitError_str;

    std::string paramConstraintValue_str;
    std::string paramConstraintError_str;

    std::string paramConstraintMode_str;

    std::string paramStackType_str;
    
    std::string paramName_str;

    std::string paramName = "";

    int paramNumber = -1;
    bool paramEnabled = false;
    bool paramEnabledP1 = false;
    bool paramEnabledP2 = false;

    double paramInitValue = 0.;
    double paramInitError = 0.;

    double paramConstraintValue = 0.;
    double paramConstraintError = 0.;

    int paramConstraintMode = MODE_CONSTRAINT_UNDEFINED;

    int paramStackType = STACK_TYPE_DEFAULT;

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
    try
    {
        ss >> paramNumber_str
           >> paramEnabled_str
           >> paramEnabledP1_str
           >> paramEnabledP2_str
           >> paramInitValue_str >> paramInitError_str
           >> paramConstraintValue_str >> paramConstraintError_str
           >> paramConstraintMode_str
           >> paramStackType_str
           >> paramName_str; 
    }
    catch(...)
    {
        std::cout << "ERROR: reading strings" << std::endl;
        ret = -1;
    }

    //std::cout << "read paramNumber_str=" << paramNumber_str << std::endl;

    // moved from below, required in below block
    try
    {
        paramNumber = std::stoi(paramNumber_str);
    }
    catch(...)
    {
        std::cout << "ERROR: paramNumber_str=" << paramNumber_str << std::endl;
        ret = -1;
    }

    if(g_pg.ParamNameToExtParamNumberMap.count(paramName_str) == 0)
    {
        g_pg.ParamNameToExtParamNumberMap[paramName_str] = paramNumber;
    }
    else
    {
        std::cout << "ERROR: g_pg.ParamNameToExtParamNumberMap already contains paramName_str=" << paramName_str << std::endl;
        ret = -1;
    }

    std::vector<std::string> unique_mc_names;
    if(process_unique_mc_names(paramName, unique_mc_names, ss, paramNumber) != 0)
    {
        ret = -1;
    }

    if(stringToBool(paramEnabled_str, paramEnabled, "enabled", "disabled") != 0)
    {
        ret = -1;
    }

    if(stringToBool(paramEnabledP1_str, paramEnabledP1, "P1_enabled", "P1_disabled") != 0)
    {
        ret = -1;
    }

    if(stringToBool(paramEnabledP2_str, paramEnabledP2, "P2_enabled", "P2_disabled") != 0)
    {
        ret = -1;
    }

    // phase initial values
    // read value from file, no other option
    try
    {
        paramInitValue = std::stod(paramInitValue_str);
    }
    catch(...)
    {
        std::cout << "ERROR: paramInitValue_str=" << paramInitValue_str << std::endl;
        ret = -1;
    }

    try
    {
        paramInitError = std::stod(paramInitError_str);
    }
    catch(...)
    {
        std::cout << "ERROR: paramInitError_str=" << paramInitError_str << std::endl;
        ret = -1;
    }

    // contraints
    // check if string is "useinit", if so use values from initial value/error
    // else read value from file
    try
    {
        function_A("useinit", paramInitValue, paramConstraintValue_str, paramConstraintValue);
    }
    catch(...)
    {
        ret = -1;
    }

    try
    {
        function_A("useinit", paramInitError, paramConstraintError_str, paramConstraintError);
    }
    catch(...)
    {
        ret = -1;
    }

    // now check if initial value strs are "useconstraint"
    // and copy values from constraints
    // also need to check for "same" again, for P2 init value and error
    // (this was an old comment - update TODO)

    // constrain mode
    if(function_paramConstraintMode(paramConstraintMode_str, paramConstraintMode) != 0)
    {
        ret = -1;
    }

    if(paramStackType_str == "ND150")
    {
        paramStackType = STACK_TYPE_ND150;
    }
    else if(paramStackType_str == "INTERNAL")
    {
        paramStackType = STACK_TYPE_INTERNAL;
    }
    else if(paramStackType_str == "EXTERNAL")
    {
        paramStackType = STACK_TYPE_EXTERNAL;
    }
    else if(paramStackType_str == "NEIGHBOUR")
    {
        paramStackType = STACK_TYPE_NEIGHBOUR;
    }
    else if(paramStackType_str == "RN222")
    {
        paramStackType = STACK_TYPE_RN222;
    }
    else if(paramStackType_str == "RN220")
    {
        paramStackType = STACK_TYPE_RN220;
    }
    else if(paramStackType_str == "NA")
    {
        paramStackType = STACK_TYPE_DEFAULT;
        // no error condition
        // this is the option for the xi_31 parameter
    }
    else
    {
        paramStackType = STACK_TYPE_DEFAULT;
        std::cout << "ERROR: STACK TYPE not recognized" << std::endl;
        ret = -1;
    }


/*
    if(paramEnabled == true)
    {
*/
        std::map<int, file_parameter>::iterator it = g_pg.file_params.find(paramNumber);
        if(it != g_pg.file_params.end())
        {
            std::cout << "ERROR: g_pg.file_params already contains paramNumber=" << paramNumber << std::endl;
            std::cout << "HALT" << std::endl;
            std::cin.get();
            ret = -1;
        }
        else
        {
            std::vector<std::string> MCNameList = unique_mc_names;
            std::string paramHumanReadableName = "UNNAMED PARAM";

            file_parameter fp;
            fp.paramNumber = paramNumber;
            fp.paramEnabled = paramEnabled;
            fp.paramEnabledP1 = paramEnabledP1;
            fp.paramEnabledP2 = paramEnabledP2;
            fp.paramInitValue = paramInitValue;
            fp.paramInitError = paramInitError;
            fp.paramConstraintValue = paramConstraintValue;
            fp.paramConstraintError = paramConstraintError;
            fp.paramConstraintMode = paramConstraintMode;
            fp.paramStackType = paramStackType;
            //fp.paramName = paramName;
            fp.paramName = paramName_str;
            fp.paramHumanReadableName = paramHumanReadableName;
            fp.MCNameList = MCNameList;

            if(g_pg.file_params.count(paramNumber) == 0)
            {
                g_pg.file_params[paramNumber] = fp;
            }
            else
            {
                std::cout << "ERROR: paramNumber=" << paramNumber << " already exists in g_pg.file_params" << std::endl;
                ret = -1;
            }

            if(paramEnabled == true)
            {
                ++ g_pg.numberEnabledParams;
            }
        }

        // this adds the parameter to the parameter group regardless of whether
        // it is enabled or not
        // so need to ignore disabled parameters later in the code

/*
    }
    else if(paramEnabled == false)
    {
        // do nothing
    }
*/
    //if(paramEnabled_str == "enabled")
    //{

#if 0
    // TODO: consider removing later, this was added after
    // numberParams was changed to 1
    if(paramNumber < g_pg.numberParams())
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

        paramConstraintModeP1Map[paramNumber] = paramConstraintModeP1;
        paramConstraintModeP2Map[paramNumber] = paramConstraintModeP2;

        paramLockModeMap[paramNumber] = paramLockMode;

        // fixed_params set here
        if(is_enabled)
        {
            /*
            if(thePhase == 0)
            {
                if(paramConstraintModeP1 == MODE_PARAM_HARD)
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
                if(paramConstraintModeP2 == MODE_PARAM_HARD)
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
            */

            // fixed params set here
            // behaviour is different depending on lock mode
            if(paramLockMode == MODE_PARAM_LOCKED)
            {
                // in locked mode, behaviour is dependent on
                // which parameter has the most severe constraint mode
                if((paramConstraintModeP1 == MODE_PARAM_HARD) ||
                   (paramConstraintModeP2 == MODE_PARAM_HARD))
                {
                    // either parameter or both are set to HARD
                    fixed_params_P1.push_back(paramNumber);
                    fixed_params_P2.push_back(paramNumber); // TODO: can I get rid of fixed_params ?
                }
                else
                {
                    // parameters are one of 4 options:
                    // soft soft
                    // soft free
                    // free soft
                    // free free
                    // in any of these cases they are considered to be free
                    free_params_P1.push_back(paramNumber);
                    free_params_P2.push_back(paramNumber);
                }
            }
            else if(paramLockMode == MODE_PARAM_UNLOCKED)
            {
                // in unlocked mode, proceed as before

                // code copied from above
                //if(thePhase == 0)
                //{
                    if(paramConstraintModeP1 == MODE_PARAM_HARD)
                    {
                        fixed_params_P1.push_back(paramNumber);
                    }
                    else
                    {
                        // mode is either SOFT or FREE
                        free_params_P1.push_back(paramNumber);
                    }
                //}
                //else if(thePhase == 1)
                //{
                    if(paramConstraintModeP2 == MODE_PARAM_HARD)
                    {
                        fixed_params_P2.push_back(paramNumber);
                    }
                    else
                    {
                        // mode is either SOFT or FREE
                        free_params_P2.push_back(paramNumber);
                    }
                //}
                //else
                //{
                //    std::cout << "ERROR: thePhase=" << thePhase << " invalid value" << std::endl;
                //}
                // end of code copied
            }
            else
            {
                std::cout << "ERROR: Unrecognized value: paramLockMode=" << paramLockMode << std::endl;
                throw "problem";
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
    std::cout << "mode: " << paramConstraintModeP1 << ", " << paramConstraintModeP2 << std::endl;
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
    //          << paramConstraintModeP1 << ", " << paramConstraintModeP2 << ", ";
    // print out the MCNameToParamNameMap
    //for(auto it = MCNameToParamNameMap.cbegin(); it != MCNameToParamNameMap.cend(); ++ it)
    //{
    //    std::cout << it->first << " -> " << it->second << std::endl;
    //}
#endif


#if 0

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
    std::cout << "         Mode: " << paramConstraintModeP1Map[paramNumber];
    if(paramConstraintModeP1Map[paramNumber] == MODE_PARAM_FREE)
    {
        std::cout << " (free)";
    }
    else if(paramConstraintModeP1Map[paramNumber] == MODE_PARAM_SOFT)
    {
        std::cout << " (soft)";
    }
    else if(paramConstraintModeP1Map[paramNumber] == MODE_PARAM_HARD)
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
    std::cout << "         Mode: " << paramConstraintModeP2Map[paramNumber];
    if(paramConstraintModeP2Map[paramNumber] == MODE_PARAM_FREE)
    {
        std::cout << " (free)";
    }
    else if(paramConstraintModeP2Map[paramNumber] == MODE_PARAM_SOFT)
    {
        std::cout << " (soft)";
    }
    else if(paramConstraintModeP2Map[paramNumber] == MODE_PARAM_HARD)
    {
        std::cout << " (hard)";
    }
    else
    {
        std::cout << " (unknown/ERROR)";
    }
    std::cout << std::endl;
    std::cout << "Parameter Lock Mode: " << paramLockModeMap[paramNumber] << std::endl;
    std::cout << std::endl;
#endif

    return ret;
}


/*
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
*/


void read_parameter_list_file()
{


    //g_pg.clear_maps();
    g_pg.reset();



    ///////////////////////////////////////////////////////////////////////////
    // parameter_names.lst
    //
    // loading of data from parameters file
    ///////////////////////////////////////////////////////////////////////////

    std::ifstream paramFile;
    std::string paramFile_filename = "parameter_names.lst";
    paramFile.open(paramFile_filename);
    if(paramFile.is_open() == false)
    {
        std::cerr << __func__ << " Error: Cannot open file" << "parameter_names.lst" << std::endl;
        throw "error";
    }
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
        //std::cout << line_count << std::endl;
        

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
            if(process_line_NAME(s, ss) != 0)
            {
                std::cout << "ERROR: Syntax error in file " << paramFile_filename << " line " << line_count << std::endl;
            }
        }
        //if(s.size() > 0)
        else
        {
            if(process_line_else(s, ss) != 0)
            {
                std::cout << "ERROR: Syntax error in file " << paramFile_filename << " line " << line_count << std::endl;
            }

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


    // TODO: I forget what these do but I seem to remember it is important

/*
    // regular parameters
    process_MCSampleNameToHumanReadableMCSampleNameMap(nExternalBkgs, ExternalBkgFiles, ExternalBkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nInternalBkgs, InternalBkgFiles, InternalBkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nRn222Bkgs, Rn222BkgFilesNew, Rn222BkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nRn220Bkgs, Rn220BkgFiles, Rn220BkgNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nNeighbours, NeighbourFiles, NeighbourNames);
    process_MCSampleNameToHumanReadableMCSampleNameMap(nNd150Samples, Nd150Files, Nd150Names);
    // add additional parameters
    MCSampleNameToHumanReadableMCSampleNameMap.insert(make_pair("axial_vector_parameter_0", "^{150}Nd 2#nu#beta#beta g_A #xi_{31}"));
*/

    // done

/*
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
*/

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
    #if 0
    if(false)
    {
        randomize_parameters();
    }
    #endif
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
#if 0
    std::cout << "List of parameters: free/fixed (P1)" << std::endl;
    for(int i = 0; i < free_params_P1.size(); ++ i)
    {
        std::cout << "free parameter (P1): " << free_params_P1[i] << std::endl;
    }
    for(int i = 0; i < fixed_params_P1.size(); ++ i)
    {
        std::cout << "fixed parameter (P1): " << fixed_params_P1[i] << std::endl;
    }
    std::cout << "List of parameters: free/fixed (P2)" << std::endl;
    for(int i = 0; i < free_params_P2.size(); ++ i)
    {
        std::cout << "free parameter (P2): " << free_params_P2[i] << std::endl;
    }
    for(int i = 0; i < fixed_params_P2.size(); ++ i)
    {
        std::cout << "fixed parameter (P2): " << fixed_params_P2[i] << std::endl;
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
#endif
    // TODO: re-enable this in parameter list class
    // TODO: this should probably be moved to just after the section which
    // reads parameter list

    // print a list of locked parameters and their modes
    //std::cout << "List of parameter lock modes: locked/unlocked" << std::endl;
    //for(int i = 0; i < 

    g_pg.init_maps();
    g_pg.set_last_parameter_values();
#if 0
    // 2020-06-17
    // set last values map
    for(int i = 0; i < numberParams; ++ i)
    {
        /*
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
        */

        paramLastValueP1Map[i] = paramInitValueP1Map[i]; 
        paramInitValueP1Map[i] = paramInitValueP1Map[i]; 

        paramLastValueP2Map[i] = paramInitValueP2Map[i]; 
        paramInitValueP2Map[i] = paramInitValueP2Map[i]; 

        int axial_vector_parameter_0_param_number = get_axial_vector_parameter_index(); 
        if(i == axial_vector_parameter_0_param_number)
        {
            //paramLastValueMap[i] = 0.0; // HSD, forces reweighting
            paramLastValueP1Map[i] = 0.0; // HSD, forces reweighting
            paramLastValueP2Map[i] = 0.0; // HSD, forces reweighting
        }
    }
#endif
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
