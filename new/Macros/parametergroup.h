#ifndef PARAMETERGROUP_H
#define PARAMETERGROUP_H



// Summers (retarded) method of accessing the MC histograms
//
// 1: allMCSamples1D[channel] = new TObjArray()
//
//    tmpHist_P1 = (TH1D*)aFile_P1->Get(name)->Clone(newname)
// 2: allMCSamples1D[channel]->Add(tmpHist_P1)
//
// 3: reweight apply:
//    for(channel)
//        for(i in channel)
//            tmpHist = allMCSamples1D[channel]->At(i)
//            tmpHist_name = tmpHist->GetName();
//            if(tmpHist_name.Contains("nd150_rot_2n2b_m4"))
//                allMCSamples[channel]->RemoveAt(i)
//                break
// this is really fucking wank, can only remove one histogram
// at a time unless the index is reset, or arbitarily reduced
// by 1, ignoring the for loop increment
//
// I'm not such a huge fan of the scaling on loading from file
// but this is probably better than having a large range of
// orders of magnitude in the Minuit system
//
// It is not practival however to move the 208Tl scaling
// to the getNumberMC function, as the scaling must occur only
// once. Could introduce another function to do this which is
// called once, but that introduces additional complexity.
//
// Point 3 above: Essentially loop through an array of all
// MC objects and check the name matches to find the one you
// want. This is utterly SHIT code.
//
//



#include "TH1.h"
#include "TH2.h"



const int MODE_CONSTRAINT_UNDEFINED = -1;
const int MODE_CONSTRAINT_FREE = 0; // completely free parameter
const int MODE_CONSTRAINT_SOFT = 1; // soft float (constrained) using penalty term
const int MODE_CONSTRAINT_HARD = 2; // hard (fixed) constrained parameter using fixed value


const int STACK_TYPE_DEFAULT = 0;
const int STACK_TYPE_RN222 = 1;
const int STACK_TYPE_RN220 = 2; // both of these get stacked into RADON anyway
const int STACK_TYPE_EXTERNAL = 3;
const int STACK_TYPE_INTERNAL = 4;
const int STACK_TYPE_NEIGHBOUR = 5;
const int STACK_TYPE_ND150 = 6;


class file_parameter
{

    public:

    file_parameter()
    {
        // default / required parameter values
        paramNumber = -1;
        paramEnabled = false;
        paramEnabledP1 = true;
        paramEnabledP1 = true;
        paramInitValue = 0.0;
        paramInitError = 0.0;
        paramLastValue = 0.0;
        paramConstraintValue = 0.0;
        paramConstraintError = 0.0;
        paramConstraintMode = MODE_CONSTRAINT_UNDEFINED;
        // required also includes MCNameList
        // however default value for this is empty vector (no MC
        // assigned)
        MCNameList.clear();

        // not required
        paramName = "";
        // automatically constructed from concatination of MC names
        
        // not required
        // parameter human readable name
        // set to default value which can be printed as something
        // visible (not blank) but that indicates the name was
        // not set
        paramHumanReadableName = "UNNAMED";

        // stacking parameters
        // not required
        // initialized to default values which put everything in the 
        // same stack with the color black
        paramStackType = STACK_TYPE_DEFAULT;
        paramStackMajorColor = kBlack; // Note: both of these currently unused
        paramStackMinorColor = kBlack; // due to fact that colors can be obtained
                                        // from stack type for major and minor is
                                        // not used

//        histogram_dimension = 0;
    }


/*
    Double_t getBinContent1D(const Int_t bin) const
    {
        if(histogram_dimension != 1)
        {
            std::cout << "histogram_dimension != 1 : histogram_dimension=" << histogram_dimension << std::endl;
            throw "histogram_dimension!=1";
        }

        Double_t total = 0.0;
        for(std::size_t index = 0; index < hist.size(); ++ index)
        {
            TObject *obj = hist[index];
            TH1D* h = (TH1D*)obj;
            Double_t content = h->GetBinContent(bin);
            total += content;
        }
        return total;
    }


    Double_t getBinContent2D(const Int_t bin_x, const Int_t bin_y) const
    {
        if(histogram_dimension != 2)
        {
            std::cout << "histogram_dimension != 2 : histogram_dimension=" << histogram_dimension << std::endl;
            throw "histogram_dimension!=2";
        }
        
        Double_t total = 0.0;
        std::vector<TObject*>::const_iterator it{hist.cbegin()};
        for(; it != hist.cend(); ++ it)
        {
            TObject *obj = *it;
            TH2D* h = (TH2D*)obj;
            Double_t content = h->GetBinContent(bin_x, bin_y);
            total += content;
        }
        return total;
    }
*/


    // required
    int paramNumber;
    bool paramEnabled;
    bool paramEnabledP1;
    bool paramEnabledP2;
    double paramInitValue;
    double paramInitError;
    double paramLastValue;
    //double parameter_last_value;
    //double parameter_last_error;
    // minuitParamCurrent
    // minuitParamLast // this is the same as above
    // minuitParamInit // TODO: can probably remove these?
    double paramConstraintValue;
    double paramConstraintError;
    int paramConstraintMode;
    std::vector<std::string> MCNameList;

//    int histogram_dimension;
//    std::vector<TObject*> hist;

    // automatic
    std::string paramName;

    // init to default value
    std::string paramHumanReadableName;

    // init to default value (all same stack)
    int paramStackType;
    // list of corresponding human readable MC names for each MC sample
    // check the file InputNameDef.h
    //std::vector<std::string> MCHumanNameList;
    Color_t paramStackMajorColor;
    Color_t paramStackMinorColor;
    


    // list of ways in which parameters must be accessed
    //
    // loop over all parameters to add them to MINUIT using
    // Minuit::Add()
    //
    // operator() provides std::vector<double> param
    // loop over this vector to access parameters,
    // need to be able to access parameter number from MC Name

};


class parameter_group
{

    public:

    parameter_group()
    {
        //numberParams = 0;
        numberEnabledParams = 0;
        //minuitParamNumberCounter = 0;

        xi_31_int_ext_param_number_set = false;
    }


    /*
    int get_number_free_params() const
    {
        int sum = 0;
        std::map<int, file_parameter>::const_iterator it{file_params.cbegin()};
        for(; it != file_params.cend(); ++ it)
        {
            if(it->second.paramConstraintMode == MODE_CONSTRAINT_FREE)
            {
                ++ sum;
            }
        }
        std::cout << "parameter_group.get_number_free_params()=" << sum << std::endl;
        return sum;
    }
    */


    int get_xi_31_int_param_number()
    {
        if(xi_31_int_ext_param_number_set == true)
        {
            return xi_31_int_param_number;
        }
        else
        {
            std::string param_name = "xi_31";

            xi_31_ext_param_number = ParamNameToExtParamNumberMap.at(param_name);
            std::map<std::string, int>::iterator it = ParamNameToExtParamNumberMap.find(param_name);
            if(it != ParamNameToExtParamNumberMap.end())
            {
                std::map<int, int>::iterator it = ExtToIntParamNumberMap.find(xi_31_ext_param_number);
                if(it != ExtToIntParamNumberMap.end())
                {
                    xi_31_int_param_number = ExtToIntParamNumberMap.at(xi_31_ext_param_number);
                }
                else
                {
                    xi_31_int_param_number = -1;
                }
                //xi_31_int_param_number = ExtToIntParamNumberMap.at(xi_31_ext_param_number);
            }
            else
            {
                xi_31_int_param_number = -1;
            }
            //xi_31_int_param_number = ExtToIntParamNumberMap.at(xi_31_ext_param_number);
            xi_31_int_ext_param_number_set = true;

            return xi_31_int_param_number;
        }
    }


    int get_xi_31_ext_param_number()
    {
        if(xi_31_int_ext_param_number_set == true)
        {
            return xi_31_ext_param_number;
        }
        else
        {
            std::string param_name = "xi_31";

            xi_31_ext_param_number = ParamNameToExtParamNumberMap.at(param_name);

            std::map<std::string, int>::iterator it = ParamNameToExtParamNumberMap.find(param_name);
            if(it != ParamNameToExtParamNumberMap.end())
            {
                std::map<int, int>::iterator it = ExtToIntParamNumberMap.find(xi_31_ext_param_number);
                if(it != ExtToIntParamNumberMap.end())
                {
                    xi_31_int_param_number = ExtToIntParamNumberMap.at(xi_31_ext_param_number);
                }
                else
                {
                    xi_31_int_param_number = -1;
                }
            }
            else
            {
                xi_31_int_param_number = -1;
            }
            xi_31_int_ext_param_number_set = true;

            return xi_31_ext_param_number;
        }
    }


    void reset()
    {
        file_params.clear();
        ParamNameToExtParamNumberMap.clear();
        MCNameToExtParamNumberMap.clear();
        ExtToIntParamNumberMap.clear();
        IntToExtParamNumberMap.clear();
        numberEnabledParams = 0;
        //minuitParamNumberCounter = 0;
    }


    void set_last_parameter_values()
    {
        std::map<int, file_parameter>::iterator it{file_params.begin()};
        for(; it != file_params.end(); ++ it)
        {
            it->second.paramLastValue = it->second.paramInitValue;
        }
        std::cout << "done setting last values" << std::endl;

        int axial_vector_parameter_0_param_number = get_xi_31_ext_param_number();
        std::cout << "axial_vector_parameter_0_param_number=" << axial_vector_parameter_0_param_number << std::endl;
        //file_params.at(axial_vector_parameter_0_param_number).paramInitValue = 0.0; // HSD, forces reweighting
        file_params.at(axial_vector_parameter_0_param_number).paramLastValue = -std::numeric_limits<double>::infinity(); // HSD, forces reweighting to initial value set in parameter_names.lst
    }


    void init_maps()
    {
        const int debug_level = 1;

        if(debug_level >= 2)
        {
            std::cout << __func__ << std::endl;
            print_map(ExtToIntParamNumberMap, "ExtToIntParamNumberMap");
            print_map(IntToExtParamNumberMap, "IntToExtParamNumberMap");
        }

        std::map<int, file_parameter>::iterator it{file_params.begin()};
        int minuit_param_number_counter = 0;
        for(; it != file_params.end(); ++ it)
        {
            if(debug_level >= 3)
            {
                std::cout << "minuit_param_number_counter=" << minuit_param_number_counter << std::endl;
            }

            if(it->second.paramEnabled)
            {
                if(it->second.paramEnabledP1 || it->second.paramEnabledP2)
                {
                    int paramNumber = it->second.paramNumber;

                    // set ExtToInt map
                    if(ExtToIntParamNumberMap.count(paramNumber) == 0)
                    {
                        ExtToIntParamNumberMap[paramNumber] = minuit_param_number_counter;
                    }
                    else
                    {
                        std::cout << "ERROR: ExtToIntParamNumberMap already contains paramNumber=" << paramNumber << std::endl;
                    }

                    // set IntToExt map
                    if(IntToExtParamNumberMap.count(minuit_param_number_counter) == 0)
                    {
                        IntToExtParamNumberMap[minuit_param_number_counter] = paramNumber;
                    }
                    else
                    {
                        std::cout << "ERROR: IntToExtParamNumberMap already contains minuit_param_number_counter=" << minuit_param_number_counter << std::endl;
                    }
                    
                    // set MCNameToExtParamNumberMap
                    // this is done in read param list function
                    // NOTE: a bit weird that some of the setting is done in this class
                    // and some done by an external function
                    // move that external function to be an internal function of this
                    // class

                    // TODO: introduce a paramName variable to the parameter_names.lst file
                    // at the moment this is set automatically using a concatinated list of
                    // the MC values, which does not make sense in the case of xi_31
                    // parameter
                    // this is a hack: using the fact that xi_31 has no associated MC
                    // to set the paramName to "NOMC"
                    // which is a bit dodgy
                    if(it->second.paramName == "xi_31")
                    {
                        if(xi_31_int_ext_param_number_set == false)
                        {
                            xi_31_int_param_number = minuit_param_number_counter;
                            xi_31_ext_param_number = paramNumber;
                            xi_31_int_ext_param_number_set = true;
                        }
                        else
                        {
                            std::cout << "ERROR: xi_31_int_ext_param_number is already set" << std::endl;
                        }
                    }

                    ++ minuit_param_number_counter;
                }
                else
                {
                    std::cout << "param: " << it->first << " is enabled but not enabled for P1 or P2; ignoring" << std::endl;
                }
            }
            else
            {
                std::cout << "param: " << it->first << " is not enabled; ignore" << std::endl;
            }
        }
        std::cout << "done initializing maps" << std::endl;

        if(debug_level >= 2)
        {
            print_map(ExtToIntParamNumberMap, "ExtToIntParamNumberMap");
            print_map(IntToExtParamNumberMap, "IntToExtParamNumberMap");
            print_map(MCNameToExtParamNumberMap, "MCNameToExtParamNumberMap");
            print_map(ParamNameToExtParamNumberMap, "ParamNameToExtParamNumberMap");
        }
    }


    // mc_name is not unique!
    // throw error in case of duplicate mc_names for cases where both params
    // are enabled (P1 and P2)
    // otherwise return P1 and P2 param numbers independently
    bool convert_MC_name_to_param_number(
        const std::string& mc_name,
        int &param_number_P1,
        int &param_number_P2)
    {

        int debuglevel = 1;

        //std::cout << __func__ << " mc_name=" << mc_name << std::endl;
        //bool success = false;
        int count_P1 = 0;
        int count_P2 = 0;
        int param_number_P1_ret = -1;
        int param_number_P2_ret = -1;
        std::map<int, file_parameter>::iterator it{file_params.begin()};
        //int param_index = 0;
        //for(; (success == false) && (it != file_params.end()); ++ it)
        for(; it != file_params.end(); ++ it)
        {
            //for(std::size_t j{0}; (success == false) && (j < it->second.MCNameList.size()); ++ j)
            for(std::size_t j{0}; j < it->second.MCNameList.size(); ++ j)
            {
                //std::cout << "param_index=" << param_index << " j=" << j << " mc_name: " << it->second.MCNameList.at(j) << std::endl;

                if(it->second.MCNameList.at(j) == mc_name)
                {
                    //success = true;
                    //param_number = param_index;
                    if(it->second.paramEnabledP1)
                    {
                        if(count_P1 == 0)
                        {
                            param_number_P1_ret = it->second.paramNumber;
                        }
                        else
                        {
                            param_number_P1_ret = -1;
                        }
                        ++ count_P1;
                    }
                    if(it->second.paramEnabledP2)
                    {
                        if(count_P2 == 0)
                        {   
                            param_number_P2_ret = it->second.paramNumber;
                        }
                        else
                        {
                            param_number_P2_ret = -1;
                        }
                        ++ count_P2;
                    }
                    if(debuglevel >= 2)
                    {
                        //std::cout << mc_name << " -> " << param_index << std::endl;
                        std::cout << mc_name << " -> " << param_number_P1_ret << " " << param_number_P2_ret << std::endl;
                    }
                }
            }

            //++ param_index;
        }

        //if(success == false)
        //{
        //    std::cout << "ERROR in " << __func__ << " could not find object with mc_name=" << mc_name << std::endl;
        //}
    
        if(count_P1 > 1)
        {
            std::cout << "Error: count_P1=" << count_P1 << std::endl;
        }
        if(count_P2 > 1)
        {
            std::cout << "Error: count_P2=" << count_P2 << std::endl;
        }

        if((count_P1 == 1) || (count_P2 == 1))
        {
            if(debuglevel >= 2)
            {
                std::cout << "success: count_P1=" << count_P1 << " count_P2=" << count_P2 << std::endl;
            }

            param_number_P1 = param_number_P1_ret;
            param_number_P2 = param_number_P2_ret;

            return true;
        }

        /*
        if(debuglevel >= 3)
        {
            std::cout << "return: success=" << success << " param_number=" << param_number << std::endl;
        }
        */
        //return success;
        return false;
    }


    void print()
    {
        std::ostream &os{std::cout};

        std::cout << std::endl;
        std::cout << "--------------------------------------------------------------------------------" << std::endl;
        std::cout << __func__ << std::endl;
        std::cout << "numberParams=" << numberParams() << std::endl;
        std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
        //std::cout << "minuitParamNumberCounter=" << minuitParamNumberCounter << std::endl;
        std::cout << "--------------------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
        
        std::map<int, file_parameter>::iterator it{file_params.begin()};
        int param_index = 0;
        for(; it != file_params.end(); ++ it)
        {
            std::cout << "param_index=" << param_index << std::endl;
            int paramNumber = param_index;
            if(it->second.paramEnabled == true)
            {
            int intParamNumber = ExtToIntParamNumberMap.at(paramNumber);
            int extParamNumber = IntToExtParamNumberMap.at(intParamNumber);
            std::cout << "intParamNumber=" << intParamNumber << std::endl;
            std::cout << "extParamNumber=" << extParamNumber << std::endl;
            std::cout << "MCNameList" << std::endl;
            for(int j = 0; j < it->second.MCNameList.size(); ++ j)
            {
                std::string mcname = it->second.MCNameList.at(j);
                int mcnameparamnumber = MCNameToExtParamNumberMap.at(it->second.MCNameList.at(j));
                std::cout << mcname << " -> " << mcnameparamnumber << std::endl;
            }
            std::cout << "paramNumber=" << it->second.paramNumber << std::endl;
            std::cout << "paramEnabled=" << it->second.paramEnabled << std::endl;
            std::cout << "paramEnabledP1=" << it->second.paramEnabledP1 << std::endl;
            std::cout << "paramEnabledP2=" << it->second.paramEnabledP2 << std::endl;
            std::cout << "paramInitValue=" << it->second.paramInitValue << std::endl;
            std::cout << "paramInitError=" << it->second.paramInitError << std::endl;
            std::cout << "paramConstraintValue=" << it->second.paramConstraintValue << std::endl;
            std::cout << "paramConstraintError=" << it->second.paramConstraintError << std::endl;
            std::cout << "paramConstraintMode=" << it->second.paramConstraintMode << std::endl;
            std::cout << "paramName=" << it->second.paramName << std::endl;
            std::cout << "paramHumanReadableName=" << it->second.paramHumanReadableName << std::endl;
            std::cout << "*****" << std::endl;
            }
            else
            {
            std::cout << "MCNameList" << std::endl;
            for(int j = 0; j < it->second.MCNameList.size(); ++ j)
            {
                std::string mcname = it->second.MCNameList.at(j);
                int mcnameparamnumber = MCNameToExtParamNumberMap.at(it->second.MCNameList.at(j));
                std::cout << mcname << " -> " << mcnameparamnumber << std::endl;
            }
            std::cout << "paramNumber=" << it->second.paramNumber << std::endl;
            std::cout << "paramEnabled=" << it->second.paramEnabled << std::endl;
            std::cout << "paramEnabledP1=" << it->second.paramEnabledP1 << std::endl;
            std::cout << "paramEnabledP2=" << it->second.paramEnabledP2 << std::endl;
            std::cout << "paramInitValue=" << it->second.paramInitValue << std::endl;
            std::cout << "paramInitError=" << it->second.paramInitError << std::endl;
            std::cout << "paramConstraintValue=" << it->second.paramConstraintValue << std::endl;
            std::cout << "paramConstraintError=" << it->second.paramConstraintError << std::endl;
            std::cout << "paramConstraintMode=" << it->second.paramConstraintMode << std::endl;
            std::cout << "paramName=" << it->second.paramName << std::endl;
            std::cout << "paramHumanReadableName=" << it->second.paramHumanReadableName << std::endl;
            std::cout << "*****" << std::endl;
            }
            ++ param_index;
        }
        std::cout << std::endl;
        std::cout << "Is the xi_31 parameter number set?" << std::endl;
        if(xi_31_int_ext_param_number_set == true)
        {
            std::cout << "yes" << std::endl;
            std::cout << "xi_31_int_param_number=" << xi_31_int_param_number << std::endl;
            std::cout << "xi_31_ext_param_number=" << xi_31_ext_param_number << std::endl;
        }
        else
        {
            std::cout << "no" << std::endl;
        }
        std::cout << "--------------------------------------------------------------------------------" << std::endl;
    }


    // map the parameters by param number (external)
    std::map<int, file_parameter> file_params;

    // convert MC name to external parameter number
    // there can be more than 1 MC name for a given external parameter
    // number
    std::map<std::string, int> MCNameToExtParamNumberMap; // TODO: probably no longer need this map

    // convert param name to external parameter number
    // there can only be 1 parmeter name for a given external parameter
    // number
    std::map<std::string, int> ParamNameToExtParamNumberMap; // TODO: I don't know if I need this map

    // to be able to convert minuit parameter numbers to external parameter
    // numbers, this is required when getting the number of MC events for
    // example
    std::map<int, int> ExtToIntParamNumberMap;

    // to be able to access parameter propeties such as whether a parameter
    // is enabled in P1 or P2, initial value, constraint, whether param is
    // fixed of free (constraint mode)
    std::map<int, int> IntToExtParamNumberMap;

    int numberParams() const
    {
        return file_params.size();
    }

    //int numberParams; // this is the same as file_params.size();
    int numberEnabledParams; // is this the same as file_params.size(). no it is not
    //int minuitParamNumberCounter;

    int xi_31_int_param_number;
    int xi_31_ext_param_number;
    bool xi_31_int_ext_param_number_set;
};

#endif // PARAMETERGROUP_H
