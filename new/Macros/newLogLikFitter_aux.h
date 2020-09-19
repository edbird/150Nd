#ifndef NEWLOGLIKFITTER_AUX_H
#define NEWLOGLIKFITTER_AUX_H


void cin_wait()
{
    std::cout << "wait for user input" << std::endl;
    std::cin.get();
}







///////////////////////////////////////////////////////////////////////////////
//
// Chi-Square Metric Calculations for Drawing Code
//
// These are used in drawchannel code. Different metrics. Each calculates
// the value of a metric for a data/MC pair. (Single Channel)
///////////////////////////////////////////////////////////////////////////////


// TODO: should this include the parameters?
double calc_chi2_draw_metric_LOGPOISSON(TH1D *data1D, TH1D *hAllMC1D)
{
    double sum = 0.0;
    for(Int_t i = 1; i < hAllMC1D->GetNbinsX(); ++ i)
    {

        double c_data = data1D->GetBinContent(i);
        double c_MC = hAllMC1D->GetBinContent(i);
        double lp = 0.0;

        if(c_MC >= 0.0)
        {
            if(c_MC == 0.0)
            {
                continue;
            }
            lp = logpoisson(c_data, c_MC);
        }

        /*
        double e_data = data1D->GetBinError(i);
        double e_MC = hAllMC1D->GetBinError(i);

        if(c_data == 0.0 && e_data == 0.0) continue;
        if(c_MC == 0.0 && e_MC == 0.0) continue;

        //e_MC = 0.0;

        double num = (c_data - c_MC);
        num = num * num;
        double den = e_MC * e_MC + e_data * e_data;
        sum += num / den;
        */

        //sum += logpoisson(c_data, c_MC);
        sum += lp;
    }
    //return sum;
    return -2.0 * sum;
}


double calc_chi2_draw_metric_CHISQ_ERR_SQRT_MC(TH1D *data1D, TH1D *hAllMC1D)
{
    double sum = 0.0;
    for(Int_t i = 1; i < hAllMC1D->GetNbinsX(); ++ i)
    {

        double c_data = data1D->GetBinContent(i);
        double c_MC = hAllMC1D->GetBinContent(i);
        //double e_MC = std::sqrt(c_MC);
        double lp = 0.0;

        if(c_MC >= 0.0)
        {
            if(c_MC == 0.0)
            {
                continue;
            }
            lp = std::pow(c_data - c_MC, 2.0) / std::abs(c_MC);
        }

        sum += lp;
    }
    return sum;
}



double calc_chi2_draw_metric_CHISQ_ERR_SQRT_DATA(TH1D *data1D, TH1D *hAllMC1D)
{
    double sum = 0.0;
    for(Int_t i = 1; i < hAllMC1D->GetNbinsX(); ++ i)
    {

        double c_data = data1D->GetBinContent(i);
        double c_MC = hAllMC1D->GetBinContent(i);
        //double e_MC = std::sqrt(c_MC);
        double lp = 0.0;

        if(c_MC >= 0.0)
        {
            if(c_MC == 0.0)
            {
                continue;
            }
            lp = std::pow(c_data - c_MC, 2.0) / std::abs(c_data);
        }

        sum += lp;
    }
    return sum;
}

///////////////////////////////////////////////////////////////////////////////
// calculate number of degrees of freedom for chisquare calculation
///////////////////////////////////////////////////////////////////////////////

// NOTE: this is not really number of degrees of freedom but number of channels
// (channels = bins)
int get_ndf_1D(const TH1D *const hist_MC, const TH1D *const hist_data)
{
    if(hist_MC->GetNbinsX() != hist_data->GetNbinsX())
    {
        std::cout << __func__ << " Error, incompatiable histograms" << std::endl;
        return -1;
    }

    //std::cout << __func__ << " " << hist_MC->GetName() << std::endl;

    int ndf = 0;
    for(Int_t bin_ix = 0; bin_ix <= hist_MC->GetNbinsX(); ++ bin_ix)
    {
        if(hist_MC->GetBinContent(bin_ix) > 0.0)
        {
            //if(hist_data->GetBinContent(bin_ix) > 0.0)
            //{
                //std::cout << "bin_ix: " << hist_MC->GetBinContent(bin_ix) << " " << hist_data->GetBinContent(bin_ix) << std::endl;
                ++ ndf;
            //}
        }
    }
    return ndf;
}


// calculate number of channels (1D)
// channels = bins
int get_nch_1D(const TH1D *const hist_MC, const TH1D *const hist_data)
{
    if(hist_MC->GetNbinsX() != hist_data->GetNbinsX())
    {
        std::cout << __func__ << " Error, incompatiable histograms" << std::endl;
        return -1;
    }

    //std::cout << __func__ << " " << hist_MC->GetName() << std::endl;

    int nch = 0;
    for(Int_t bin_ix = 0; bin_ix <= hist_MC->GetNbinsX(); ++ bin_ix)
    {
        if(hist_MC->GetBinContent(bin_ix) > 0.0)
        {
            //if(hist_data->GetBinContent(bin_ix) > 0.0)
            //{
                //std::cout << "bin_ix: " << hist_MC->GetBinContent(bin_ix) << " " << hist_data->GetBinContent(bin_ix) << std::endl;
            ++ nch;
            //}
        }
    }
    return nch;
}



///////////////////////////////////////////////////////////////////////////////
// internal / external index manipulation
///////////////////////////////////////////////////////////////////////////////

/*
int get_axial_vector_parameter_index()
{

    std::string mc_name = "axial_vector_parameter_0";
    std::string search_object = MCNameToParamNameMap.at(mc_name);
    int axial_vector_parameter_0_param_number = -1;
    if(paramNameToNumberMap.count(search_object) > 0)
    {
        int param_number = paramNameToNumberMap.at(search_object);
        axial_vector_parameter_0_param_number = param_number;
   
        if(param_number != 1)
        {
            throw "param_number != 1";
        }

    }
    else
    {
        throw "mc_name not found in paramNameToNumberMap";
    }
    
    return axial_vector_parameter_0_param_number;
}
*/



///////////////////////////////////////////////////////////////////////////////
// Get initial value and error of parameter depending on phase
///////////////////////////////////////////////////////////////////////////////

/*
void get_paramInitValueError(const Int_t thePhase, const int param_number, double &param_init_value, double& param_init_error)
{
    const int j = param_number;

    if(thePhase == 0)
    {
        param_init_value = paramInitValueP1Map[j];
        param_init_error = paramInitErrorP1Map[j];
    }
    else if(thePhase == 1)
    {
        param_init_value = paramInitValueP2Map[j];
        param_init_error = paramInitErrorP2Map[j];
    }
    else
    {
        std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << " in function " << __func__ << std::endl;
    }
}

void get_paramConstraintValueError(const Int_t thePhase, const int param_number, double &param_constraint_value, double &param_constraint_error)
{
    const int i = param_number;

    if(thePhase == 0)
    {
        param_constraint_value = paramConstraintValueP1Map[i];
        param_constraint_error = paramConstraintErrorP1Map[i];
    }
    else if(thePhase == 1)
    {
        param_constraint_value = paramConstraintValueP2Map[i];
        param_constraint_error = paramConstraintErrorP2Map[i];
    }
    else
    {
        std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << " in function " << __func__ << std::endl;
    }
}
*/



///////////////////////////////////////////////////////////////////////////////
// convert name of histogram to minuit (internal) parameter number
///////////////////////////////////////////////////////////////////////////////

#if 0
bool fit_histogram_name_to_param_number(const TString &name, int &which_param)
{

    // search through parameters to find right one
    // the histogram names are formatted like:
    // hTotalE_bi214_mylar_fit
    // histogram_name + "_" + mc_sample_name + "_fit"


    // used later
    //double activity_scale_branching_ratio = 1.0;

    which_param = -1;
    bool found_param = false;

    // get index for parameter

    std::string tmp_hist_name(name);
    auto i_start = tmp_hist_name.find('_') + 1;
    auto i_end = tmp_hist_name.rfind('_');
    if(i_end - i_start > 0)
    {
        std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);

        // set branching ratio fraction
        /*
        if(tmp_sample_name == std::string("tl208_int_rot") ||
           tmp_sample_name == std::string("tl208_feShield") ||
           tmp_sample_name == std::string("tl208_pmt"))
        {
            activity_scale_branching_ratio = 0.36;
        }
        */

        if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
        {
            int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
            //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> name=" << name << std::endl;
            which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
            found_param = true;
        }
        else
        {
           std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
        }
    }

    return found_param;

}
#endif


///////////////////////////////////////////////////////////////////////////////
// convert MC name (background file name) to scale factor
// helper function for book1DHistograms_helper and book2DHistograms_helper
///////////////////////////////////////////////////////////////////////////////

#if 0
bool convert_MC_name_to_scale_factor(
    const std::string &mc_name,
    int &param_number,
    Double_t &scale_factor)
{

    bool ret = false;

    // example: "bi214_int_rot" -> "bi214_int_rot,pb214_int_rot"
    std::string search_object = MCNameToParamNameMap.at(mc_name);

    // example: "bi214_int_rot,pb214_int_rot" -> 4
    // check if parameter number exists
    // (was defined by parameter_names.lst)
    if(paramNameToNumberMap.count(search_object) > 0)
    {
        // convert from mc sample name to param number
        /*int*/ param_number = paramNameToNumberMap.at(search_object);

        // check if this parameter number is enabled
        if(std::find(enabled_params.begin(), enabled_params.end(), param_number) != enabled_params.end())
        {

            // scale by activity

            // convert parameter number to minuit parameter number
            //minuit_param_number = paramNumberToMinuitParamNumberMap.at(param_number);

            // TODO: change such that samples are pre-scaled by activity input value
            // get initial parameter values and error
            Double_t param_init_value = 0.;
            Double_t param_init_error = 0.;
            get_paramInitValueError(thePhase, param_number, param_init_value, param_init_error);
            /*Double_t*/ scale_factor = param_init_value;

            // account for 208 Tl branching ratio of 36 %
            // TODO: should I move this into fit_2e code
            // and apply using ->Fill() function call with
            // weight = 0.36
            if(mc_name == std::string("tl208_int_rot") ||
               mc_name == std::string("tl208_feShield") ||
               mc_name == std::string("tl208_pmt"))
               // TODO: do not apply to tl208_air ?
            {
                //std::cout << "mc_name=" << mc_name << " applying additional scaling factor of 0.36" << std::endl;
                //std::cin.get();
                scale_factor *= 0.36;
            }

            ret = true;
        }
        else
        {
            // paramter not enabled, do not load histogram/sample
            std::cout << "parameter number " << param_number << " is not enabled (not found in vector)" << std::endl;

            ret = false;
        }
    }
    else
    {
        std::cout << "!!!!! ERROR: search_object=" << search_object << " not found in paramNameToNumberMap" << std::endl;
        std::cout << "mc_name=" << mc_name << std::endl;

        print_map(paramNameToNumberMap, "paramNameToNumberMap");
        print_map(MCNameToParamNameMap, "MCNameToParamNameMap");

        ret = false;
    }

    return ret;
}
#endif


///////////////////////////////////////////////////////////////////////////////
// convert MC name (background file name) to (external) parameter number
///////////////////////////////////////////////////////////////////////////////

//bool MC_name_to_param_number(const std::string& mc_name, int





// split filename into base name and extension
// split string based on rfind '.' char
// return 0 if failed to split
// return 1 if success
// if '.' is found, it is contained and returned in output_extension
int filename_split_extension(const std::string &input, std::string &output_base, std::string &output_extension)
{
    std::size_t pos = input.rfind('.');
    if(pos != std::string::npos)
    {
        output_base = input.substr(0, pos);
        output_extension = input.substr(pos);
        return 1;
    }
    else
    {
        output_base = input;
        output_extension = "";
        return 0;
    }
}


void canvas_saveas_helper(const std::string &directory, const std::string &saveas_filename, TCanvas *canvas)
{
    if(saveas_filename.size() > 0)
    {
        std::size_t length = saveas_filename.size();
        if(
            (length >= 2) &&
            ((saveas_filename[length - 1] == '*') && (saveas_filename[length - 2] == '.'))
          )
        {

            // TODO: get the name based on the histogram type
            // from data accessable via loop index i
            std::string base_name = saveas_filename.substr(0, length - 2);
            std::vector<std::string> extensions;
            extensions.push_back("C");
            extensions.push_back("png");
            extensions.push_back("eps");
            extensions.push_back("pdf");
            for(auto it{extensions.begin()}; it != extensions.end(); ++ it)
            {
                //std::string name = base_name + "_" + std::string(i_str) + "." + *it;
                std::string name = base_name + "." + *it;
                std::string fullname;
                if(directory.size() > 0)
                {
                    fullname = directory + "/";
                }
                fullname += name;
                std::cout << "extension=\".*\": multisave: saving as " << fullname << std::endl;
                canvas->SaveAs(fullname.c_str());
            }
        }
        else
        {
            //std::cout << "saveas_filename=" << saveas_filename << std::endl;

            
            std::string base_name;
            std::string extension;
            if(filename_split_extension(saveas_filename, base_name, extension) == 1)
            {
                // contains an extension
            }
            else
            {
                // does not contain an extension
                //extension = ".png";
            }
            /*
            std::string base_name;
            std::string extension;
            std::size_t pos = saveas_filename.rfind('.');
            if(pos != std::string::npos)
            {
                base_name = saveas_filename.substr(0, pos);
                extension = saveas_filename.substr(pos);

                //std::cout << "base_name=" << base_name << std::endl;
                //std::cout << "extension=" << extension << std::endl;
            }
            else
            {
                base_name = saveas_filename;
            }
            */

            //std::string name = base_name + "_c" + "_" + std::string(i_str);
            std::string name = base_name;
            if(extension.size() > 0)
            {
                name += extension;
            }
            std::string fullname;
            if(directory.size() > 0)
            {
                fullname = directory + "/";
            }
            fullname += name;
            std::cout << "saving as " << fullname << std::endl;
            canvas->SaveAs(fullname.c_str());
        }
    }
    else
    {
        std::cout << "Error: saveas_filename.size() == 0" << std::endl;
    }
}

#endif // NEWLOGLIKFITTER_AUX_H
