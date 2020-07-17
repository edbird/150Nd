#ifndef NEWLOGLIKFITTER_AUX_H
#define NEWLOGLIKFITTER_AUX_H


///////////////////////////////////////////////////////////////////////////////
// internal / external index manipulation
///////////////////////////////////////////////////////////////////////////////

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


///////////////////////////////////////////////////////////////////////////////
// Get initial value and error of parameter depending on phase
///////////////////////////////////////////////////////////////////////////////

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
                std::cout << "saving as " << fullname << std::endl;
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
