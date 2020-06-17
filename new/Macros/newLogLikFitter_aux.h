#ifndef NEWLOGLIKFITTER_AUX_H
#define NEWLOGLIKFITTER_AUX_H


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

#endif // NEWLOGLIKFITTER_AUX_H
