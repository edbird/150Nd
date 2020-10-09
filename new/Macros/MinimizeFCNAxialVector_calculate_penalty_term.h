#ifndef MINIMIZEFCNAXIALVECTOR_CALCULATE_PENALTY_TERM_H
#define MINIMIZEFCNAXIALVECTOR_CALCULATE_PENALTY_TERM_H


void
MinimizeFCNAxialVector::calculate_penalty_term(double& penalty_term_ret, const std::vector<double> &param) const
{
    // penalty terms section
    double penalty_term = 0.0;

    // loop over all the parameters
    std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
    for(; it != g_pg.file_params.end(); ++ it)
    {
        int paramNumberInt = -1;

        int paramNumber = it->second.paramNumber;
        bool paramEnabled = it->second.paramEnabled;
        bool paramEnabledP1 = it->second.paramEnabledP1;
        bool paramEnabledP2 = it->second.paramEnabledP2;
        double paramInitValue = it->second.paramInitValue;
        double paramInitError = it->second.paramInitError;
        double paramConstraintValue = it->second.paramConstraintValue;
        double paramConstraintError = it->second.paramConstraintError;
        int paramConstraintMode = it->second.paramConstraintMode;

        // stop it crashing
        // (below)
        if(paramEnabled == false)
        {
            continue;
        }

        // ignore parameter if it is of background type and FORCE_BKG_HARD
        // is set
        if(FORCE_BKG_HARD == true)
        {
            if(paramNumber > 1)
            {
                continue;
            }
        }

        paramNumberInt = g_pg.ExtToIntParamNumberMap.at(paramNumber);

        if(debuglevel >= 5)
        {
            std::cout << "paramNumber=" << paramNumber << " -> " << paramNumberInt << std::endl;
        }

        if(paramEnabled == true)
        {
            if((paramEnabledP1 == true) || (paramEnabledP2 == true))
            {
                // do nothing
            }
            else
            {
                // not enabled
                continue;
            }
        }
        else
        {
            // not enabled
            continue;
        }

        if(paramConstraintMode == MODE_CONSTRAINT_SOFT)
        {
            if(debuglevel >= 6)
            {
                std::cout << "soft" << std::endl;
            }
            // do nothing
        }
        else if(paramConstraintMode == MODE_CONSTRAINT_FREE)
        {
            //-- ndf;
            if(debuglevel >= 6)
            {
                std::cout << "free" << std::endl;
            }
            continue;
        }
        else if(paramConstraintMode == MODE_CONSTRAINT_HARD)
        {
            if(debuglevel >= 6)
            {
                std::cout << "hard" << std::endl;
            }
            continue;
        }
        else
        {
            std::cout << "ERROR: Invalid value for paramNumber=" << paramNumber << ", paramConstraintMode=" << paramConstraintMode << std::endl;
        }

        // this parameter is from minuit internal and thus is in minuit
        // internal units (not Bq)
        // have to convert to Bq units
        double param_value = param.at(paramNumberInt); // TODO this might break if we try to fit P1 seperatly and P2 seperatly
        double activity_value_Bq = paramInitValue;
    
        // convert to Bq
        // multiply by the initial value
        // activity_value_Bq should really be called param_initial_value
        //double activity_value_Bq = 0.0;
        //double tmp_err;
        //get_paramInitValueError(thePhase, i, activity_value_Bq, tmp_err);

        double value = param_value * activity_value_Bq;
        //double penalty = std::pow((param_value - constraint) / error, 2.0);
        double penalty = 0.0;

        // TODO
        /*
        if(EMODE == 1)
        {
            // data
        }
        else if(EMODE == 2)
        {
            // MC
        }   
        else if(EMODE == 3)
        {
            // quadrature
        }
        */
        
        penalty = std::pow((value - paramConstraintValue) / paramConstraintError, 2.0);
        if(debuglevel >= 5)
        {
            //std::cout << "j=" << j << std::endl;
            std::cout << "paramNumber=" << paramNumber
                      << " value=" << value
                      << " param_value=" << param_value
                      << " paramConstraintValue=" << paramConstraintValue
                      << " paramConstraintError=" << paramConstraintError
                      << " penalty=" << penalty
                      << std::endl;
        }


        // TODO: is this the correct error term?
        // error on constraint rather than error on current fit value?

        penalty_term += penalty;
    }

    penalty_term_ret = penalty_term;
}

#endif // MINIMIZEFCNAXIALVECTOR_CALCULATE_PENALTY_TERM_H
