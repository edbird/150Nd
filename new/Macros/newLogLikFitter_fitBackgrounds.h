#ifndef NEWLOGLIKFITTER_FITBACKGROUNDS_H
#define NEWLOGLIKFITTER_FITBACKGROUNDS_H


///////////////////////////////////////////////////////////////////////////////
// fitBackgrounds
///////////////////////////////////////////////////////////////////////////////



void fitBackgrounds_init(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
//    double *AdjustActs,
//    double *AdjustActs_Err
    const double xi_31_value, const double xi_31_error
    )
{
//    std::cout << " 1 " << AdjustActs[1] << std::endl;
//    std::cin.get();

    int debuglevel = 1;
    if(debuglevel >= 3)
    {
        std::cout << ">>>>> fitBackgrounds_init()" << std::endl;
    }

    //TVirtualFitter::SetDefaultFitter("Minuit2");
    
    // TODO:
    // want to change Minuit so that it does not store disabled parameters
    // do this using a map to translate from parameter number to minuit
    // parameter numbers
    // and change numberParams to numberEnabled params
    //TMinuit *minuit = new TMinuit(numberParams);
    //std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
    //TMinuit *minuit = new TMinuit(numberEnabledParams);
    // TODO working here need to check all instances of numberParams


    //std::cout << "Fit created" << std::endl;

    //std::cout << "minuit tests" << std::endl;
    //minuit->DefineParameter(1, "test1", 100.0, 10.0, -500.0, 500.0);
    //std::cin.get();

    // set debug level
    //minuit->SetPrintLevel(3);

    // moved reading of parameter list file to before book histogram function
    // calls
  
    // TODO: printing sample_names array, what does this do, is it still used?
    if(debuglevel >= 4)
    {
        //print_paramNameMap();
        //g_pg.print();
    }

    if(debuglevel >= 5)
    {
        std::cout << "set errors" << std::endl;
    }

    // TODO: de-interlace
    // NOTE: done

    ///////////////////////////////////////////////////////////////////////////
    // Phase 1 and 2: MINUIT PARAMETERS
    ///////////////////////////////////////////////////////////////////////////

    // loop through all parameters in g_pg.fileparams and initialize them
    //for(int i = 0; i < numberParams; i++)
    //int i_minuit = 0; // internal minuit counter
    //for(int i = 0; i < g_pg.numberParams(); ++ i)
    std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
    for(; it != g_pg.file_params.end(); ++ it)
    {
        int minuit_param_number = -1;
        
        int paramNumber = it->second.paramNumber;
        bool paramEnabled = it->second.paramEnabled;
        bool paramEnabledP1 = it->second.paramEnabledP1;
        bool paramEnabledP2 = it->second.paramEnabledP2;
        double paramInitValue = it->second.paramInitValue;
        double paramInitError = it->second.paramInitError;
        int paramConstraintMode = it->second.paramConstraintMode;

        if(debuglevel >= 5)
        {
            std::cout << "paramNumber=" << paramNumber << std::endl;
        }

        // decide what to do depending on whether parameter is enabled
        // and for which phases

        bool ok = false;
        if(paramEnabled == true)
        {
            if(gEnablePhase1 == true)
            {
                if(paramEnabledP1 == true)
                {
                    ok = true;
                }
            }

            if(gEnablePhase2 == true)
            {
                if(paramEnabledP2 == true)
                {
                    ok = true;
                }
            }
        }
        else
        {
            ok = true;
            continue;
        }
        if(ok == false)
        {
            std::cout << __func__ << " ok == false" << std::endl;
            std::cin.get();
            continue;
        }
        // TODO: this correctly ignores any parameter which is disabled such that
        // paramEnabled == false
        // however, it also ignores parameters which are set as disabled for P1
        // and P2, and when these are irrelevent due to the value of gEnablePhaseX
        // so... the internal and external index will not match
        // need to add some code to fix this when the parameters are read from
        // file, (probably)
        // unless I just ignore that here... perhaps paramEnabled dictates
        // whether parameter is drawn and the phase1/phase2 enable flag
        // is to decide whether minuit does the fit or not (in which case
        // the param may still contribute to chisquare but may not be minimized
        // by minuit)


        TString paramNumber_str;
        paramNumber_str.Form("%i", paramNumber);
        minuit_param_number = g_pg.ExtToIntParamNumberMap.at(paramNumber);
        TString minuit_param_number_str;
        minuit_param_number_str.Form("%i", minuit_param_number);

        if(paramConstraintMode == MODE_CONSTRAINT_HARD)
        {
            // this is a fixed parameter
            // define parameter using constrained value if hard constrained
            if(debuglevel >= 6)
            {
                std::cout << "HARD" << std::endl;
            }

            TString minuit_param_name = "_" + paramNumber_str + "_" + minuit_param_number_str + "_FIXED";
            
            if(paramNumber == g_pg.get_xi_31_ext_param_number())
            {
                // xi_31 parameter
                theParameterState.Add(std::string(minuit_param_name), xi_31_value, xi_31_error); // instead of _Err was 0.5
                //theParameterState.SetLowerLimit(i, 0.0); // no limit for xi
            }
            else
            {
                // MC sample amplitude parameter
                theParameterState.Add(std::string(minuit_param_name), 1.0, 0.5);
                //theParameterState.SetLowerLimit(minuit_param_number, 0.0);
            }
            theParameterState.Fix(std::string(minuit_param_name));
        }
        else
        {
            // param is either free or soft constrained
            if(debuglevel >= 6)
            {
                std::cout << "SOFT / FREE" << std::endl;
            }

            TString minuit_param_name = "_" + paramNumber_str + "_" + minuit_param_number_str + "_";

            if(g_pg.get_xi_31_int_param_number() != -1)
            {
                // TODO: this does not work if xi_31 param is disabled
                // 2020-09-23: FIXED
                if(paramNumber == g_pg.get_xi_31_ext_param_number())
                {
                    // xi_31 parameter
                    theParameterState.Add(std::string(minuit_param_name), xi_31_value, xi_31_error); // instead of _Err was 0.5
                    theParameterState.SetLowerLimit(minuit_param_number, -1.0); // was -0.4
                }
                else
                {
                    // MC sample amplitude parameter
                    theParameterState.Add(std::string(minuit_param_name), 1.0, 0.5);
                    theParameterState.SetLowerLimit(minuit_param_number, 0.0);
                }
            }
            else
            {
                    // MC sample amplitude parameter
                    theParameterState.Add(std::string(minuit_param_name), 1.0, 0.5);
                    theParameterState.SetLowerLimit(minuit_param_number, 0.0);
            }

            // TODO: set initial error using initError/initValue


            if(FORCE_BKG_HARD == true)
            {
                if(minuit_param_number > 1)
                {
                    std::cout << "Fixing parameter number: " << minuit_param_number << std::endl;
                    theParameterState.Fix(std::string(minuit_param_name));
                }
            }
        }


        // note, only enabled parameters
        // get through to this execution block
        // if parameter is enabled, then it is either fixed or free
        // it has been copied and edited above

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


    if(debuglevel >= 5)
    {
        std::cout << "all parameters fixed" << std::endl;
    }

    // TODO: check what parameters Summer had fixed and check the MC cuts
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




    /*
    if(debuglevel >= 5)
    {
        std::cout << "Ready to exec fix" << std::endl;
    }
    */


    //minuit->SetErrorDef(0.5);
    // TODO ? required ?
    // 1.0 = chisquare
    // 0.5 = negative log likelihood
    

    // mnsimp()?

    
    //minuit->SetMaxIterations(50000);
    //minuit->SetMaxIterations(1000); TODO
    //minuit->mnexcm("SET EPS", 0.01);
    //minuit->SetEPS(1.0e-3); // TODO
    //give it the function to use
    //minuit->SetFCN(logLikelihood); // done elsewhere
    //std::cout << "calling: minuit->mnsimp()" << std::endl;
    //minuit->mnsimp();
    //std::cout << "calling: minuit->Migrad()" << std::endl;


// draw 1D histograms, write chisquare values to file for range of
// different gA values
//    newloglikfitter_gA_chisquaretest(minuit, AdjustActs, AdjustActs_Err);

    /*
    Int_t npar = 0;
    Double_t fval;
    logLikelihood(npar, nullptr, fval, AdjustActs, 0);
    std::cout << "fval=" << fval << std::endl;
    draw(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);
    std::cin.get();

    AdjustActs[1] = 0.296 -0.1;
    logLikelihood(npar, nullptr, fval, AdjustActs, 0);
    std::cout << "fval=" << fval << std::endl;
    draw(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);
    std::cin.get();

    AdjustActs[1] = 0.296 + 0.1;
    logLikelihood(npar, nullptr, fval, AdjustActs, 0);
    std::cout << "fval=" << fval << std::endl;
    draw(AdjustActs, AdjustActs_Err, CovMatrix, number_free_params);
    std::cin.get();
    */


    // TODO: something should go here?
    // to set some values in the header file?
    // they might be set elsewhere, if so remove them
    #if 0
    {
        /*
        int num_params = minuit->GetNumFreePars();
        minuitParamCurrent = new double[num_params];
        minuitParamInit = new double[num_params];
        minuitParamLast = new double[num_params];
        for(int i = 0; i < num_params; ++ i)
        {
            int free_param_index_ext = free_params.at(i);
            double value, error;
            minuit->GetParameter(free_param_index_ext, value, error);
            int free_param_index_int = paramNumberToMinuitParamNumberMap.at(i);
            minuitParamCurrent[free_param_index_int] = value;
            minuitParamInit[free_param_index_int] = value;
            minuitParamLast[free_param_index_int] = value;
        }
        */

        minuitParamCurrent = new double[numberEnabledParams];
        minuitParamInit = new double[numberEnabledParams];
        minuitParamLast = new double[numberEnabledParams];
        for(int i = 0; i < numberEnabledParams; ++ i)
        {
            double value, error;
            //minuit->GetParameter(i, value, error);
            value = theParameterState.Value(i);
            error = theParameterState.Error(i);
            //int j = minuitParamNumberTo
            minuitParamCurrent[i] = value;
            minuitParamInit[i] = value;
            minuitParamLast[i] = value;

        }
    }
    #endif


    //return minuit;
}




void fitBackgrounds_phasespace_init(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
    //double *AdjustActs,
    //double *AdjustActs_Err,
    const double Nd150_A_value,
    const double Nd150_A_error,
    const double xi_31_value,
    const double xi_31_error
    )
{
    //std::cout << ">>>>> fitBackgrounds_phasespace_init()" << std::endl;
    //std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;

    const int debuglevel = 1;


    ///////////////////////////////////////////////////////////////////////////
    // Phase 1: MINUIT PARAMETERS
    ///////////////////////////////////////////////////////////////////////////

    
    std::map<int, file_parameter>::iterator it{g_pg.file_params.begin()};
    for(; it != g_pg.file_params.end(); ++ it)
    {
        int minuit_param_number = -1;
        
        int paramNumber = it->second.paramNumber;
        bool paramEnabled = it->second.paramEnabled;
        bool paramEnabledP1 = it->second.paramEnabledP1;
        bool paramEnabledP2 = it->second.paramEnabledP2;
        double paramInitValue = it->second.paramInitValue;
        double paramInitError = it->second.paramInitError;
        int paramConstraintMode = it->second.paramConstraintMode;

        if(debuglevel >= 3)
        {
            std::cout << "paramNumber=" << paramNumber << std::endl;
        }

        // decide what to do depending on whether parameter is enabled
        // and for which phases

        bool ok = false;
        if(paramEnabled == true)
        {
            if(gEnablePhase1 == true)
            {
                if(paramEnabledP1 == true)
                {
                    ok = true;
                }
            }

            if(gEnablePhase2 == true)
            {
                if(paramEnabledP2 == true)
                {
                    ok = true;
                }
            }
        }
        else
        {
            ok = true;
            continue;
        }
        if(ok == false)
        {
            std::cout << __func__ << " ok == false" << std::endl;
            std::cin.get();
            continue;
        }
        // TODO: this correctly ignores any parameter which is disabled such that
        // paramEnabled == false
        // however, it also ignores parameters which are set as disabled for P1
        // and P2, and when these are irrelevent due to the value of gEnablePhaseX
        // so... the internal and external index will not match
        // need to add some code to fix this when the parameters are read from
        // file, (probably)
        // unless I just ignore that here... perhaps paramEnabled dictates
        // whether parameter is drawn and the phase1/phase2 enable flag
        // is to decide whether minuit does the fit or not (in which case
        // the param may still contribute to chisquare but may not be minimized
        // by minuit)


        TString paramNumber_str;
        paramNumber_str.Form("%i", paramNumber);
        minuit_param_number = g_pg.ExtToIntParamNumberMap.at(paramNumber);
        TString minuit_param_number_str;
        minuit_param_number_str.Form("%i", minuit_param_number);

        if(paramConstraintMode == MODE_CONSTRAINT_HARD)
        {
            // this is a fixed parameter
            // define parameter using constrained value if hard constrained
            if(debuglevel >= 3)
            {
                std::cout << "HARD" << std::endl;
            }

            TString minuit_param_name = "_" + paramNumber_str + "_" + minuit_param_number_str + "_FIXED";
            
            if(paramNumber == 0)
            {
                theParameterState.Add(std::string(minuit_param_name), Nd150_A_value, Nd150_A_error);
                //theParameterState.Fix(std::string(minuit_param_name));
            }
            //else if(paramNumber == g_pg.get_xi_31_ext_param_number())
            else if(paramNumber == 1)
            {
                // xi_31 parameter
                //theParameterState.Add(std::string(minuit_param_name), xi_31_value, xi_31_error); // instead of _Err was 0.5
                //theParameterState.SetLowerLimit(i, 0.0); // no limit for xi
                theParameterState.Add(std::string(minuit_param_name), xi_31_value, xi_31_error);
                //theParameterState.Fix(std::string(minuit_param_name));
            }
            else
            {
                // MC sample amplitude parameter
                theParameterState.Add(std::string(minuit_param_name), 1.0, 0.1);
                //theParameterState.SetLowerLimit(minuit_param_number, 0.0);
            }
            theParameterState.Fix(std::string(minuit_param_name));
        }
        else
        {
            // param is either free or soft constrained
            if(debuglevel >= 3)
            {
                std::cout << "SOFT / FREE" << std::endl;
            }

            TString minuit_param_name = "_" + paramNumber_str + "_" + minuit_param_number_str + "_";

            if(paramNumber == 0)
            {
                theParameterState.Add(std::string(minuit_param_name), Nd150_A_value, Nd150_A_error);
                theParameterState.Fix(std::string(minuit_param_name));
            }
            //if(paramNumber == g_pg.get_xi_31_ext_param_number())
            else if(paramNumber == 1)
            {
                // xi_31 parameter
                theParameterState.Add(std::string(minuit_param_name), xi_31_value, xi_31_error); // instead of _Err was 0.5
                theParameterState.Fix(std::string(minuit_param_name));
                //theParameterState.SetLowerLimit(minuit_param_number, -1.0); // was -0.4
            }
            else
            {
                // MC sample amplitude parameter
                theParameterState.Add(std::string(minuit_param_name), 1.0, 0.1);
                theParameterState.SetLowerLimit(minuit_param_number, 0.0);
            }

           // TODO: set initial error using initError/initValue


            if(FORCE_BKG_HARD == true)
            {
                if(minuit_param_number > 1)
                {
                    std::cout << "Fixing parameter number: " << minuit_param_number << std::endl;
                    theParameterState.Fix(std::string(minuit_param_name));
                }
            }
        }

    }

}


void fitBackgrounds_setparams(TMinuit *minuit, double *AdjustActs, double *AdjustActs_Err)
{

}


ROOT::Minuit2::FunctionMinimum fitBackgrounds_exec(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
    MinimizeFCNAxialVector &theFCN
    )
{

    ll_walk.clear();

    ROOT::Minuit2::MnStrategy theStrategy(1);
    ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, theParameterState, theStrategy);
    ll_walk_save = ll_walk;
    std::cout << "walk length: " << ll_walk_save.size() << std::endl;
    ll_walk.clear();

    #if 0
    // Then get results
    //for(int i = 0; i < numberParams; i++)
    for(int i = 0; i < numberEnabledParams; i++)
    {
        //minuit->GetParameter(i, AdjustActs[i], AdjustActs_Err[i]);
        double value = theParameterState.Value(i);
        double error = theParameterState.Error(i);
        AdjustActs[i] = value;
        AdjustActs_Err[i] = error;
    }
    #endif

    return FCN_min;

}


ROOT::Minuit2::FunctionMinimum fitBackgrounds_phasespace_exec(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
    MinimizeFCNAxialVector &theFCN
    )
{
    ROOT::Minuit2::MnStrategy theStrategy(1);
    ROOT::Minuit2::FunctionMinimum FCN_min = theMinimizer.Minimize(theFCN, theParameterState, theStrategy);

    return FCN_min;

}



void fitBackgrounds_postexectest(TMinuit *minuit, double *AdjustActs, double *AdjustActs_Err)
{
#if 0

    // TODO: remove AdjustActs, AdjustActs_Err arguments?
    // Do I want these to be copy of original values

    if(0)
    {
        //newloglikfitter_100Mo_chisquaretest(minuit, AdjustActs, AdjustActs_Err);
    }

    // Some people have seen further by standing on the shoulders of giants.
    // In my case, my vision has been obscured by floundering hopeless idiots
    // doing about as much as it was possible to do to inhibit my ability to
    // see anything, either simply by being incompetent, or by actively
    // doing everything possible to throw as many spanners into the works as
    // time would allow for.

    if(0)
    {
        //newloglikfitter_testmyphasespace(minuit, AdjustActs, AdjustActs_Err);
    }

#endif

}


void fitBackgrounds_getcovmatrix(TMinuit* minuit, double *&CovMatrix, int& number_free_params)
{
#if 0
    // TODO: this no longer works, or does it?
    // needs to take into account the number of ENABLED free params
    // NOTE: 2020-04-16 fixed
    number_free_params = minuit->GetNumFreePars();
    CovMatrix = new Double_t[number_free_params * number_free_params];
    for(int ix{0}; ix < number_free_params * number_free_params; ++ ix)
    {
        CovMatrix[ix] = 0.;
    }
    minuit->mnemat(CovMatrix, number_free_params);  
#endif

}




//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
//TMinuit* fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
#if 0
ROOT::Minuit2::FunctionMinimum fitBackgrounds(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
    MinimizeFCNAxialVector &theFCN,
    //double *AdjustActs,
    //double *AdjustActs_Err,
    //double *&CovMatrix,
    //int &number_free_params
    const double xi_31_value, const double xi_31_error
    )
{
#if 0
    std::cout << ">>>>> fitBackgrounds()" << std::endl;
    
    //TMinuit* minuit = fitBackgrounds_init(AdjustActs, AdjustActs_Err);
    //fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
    fitBackgrounds_init(theParameterState, theMinimizer, xi_31_value, xi_31_error);
    /*ROOT::Minuit2::FunctionMinimum FCN_min*/
    if(1)
    {
        //fitBackgrounds_exec(minuit, AdjustActs, AdjustActs_Err);
        ROOT::Minuit2::FunctionMinimum FCN_min = fitBackgrounds_exec(theParameterState, theMinimizer, theFCN);
        //fitBackgrounds_postexectest(minuit, AdjustActs, AdjustActs_Err);
        //fitBackgrounds_getcovmatrix(minuit, CovMatrix, number_free_params);
        // TODO

    
    return FCN_min;
    }


#endif
}




ROOT::Minuit2::FunctionMinimum fitBackgrounds_phasespace(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
    MinimizeFCNAxialVector &theFCN,
    //double *AdjustActs,
    //double *AdjustActs_Err,
    //double *&CovMatrix,
    const double Nd150_A,
    const double xi_31
    )
{
#if 0
    std::cout << ">>>>> fitBackgrounds_phasespace()" << std::endl;
    
    //TMinuit* minuit = fitBackgrounds_init(AdjustActs, AdjustActs_Err);
    //fitBackgrounds_phasespace_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err, Nd150_A, xi_31);
    fitBackgrounds_phasespace_init(theParameterState, theMinimizer, Nd150_A, xi_31);
    /*ROOT::Minuit2::FunctionMinimum FCN_min*/
    //fitBackgrounds_exec(minuit, AdjustActs, AdjustActs_Err);
    ROOT::Minuit2::FunctionMinimum FCN_min = fitBackgrounds_phasespace_exec(theParameterState, theMinimizer, theFCN); // TODO: same function?
    //fitBackgrounds_postexectest(minuit, AdjustActs, AdjustActs_Err);
    //fitBackgrounds_getcovmatrix(minuit, CovMatrix, number_free_params);
    // TODO

    return FCN_min;
#endif
}
#endif



#endif // NEWLOGLIKFITTER_FITBACKGROUNDS_H
