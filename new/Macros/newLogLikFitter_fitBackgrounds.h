#ifndef NEWLOGLIKFITTER_FITBACKGROUNDS_H
#define NEWLOGLIKFITTER_FITBACKGROUNDS_H


///////////////////////////////////////////////////////////////////////////////
// fitBackgrounds
///////////////////////////////////////////////////////////////////////////////



void fitBackgrounds_init(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
    double *AdjustActs,
    double *AdjustActs_Err
    )
{
//    std::cout << " 1 " << AdjustActs[1] << std::endl;
//    std::cin.get();

    std::cout << ">>>>> fitBackgrounds_init()" << std::endl;

    //TVirtualFitter::SetDefaultFitter("Minuit2");
    
    // TODO:
    // want to change Minuit so that it does not store disabled parameters
    // do this using a map to translate from parameter number to minuit
    // parameter numbers
    // and change numberParams to numberEnabled params
    //TMinuit *minuit = new TMinuit(numberParams);
    std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
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
    print_paramNameMap();

    std::cout << "set errors" << std::endl;

    // Set the parameters for the fit, and give them an arbitrary 10% error to start.
    //for(int i = 0; i < numberParams; i++)
    //int i_minuit = 0; // internal minuit counter
    for(int i = 0; i < numberParams; i++)
    {
        // internal (minuit) parameter number
        int minuit_param_number = -1;

        // check if parameter enabled
        if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
        {
            // NOT enabled
            // ignore
            continue;
        }
        else
        {
            // is enabled
            // do nothing (exec code in following block)

            // set internal parameter number
            minuit_param_number = paramNumberToMinuitParamNumberMap.at(i);

            /*
            std::cout << "contents of paramNumberToMinuitParamNumberMap" << std::endl;
            for(auto it = paramNumberToMinuitParamNumberMap.cbegin(); it != paramNumberToMinuitParamNumberMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
            std::cin.get();
            */
        }
            

        // AdjustActs[i] = 1.;
        //AdjustActs_Err[i] = 1.;
        TString i_str;
        i_str.Form("%i", i);
        TString minuit_param_number_str;
        minuit_param_number_str.Form("%i", minuit_param_number);
        std::cout << "DefineParameter: i=" << i << " -> minuit_param_number=" << minuit_param_number << std::endl;
        //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 1000.0);


        if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            // define parameter using constrained value if hard constrained

            std::cout << "minuit: fixed parameter i=" << i << std::endl;
            TString minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_FIXED";
            
            //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 2.0);
            //minuit->DefineParameter(i, "_" + i_str + "_FIXED", param_init_value, param_init_error, -2.0, 1.0e+06);
            //minuit->FixParameter(i);
            
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", param_init_value, param_init_error);
            // TODO: change such that samples are pre-scaled by activity input value
            
            // 2020-06-17
            if(i == 1)
            {
//                std::cout << "i 1 fixed " << AdjustActs[i] << std::endl;
//                std::cin.get();
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", AdjustActs[i], 0.5, 0.0, 50.0);
                //TString minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_FIXED";
                theParameterState.Add(std::string(minuit_param_name), AdjustActs[i], AdjustActs_Err[i]); // instead of _Err was 0.5
                //theParameterState.SetLowerLimit(i, 0.0); // no limit for xi
            }
            else
            {
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
                theParameterState.Add(std::string(minuit_param_name), 1.0, 0.5);
                theParameterState.SetLowerLimit(i, 0.0);
            }
            /*
            if(i == 1)
            {
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", AdjustActs[i], 0.5, 0.0, 50.0);
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
            }
            */

            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_FIXED", 1.0, 0.5, 0.0, 50.0);
            //minuit->FixParameter(minuit_param_number);
            theParameterState.Fix(std::string(minuit_param_name));
        }
        else
        {
            // define parameter using initial value if free/soft constrained
            
            std::cout << "minuit: parameter i=" << i << " is enabled and not fixed, leaving free" << std::endl;
            TString minuit_param_name = "_" + i_str + "_" + minuit_param_number_str + "_";

            //minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 1.0e+06);
            
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", param_init_value, param_init_error, 0.0, 1.0e+05);
            // TODO: change such that samples are pre-scaled by activity input value
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 50.0);
            //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);

            // 2020-06-17
            
            if(i == 1)
            {
//                std::cout << "i 1 not fixed " << AdjustActs[i] << std::endl;
//                std::cin.get();
                // TODO: fix this
                // does not work if xi_31 paramter is not number 1
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", xi_31_init, 0.5, 0.0, 1000.0);
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, 0.0, 1000.0);
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, -1.0, 5.0);
                theParameterState.Add(std::string(minuit_param_name), AdjustActs[i], AdjustActs_Err[i]); // instead of _Err was 0.5
                theParameterState.SetLowerLimit(i, -0.4);
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, -0.4, 5.0);
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, -1.0, 5.0);
                //std::cout << "define parameter" << AdjustActs[i] << std::endl;
                //std::cin.get();
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 2.0 * AdjustActs[i], 0.5, -1.0, 5.0);
            }
            else
            {

            //    if(i == 0)
            //    {
            //        //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
            //        minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.705, 0.005, 0.0, 1000.0);
            //
            //        //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 2.5, 0.5, 0.0, 1000.0);
            //    }
            //    else
            //    {
                    //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
                    theParameterState.Add(std::string(minuit_param_name), 1.0, 0.5);
                    theParameterState.SetLowerLimit(i, 0.0);
            //    }
                //if(i == 0)
                //{
                //    //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
                //    minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.705, 0.005, 0.0, 1000.0);
                //    //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 2.5, 0.5, 0.0, 1000.0);
                //}
                //else
                //{
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
                //}
            }
            
            /*
            if(i == 1)
            {
                // TODO: fix this
                // does not work if xi_31 paramter is not number 1
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", xi_31_init, 0.5, 0.0, 1000.0);
                //minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", AdjustActs[i], 0.5, 0.0, 1000.0);
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, -1.0, 5.0);
            }
            else
            {
                minuit->DefineParameter(minuit_param_number, "_" + i_str + "_" + minuit_param_number_str + "_", 1.0, 0.5, 0.0, 1000.0);
            }
            */


            // TODO: set initial error using initError/initValue
            // TODO: limits were set to 50. minuit trying to exceed 50 for some
            // backgrounds when only using nd, mo bb, zr bb, ca bb an Ca Y90
            // allow larger range of values, check fit result
            // then put back and include more backgrounds, does fit result
            // still try and exceed 50x
            // TODO: there is still the issue of the discontinuity in
            // the chisquare plots. MPS. is this due to Poisson function
            // failing with negative events?
            
        }
        // TODO: I have changed the input, now I have to change the output
        // and also the initial scaling given to histograms

        // note, no longer need this section, because only enabled parameters
        // get through to this execution block
        // if parameter is enabled, then it is either fixed or free
        // it has been copied and edited above
    /*
        // TODO: not sure if this is the correct thing to do
        // I think that disabled params do not necessarily appear in the
        // fixed params list, and therefore I have to fix them (not sure
        // about that) here if they are disabled
        if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
        {
            // NOT enabled, therefore fix
            
            std::cout << "minuit: disabled parameter i=" << i << std::endl;
         
            // leave parameter (amplitude) at default value of 1.0, which
            // is more likely to show up any errors in code
            minuit->DefineParameter(i, "_" + i_str + "_DISABLED", 1.0, 0.1, 0.0, 2.0);
            //minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 2.0);
            minuit->FixParameter(i);

            // TODO: this code never executed, because disabled parameters are
            // always fixed? or always set as "hard" in parameter list file?
        }
        else if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            // define parameter using constrained value if hard constrained

            std::cout << "minuit: fixed parameter i=" << i << std::endl;
            
            //minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 2.0);
            minuit->DefineParameter(i, "_" + i_str + "_FIXED", param_init_value, param_init_error, -2.0, 50000.0);
            minuit->FixParameter(i);
        }
        else
        {
            // define parameter using initial value if free/soft constrained
            
            std::cout << "minuit: parameter i=" << i << " is enabled and not fixed, leaving free" << std::endl;

            minuit->DefineParameter(i, "_" + i_str + "_", param_init_value, param_init_error, 0.0, 10000.0);
            
        }
        // TODO: what about disabled parameters
    */

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




    // Array to pass arguments to fitter
    //double arglist[100];
    //arglist[0] = 0;

    // Set print level
    //  minuit->ExecuteCommand("SET PRINT",arglist,2);
    std::cout << "Ready to exec fix" << std::endl;

    // Do minimisation
    // arglist[0] = 50000;  // number of function calls
    // arglist[1] = 0.1;  // tolerance

    //minuit->SetErrorDef(0.5);
    // TODO ? required ?
    // 1.0 = chisquare
    // 0.5 = negative log likelihood
    

    // mnsimp()?

    
    // MARKER
    // disable the 150 Nd gA parameter
    //minuit->FixParameter(1);




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

    std::cout << "return" << std::endl;

    //return minuit;
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
//                std::cout << "i 1 exec: " << AdjustActs[i] << std::endl;
//                std::cin.get();

    // TODO: re-enable
    ll_walk.clear();
    //minuit->Migrad();
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



void fitBackgrounds_postexectest(TMinuit *minuit, double *AdjustActs, double *AdjustActs_Err)
{


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


}


void fitBackgrounds_getcovmatrix(TMinuit* minuit, double *&CovMatrix, int& number_free_params)
{

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

}




//void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
//TMinuit* fitBackgrounds(double *AdjustActs, double *AdjustActs_Err, double *&CovMatrix, int& number_free_params, Int_t thePhase)
ROOT::Minuit2::FunctionMinimum fitBackgrounds(
    ROOT::Minuit2::MnUserParameterState& theParameterState,
    ROOT::Minuit2::VariableMetricMinimizer& theMinimizer,
    MinimizeFCNAxialVector &theFCN,
    double *AdjustActs,
    double *AdjustActs_Err,
    double *&CovMatrix,
    int &number_free_params,
    Int_t thePhase
    )
{

    std::cout << ">>>>> fitBackgrounds()" << std::endl;
    
    //TMinuit* minuit = fitBackgrounds_init(AdjustActs, AdjustActs_Err);
    fitBackgrounds_init(theParameterState, theMinimizer, AdjustActs, AdjustActs_Err);
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


}



#endif // NEWLOGLIKFITTER_FITBACKGROUNDS_H
