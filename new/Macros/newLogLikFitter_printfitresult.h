#ifndef NEWLOGLIKFITTER_PRINTFITRESULT_H
#define NEWLOGLIKFITTER_PRINTFITRESULT_H


void print_adjustacts(std::ofstream &os)
{


	os << "The following adjustments (in minuit parameter units) should be made:" << std::endl;
 
	os << "Note that gA (1) is a special parameter" << std::endl;

    // TODO: after homogenizing gA is no longer a special parameter
    //for(int i = 0; i < numberParams; i++)
    //std::ofstream myFileFitResults("fit_results.txt", std::ios::out | std::ios::app);
    timestamp(myFileFitResults);
    for(int i = 0; i < numberEnabledParams; i++)
    {
    	os << i << " :\t" << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;
       	//myFileFitResults << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;
    }
	//myFileFitResults.close();

    os << "The following adjustments (in units of Bq) should be made:" << std::endl;
    //for(int i = 0; i < numberParams; i++)
    //std::cout << "numberEnabledParams=" << numberEnabledParams << std::endl;
    for(int i = 0; i < numberEnabledParams; i++)
    {
        //if(i == 1)
        //{
        //    std::cout << "skipping gA parameter TODO FIX" << std::endl;
        //}

        int j = minuitParamNumberToParamNumberMap.at(i);
        //std::cout << "i=" << i << " j=" << j << std::endl;
        double param_init_value = 0.;
        double param_init_error = 0.; 
        get_paramInitValueError(thePhase, j, param_init_value, param_init_error);
        //std::cout << "value=" << param_init_value << " err=" << param_init_error << " AdjustActs[i]=" << AdjustActs[i] << std::endl;
        
        // 2020-06-17
        
        if(i != 1)
        {
            os << i << " :\t" << AdjustActs[i] * param_init_value
                           << " +- " << AdjustActs_Err[i] * param_init_value;
            // TODO: put the mutiplication by xi_31_init INSIDE the reweight/fit functions,
            // to restore uniformity in minuit parameters
        }
        else
        {
            os << i << " :\t" << AdjustActs[i]
                           << " +- " << AdjustActs_Err[i];
        }
        
        /*
        std::cout << i << " :\t" << AdjustActs[i] * param_init_value
                       << " +- " << AdjustActs_Err[i] * param_init_value;
        */
        // TODO: put the mutiplication by xi_31_init INSIDE the reweight/fit functions,
        // to restore uniformity in minuit parameters
        // xi works as an additive parameter so cannot do this, in case where HSD
        // is zero does not work


        Double_t change = 0.0;
        // 2020-06-17
        if(i != 1)
        {
            change = 100.0 * (AdjustActs[i] - 1.0);
        }
        else
        {
            // 2020-06-17
            //change = 100.0 * (AdjustActs[i] - xi_31_init);
            change = 100.0 * (AdjustActs[i] - xi_31_init_value);
            //change = 100.0 * (AdjustActs[i] - 1.0);
        }
        /*change = 100.0 * (AdjustActs[i] - 1.0);*/
        if(change >= 0)
        {
            os << " -> +";
        }
        else
        {
            os << " -> -";
        }
        os << std::abs(change) << " %" << std::endl;  
    }

}


#endif // NEWLOGLIKFITTER_PRINTFITRESULT_H
