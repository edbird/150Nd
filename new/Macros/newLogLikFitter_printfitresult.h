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

}


#endif // NEWLOGLIKFITTER_PRINTFITRESULT_H
