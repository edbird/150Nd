#ifndef MINIMIZEFCNAXIALVECTOR_H
#define MINIMIZEFCNAXIALVECTOR_H



// Root headers, minuit
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"



Double_t getNumberMC1D(const Int_t channel, const Int_t bin_ix, const std::vector<double> &p);
Double_t getNumberMC2D(const Int_t channel, const Int_t bin_ix, const Int_t bin_iy, const std::vector<double> &p);



// TODO: may need to change precision to float not double using
// MnMachinePrecision::setPrecision

class MinimizeFCNAxialVector : public ROOT::Minuit2::FCNBase
{


    public:


    MinimizeFCNAxialVector()
        : error_def{0.5} // LL // TODO: check as I am actually minimizing chisquare
    {
    }

    virtual
    ~MinimizeFCNAxialVector()
    {
    }

    virtual
    double Up() const
    {
        return error_def;
    }

    virtual
    double operator()(const std::vector<double> &param) const
    {
       
        Int_t nPar = param.size();
        Double_t fval = 0.0;

        if(allFakeDataSamples1D == nullptr)
        {
            build_fake_data();
        }


        // save the walk
        ll_walk.push_back(std::make_pair(param.at(1), param.at(0)));


        bool debugprint = false;
        bool mode_fake_data = false; //true;


        // draw the output
        TString fname;
        fname.Form("lliter_%d", counter);
        //draw_channel(1, p, std::string(fname));



        // error mode
        // 1 = data
        // 2 = MC
        // 3 = both in quadrature
        const int EMODE = 2;


        int axial_vector_parameter_0_param_number = get_axial_vector_parameter_index(); 


        if(debugprint)
        {
            std::cout << std::scientific;
            std::cout << "logLikelihood" << std::endl;
            std::cout << "param[0]=" << param[0] << " param[" << axial_vector_parameter_0_param_number << "]="
                      << param[axial_vector_parameter_0_param_number] << std::endl;
        }


        // TODO: will not work if parameter number changes
        // 2020-06-17
        //if(p[1] != last_xi_31_parameter_value)
        //if(p[1] != paramLastValueMap[1])
        if(param[axial_vector_parameter_0_param_number] != paramLastValueMap[axial_vector_parameter_0_param_number])
        {

            // TODO: rebuild nd150 xi_31 paramter histogram here
            //std::cout << "rebuilding 150 Nd MC" << std::endl;

            const double xi_31{param[axial_vector_parameter_0_param_number]};
            rebuild_150Nd_MC(xi_31, xi_31_baseline);


    /*
            TCanvas *ctmp = new TCanvas("ctmp", "ctmp");
            hSingleEnergy->Draw();
            TString fname;
            fname.Form("ctmp_%d.png", counter);
            ctmp->SaveAs(fname);
            ++ counter;
    */


            ///////////////////////////////////////////////////////////////////////
            // set last value
            ///////////////////////////////////////////////////////////////////////
            

            last_xi_31_parameter_value = xi_31;
        }


        // TODO: add check here to see if any disabled parameters are accessed


        ///////////////////////////////////////////////////////////////////////////
        // loglikelihood, 1D channels
        ///////////////////////////////////////////////////////////////////////////


        double loglik = 0.0;

    //   std::cout << "getting 1D histograms" << std::endl;

        TH1F *tmpData1D;
        TH1F *tmpFakeData1D;
        // std::cout << allDataSamples1D->GetEntries()  << std::endl;

        // there are i samples for each channel
        for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
        {

            if(channel_enable_1D[channel] == 0)
            {
                if(debugprint)
                {
                    std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
                }
                continue;
            }
            

            double ll_channel = 0.0;

            TString i_str;
            i_str.Form("%i", channel);
            //std::cout << i << std::endl;

            tmpData1D = (TH1F*)allDataSamples1D->At(channel);
            tmpFakeData1D = (TH1F*)allFakeDataSamples1D->At(channel);

            // std::cout << tmpData1D->Integral() << std::endl;

            int nBinsX = tmpData1D->GetNbinsX();
            for(int bin_ix = 1; bin_ix <= nBinsX; ++ bin_ix)
            {
                //Int_t nData = 0;
                Double_t nData = 0;
                if(mode_fake_data == false)
                {
                    nData = tmpData1D->GetBinContent(bin_ix);
                }
                if(mode_fake_data == true)
                {
                    nData = tmpFakeData1D->GetBinContent(bin_ix);
                }
                // i is the index of the sample ?
                // ix is the bin index
                // p is a pointer to an array of parameter values

                // fake data
                //double nFakeData = getNumberMC1D(channel, bin_ix, paramInitValueMap);
                //for(int k = 0; k < numberEnabledParams; ++ k)
                //{
                //    std::cout << "minuitParamInit[" << k << "]=" << minuitParamInit[k] << " p[" << k << "]=" << p[k] << std::endl;
                //}
                //double nFakeData = getNumberMC1D(channel, bin_ix, minuitParamInit);
                // TODO: double or int
                // TODO: I don't know if this works. param index might be different
                // internal vs external issue
                // the index for paramInitValueMap must be an external param index
                // this is because paramInitValueP1Map and paramInitValueP2Map
                // are indexed using an external param index and therefore
                // paramInitValueMap is also
                // The function getNumberMC1D uses internal index format

                double nMC = getNumberMC1D(channel, bin_ix, param);

                //std::cout << "bin_ix=" << bin_ix << " nFakeData=" << nFakeData << " nMC=" << nMC << std::endl;
                //std::cout << "for bin_ix=" << bin_ix << " nMC=" << nMC << " nData=" << nData << std::endl;
                

                // poisson is a probability, so values are between 0 and 1 which means that
                // log of this value is always negative
                // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                if(nMC >= 0.0)
                {
                    double lp = logpoisson(nData, nMC);
                    //double lp = logpoisson(nFakeData, nMC);
                    //double lp = logpoisson_sterling(nFakeData, nMC);
                    ll_channel += lp;

                    if(debugprint)
                    {
                        std::cout << "bin_ix=" << bin_ix << " lp=" << lp << " nData=" << nData << " nMC=" << nMC << " (1)" << std::endl;
                    }
                }
                // not sure we are dealing with zero bins correctly, should
                // ignore?
                // this appears to happen a lot
                else
                {
                    std::cout << "MC WENT LOWER THAN ZERO" << std::endl;

                    double lp = logpoisson(nData, 0.0);
                    //double lp = logpoisson(nFakeData, 1.0e-05);
                    //double lp = logpoisson_sterling(nFakeData, 1.0e-05);
                    ll_channel += lp;

                    if(debugprint)
                    {
                        std::cout << "bin_ix=" << bin_ix << " lp=" << lp << " nData=" << nData << " nMC=" << nMC << " (2)" << std::endl;
                    }
                }
            
            } //~bins

            if(debugprint)
            {
                std::cout << "1D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
            }
            loglik += ll_channel;

        } //channels
        

        // add contribution from g_A related histograms (single electron energy)
        // either as single electron
        // or as high+low energy 1D
        // or as high+low energy 2D


        ///////////////////////////////////////////////////////////////////////////
        // loglikelihood, 2D channels
        ///////////////////////////////////////////////////////////////////////////

        
        //  std::cout << "getting 2D histograms" << std::endl;

        TH2F *tmpData2D;
        TH2F *tmpFakeData2D;
        // std::cout << allDataSamples2D->GetEntries()  << std::endl;

        // there are i samples for each channel
        for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
        {

            if(channel_enable_2D[channel] == 0)
            {
                if(debugprint)
                {
                    std::cout << "2D: channel " << channel << " disabled, skip" << std::endl;
                }
                continue;
            }

            double ll_channel = 0.0;

            TString i_str;
            i_str.Form("%i", channel);
            //std::cout << i << std::endl;

            // TODO: can I remove this Clone() call safely to improve speed?
            //tmpData1D = (TH1F*)allDataSamples1D->At(i)->Clone("tmpData1D" + i_str + "_");
            tmpData2D = (TH2F*)allDataSamples2D->At(channel);
            tmpFakeData2D = (TH2F*)allFakeDataSamples2D->At(channel);

            // std::cout << tmpData2D->Integral() << std::endl;

            int nBinsX = tmpData2D->GetNbinsX();
            int nBinsY = tmpData2D->GetNbinsY();
            for(int bin_ix = 1; bin_ix <= nBinsX; ++ bin_ix)
            {
                for(int bin_iy = 1; bin_iy <= nBinsY; ++ bin_iy)
                {
                    //Int_t nData = (Int_t)tmpData2D->GetBinContent(bin_ix, bin_iy);
                    // i is the index of the sample (actually histogram type/channel)
                    // ix is the bin index
                    // p is a pointer to an array of parameter values
                    double nData = 0.0;
                    if(mode_fake_data == false)
                    {
                        nData = tmpData2D->GetBinContent(bin_ix, bin_iy);
                    }
                    if(mode_fake_data == true)
                    {
                        nData = tmpFakeData2D->GetBinContent(bin_ix, bin_iy);
                    }

                    
                    double nMC = getNumberMC2D(channel, bin_ix, bin_iy, param);

                    if(nMC >= 0.0)
                    {
                        double lp = logpoisson(nData, nMC);
                        ll_channel += lp;
                    }
                    else
                    {
                        std::cout << "2D MC WENT LOWER THAN ZERO" << std::endl;

                        double lp = logpoisson(nData, 0.0);
                        ll_channel += lp;
                    }

                } // binsY
            } // binX

            if(debugprint)
            {
                std::cout << "2D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
            }
            loglik += ll_channel;

        } // channels

     
        // add constraints to improve likelihood
        //will eventually add gaussian constraint
        // disabled
        /*
        for ( int i = 0; i < 14; i++ ) {
            double constraint = 0.;
            //constraint = pow( (p[i] - 1) / (0.01), 2);
            // std::cout << i << "   " << p[i];
            bool fixed = false;
            for ( int j = 0; j < fixed_params.size(); j++ ) {
                if ( i == fixed_params.at(j) ) fixed = true;
            }

            //if ( (i <2) || (i == 7 ) || (i==10) || (i == 14) || (i==15) || (i==16) || (i == 18) || (i == 19) || (i == 21) || 9i==22) ) {
            // constraint = pow( (p[i] - 1.) / (paramActErrMap[i].front() / paramActMap[i].front()), 2);
            constraint = TMath::Gaus(p[i],1.,(paramActErrMap[i].front() / paramActMap[i].front()),true);
          

            // std::cout << "   " << constraint << std::endl;
            //  loglik += log(constraint);  
            // }
        }
        */

        int mode = MODE_PARAM_UNDEFINED;

        // penalty terms section
        double penalty_term = 0.0;

        // TODO: I don't like this - should loop over the enabled params
        // however, this should still work as it is
        for(int i = 0; i < numberParams; ++ i)
        {

            if(std::find(enabled_params.begin(), enabled_params.end(), i) == enabled_params.end())
            {
                //std::cout << "parameter number " << param_number << " is disabled" << std::endl;
                //std::cout << "ERROR: i=" << i << " - parameter is DISABLED" << std::endl;
                //std::cin.get();
                continue;
            }


            if(thePhase == 0)
            {
                mode = paramConstrainModeP1Map[i];

                if(mode == MODE_PARAM_SOFT)
                {
                    // do nothing, soft constraint will be applied below
                }
                else if(mode == MODE_PARAM_HARD)
                {
                    // parameter fixed by minuit, continue to next param
                    
                    // NOTE: changed to ignore HARD
                    continue;
                }
                else if(mode == MODE_PARAM_FREE)
                {
                    // no constraint to apply, continue to next param
                    continue;
                }
                else
                {
                    std::cout << "ERROR: Invalid value in paramConstrainModeP1Map: paramConstrainModeP1Map[" << i << "]=" << paramConstrainModeP1Map[i] << std::endl;
                }
            }
            else if(thePhase == 1)
            {
                mode = paramConstrainModeP2Map[i];

                if(mode == MODE_PARAM_SOFT)
                {
                    // do nothing, soft constraint will be applied below
                }
                else if(mode == MODE_PARAM_HARD)
                {
                    // parameter fixed by minuit, continue to next param
                    
                    // NOTE: changed to ignore HARD
                    continue;
                }
                else if(mode == MODE_PARAM_FREE)
                {
                    // no constraint to apply, continue to next param
                    continue;
                }
                else
                {
                    std::cout << "ERROR: Invalid value in paramConstrainModeP2Map: paramConstrainModeP2Map[" << i << "]=" << paramConstrainModeP2Map[i] << std::endl;
                }
            }
            else
            {
                std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
            }
            
            // soft constraint is applied here
            double constraint = 0.;
            double error = 0.;
            get_paramConstraintValueError(thePhase, i, constraint, error);
            // NOTE: these values read from parameter list file and thus are in
            // units of activity (Bq)

            // TODO: can optimize this code

            if(error < 0.0)
            {
                std::cout << "ERROR: Invalid error value: error=" << error << std::endl;
            }

            // check if hard parameter and error == 0.0
            /*
            if(mode == MODE_PARAM_HARD)
            {
                if(error == 0.0)
                {
                    // this parameter is a "constant" (according to minuit)
                    // so ignore
                    continue;
                }
                else
                {
                    // do nothing, add constraint for hard param
                }
            }
            */

            //double param_value = p[i];
            // i is index of the parameter number (external / non minuit)
            // convert to internal parameter number (minuit param number)
            int j = paramNumberToMinuitParamNumberMap.at(i);
            if(j < numberEnabledParams)
            {
                // ok
            }
            else
            {
                throw std::runtime_error("error: invalid value of j (internal param number)");
            }
            double param_value = param[j];
            // this parameter is from minuit internal and thus is in minuit
            // internal units (not Bq)
            // have to convert to Bq units
        
            // convert to Bq
            // multiply by the initial value
            // activity_value_Bq should really be called param_initial_value
            double activity_value_Bq = 0.0;
            double tmp_err;
            get_paramInitValueError(thePhase, i, activity_value_Bq, tmp_err);

            double value = param_value * activity_value_Bq;
            //double penalty = std::pow((param_value - constraint) / error, 2.0);
            double penalty = 0.0;

            // TODO
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
            
            penalty = std::pow((value - constraint) / error, 2.0);


            // TODO: is this the correct error term?
            // error on constraint rather than error on current fit value?

            penalty_term += penalty;
        }


        if(debugprint)
        {
            std::cout << "penalty_term=" << penalty_term << std::endl;
        }
        if(penalty_term != 0.0)
        {
            std::cout << "penalty_term=" << penalty_term << std::endl;
        }
      
        fval = -2.0 * loglik + penalty_term;

        //tmpData->Delete();

        // hook
        global_chisquare = fval;


        // set last parameter values
        // could also loop over nPar?
        for(int i = 0; i < param.size(); ++ i)
        {
            double value = param.at(i);
            paramLastValueMap[i] = value;
            minuitParamCurrent[i] = value;
            minuitParamLast[i] = value;
        }

        return fval;
    }



    protected:

    double error_def;


};




///////////////////////////////////////////////////////////////////////////////
// get number MC 1D
///////////////////////////////////////////////////////////////////////////////

Double_t getNumberMC1D(
    const Int_t channel,
    const Int_t bin_ix,
    const std::vector<double> &p
    //const Double_t *const p
    )
{


    //std::cout << "getNumberMC1D() called with channel=" << channel << " binx=" << binx << std::endl;
    //std::cout << "printing contents of parameter array" << std::endl;
    //for(int i = 0; i < numberParams; ++ i)
    //{
    //    std::cout << "i=" << i << " p[i]=" << p[i] << std::endl;
    //}


    // std::cout <<"getting number of1D MC... "  <<channel << std::endl;

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1F *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;

    //std::cout << "channel=" << channel << std::endl;
    //std::cout << allMCSamples1D[channel]->GetEntries() << std::endl;


    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); j++)
    {

        // disabled params do not appear in this array

        tmpHist = (TH1F*)allMCSamples1D[channel]->At(j);

        TString tmpName = tmpHist->GetName();

        int which_param = -1;
        bool found_param = false;

        //std::cout << "NEW CODE" << std::endl;
        //try
        //{
        // TODO: remove TString
        {
            std::string tmp_hist_name(tmpName);
            auto i_start = tmp_hist_name.find('_') + 1;
            auto i_end = tmp_hist_name.rfind('_');
            if(i_end - i_start > 0)
            {
                std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);
                //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                {
                    int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                    // TODO: removed std::string, change tmpName type to be std::string from TString
                
                    //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                    //which_param = paramNumber;
                    which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                    found_param = true;

                    //std::cout << tmp_sample_name << " -> " << paramNumber << std::endl;

                    //std::cout << "DEBUG: found parameter with minuit (internal) number: " << which_param << std::endl;
                    //std::cin.get();
                }
                else
                {
                   std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                }
            }
        }
        /*
        }
        catch(std::exception &e)
        {
            std::cout << "e.what(): " << e.what() << std::endl;
            std::cout << "tmpName=" << tmpName << std::endl;
            std::cout << "contents of map" << std::endl;
            for(auto it = MCNameToParamNumberMap.cbegin(); it != MCNameToParamNumberMap.cend(); ++ it)
            {
                std::cout << it->first << " -> " << it->second << std::endl;
            }
        }
        std::cin.get();
        */

        //std::cout << "bin_ix=" << bin_ix << " tmpHist->GetName()=" << tmpHist->GetName() << " which_param=" << which_param << std::endl;

        /*
        if(foundParam)
        {
            //std::cout << "adding to nMC with index of which_param=" << which_param << std::endl;
            // TODO: think this is collecting the wrong parameter? or is it?
            nMC += p[which_param] * tmpHist->GetBinContent(binx);
        }
        else
        {
            std::cout << "error could not find histogram: " << tmpName << std::endl;
        }
        */

        if(found_param == true)
        {

            // check here to see if param is disabled
            // TODO
            //
            // since which_param must be an enabled parameter, no longer need this
            // leave for now as a check?
            //if(std::find(enabled_params.begin(), enabled_params.end(), which_param) == enabled_params.end())
            if(std::find(enabled_params.begin(), enabled_params.end(), minuitParamNumberToParamNumberMap.at(which_param)) == enabled_params.end())
            {
                std::cout << "ERROR: which_param=" << which_param << " - parameter is DISABLED" << std::endl;
                std::cin.get();
            }

            if(which_param == 1)
            {
                std::cout << "Error: which_param == 1 !" << std::endl;
                std::cin.get();
            }

            nMC += p[which_param] * (double)tmpHist->GetBinContent(bin_ix);

        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        } 

    }

    //std::cin.get();

    return nMC;

}




Double_t getNumberMC2D(
    const Int_t channel,
    const Int_t bin_ix,
    const Int_t bin_iy,
    //const Double_t *const p
    const std::vector<double> &p
    )
{

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1F *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;
    //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples2D[channel]->GetEntries(); j++)
    {

        // get histogram name
        tmpHist = (TH1F*)allMCSamples2D[channel]->At(j);
        TString tmpName = tmpHist->GetName();

        // set paramter to default value
        int which_param = -1;
        bool found_param = false;

        {
            std::string tmp_hist_name(tmpName);
            auto i_start = tmp_hist_name.find('_') + 1;
            auto i_end = tmp_hist_name.rfind('_');
            if(i_end - i_start > 0)
            {
                // strip off WHAT - TODO comment
                std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);
                //std::cout << "tmp_sample_name=" << tmp_sample_name << std::endl;
                //std::cin.get(); // check still working

                if(MCNameToParamNumberMap.count(tmp_sample_name) > 0)
                {
                    int paramNumber = MCNameToParamNumberMap.at(tmp_sample_name);
                
                    // convert to minuit (internal) parameter number
                    which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                    found_param = true;

                    //std::cout << "DEBUG: found parameter with minuit (internal) number: " << which_param << std::endl;
                    //std::cin.get();
                }
                else
                {
                   std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                }
            }
        }

        if(found_param == true)
        {

            // check here to see if param is disabled
            // TODO
            //
            // since which_param must be an enabled parameter, no longer need this
            // leave for now as a check?
            //if(std::find(enabled_params.begin(), enabled_params.end(), which_param) == enabled_params.end())
            if(std::find(enabled_params.begin(), enabled_params.end(), minuitParamNumberToParamNumberMap.at(which_param)) == enabled_params.end())
            {
                //std::cout << "parameter number " << param_number << " is disabled" << std::endl;
                std::cout << "ERROR: which_param=" << which_param << " - parameter is DISABLED" << std::endl;
                std::cin.get();
            }

            if(which_param == 1)
            {
                std::cout << "Error: which_param == 1 !" << std::endl;
                std::cin.get();
            }

            // add to number of MC
            nMC += p[which_param] * (double)tmpHist->GetBinContent(bin_ix, bin_iy);
        }
        else
        {
            std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
        }

    } // j

    return nMC;


}





#endif // MINIMIZEFCNAXIALVECTOR_H
