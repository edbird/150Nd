#ifndef NEWLOGLIKFITTER_LOGLIKELIHOOD_H
#define NEWLOGLIKFITTER_LOGLIKELIHOOD_H


///////////////////////////////////////////////////////////////////////////////
// function declarations
///////////////////////////////////////////////////////////////////////////////


void logLikelihood(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */); 
Double_t getNumberMC1D(const Int_t channel, const Int_t bin_ix, const Double_t *const p);
Double_t getNumberMC2D(const Int_t channel, const Int_t bin_ix, const Int_t bin_iy, const Double_t *const p);





///////////////////////////////////////////////////////////////////////////////
// logLikelihood
///////////////////////////////////////////////////////////////////////////////


// TODO don't appear to work with parameters with more than one MC
void logLikelihood(Int_t & nPar, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{




    // 2020-06-17
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





    std::cout << std::scientific;
    std::cout << "logLikelihood" << std::endl;
    //std::cout << "p[0]=" << p[0] << " p[1]=" << p[1] << std::endl;
    std::cout << "p[0]=" << p[0] << " p[" << axial_vector_parameter_0_param_number << "]="
              << p[axial_vector_parameter_0_param_number] << std::endl;
    // TODO: use non fixed parameter number index


    // TODO: will not work if parameter number changes
    // 2020-06-17
    //if(p[1] != last_xi_31_parameter_value)
    //if(p[1] != paramLastValueMap[1])
    if(p[axial_vector_parameter_0_param_number] != paramLastValueMap[axial_vector_parameter_0_param_number])
    // TODO: will not work if xi parameter is not param number 1
    {

        // TODO: rebuild nd150 xi_31 paramter histogram here


        ///////////////////////////////////////////////////////////////////////
        // reweight hMinMaxEnergy_
        // reweight all
        ///////////////////////////////////////////////////////////////////////


        // new code to reweight 150Nd by xi_{31} parameter

        // list of MC index to reweight, 1D
        //std::vector<int> reweight_index_1D[number1DHists];

        // list of MC index to reweight, 2D
        //std::vector<int> reweight_index_2D[number2DHists];

        // pointers of histograms to pass to reweight function
        TH1F *hTotalE = nullptr;
        TH1F *hSingleEnergy = nullptr;
        TH1F *hHighEnergy = nullptr;
        TH1F *hLowEnergy = nullptr;
        TH1F *hEnergySum = nullptr;
        TH1F *hEnergyDiff = nullptr;
        TH2F *hHighLowEnergy = nullptr;

        //std::map<int, TH1F*> channel_to_pointer_map_1D;
        //std::map<int, TH2F*> channel_to_pointer_map_2D;

        // search through each channel for 150nd samples
        //for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
        for(int channel = 0; channel < number1DHists; ++ channel)
        {

            // search through each sample for this channel
            for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
            {
                TH1F *tmpHist = (TH1F*)allMCSamples1D[channel]->At(i);
                TString tmpHist_name = tmpHist->GetName();
                // TODO: had to add "_fit" here - might not work after 1 iteration
                //if(tmpHist_name.CompareTo("hTotalE_nd150_rot_2n2b_m4_fit") == 0 ||
                //   tmpHist_name.CompareTo("hTotalE_nd150_rot_2b2n_m4_fit") == 0)

                // if name contains this string, it needs to be reweighted
                if(tmpHist_name.Contains("nd150_rot_2n2b_m4") ||
                   tmpHist_name.Contains("nd150_rot_2b2n_m4"))
                {
                    // select correct histogram
                    /*if(tmpHist_name.Contains("hTotalE"))
                    {
                        //hTotalE = tmpHist;
                        //channel_to_pointer_map[channel] = tmpHist;
                    }
                    else if(tmpHist_name.Contains("hSingleEnergy"))
                    {
                        //hSingleEnergy = tmpHist;
                    }
                    else if(tmpHist_name.Contains("hLowEnergy"))
                    {
                        //hLowEnergy = tmpHist;
                    }
                    else if(tmpHist_name.Contains("hHighEnergy"))
                    {
                        //hHighEnergy = tmpHist;
                    }

                    // signal that this index needs to be removed
                    reweight_index_1D[channel].push_back(i);
                    */

                    // remove histogram
                    allMCSamples1D[channel]->RemoveAt(i);
                    break;
                    // TODO: will only work for a maximum of 1 histograms
                    // removed because index i will shift due to removing
                    // histogram at i

                }

            } // i

        } // channel


        // search through each channel for 150nd samples
        //for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
        for(int channel = 0; channel < number2DHists; ++ channel)
        {

            // search through each sample for this channel
            for(int i = 0; i < allMCSamples2D[channel]->GetEntries(); ++ i)
            {
                TH1F *tmpHist = (TH1F*)allMCSamples2D[channel]->At(i);
                TString tmpHist_name = tmpHist->GetName();

                // if name contains this string, it needs to be reweighted
                if(tmpHist_name.Contains("nd150_rot_2n2b_m4") ||
                   tmpHist_name.Contains("nd150_rot_2b2n_m4"))
                {
                    // select correct histogram
                    //if(tmpHist_name.Contains("hHighLowEnergy"))
                    //{
                    //    hHighLowEnergy = tmpHist;
                    //}

                    // signal that this index needs to be removed
                    //reweight_index_2D[channel].push_back(i);

                    // remove histogram
                    allMCSamples2D[channel]->RemoveAt(i);
                    break;
                    // TODO: will only work for a maximum of 1 histograms
                    // removed because index i will shift due to removing
                    // histogram at i

                }

            } // i

        } // channel


        // TODO: this will not work if parameter number changes
        //const double xi_31{xi_31_init * p[1]};
        /*
        const double xi_31_init_P1 = paramInitValueP1Map[axial_vector_parameter_0_param_number];
        const double xi_31_init_P2 = paramInitValueP2Map[axial_vector_parameter_0_param_number];
        //const double xi_31_init_P1 = paramInitValueP1Map[1];
        //const double xi_31_init_P2 = paramInitValueP2Map[1];
        double xi_31_init = 0.0;
        if(thePhase == 0)
        {
            xi_31_init = xi_31_init_P1;
        }
        else if(thePhase == 1)
        {
            xi_31_init = xi_31_init_P2;
        }
        else
        {
            std::cout << "invalid value for thePhase" << std::endl;
            throw "invalid value for thePhase";
        }
        */
        ////double param_init_value = 0.0;
        ////double param_init_error = 0.0;
        ////get_paramInitValueError(thePhase, axial_vector_parameter_0_param_number, param_init_value, param_init_error);
        //const double xi_31_init{paramInitValueMap[axial_vector_parameter_0_param_number]};
        ////const double xi_31_init{param_init_value};
        ////const double xi_31{xi_31_init * p[axial_vector_parameter_0_param_number]};
        const double xi_31{p[axial_vector_parameter_0_param_number]};
        std::cout << "xi_31=" << xi_31 << " xi_31_baseline=" << xi_31_baseline << std::endl;
        //const double xi_31_baseline{0.296};
        // NOTE: 2020-06-17 this was a bug, removed

        TH1F *hWeight = nullptr;
        reweight_apply(
            hTotalE,
            hSingleEnergy,
            hHighEnergy,
            hLowEnergy,
            hEnergySum,
            hEnergyDiff,
            hHighLowEnergy,
            hWeight,
            "Nd150_2eNg_output_truth_postprocessed_small.root",
            xi_31,
            xi_31_baseline,
            h_nEqNull,
            h_nEqTwo,
            psiN0,
            psiN2,
            bb_Q);

        // TODO: just another example of manual code edits
        // make a file describing the channels to fit as well as the parameters
        allMCSamples1D[0]->Add(hTotalE);
        allMCSamples1D[1]->Add(hSingleEnergy);
        allMCSamples1D[2]->Add(hHighEnergy);
        allMCSamples1D[3]->Add(hLowEnergy);
        allMCSamples1D[4]->Add(hEnergySum);
        allMCSamples1D[5]->Add(hEnergyDiff);
        allMCSamples2D[0]->Add(hHighLowEnergy);

        /*
        for(int channel{0}; channel < number1DHists; ++ channel)
        {
            for(std::vector<int>::const_iterator it{reweight_index_1D[channel].cbegin()};
                it != reweight_index_1D[channel].cend();
                ++ it)
            {
                allMCSamples1D[channel]->RemoveAt(i);
                allMCSamples1D[channel]->Add();
            }

        }
        allMCSamples1D[channel]->RemoveAt(i);
        allMCSamples1D[channel]->Add(tmpHist_reweight);
        */



        ///////////////////////////////////////////////////////////////////////
        // set last value
        ///////////////////////////////////////////////////////////////////////
        

        //last_xi_31_parameter_value = p[1];
        // 2020-06-17
        /*
        std::string mc_name = "axial_vector_parameter_0";
        std::string search_object = MCNameToParamNameMap.at(mc_name);
        if(paramNameToNumberMap.count(search_object) > 0)
        {
            int param_number = paramNameToNumberMap.at(search_object);
       
            if(param_number != 1)
            {
                throw "param_number != 1";
            }
            const double xi_31_init_P1 = paramInitValueP1Map[param_number];
            const double xi_31_init_P2 = paramInitValueP2Map[param_number];
            //const double xi_31_init_P1 = paramInitValueP1Map[1];
            //const double xi_31_init_P2 = paramInitValueP2Map[1];
            double xi_31_init = 0.0;
            if(thePhase == 0)
            {
                xi_31_init = xi_31_init_P1;
            }
            else if(thePhase == 1)
            {
                xi_31_init = xi_31_init_P2;
            }
            else
            {
                std::cout << "invalid value for thePhase" << std::endl;
                throw "invalid value for thePhase";
            }
            last_xi_31_parameter_value = p[param_number] * xi_31_init;
            //last_xi_31_parameter_value = p[1] * xi_31_init;
        }
        else
        {
            throw "mc_name not found in paramNameToNumberMap";
        }
        */
        /*
        const double xi_31_init_P1 = paramInitValueP1Map[axial_vector_parameter_0_param_number];
        const double xi_31_init_P2 = paramInitValueP2Map[axial_vector_parameter_0_param_number];
        //const double xi_31_init_P1 = paramInitValueP1Map[1];
        //const double xi_31_init_P2 = paramInitValueP2Map[1];
        double xi_31_init = 0.0;
        if(thePhase == 0)
        {
            xi_31_init = xi_31_init_P1;
        }
        else if(thePhase == 1)
        {
            xi_31_init = xi_31_init_P2;
        }
        else
        {
            std::cout << "invalid value for thePhase" << std::endl;
            throw "invalid value for thePhase";
        }
        */
        //last_xi_31_parameter_value = p[axial_vector_parameter_0_param_number] * xi_31_init;
        last_xi_31_parameter_value = xi_31;
        //last_xi_31_parameter_value = p[1] * xi_31_init;
        // TODO: will not work if parameter number changes
        // NOTE: fixed 2020-06-17
    }



    // TODO: add check here to see if any disabled parameters are accessed



    ///////////////////////////////////////////////////////////////////////////
    // loglikelihood, 1D channels
    ///////////////////////////////////////////////////////////////////////////


    double loglik = 0.; 
    //double tmp;


//   std::cout << "getting 1D histograms" << std::endl;

    TH1F *tmpData1D;
    // std::cout << allDataSamples1D->GetEntries()  << std::endl;

    // there are i samples for each channel
    for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
    {

        if(channel_enable_1D[channel] == 0)
        {
            std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
            continue;
        }
        

        double ll_channel = 0.0;

        // allDataSamples1D only contains one object
        
        TString i_str;
        i_str.Form("%i", channel);
        //std::cout << i << std::endl;

        // TODO: can I remove this Clone() call safely to improve speed?
        //tmpData1D = (TH1F*)allDataSamples1D->At(i)->Clone("tmpData1D" + i_str + "_");
        tmpData1D = (TH1F*)allDataSamples1D->At(channel);

        // std::cout << tmpData1D->Integral() << std::endl;

        int nBinsX = tmpData1D->GetNbinsX();
        for(int bin_ix = 1; bin_ix <= nBinsX; ++ bin_ix)
        {
            Int_t nData = (Int_t)tmpData1D->GetBinContent(bin_ix);
            // i is the index of the sample ?
            // ix is the bin index
            // p is a pointer to an array of parameter values

            // TODO:
            // think there is a bug here
            // i appears to be the index of the data sample not the MC sample?
            // it becomes channel index
            double nMC = getNumberMC1D(channel, bin_ix, p);

            //std::cout << "for bin_ix=" << bin_ix << " nMC=" << nMC << " nData=" << nData << std::endl;

            //Int_t new_i = -1;
            //TString name = names.at(i);
            //new_i = paramNameToNumberMap[name];
            //std::cout << "the new i value is new_i=" << new_i << std::endl;
            //double nMC = getNumberMC1D(new_i, ix, p);

            //std::cin.get();

            
            // 2020-04-17 Notes:
            // output histograms do not look right (2d MPS plots)
            //
            // range of values for parameter include values such as
            // -8 to 10
            // -15 to 15
            // these ranges seem too large? / are indicating very large uncertainty
            // update: caused by a bug in the min/max parameter settings, now fixed
            // (not fixed in h_mps)
            //
            // there is a white square in the center. what value does this have?
            // is it negative, or zero?
            //
            // h_mps looks different to h_mps_10_0, they should be the same!
            // update: they are the same if ranges set the same and fval_min
            // is subtracted from fval before filling
            //
            // my guess was that the change in parameter values and thus n_mc
            // is having a much weaker effect compared to the penalty term
            // why is this?
            // this may not be correct, since all plots appear identical
            // indicating that something is not being computed correctly

            if(nMC > 0.)
            {
                Double_t poisson = TMath::Poisson(nData, nMC);
                if(poisson > 0.)
                {
                    //std::cout << "adding loglik value : " << TMath::Log(poisson) << " bin_ix=" << bin_ix << " poisson=" << poisson << std::endl;
	                ll_channel += TMath::Log(poisson);
                    // adding positive number makes fval go down
                    // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                    // log is taken here, should I take log of the penalty term?
                    // TODO: answer above question
                    // TODO: are there any conditions for which this can be negative?
                    // poisson is a probability, so values are between 0 and 1 which means that
                    // log of this value is always negative
                }
                else
                {
                    // MARK have not tested this yet
                    // can this ever happen? is this the correct way to deal
                    // with the problem?
                    //std::cout << "catch: poisson" << std::endl;
                    // this does not appear to happen

                    //std::cout << "adding penalty of -10. to loglik bin_ix=" << bin_ix << " nMC=" << nMC << " poisson=" << poisson << std::endl;
                    // TODO: should this be removed? check for bins where ndata = 0?

                    // TODO: there were a lot of failures here
                    //std::cout << "ERROR: failed to evaluate TMath::Poisson()=" << poisson << " -> nData=" << nData << ", nMC=" << nMC << "; bin_ix=" << bin_ix << std::endl;
	                //loglik -= 10.;
                    // subtracting positive number makes fval go up
                    // TODO: this may not be a large enough penalty

                    // TODO: removed this
                }
            }
            else
            {
                // MARK have not tested this yet
                // not sure we are dealing with zero bins correctly, should
                // ignore?
                //std::cout << "catch2: poisson" << std::endl;
                // this appears to happen a lot

                // if nMC <= 0., just add penalty and cout nothing
                // subtracting positive number makes fval go up
                //loglik -= 10.;
                // 2020-04-17: removed, should I have something here?
            }
        
        } //~bins

        std::cout << "1D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
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
    // std::cout << allDataSamples2D->GetEntries()  << std::endl;

    // there are i samples for each channel
    for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
    {

        if(channel_enable_2D[channel] == 0)
        {
            std::cout << "2D: channel " << channel << " disabled, skip" << std::endl;
            continue;
        }

        double ll_channel = 0.;

        // allDataSamples2D only contains one object
        
        TString i_str;
        i_str.Form("%i", channel);
        //std::cout << i << std::endl;

        // TODO: can I remove this Clone() call safely to improve speed?
        //tmpData1D = (TH1F*)allDataSamples1D->At(i)->Clone("tmpData1D" + i_str + "_");
        tmpData2D = (TH2F*)allDataSamples2D->At(channel);

        // std::cout << tmpData2D->Integral() << std::endl;

        int nBinsX = tmpData2D->GetNbinsX();
        int nBinsY = tmpData2D->GetNbinsY();
        for(int bin_ix = 1; bin_ix <= nBinsX; ++ bin_ix)
        {
            for(int bin_iy = 1; bin_iy <= nBinsY; ++ bin_iy)
            {
                Int_t nData = (Int_t)tmpData2D->GetBinContent(bin_ix, bin_iy);
                // i is the index of the sample (actually histogram type/channel)
                // ix is the bin index
                // p is a pointer to an array of parameter values

                double nMC = getNumberMC2D(channel, bin_ix, bin_iy, p);

                if(nMC > 0.)
                {
                    Double_t poisson = TMath::Poisson(nData, nMC);
                    if(poisson > 0.)
                    {
                        ll_channel += TMath::Log(poisson);
                        // adding positive number makes fval go down
                        // NOTE: Log(poisson) is always negative! so fval goes UP NOT DOWN
                        // log is taken here, should I take log of the penalty term?
                    }
                }

            } // binsY
        } // binX

        std::cout << "2D: channel " << channel << " enabled, ll=" << ll_channel << std::endl;
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

            //if(paramConstrainModeP1Map[i] == MODE_PARAM_SOFT)
            if(mode == MODE_PARAM_SOFT)
            {
                // do nothing, soft constraint will be applied below
            }
            //else if(paramConstrainModeP1Map[i] == MODE_PARAM_HARD)
            else if(mode == MODE_PARAM_HARD)
            {
                // parameter fixed by minuit, continue to next param
                //continue;
                // although param is fixed by minuit, we still want to add
                // the penalty term, if available (TODO: will it always
                // be available?)
                // can check if error == 0.0 first
                
                // NOTE: changed to ignore HARD
                continue;
            }
            //else if(paramConstrainModeP1Map[i] == MODE_PARAM_FREE)
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

            //if(paramConstrainModeP2Map[i] == MODE_PARAM_SOFT)
            if(mode == MODE_PARAM_SOFT)
            {
                // do nothing, soft constraint will be applied below
            }
            //else if(paramConstrainModeP2Map[i] == MODE_PARAM_HARD)
            else if(mode == MODE_PARAM_HARD)
            {
                // parameter fixed by minuit, continue to next param
                //continue;
                // although param is fixed by minuit, we still want to add
                // the penalty term, if available (TODO: will it always
                // be available?)
                // can check if error == 0.0 first
                
                // NOTE: changed to ignore HARD
                continue;
            }
            //else if(paramConstrainModeP2Map[i] == MODE_PARAM_FREE)
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

        /*
        if(thePhase == 0)
        {
            constraint = paramConstraintValueP1Map[i];
            error = paramConstraintErrorP1Map[i];
        }
        else if(thePhase == 1)
        {
            constraint = paramConstraintValueP2Map[i];
            error = paramConstraintErrorP2Map[i];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        */
        
        // TODO: can optimize this code

        if(error < 0.0)
        {
            std::cout << "ERROR: Invalid error value: error=" << error << std::endl;
        }

        // check if hard parameter and error == 0.0
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
        double param_value = p[j];
        // this parameter is from minuit internal and thus is in minuit
        // internal units (not Bq)
        // have to convert to Bq units
    
        // convert to Bq
        // multiply by the initial value
        // activity_value_Bq should really be called param_initial_value
        double activity_value_Bq = 0.0;
        double tmp_err;
        get_paramInitValueError(thePhase, i, activity_value_Bq, tmp_err);
        /*
        if(thePhase == 0)
        {
            activity_value_Bq = paramInitValueP1Map[i];
        }
        else if(thePhase == 1)
        {
            activity_value_Bq = paramInitValueP2Map[i];
        }
        else
        {
            std::cout << "ERROR: Invalid value for thePhase: thePhase=" << thePhase << std::endl;
        }
        */
        //param_value *= activity_value_Bq;

        double value = param_value * activity_value_Bq;
        //double penalty = std::pow((param_value - constraint) / error, 2.0);
        double penalty = std::pow((value - constraint) / error, 2.0);
        // TODO: penalty term should be a Gaussian constraint?
        // NOTE: gaussian constraint, after taking log, is the same as
        // quadratic constraint - however there is the issue of a remaining
        // constant term which I do not yet fully understand
        // TODO: is this the correct error term?
        // error on constraint rather than error on current fit value?
        // TODO: is the value correct?
        // NOTE: I think it's all correct

        //double value = param_value * activity_value_Bq;
        //double penalty = std::pow((value - constraint) / error, 2.0);

        //std::cout << "adding (but it has to be subtracting!) penalty (i=" << i << ") : " << penalty << std::endl;

        // subtracting positive number makes fval go up
        //loglik -= penalty;
        penalty_term += penalty;
    }


    std::cout << "penalty_term=" << penalty_term << std::endl;
    //std::cin.get();
  
    //fval = -2.0 * loglik; 
    // equivalent to
    //fval = -2.0 * (loglik_no_penalty_terms + penalty_terms); 
    //fval = -2.0 * loglik_no_penalty_terms - 2.0 * penalty_terms; 
    //fval = -2.0 * loglik_no_penalty_terms + 2.0 * penalty_terms_positive_sign;
    // then fix factor of 2.0 bug to get
    //fval = -2.0 * loglik_no_penalty_terms + penalty_terms_positive_sign;
    fval = -2.0 * loglik + penalty_term;
    //tmpData->Delete();

    // hook
    global_chisquare = fval;


    // set last parameter values
    // could also loop over nPar?
    for(int i = 0; i < numberParams; ++ i)
    {
        paramLastValueMap[i] = p[i];
    }

    return;

}



///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

Double_t getNumberMC1D(
    const Int_t channel,
    const Int_t bin_ix,
    const Double_t *const p
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
    const Double_t *const p)
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


#endif // NEWLOGLIKFITTER_LOGLIKELIHOOD_H
