#ifndef NEWLOGLIKFITTER_LOGLIKELIHOOD_H
#define NEWLOGLIKFITTER_LOGLIKELIHOOD_H


///////////////////////////////////////////////////////////////////////////////
// function declarations
///////////////////////////////////////////////////////////////////////////////


void logLikelihood(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */); 
Double_t getNumberMC1D(const Int_t channel, const Int_t bin_ix, const Double_t *const p);
Double_t getNumberMC2D(const Int_t channel, const Int_t bin_ix, const Int_t bin_iy, const Double_t *const p);



/*
double logpoisson(const double nData, const double nMC)
{
    double x = nData;
    double par = nMC;
    if(par < 0.0)
    {
        par = 0.0;
    }
    double ret = 0.0;
    if(x > 0.0)
    {
        ret = x * std::log(par) - par - TMath::LnGamma(x + 1.0);
    }
    else
    {
        ret = -par;
    }
    return ret;
}

*/


//double logpoisson_sterling(const double nData, const double nMC)
double logpoisson(const double nData, const double nMC)
{
    //std::cout << __func__ << "(nData=" << nData << ", nMC=" << nMC << ")" << std::endl;

    if(nData < 0.0) throw "error";

    double mnu = nMC;
    double dnu = nData;
    //std::cout << __func__ << "(dnu=" << dnu << ", mnu=" << mnu << ")" << std::endl;

    if(mnu < 0.0)
    {
        std::cout << "settng mnu=0" << std::endl;
        mnu = 0.0;
    }
    double ret = 0.0;
    if(dnu > 0.0)
    {   
        //std::cout << "dnu > 0.0" << std::endl;
        ret = -1.0 * (mnu - dnu + dnu * std::log(dnu / mnu));
        //std::cout << "ret=" << ret << std::endl;
    }
    else
    {
        //std::cout << "dnu ! > 0.0" << std::endl;
        //ret = -1.0 * (mnu - dnu);
        ret = -mnu;
        //std::cout << "ret=" << ret << std::endl;
    }
    //std::cout << "ll: nData=" << dnu << " nMC=" << mnu << " ll=" << ret << std::endl;
    return ret;
}


/*
double logpoisson(const double nData, const double nMC)
{
    double mnu = nMC;
    double dnu = nData;

    if(mnu < 0.0)
    {
        mnu = 0.0;
    }
    double ret = 0.0;
    if(dnu > 0.0)
    {   
        ret = TMath::Log(TMath::Poisson(nData, nMC));
    }
    else
    {
        //ret = -1.0 * (mnu - dnu);
        ret = TMath::Log(TMath::Poisson(nData, nMC));
    }
    //std::cout << "ll: nData=" << dnu << " nMC=" << mnu << " ll=" << ret << std::endl;
    return ret;
}
*/


// changing resolution has moved fit point
// does changing x maximum move the contours? (it used to)
// does not appear to do so by changing x upper limit
// or x lower limit
// or y axis upper limit
// disable other backgrounds and check result


void rebuild_150Nd_MC(const double, const double);


// NOTE:
// this function does not work if reweight is called before this function is
// called
// for obvious reasons
void build_fake_data()
{

    std::cout << "build_fake_data()" << std::endl;


    std::cout << "rebuild: xi_31_baseline=" << xi_31_baseline << std::endl;
    //rebuild_150Nd_MC(/*xi_31_baseline*/ 0.296, xi_31_baseline);
    rebuild_150Nd_MC(xi_31_baseline, xi_31_baseline);


    TH1D *hAllMC1D[number1DHists];
    TH2D *hAllMC2D[number2DHists]; // TODO

    std::cout << "debug: number of data samples: " << allDataSamples1D->GetEntries() << std::endl;
    std::cout << "debug: number of MC samples: " << allMCSamples1D[0]->GetEntries() << std::endl;
        
    allFakeDataSamples1D = new TObjArray();
    allFakeDataSamples2D = new TObjArray();


    // each channel 1D hists
    // this is for(i = 0; i < 1; ++ i)
    // TODO: this isn't right. should this be iterating over the "channel" ?
    for(int i = 0; i < allDataSamples1D->GetEntries(); i++)
    {
//        std::cout << "i=" << i << std::endl;

        // because this isn't right TODO
        // uses at(i), but i should always be zero and there should be an
        // additional array index
        //data1D[i] = (TH1D*)allDataSamples1D->At(i)->Clone();
        

        TH1D *tmpHist;

        // j list MC samples for this channel i
        //std::cout << "debug: number of MC samples (i=" << i << "): " << allMCSamples1D[i]->GetEntries() << std::endl;


        // allMCSamples1D[0] contains objects such as: "zr96_rot_k40_2e_P2"

        // TODO i should be channel here?
        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {

            //std::cout << "j=" << j << std::endl;

            TString j_str;
            j_str.Form("%i", j);

            tmpHist = (TH1D*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist->GetName();

            //std::cout << "(1) tmpName=" << tmpName << std::endl;

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param = -1;
            bool found_param = false;

            // search through parameters to find right one
            // the histogram names are formatted like:
            // hTotalE_bi214_mylar_fit
            // histogram_name + "_" + mc_sample_name + "_fit"
            
            // used later
            double activity_scale_branching_ratio = 1.0;

            {
                std::string tmp_hist_name(tmpName);
                auto i_start = tmp_hist_name.find('_') + 1;
                auto i_end = tmp_hist_name.rfind('_');
                if(i_end - i_start > 0)
                {
                    std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);

                    // TODO
                    // do not have to scale by 0.36 after reading from file
                    // as scaling is done when reading
                    /*
                    // set branching ratio fraction
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
                        // TODO: removed std::string, change tmpName type to be std::string from TString
                    
                        //which_param = paramNumber;
                        which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                        found_param = true;
//                        std::cout << "j=" << j << ": paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << " which_param=" << which_param << std::endl;
                    }
                    else
                    {
                       std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                    }
                }
            }

//            std::cin.get();

            if(found_param == true)
            {
                //std::cout << "found histogram: tmpName=" << tmpName << " which_param=" << which_param << std::endl;

                // scale histogram to correct size using output parameter
                // from fit
                if(which_param >= numberEnabledParams)
                {
                    std::cout << "throwing exception, which_param=" << which_param << std::endl;
                    throw std::runtime_error("which_param invalid value");
                }



                // no error thrown, which_param is presumably the correct index
                //Double_t activity_scale = AdjustActs[which_param] * activity_scale_branching_ratio;
                Double_t activity_scale = paramInitValueMap[which_param]; // * activity_scale_branching_ratio;
                if(which_param == 0)
                {
                    std::cout << "in build_fake_data(): activity_scale=" << activity_scale << std::endl;

                    // SSD
                    //activity_scale *= 1.2;
                    //tmpHist->Scale(1.5);
                }
                else if(which_param == 1)
                {
                    std::cout << "ERROR: which_param == 1 !" << std::endl;
                    std::cin.get();
                }
                //std::cout << "activity_scale=" << activity_scale << std::endl;
//                tmpHist->Scale(activity_scale);
// TODO: has already been scaled by this activity when read in

                if(tmpHist->Integral() > 0)
                {

                    TString hname = tmpHist->GetName();

                    if(j == 0)
                    {
                        //std::cout << "Clone() done" << "j=" << j << std::endl;
                        // TODO: bug here if Integral() for j == 0 is zero
                        
                        hAllMC1D[i] = (TH1D*)tmpHist->Clone("Total MC Fake Data");
                        //hAllMC1D[i] = (TH1D*)tmpHist_drawpointer->Clone("Total MC");
                        
                        /*
                        std::cout << "j=" << j << std::endl;
                        for(int k = 0; k < tmpHist->GetNbinsX(); ++ k)
                        {
                            std::cout << "k=" << k << " " << tmpHist->GetBinContent(k) << std::endl;
                        }
                        */
                    }
                    else
                    {
                        hAllMC1D[i]->Add((TH1D*)tmpHist);
                        //hAllMC1D[i]->Add((TH1D*)tmpHist_drawpointer);
                        
                        /*
                        std::cout << "j=" << j << std::endl;
                        for(int k = 0; k < tmpHist->GetNbinsX(); ++ k)
                        {
                            std::cout << "k=" << k << " " << tmpHist->GetBinContent(k) << std::endl;
                        }
                        */
                    }
	            }
                else
                {
                    //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist->GetName() << std::endl;
                }
            }
            else
            {
                std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }

        /*
        std::cout << "integral for fakedata sample " << i << " is " << hAllMC1D[i]->Integral() << std::endl;
        for(int k = 0; k < hAllMC1D[i]->GetNbinsX(); ++ k)
        {
            std::cout << "k=" << k << " " << hAllMC1D[i]->GetBinContent(k) << std::endl;
        }
        */
        allFakeDataSamples1D->Add((TH1D*)hAllMC1D[i]);


        std::cout << "The integral for fake data is " << hAllMC1D[i]->Integral() << std::endl;
    }







    // each channel 2D hists
    for(int i = 0; i < allDataSamples2D->GetEntries(); i++)
    {

        TH2D *tmpHist;

        // allMCSamples2D[0] contains objects such as: "zr96_rot_k40_2e_P2"

        for(int channel = 0; channel < allMCSamples2D[i]->GetEntries(); ++ channel)
        {

            TString channel_str;
            channel_str.Form("%i", channel);

            tmpHist = (TH2D*)allMCSamples2D[i]->At(channel)->Clone();
            TString tmpName = tmpHist->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param = -1;
            bool found_param = false;

            // search through parameters to find right one
            // the histogram names are formatted like:
            // hTotalE_bi214_mylar_fit
            // histogram_name + "_" + mc_sample_name + "_fit"
            
            // used later
            double activity_scale_branching_ratio = 1.0;

            {
                std::string tmp_hist_name(tmpName);
                auto i_start = tmp_hist_name.find('_') + 1;
                auto i_end = tmp_hist_name.rfind('_');
                if(i_end - i_start > 0)
                {
                    std::string tmp_sample_name = tmp_hist_name.substr(i_start, i_end - i_start);

                    // TODO
                    // do not have to scale by 0.36 after reading from file
                    // as scaling is done when reading
                    /*
                    // set branching ratio fraction
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
                        // TODO: removed std::string, change tmpName type to be std::string from TString
                    
                        //which_param = paramNumber;
                        which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                        found_param = true;
//                        std::cout << "j=" << j << ": paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << " which_param=" << which_param << std::endl;
                    }
                    else
                    {
                       std::cout << "ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
                    }
                }
            }

            if(found_param == true)
            {
                //std::cout << "found histogram: tmpName=" << tmpName << " which_param=" << which_param << std::endl;

                // scale histogram to correct size using output parameter
                // from fit
                if(which_param >= numberEnabledParams)
                {
                    std::cout << "throwing exception, which_param=" << which_param << std::endl;
                    throw std::runtime_error("which_param invalid value");
                }

                // no error thrown, which_param is presumably the correct index
                //Double_t activity_scale = AdjustActs[which_param] * activity_scale_branching_ratio;
                Double_t activity_scale = paramInitValueMap[which_param]; // * activity_scale_branching_ratio;
                if(which_param == 0)
                {
                    std::cout << "in build_fake_data(): activity_scale=" << activity_scale << std::endl;

                    // SSD
                    //activity_scale *= 1.2;
                    //tmpHist->Scale(1.5);
                }
                else if(which_param == 1)
                {
                    std::cout << "ERROR: which_param == 1 !" << std::endl;
                    std::cin.get();
                }
                //std::cout << "activity_scale=" << activity_scale << std::endl;
//                tmpHist->Scale(activity_scale);
// TODO: has already been scaled by this activity when read in

                if(tmpHist->Integral() > 0)
                {

                    TString hname = tmpHist->GetName();

                    if(channel == 0)
                    {
                        //std::cout << "Clone() done" << "j=" << j << std::endl;
                        // TODO: bug here if Integral() for j == 0 is zero
                        
                        hAllMC2D[i] = (TH2D*)tmpHist->Clone("Total MC Fake Data");
                    }
                    else
                    {
                        hAllMC2D[i]->Add((TH2D*)tmpHist);
                    }
	            }
                else
                {
                    //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist->GetName() << std::endl;
                }
            }
            else
            {
                std::cout << "error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }

        allFakeDataSamples2D->Add((TH2D*)hAllMC2D[i]);


        std::cout << "The integral for fake data is " << hAllMC2D[i]->Integral() << std::endl;
    }






}


















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





                //
                // if n is large, then exp(-n) may fail in Poisson calc
                // n! may also fail, so use Stirling
                //
                // so calculate LL from expansion of log
                //
                //if(nMC > 10.0)
                //{
                //    // nMC > 1.0e-05 here already
                //    const double l = nMC;
                //    const double n = nData;
                //    const double stirling = n*TMath::Log(n) - n;
                //    ll_channel += -l + n*TMath::Log(l) - stirling;
                //}




void rebuild_150Nd_MC(const double xi_31, const double xi_31_baseline)
{

    bool debugprint = false;


    ///////////////////////////////////////////////////////////////////////
    // reweight hMinMaxEnergy_
    // reweight all
    ///////////////////////////////////////////////////////////////////////


    // new code to reweight 150Nd by xi_{31} parameter

    // pointers of histograms to pass to reweight function
    TH1D *hTotalE = nullptr;
    TH1D *hSingleEnergy = nullptr;
    TH1D *hHighEnergy = nullptr;
    TH1D *hLowEnergy = nullptr;
    TH1D *hEnergySum = nullptr;
    TH1D *hEnergyDiff = nullptr;
    TH2D *hHighLowEnergy = nullptr;

    /*
    TH1D *hSingleEnergyClone = nullptr;
    */

    // search through each channel for 150nd samples
    //for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        // search through each sample for this channel
        for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
        {
            TH1D *tmpHist = (TH1D*)allMCSamples1D[channel]->At(i);
            TString tmpHist_name = tmpHist->GetName();
            // TODO: had to add "_fit" here - might not work after 1 iteration
            //if(tmpHist_name.CompareTo("hTotalE_nd150_rot_2n2b_m4_fit") == 0 ||
            //   tmpHist_name.CompareTo("hTotalE_nd150_rot_2b2n_m4_fit") == 0)

            // if name contains this string, it needs to be reweighted
            if(tmpHist_name.Contains("nd150_rot_2n2b_m4") ||
               tmpHist_name.Contains("nd150_rot_2b2n_m4"))
            {
                /*
                if(channel == 1)
                {
                    hSingleEnergyClone = (TH1D*)allMCSamples1D[channel]->At(i)->Clone();
                }
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
            TH1D *tmpHist = (TH1D*)allMCSamples2D[channel]->At(i);
            TString tmpHist_name = tmpHist->GetName();

            // if name contains this string, it needs to be reweighted
            if(tmpHist_name.Contains("nd150_rot_2n2b_m4") ||
               tmpHist_name.Contains("nd150_rot_2b2n_m4"))
            {
                // remove histogram
                allMCSamples2D[channel]->RemoveAt(i);
                break;
                // TODO: will only work for a maximum of 1 histograms
                // removed because index i will shift due to removing
                // histogram at i
            }

        } // i
    } // channel


    if(debugprint)
    {
        std::cout << "xi_31=" << xi_31 << " xi_31_baseline=" << xi_31_baseline << std::endl;
    }
    //const double xi_31_baseline{0.296};
    // NOTE: 2020-06-17 this was a bug, removed

    TH1D *hWeight = nullptr;
    if(debugprint || false)
    {
        std::cout << "before reweight_apply()" << std::endl;
    }
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
    if(debugprint || false)
    {
        std::cout << "after reweight_apply()" << std::endl;
    }

    /*
    hSingleEnergyClone->Divide(hSingleEnergy);
    for(Int_t ii = 1; ii <= hSingleEnergyClone->GetNbinsX(); ++ ii)
    {
        float content = hSingleEnergyClone->GetBinContent(ii);
        std::cout << "content: " << ii << " " << content << std::endl;
    }
    std::cin.get();
    */



    // TODO: just another example of manual code edits
    // make a file describing the channels to fit as well as the parameters
    allMCSamples1D[0]->Add(hTotalE);
    allMCSamples1D[1]->Add(hSingleEnergy);
    allMCSamples1D[2]->Add(hHighEnergy);
    allMCSamples1D[3]->Add(hLowEnergy);
    allMCSamples1D[4]->Add(hEnergySum);
    allMCSamples1D[5]->Add(hEnergyDiff);
    allMCSamples2D[0]->Add(hHighLowEnergy);
}



















///////////////////////////////////////////////////////////////////////////////
// logLikelihood
///////////////////////////////////////////////////////////////////////////////

static int counter = 0;

// TODO don't appear to work with parameters with more than one MC
void logLikelihood(Int_t & nPar, Double_t* /*grad*/, Double_t &fval, Double_t *p, Int_t /*iflag */)
{


    /*
    if(p[1] < -0.4)
    {
        std::cout << __func__ << " p[1]=" << p[1] << std::endl;
        
        for(Int_t pix = 0; pix < nPar; ++ pix)
        {
            std::cout << p[pix] << std::endl;
        }

        std::cin.get();

    }
    */


    if(allFakeDataSamples1D == nullptr)
    {
        build_fake_data();
    }


    // save the walk
    ll_walk.push_back(std::make_pair(p[1], p[0]));
    //std::cout << __func__ << " ll_walk: " << "p[1]=" << p[1] << " p[0]=" << p[0] << std::endl;


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
        std::cout << "p[0]=" << p[0] << " p[" << axial_vector_parameter_0_param_number << "]="
                  << p[axial_vector_parameter_0_param_number] << std::endl;
    }


    // TODO: will not work if parameter number changes
    // 2020-06-17
    //if(p[1] != last_xi_31_parameter_value)
    //if(p[1] != paramLastValueMap[1])
    if(p[axial_vector_parameter_0_param_number] != paramLastValueMap[axial_vector_parameter_0_param_number])
    {

        // TODO: rebuild nd150 xi_31 paramter histogram here
        //std::cout << "rebuilding 150 Nd MC" << std::endl;

        const double xi_31{p[axial_vector_parameter_0_param_number]};
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

    TH1D *tmpData1D;
    TH1D *tmpFakeData1D;
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

        tmpData1D = (TH1D*)allDataSamples1D->At(channel);
        tmpFakeData1D = (TH1D*)allFakeDataSamples1D->At(channel);

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

            double nMC = getNumberMC1D(channel, bin_ix, p);

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

    TH2D *tmpData2D;
    TH2D *tmpFakeData2D;
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
        //tmpData1D = (TH1D*)allDataSamples1D->At(i)->Clone("tmpData1D" + i_str + "_");
        tmpData2D = (TH2D*)allDataSamples2D->At(channel);
        tmpFakeData2D = (TH2D*)allFakeDataSamples2D->At(channel);

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

                
                double nMC = getNumberMC2D(channel, bin_ix, bin_iy, p);

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
    for(int i = 0; i < nPar; ++ i)
    {
        double value = p[i];
        paramLastValueMap[i] = value;
        minuitParamCurrent[i] = value;
        minuitParamLast[i] = value;
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
    TH1D *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;

    //std::cout << "channel=" << channel << std::endl;
    //std::cout << allMCSamples1D[channel]->GetEntries() << std::endl;


    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples1D[channel]->GetEntries(); j++)
    {

        // disabled params do not appear in this array

        tmpHist = (TH1D*)allMCSamples1D[channel]->At(j);

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
    const Double_t *const p)
{

    double nMC = 0.;

    // (1) grab a hist from the sample list of this channel
    // (2) figure out which parameter it corresponds to
    TH1D *tmpHist;
    //int which_param;

    //std::cout << "getting number of MC... "  << channel << std::endl;
    //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

    // copied from above
    //for(int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++)
    for(int j = 0; j < allMCSamples2D[channel]->GetEntries(); j++)
    {

        // get histogram name
        tmpHist = (TH1D*)allMCSamples2D[channel]->At(j);
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


#endif // NEWLOGLIKFITTER_LOGLIKELIHOOD_H
