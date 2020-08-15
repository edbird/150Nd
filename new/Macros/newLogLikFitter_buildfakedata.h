#ifndef NEWLOGLIKFITTER_BUILDFAKEDATA_H
#define NEWLOGLIKFITTER_BUILDFAKEDATA_H


///////////////////////////////////////////////////////////////////////////////
// function declarations
///////////////////////////////////////////////////////////////////////////////



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
#if 0
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
                       std::cout << __func__ << " ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
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
                std::cout << __func__ << " error could not find histogram: tmpName=" << tmpName << std::endl;
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
                       std::cout << __func__ << " ERROR: could not find " << tmp_sample_name << " in MCNameToParamNumberMap" << std::endl;
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
                std::cout << __func__ << " error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }

        allFakeDataSamples2D->Add((TH2D*)hAllMC2D[i]);


        std::cout << "The integral for fake data is " << hAllMC2D[i]->Integral() << std::endl;
    }






}
#endif


















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

    bool debugprint = true;


    ///////////////////////////////////////////////////////////////////////
    // reweight hMinMaxEnergy_
    // reweight all
    ///////////////////////////////////////////////////////////////////////


    // new code to reweight 150Nd by xi_{31} parameter

    // pointers of histograms to pass to reweight function
    TH1D *hTotalE_P1 = nullptr;
    TH1D *hSingleEnergy_P1 = nullptr;
    TH1D *hHighEnergy_P1 = nullptr;
    TH1D *hLowEnergy_P1 = nullptr;
    TH1D *hEnergySum_P1 = nullptr;
    TH1D *hEnergyDiff_P1 = nullptr;
    TH2D *hHighLowEnergy_P1 = nullptr;

    TH1D *hTotalE_P2 = nullptr;
    TH1D *hSingleEnergy_P2 = nullptr;
    TH1D *hHighEnergy_P2 = nullptr;
    TH1D *hLowEnergy_P2 = nullptr;
    TH1D *hEnergySum_P2 = nullptr;
    TH1D *hEnergyDiff_P2 = nullptr;
    TH2D *hHighLowEnergy_P2 = nullptr;

    /*
    TH1D *hSingleEnergyClone = nullptr;
    */

    // search through each channel for 150nd samples
    //for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        std::cout << "channel=" << channel << std::endl;
 
        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1 = histname + std::string(Nd150Files[0]) + "_P1";
        std::string search_object_P2 = histname + std::string(Nd150Files[0]) + "_P2";

        TObject *tmpobj_P1 = allDataSamples1D->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allDataSamples1D->FindObject(search_object_P2.c_str());
        allDataSamples1D->Remove(tmpobj_P1);
        allDataSamples1D->Remove(tmpobj_P2);
    } // channel

#if 0
    // search through each channel for 150nd samples
    //for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        std::cout << "channel=" << channel << std::endl;

        // search through each sample for this channel
        for(int i = 0; i < allMCSamples2D[channel]->GetEntries(); ++ i)
        {
            TH1D *tmpHist = (TH1D*)allMCSamples2D[channel]->At(i);
            TString tmpHist_name = tmpHist->GetName();
            std::cout << "tmpHist_name=" << tmpHist_name << std::endl;

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
#endif

    if(debugprint)
    {
        std::cout << "xi_31=" << xi_31 << " xi_31_baseline=" << xi_31_baseline << std::endl;
    }
    //const double xi_31_baseline{0.296};
    // NOTE: 2020-06-17 this was a bug, removed

    TH1D *hWeight_P1 = nullptr;
    TH1D *hWeight_P2 = nullptr; // TODO: move?
    if(debugprint || false)
    {
        std::cout << "before reweight_apply()" << std::endl;
    }
    reweight_apply(
        hTotalE_P1,
        hSingleEnergy_P1,
        hHighEnergy_P1,
        hLowEnergy_P1,
        hEnergySum_P1,
        hEnergyDiff_P1,
        hHighLowEnergy_P1,
        hWeight_P1,
        hTotalE_P2,
        hSingleEnergy_P2,
        hHighEnergy_P2,
        hLowEnergy_P2,
        hEnergySum_P2,
        hEnergyDiff_P2,
        hHighLowEnergy_P2,
        hWeight_P2,
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
    std::cout << "adding P1 histograms" << std::endl;
    allMCSamples1D[0]->Add(hTotalE_P1);
    allMCSamples1D[1]->Add(hSingleEnergy_P1);
    allMCSamples1D[2]->Add(hHighEnergy_P1);
    allMCSamples1D[3]->Add(hLowEnergy_P1);
    allMCSamples1D[4]->Add(hEnergySum_P1);
    allMCSamples1D[5]->Add(hEnergyDiff_P1);

    std::cout << "adding P2 histograms" << std::endl;
    allMCSamples1D[0]->Add(hTotalE_P2);
    allMCSamples1D[1]->Add(hSingleEnergy_P2);
    allMCSamples1D[2]->Add(hHighEnergy_P2);
    allMCSamples1D[3]->Add(hLowEnergy_P2);
    allMCSamples1D[4]->Add(hEnergySum_P2);
    allMCSamples1D[5]->Add(hEnergyDiff_P2);


    allMCSamples2D[0]->Add(hHighLowEnergy_P1);
    allMCSamples2D[0]->Add(hHighLowEnergy_P2);
}










#endif // NEWLOGLIKFITTER_BUILDFAKEDATA_H
