#ifndef NEWLOGLIKFITTER_REBUILD_150ND_MC_H
#define NEWLOGLIKFITTER_REBUILD_150ND_MC_H



void rebuild_150Nd_MC(const double xi_31, const double xi_31_baseline)
{

    int debuglevel = 1;

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
        if(debuglevel >= 3)
        {
            std::cout << "channel=" << channel << std::endl;
        }
 
        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1 = histname + std::string(Nd150Files[0]) + "_P1_fit";
        std::string search_object_P2 = histname + std::string(Nd150Files[0]) + "_P2_fit";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allMCSamples1D[channel]->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allMCSamples1D[channel]->FindObject(search_object_P2.c_str());
        if(tmpobj_P1 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P1 << std::endl;
            /*for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples1D[channel]->At(i)->GetName() << std::endl;
            }*/
            throw "problem";
        }
        if(tmpobj_P2 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P2 << std::endl;
            throw "problem";
        }

        if(debuglevel >= 5)
        {
            std::cout << "before ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples1D[channel]->At(i)->GetName() << std::endl;
            }
        }
        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allMCSamples1D[channel]->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allMCSamples1D[channel]->Compress();
        if(allMCSamples1D[channel]->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allMCSamples1D[channel]->Compress();
        if(debuglevel >= 5)
        {
            std::cout << "after ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples1D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples1D[channel]->At(i)->GetName() << std::endl;
            }
        }
    } // channel

    // search through each channel for 150nd samples
    //for(int channel = 0; channel < allDataSamples2D->GetEntries(); ++ channel)
    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        if(debuglevel >= 3)
        {
            std::cout << "channel=" << channel << std::endl;
        }

        std::string histname = std::string(channel_histname_2D[channel]);
        std::string search_object_P1 = histname + std::string(Nd150Files[0]) + "_P1_fit";
        std::string search_object_P2 = histname + std::string(Nd150Files[0]) + "_P2_fit";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allMCSamples2D[channel]->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allMCSamples2D[channel]->FindObject(search_object_P2.c_str());
        if(tmpobj_P1 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P1 << std::endl;
            throw "problem";
        }
        if(tmpobj_P2 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P2 << std::endl;
            throw "problem";
        }

        if(debuglevel >= 5)
        {
            std::cout << "before ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples2D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples2D[channel]->At(i)->GetName() << std::endl;
            }
        }
        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allMCSamples2D[channel]->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allMCSamples2D[channel]->Compress();
        if(allMCSamples2D[channel]->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allMCSamples2D[channel]->Compress();
        if(debuglevel >= 5)
        {
            std::cout << "after ->Remove():" << std::endl;
            for(int i = 0; i < allMCSamples2D[channel]->GetEntries(); ++ i)
            {
                std::cout << allMCSamples2D[channel]->At(i)->GetName() << std::endl;
            }
        }
    } // channel

    if(debuglevel >= 2)
    {
        std::cout << "xi_31=" << xi_31 << " xi_31_baseline=" << xi_31_baseline << std::endl;
    }
    //const double xi_31_baseline{0.296};
    // NOTE: 2020-06-17 this was a bug, removed

    TH1D *hWeight_P1 = nullptr;
    TH1D *hWeight_P2 = nullptr; // TODO: move?
    if(debuglevel >= 5)
    {
        std::cout << "before reweight_apply()" << std::endl;
    }
    reweight_apply_MC(
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
    if(debuglevel >= 5)
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
    if(debuglevel >= 4)
    {
        std::cout << "adding P1 histograms" << std::endl;
    }
    allMCSamples1D[0]->Add(hTotalE_P1);
    allMCSamples1D[1]->Add(hSingleEnergy_P1);
    allMCSamples1D[2]->Add(hHighEnergy_P1);
    allMCSamples1D[3]->Add(hLowEnergy_P1);
    allMCSamples1D[4]->Add(hEnergySum_P1);
    allMCSamples1D[5]->Add(hEnergyDiff_P1);

    allMCSamples2D[0]->Add(hHighLowEnergy_P1);
    
    if(debuglevel >= 4)
    {
        std::cout << "adding P2 histograms" << std::endl;
    }
    allMCSamples1D[0]->Add(hTotalE_P2);
    allMCSamples1D[1]->Add(hSingleEnergy_P2);
    allMCSamples1D[2]->Add(hHighEnergy_P2);
    allMCSamples1D[3]->Add(hLowEnergy_P2);
    allMCSamples1D[4]->Add(hEnergySum_P2);
    allMCSamples1D[5]->Add(hEnergyDiff_P2);


    allMCSamples2D[0]->Add(hHighLowEnergy_P2);
}


#endif // NEWLOGLIKFITTER_REBUILD_150ND_MC_H

