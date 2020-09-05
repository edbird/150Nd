#ifndef NEWLOGLIKFITTER_REBUILD_150ND_DATA_H
#define NEWLOGLIKFITTER_REBUILD_150ND_DATA_H



void rebuild_150Nd_data()
{

    int debuglevel = 1;

    ///////////////////////////////////////////////////////////////////////
    // apply systematics to data sample
    ///////////////////////////////////////////////////////////////////////


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
        std::string search_object_P1 = histname + std::string(DataFile) + "_P1";
        std::string search_object_P2 = histname + std::string(DataFile) + "_P2";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allDataSamples1D->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allDataSamples1D->FindObject(search_object_P2.c_str());
        if(tmpobj_P1 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P1 << std::endl;
            /*for(int i = 0; i < allDataSamples1D->GetEntries(); ++ i)
            {
                std::cout << allDataSamples1D->At(i)->GetName() << std::endl;
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
            for(int i = 0; i < allDataSamples1D->GetEntries(); ++ i)
            {
                std::cout << "allDataSamples1D->At(" << i << ")=" << allDataSamples1D->At(i) << std::endl;
                //std::cout << allDataSamples1D->At(i)->GetName() << std::endl;
            }
        }
        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allDataSamples1D->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allDataSamples1D->Compress();
        if(allDataSamples1D->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allDataSamples1D->Compress();
        if(debuglevel >= 5)
        {
            std::cout << "after ->Remove():" << std::endl;
            std::cout << "there are " << allDataSamples1D->GetEntries() << " samples to search" << std::endl;
            for(int i = 0; i < allDataSamples1D->GetEntries(); ++ i)
            {
                std::cout << "allDataSamples1D->At(" << i << ")=" << allDataSamples1D->At(i) << std::endl;
                //std::cout << allDataSamples1D->At(i)->GetName() << std::endl;
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
        std::string search_object_P1 = histname + std::string(DataFile) + "_P1";
        std::string search_object_P2 = histname + std::string(DataFile) + "_P2";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allDataSamples2D->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allDataSamples2D->FindObject(search_object_P2.c_str());
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
            for(int i = 0; i < allDataSamples2D->GetEntries(); ++ i)
            {
                std::cout << allDataSamples2D->At(i)->GetName() << std::endl;
            }
        }
        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allDataSamples2D->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allDataSamples2D->Compress();
        if(allDataSamples2D->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            throw "failed to remove";
        }
        allDataSamples2D->Compress();
        if(debuglevel >= 5)
        {
            std::cout << "after ->Remove():" << std::endl;
            for(int i = 0; i < allDataSamples2D->GetEntries(); ++ i)
            {
                std::cout << allDataSamples2D->At(i)->GetName() << std::endl;
            }
        }
    } // channel


    TH1D *hWeight_P1 = nullptr;
    TH1D *hWeight_P2 = nullptr; // TODO: move?
    if(debuglevel >= 5)
    {
        std::cout << "before reweight_apply()" << std::endl;
    }

    reweight_apply_data(
        hTotalE_P1,
        hSingleEnergy_P1,
        hHighEnergy_P1,
        hLowEnergy_P1,
        hEnergySum_P1,
        hEnergyDiff_P1,
        hHighLowEnergy_P1,
        hTotalE_P2,
        hSingleEnergy_P2,
        hHighEnergy_P2,
        hLowEnergy_P2,
        hEnergySum_P2,
        hEnergyDiff_P2,
        hHighLowEnergy_P2);

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
    allDataSamples1D->Add(hTotalE_P1);
    allDataSamples1D->Add(hSingleEnergy_P1);
    allDataSamples1D->Add(hHighEnergy_P1);
    allDataSamples1D->Add(hLowEnergy_P1);
    allDataSamples1D->Add(hEnergySum_P1);
    allDataSamples1D->Add(hEnergyDiff_P1);

    allDataSamples2D->Add(hHighLowEnergy_P1);
    
    if(debuglevel >= 4)
    {
        std::cout << "adding P2 histograms" << std::endl;
    }
    allDataSamples1D->Add(hTotalE_P2);
    allDataSamples1D->Add(hSingleEnergy_P2);
    allDataSamples1D->Add(hHighEnergy_P2);
    allDataSamples1D->Add(hLowEnergy_P2);
    allDataSamples1D->Add(hEnergySum_P2);
    allDataSamples1D->Add(hEnergyDiff_P2);


    allDataSamples2D->Add(hHighLowEnergy_P2);


}


#endif // NEWLOGLIKFITTER_REBUILD_150ND_DATA_H
