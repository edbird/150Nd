#ifndef NEWLOGLIKFITTER_REBUILD_150ND_FAKEDATA
#define NEWLOGLIKFITTER_REBUILD_150ND_FAKEDATA



///////////////////////////////////////////////////////////////////////////////
// build_fake_data() function
// with systematics
///////////////////////////////////////////////////////////////////////////////


void rebuild_fake_data_systematics(const double xi_31, const double xi_31_baseline)
{

    std::cout << "calling " << __func__ << std::endl;
    std::cout << "--- systematics ---" << std::endl;
    std::cout << gSystematics.systematic_energy_offset << std::endl;
    std::cout << "xi_31=" << xi_31 << " xi_31_baseline=" << xi_31_baseline << std::endl;
    //std::cin.get();


    int debuglevel = 9;

    ///////////////////////////////////////////////////////////////////////////
    // remove existing fake data from memory
    ///////////////////////////////////////////////////////////////////////////

    // pointers to new histograms
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

    // search through each channel for fake data samples (1D)
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(debuglevel >= 3)
        {
            std::cout << "channel=" << channel << std::endl;
            if(debuglevel >= 10)
            {
                std::cin.get();
            }
        }
 
        std::string histname = std::string(channel_histname_1D[channel]);
        //std::string search_object_P1 = histname + std::string("fake_") + std::string(DataFile) + "_P1";
        //std::string search_object_P2 = histname + std::string("fake_") + std::string(DataFile) + "_P2";
        std::string search_object_P1 = histname + std::string("fakedata") + "_P1";
        std::string search_object_P2 = histname + std::string("fakedata") + "_P2";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allFakeDataSamples1D->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allFakeDataSamples1D->FindObject(search_object_P2.c_str());
        if(tmpobj_P1 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P1 << std::endl;
            //throw "problem";
            // TODO: disabled, as this will always trigger on first iteration
        }
        if(tmpobj_P2 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P2 << std::endl;
            //throw "problem";
            // TODO: disabled, as this will always trigger on first iteration
        }

        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allFakeDataSamples1D->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            //throw "failed to remove";
            // TODO: disabled, as this will always trigger on first iteration
        }
        allFakeDataSamples1D->Compress();

        if(allFakeDataSamples1D->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            //throw "failed to remove";
            // TODO: disabled, as this will always trigger on first iteration
        }
        allFakeDataSamples1D->Compress();
        
    } // channel
    

    // search through each channel for fake data samples (2D)
    for(int channel = 0; channel < number2DHists; ++ channel)
    {
        if(debuglevel >= 3)
        {
            std::cout << "channel=" << channel << std::endl;
        }

        std::string histname = std::string(channel_histname_2D[channel]);
        std::string search_object_P1 = histname + std::string("fakedata") + "_P1";
        std::string search_object_P2 = histname + std::string("fakedata") + "_P2";

        if(debuglevel >= 3)
        {
            std::cout << "search_object_P1=" << search_object_P1 << std::endl;
            std::cout << "search_object_P2=" << search_object_P2 << std::endl;
        }

        TObject *tmpobj_P1 = allFakeDataSamples2D->FindObject(search_object_P1.c_str());
        TObject *tmpobj_P2 = allFakeDataSamples2D->FindObject(search_object_P2.c_str());
        if(tmpobj_P1 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P1 << std::endl;
            //throw "problem";
            // TODO: disabled, as this will always trigger on first iteration
        }
        if(tmpobj_P2 == nullptr)
        {
            std::cout << "could not find object: " << search_object_P2 << std::endl;
            //throw "problem";
            // TODO: disabled, as this will always trigger on first iteration
        }

        // NOTE: this just removes the object from the TObjArray
        // it does not delete the object from gROOT memory
        if(allFakeDataSamples2D->Remove(tmpobj_P1) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            //throw "failed to remove";
            // TODO: disabled, as this will always trigger on first iteration
        }
        allFakeDataSamples2D->Compress();
        
        if(allFakeDataSamples2D->Remove(tmpobj_P2) == nullptr)
        {
            std::cout << "FAILED TO REMOVE" << std::endl;
            //throw "failed to remove";
            // TODO: disabled, as this will always trigger on first iteration
        }
        allFakeDataSamples2D->Compress();

    } // channel
    

    ///////////////////////////////////////////////////////////////////////////
    // re-create fake data samples by reading in from files
    ///////////////////////////////////////////////////////////////////////////

    if(debuglevel >= 2)
    {
        std::cout << "recreating fake data by reading from file and applying systematics" << std::endl;
    }

    // look at 150 Nd xi_31 reweighting code to find closest example of what to do here
    // here we need to apply to all samples not just one
    reweight_apply_fakedata(
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
        hHighLowEnergy_P2,
        xi_31,
        xi_31_baseline,
        h_nEqNull,
        h_nEqTwo,
        psiN0,
        psiN2,
        bb_Q);


    if(debuglevel >= 4)
    {
        std::cout << "adding P1 histograms" << std::endl;
    }
    allFakeDataSamples1D->Add(hTotalE_P1);
    allFakeDataSamples1D->Add(hSingleEnergy_P1);
    allFakeDataSamples1D->Add(hHighEnergy_P1);
    allFakeDataSamples1D->Add(hLowEnergy_P1);
    allFakeDataSamples1D->Add(hEnergySum_P1);
    allFakeDataSamples1D->Add(hEnergyDiff_P1);

    allFakeDataSamples2D->Add(hHighLowEnergy_P1);
    
    if(debuglevel >= 4)
    {
        std::cout << "adding P2 histograms" << std::endl;
    }
    allFakeDataSamples1D->Add(hTotalE_P2);
    allFakeDataSamples1D->Add(hSingleEnergy_P2);
    allFakeDataSamples1D->Add(hHighEnergy_P2);
    allFakeDataSamples1D->Add(hLowEnergy_P2);
    allFakeDataSamples1D->Add(hEnergySum_P2);
    allFakeDataSamples1D->Add(hEnergyDiff_P2);

    allFakeDataSamples2D->Add(hHighLowEnergy_P2);

}



#endif // NEWLOGLIKFITTER_REBUILD_150ND_FAKEDATA
