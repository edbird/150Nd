#ifndef MINIMIZEFCNAXIALVECTOR_SET_D_H
#define MINIMIZEFCNAXIALVECTOR_SET_D_H


void
MinimizeFCNAxialVector::set_D() const
{

    // loop over all channels
    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(debuglevel >= 4)
        {
            std::cout << "channel=" << channel << std::endl;
        }

        // check channel enabled
        if(channel_enable_1D[channel] == 0)
        {
            if(debuglevel >= 5)
            {
                std::cout << "1D: channel " << channel << " disabled, skip" << std::endl;
            }
            continue;
        }

        // TODO: this no longer works because there are multiple phases
        // in the tmpData1D arrays
        // copy code from below to find correct object
        // TODO: phase 1 and phase 2 data objects
        //Double_t *tmpData1D_P1 = nullptr; //(TH1D*)allDataSamples1D->At(channel);
        //Double_t *tmpData1D_P2 = nullptr;
        //Double_t *tmpFakeData1D_P1 = nullptr; //(TH1D*)allFakeDataSamples1D->At(channel);
        //Double_t *tmpFakeData1D_P2 = nullptr;

        std::string histname = std::string(channel_histname_1D[channel]);
        std::string search_object_P1;
        std::string search_object_P2;
        if(g_mode_fake_data == false)
        {
            search_object_P1 = histname + std::string(DataFile) + "_P1";
            search_object_P2 = histname + std::string(DataFile) + "_P2";
        }
        else if(g_mode_fake_data == true)
        {
            search_object_P1 = histname + std::string("fakedata") + "_P1";
            search_object_P2 = histname + std::string("fakedata") + "_P2";
        }
        TH1D *tmpDataHist1D_P1 = nullptr;
        TH1D *tmpDataHist1D_P2 = nullptr;
        
        if(debuglevel >= 6)
        {
            std::cout << "search_object_P1=" << search_object_P1
                      << " search_object_P2=" << search_object_P2 << std::endl;
        }

        if(g_mode_fake_data == false)
        {
            tmpDataHist1D_P1 = (TH1D*)allDataSamples1D->FindObject(search_object_P1.c_str());
            tmpDataHist1D_P2 = (TH1D*)allDataSamples1D->FindObject(search_object_P2.c_str());
        }
        else if(g_mode_fake_data == true)
        {
            tmpDataHist1D_P1 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P1.c_str());
            tmpDataHist1D_P2 = (TH1D*)allFakeDataSamples1D->FindObject(search_object_P2.c_str());
        }

        if(tmpDataHist1D_P1 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }
        if(tmpDataHist1D_P2 == nullptr)
        {
            std::cout << "ERROR: Could not find object " << search_object_P1 << std::endl;
            throw "problem";
        }

        //tmpData1D_P1 = new Double_t[tmpDataHist1D_P1->GetNbinsX()]; //tmpHistP1->Clone("tmpData1D_P1");
        //tmpData1D_P2 = new Double_t[tmpDataHist1D_P2->GetNbinsX()]; //tmpHistP2->Clone("tmpData1D_P2");
        for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P1->GetNbinsX(); ++ bin_x)
        {
            //tmpData1D_P1[bin_x] = tmpDataHist1D_P1->GetBinContent(bin_x);
            //Int_t super_index = channel * 2 * 50 + bin_x;
            //Int_t super_index = bin_x;
            Double_t content = tmpDataHist1D_P1->GetBinContent(bin_x + 1);
            //D_1D_P1[channel]->SetBinContent(super_index + 1, 1, content_output);
            #if VECTOR_RANGE_CHECK
            D_1D_P1_data[channel]->at(bin_x) = content;
            #else
            D_1D_P1_data[channel]->operator[](bin_x) =  content;
            #endif
        }
        
        //std::cout << "LIST OF P2 DATA channel=" << channel << std::endl;
        //std::cin.get();
        for(Int_t bin_x{0}; bin_x < tmpDataHist1D_P2->GetNbinsX(); ++ bin_x)
        {
            //tmpData1D_P2[bin_x] = tmpDataHist1D_P2->GetBinContent(bin_x);
            //Int_t super_index = channel * 2 * 50 + 50 + bin_x;
            //Int_t super_index = bin_x;
            Double_t content = tmpDataHist1D_P2->GetBinContent(bin_x + 1);
            //D_1D_P2[channel]->SetBinContent(super_index + 1, 1, content_output);
            #if VECTOR_RANGE_CHECK
            D_1D_P2_data[channel]->at(bin_x) = content;
            #else
            D_1D_P2_data[channel]->operator[](bin_x) = content;
            #endif
            //std::cout << "bin_x=" << bin_x + 1 << " content=" << content_output << std::endl;
        }
    }
}


#endif // MINIMIZEFCNAXIALVECTOR_SET_D_H
