#ifndef NEWLOGLIKFITTER_DRAW_ALL_H
#define NEWLOGLIKFITTER_DRAW_ALL_H



void draw_all(const double* const AdjustActs, const double* const AdjustActs_Err, const double* const CovMatrix, const int number_free_params, const double fval)
{
    // TODO: load AdjustActs, AdjustActs_Err from file

    ///////////////////////////////////////////////////////////////////////////
    // read from input file
    ///////////////////////////////////////////////////////////////////////////

/*
    if(allDataSamples1D->GetEntries() == 0)
    {

        TFile *fin = new TFile("Nd150_2e_P" + Phase + "_fit_histograms.root", "UPDATE");
        std::cout << "reading histograms from \"\"" << std::endl;

        std::ifstream ifinaux("Nd150_2e_P" + Phase + "_fit_histograms.txt", std::ifstream::in);

        std::string auxname;

        for(int channel = 0; channel < number1DHists; ++ channel)
        {
            std::cout << "reading: 1D: channel=" << channel << std::endl;
            
            TString channel_str;
            channel_str.Form("%i", channel);
            //fout->cd("/1D");
            //fout->cd("/");
            //TDirectory *dir = gDirectory;
            //TDirectory *dir_histogram = dir->mkdir("channel_" + channel_str);
            //TDirectory *dir_histogram = dir->mkdir("channel_1D_" + channel_str);

            ifinaux >> auxname;

            TString hfullname = TString(auxname);
            //TString hfullname = "channel_1D_" + channel_str + "/" + hname;

            allDataSamples1D->Add((TH1F*)fin->Get(hfullname));
            allMCSamples1D[channel]->Add((TH1F*)fin->Get(hfullname));

        }

        for(int channel = 0; channel < number2DHists; ++ channel)
        {
            std::cout << "reading: 2D: channel=" << channel << std::endl;

            TString channel_str;
            channel_str.Form("%i", channel);
            //fout->cd("/2D");
            //fout->cd("/");
            //TDirectory *dir = gDirectory;
            //TDirectory *dir_histogram = dir->mkdir("channel_" + channel_str);
            //TDirectory *dir_histogram = dir->mkdir("channel_2D_" + channel_str);
            
            ifinaux >> auxname;

            TString hfullname = TString(auxname);
            //TString hfullname = "channel_2D_" + channel_str + "/" + hname;

            allDataSamples2D->Add((TH2F*)fin->Get(hfullname));
            allMCSamples2D[channel]->Add((TH2F*)fin->Get(hfullname));

        }

        ifinaux.close();

    }
*/

    TH1F *hHighEnergy_allMC = nullptr;
    TH1F *hLowEnergy_allMC = nullptr;
    TH1F *hHighEnergy_data = nullptr;
    TH1F *hLowEnergy_data = nullptr;

    draw(AdjustActs, AdjustActs_Err, fval,
         hHighEnergy_allMC, hLowEnergy_allMC,
         hHighEnergy_data, hLowEnergy_data,
         "hTotalE.*");
    
    draw_2D(AdjustActs, AdjustActs_Err, "hHighLowEnergy.*",
            hHighEnergy_allMC, hLowEnergy_allMC,
            hHighEnergy_data, hLowEnergy_data);

    draw_covariance_matrix(CovMatrix, number_free_params, "cov_matrix.*");


    //draw_outputdiff(AdjustActs, 0.296, "c_outputdiff_.png", -1);
    //draw_outputdiff(AdjustActs, 0.0, "c_outputdiff_.png", -1);
//    draw_outputdiff(AdjustActs, 0.0, "c_outputdiff_.png", 4);
// TODO: this is replacing existing histograms with same name?

}

#endif //NEWLOGLIKFITTER_DRAW_ALL_H
