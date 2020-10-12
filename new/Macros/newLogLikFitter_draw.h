#ifndef NEWLOGLIKFITTER_DRAW_H
#define NEWLOGLIKFITTER_DRAW_H




//void draw(const Double_t *const AdjustActs, const Double_t *const AdjustActs_Err, const std::string& saveas_filename = "",
//          TH1D* hHighEnergy_allMC = nullptr, TH1D* hLowEnergy_allMC = nullptr,
//          TH1D* hHighEnergy_data = nullptr, TH1D* hLowEnergy_data = nullptr);

void draw_covariance_matrix(const double * const CovMatrix, const int number_free_params, const std::string& saveas_filename);

/*
void draw_2D(
    const Double_t *const AdjustActs,
    const Double_t *const AdjustActs_Err,
    const std::string &saveas_filename,
    TH1D*, TH1D*);
*/




void untitledfunction()
{

    ///////////////////////////////////////////////////////////////////////////
    // draw results - 2d
    ///////////////////////////////////////////////////////////////////////////

    #if 0

    // each channel 1D hists
    THStack *stacks2D[number2DHists];
    TH2D *data2D[number2DHists];
    for(int i = 0; i < allDataSamples2D->GetEntries(); i++)
    {
        data2D[i] = (TH2D*)allDataSamples2D->At(i)->Clone();
        TString i_str;
        i_str.Form("%i",i);
        c = new TCanvas("c_"+i_str);
        c->SetFillColor(kWhite);
        stacks2D[i] = new THStack("stacks2D"+i_str,i_str);

        // j list MC samples for this channel i
        TH2D *tmpHist_draw2D;
        for(int j = 0; j < allMCSamples2D[i]->GetEntries(); j++)
        {
            TString j_str;
            j_str.Form("%i",j);
            tmpHist_draw2D = (TH2D*)allMCSamples2D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw2D->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            
            // k searching through array of params for the right one
            int which_param;
            bool foundParam = false;
            for(int k = 0; (k < numberParams) && !foundParam; k++)
            {
            
                // match up the isotope with parm
                for(int n = 0; (n < paramNameMap[k].size()) && !foundParam; n++)
                {
	                //std::cout <<"is it...  " <<paramNameMap[k].at(j) << std::endl;

	                if(tmpName.Contains(paramNameMap[k].at(n)))
                    {
	                    foundParam = true;
                        which_param = k;
	                }
	            }
            }

            if(foundParam)
            {
                tmpHist_draw2D->Scale(AdjustActs[which_param]);
     
                if(tmpHist_draw2D->Integral() > 0)
                {
                    stacks2D[i]->Add(tmpHist_draw2D);
                }
            }
            else
            {
                std::cout << "error could not find histograms" << std::endl;
            }

        }

        stacks2D[i]->Draw("");
        data2D[i]->SetLineWidth(2);
        data2D[i]->SetMarkerStyle(20);
        data2D[i]->SetMarkerSize(0.5);
        data2D[i]->Draw("PEsames");

        c->SaveAs("finalHisto2D_" + i_str + ".C");
        c->SaveAs("finalHisto2D_" + i_str + ".eps");
        c->SaveAs("finalHisto2D_" + i_str + ".png");
        std::cout << "saved 2D histogram to file" << std::endl;

    }

    #endif

}



void draw_inputdata()
{

    // debug
    TCanvas *c_nEqNull;
    c_nEqNull = new TCanvas("c_nEqNull", "c_nEqNull"); //, 4000, 3000);
    TH2D *h_nEqNull_clone = (TH2D*)h_nEqNull->Clone();
    h_nEqNull_clone->SetTitle("");
    h_nEqNull_clone->GetXaxis()->SetTitle("Electron Energy 2");
    h_nEqNull_clone->GetXaxis()->SetTitleFont(43);
    h_nEqNull_clone->GetXaxis()->SetTitleSize(20);
    h_nEqNull_clone->GetXaxis()->SetLabelFont(43);
    h_nEqNull_clone->GetXaxis()->SetLabelSize(20);
    h_nEqNull_clone->GetYaxis()->SetTitle("Electron Energy 1");
    h_nEqNull_clone->GetYaxis()->SetTitleFont(43);
    h_nEqNull_clone->GetYaxis()->SetTitleSize(20);
    h_nEqNull_clone->GetYaxis()->SetLabelFont(43);
    h_nEqNull_clone->GetYaxis()->SetLabelSize(20);
    h_nEqNull_clone->GetZaxis()->SetTitle("Decay Rate");
    h_nEqNull_clone->GetZaxis()->SetTitleFont(43);
    h_nEqNull_clone->GetZaxis()->SetTitleSize(20);
    h_nEqNull_clone->GetZaxis()->SetLabelFont(43);
    h_nEqNull_clone->GetZaxis()->SetLabelSize(20);
    c_nEqNull->SetRightMargin(0.15);
    h_nEqNull_clone->SetContour(1000);
    h_nEqNull_clone->Draw("colz");

    TLatex latexlabel;
    latexlabel.SetNDC();
    latexlabel.SetTextAlign(33);
    latexlabel.SetTextFont(43);
    latexlabel.SetTextSize(30);

    latexlabel.DrawLatex(0.80, 0.85, "#frac{dG_{0}}{dT_{ee}}");

    c_nEqNull->SaveAs("c_nEqNull.png");
    c_nEqNull->SaveAs("c_nEqNull.pdf");
    c_nEqNull->SaveAs("c_nEqNull.C");
    //delete c_nEqNull;

    // debug
    TCanvas *c_nEqTwo;
    c_nEqTwo = new TCanvas("c_nEqTwo", "c_nEqTwo"); //, 4000, 3000);
    TH2D *h_nEqTwo_clone = (TH2D*)h_nEqTwo->Clone();
    h_nEqTwo_clone->SetTitle("");
    h_nEqTwo_clone->GetXaxis()->SetTitle("Electron Energy 2");
    h_nEqTwo_clone->GetXaxis()->SetTitleFont(43);
    h_nEqTwo_clone->GetXaxis()->SetTitleSize(20);
    h_nEqTwo_clone->GetXaxis()->SetLabelFont(43);
    h_nEqTwo_clone->GetXaxis()->SetLabelSize(20);
    h_nEqTwo_clone->GetYaxis()->SetTitle("Electron Energy 1");
    h_nEqTwo_clone->GetYaxis()->SetTitleFont(43);
    h_nEqTwo_clone->GetYaxis()->SetTitleSize(20);
    h_nEqTwo_clone->GetYaxis()->SetLabelFont(43);
    h_nEqTwo_clone->GetYaxis()->SetLabelSize(20);
    h_nEqTwo_clone->GetZaxis()->SetTitle("Decay Rate");
    h_nEqTwo_clone->GetZaxis()->SetTitleFont(43);
    h_nEqTwo_clone->GetZaxis()->SetTitleSize(20);
    h_nEqTwo_clone->GetZaxis()->SetLabelFont(43);
    h_nEqTwo_clone->GetZaxis()->SetLabelSize(20);
    c_nEqTwo->SetRightMargin(0.15);
    //h_nEqTwo_clone->SetNumberContours(256);
    h_nEqTwo_clone->SetContour(1000);
    h_nEqTwo_clone->Draw("colz");

    latexlabel.DrawLatex(0.80, 0.85, "#frac{dG_{2}}{dT_{ee}}");

    c_nEqTwo->SaveAs("c_nEqTwo.png");
    c_nEqTwo->SaveAs("c_nEqTwo.pdf");
    c_nEqTwo->SaveAs("c_nEqTwo.C");
    //delete c_nEqTwo;
}













void draw(const draw_input_data &drawinputdata,
          const std::vector<double> &params,
          const std::vector<double> &param_errs)
{

    ///////////////////////////////////////////////////////////////////////////
    // draw result
    ///////////////////////////////////////////////////////////////////////////


    /*
    THStack *stacks1D[number1DHists];

    TH1D *h_2nubb[number1DHists];
    TH1D *h_tl208_int[number1DHists];
    TH1D *h_bi214_int[number1DHists];
    TH1D *h_bi207_int[number1DHists];
    TH1D *h_internal[number1DHists];
    TH1D *h_external[number1DHists];
    TH1D *h_radon[number1DHists];
    TH1D *h_neighbours[number1DHists];
    TH1D *h_other[number1DHists];

    for(int i = 0; i < number1DHists; ++ i)
    {
        h_2nubb[i] = nullptr;
        h_tl208_int[i] = nullptr;
        h_bi214_int[i] = nullptr;
        h_bi207_int[i] = nullptr;
        h_internal[i] = nullptr;
        h_external[i] = nullptr;
        h_radon[i] = nullptr;
        h_neighbours[i] = nullptr;
        h_other[i] = nullptr;
    }


    TCanvas *c[number1DHists];
    TPad *p0[number1DHists]; // 2020-06-12
    TPad *p1[number1DHists];
    TH1D *hRatio[number1DHists];
    TH1D *hAllMC1D[number1DHists];
    TH1D *data1D[number1DHists];
    TH1D *fakeData1D[number1DHists];
    */

    //std::cout << "debug: number of data samples: " << allDataSamples1D->GetEntries() << std::endl;
    //std::cout << "debug: number of MC samples: " << allMCSamples1D[0]->GetEntries() << std::endl;


    // each channel 1D hists
    // this is for(i = 0; i < 1; ++ i)
    // TODO: this isn't right. should this be iterating over the "channel" ?
    //for(int channel = 0; channel < allDataSamples1D->GetEntries(); ++ channel)
    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        if(channel_enable_draw_1D[channel] == false)
        {
            continue;
        }

        //if(i != 1) continue;
        // TODO: remove



        //draw_aux_data drawauxdata;
        draw_channel(channel,
                     drawinputdata,
                     params,
                     param_errs);

        /*
        stacks1D[channel] = drawauxdata.stacks1D;
        h_2nubb[channel] = drawauxdata.h_2nubb;
        h_tl208_int[channel] = drawauxdata.h_tl208_int;
        h_bi214_int[channel] = drawauxdata.h_bi214_int;
        h_bi207_int[channel] = drawauxdata.h_bi207_int;
        h_internal[channel] = drawauxdata.h_internal;
        h_external[channel] = drawauxdata.h_external;
        h_radon[channel] = drawauxdata.h_radon; 
        h_neighbours[channel] = drawauxdata.h_neighbours;
        h_other[channel] = drawauxdata.h_other;


        c[channel] = drawauxdata.c;
        p0[channel] = drawauxdata.p0;
        p1[channel] = drawauxdata.p1;
        hRatio[channel] = drawauxdata.hRatio;
        hAllMC1D[channel] = drawauxdata.hAllMC1D;
        data1D[channel] = drawauxdata.data1D;
        fakeData1D[channel] = drawauxdata.fakeData1D;
        */
    }




}


// I stare down the concrete building,
// My workplace, something that is masquerading as an institution of learning
// but has long since lost its nerve
// Stride boldly on my right foot,
// Into the den of criminals who seek to destroy me
// I collect a sphere of repulsion,
// My sheild that I carry throughout the day to repell all those who set eyes on me
// On the evening I return to my home, to my bed,
// Where I lie troubled and in deep contemplation, yet still peaceful
// All others rush into the night, they dash through their door and bolt it shut tight
// Hand quivering and fumbling the key
// They race to their beds and cling tightly to their duvet as if disturbed
// by the hand of death himself
// Bite and tear at the cloth, a duvet endowed with the scripture of their religion
// Diversity is Our Strength.
// It is of no concolace.
// A simple religion with no depth. A single word to repeat over.
// And over. And over.
// Long into the night.
// Until they open their eyes to the darkness.
// And for long they scream





#if 0
void draw_get_total_data_MC(
    const std::vector<double> &params,
    const std::vector<double> &param_errs,
    const double fval,
    TH1D *&hHighEnergy_allMC_out,
    TH1D *&hLowEnergy_allMC_out,
    TH1D *&hHighEnergy_data_out,
    TH1D *&hLowEnergy_data_out,
    bool g_mode_fake_data = false,
    int channel_index)
{

    ///////////////////////////////////////////////////////////////////////////
    // get total data and MC
    ///////////////////////////////////////////////////////////////////////////


    THStack *stacks1D[number1DHists];

    TH1D *h_2nubb[number1DHists];
    TH1D *h_tl208_int[number1DHists];
    TH1D *h_bi214_int[number1DHists];
    TH1D *h_bi207_int[number1DHists];
    TH1D *h_internal[number1DHists];
    TH1D *h_external[number1DHists];
    TH1D *h_radon[number1DHists];
    TH1D *h_neighbours[number1DHists];
    TH1D *h_other[number1DHists];

    for(int i = 0; i < number1DHists; ++ i)
    {
        h_2nubb[i] = nullptr;
        h_tl208_int[i] = nullptr;
        h_bi214_int[i] = nullptr;
        h_bi207_int[i] = nullptr;
        h_internal[i] = nullptr;
        h_external[i] = nullptr;
        h_radon[i] = nullptr;
        h_neighbours[i] = nullptr;
        h_other[i] = nullptr;
    }


    TH1D *hAllMC1D[number1DHists];
    TH1D *data1D[number1DHists];
    TH1D *fakeData1D[number1DHists];

    // each channel 1D hists
    // this is for(i = 0; i < 1; ++ i)
    // TODO: this isn't right. should this be iterating over the "channel" ?
    for(int i = 0; i < allDataSamples1D->GetEntries(); i++)
    {
        if(draw_index == -1)
        {
            // do nothing
            // draw all channels
        }
        else if((draw_index >= 0) && (draw_index < allDataSamples1D->GetEntries()))
        {
            if(i != draw_index)
            {
                continue;
            }
        }

        //if(i != 1) continue;
        // TODO: remove


        // because this isn't right TODO
        // uses at(i), but i should always be zero and there should be an
        // additional array index
        data1D[i] = (TH1D*)allDataSamples1D->At(i)->Clone();
        if(g_mode_fake_data == true)
        {
            fakeData1D[i] = (TH1D*)allFakeDataSamples1D->At(i)->Clone(); // TODO
            // TODO: will not work if logLikelihood not called before
            // because LL function calls function to construct fakedata

            /*
            for(int k = 0; k < fakeData1D[i]->GetNbinsX(); ++ k)
            {
                std::cout << "k=" << k << " " << fakeData1D[i]->GetBinContent(k) << std::endl;
            }
            */
        }

        TString i_str;
        i_str.Form("%i", i);
    
        c[i] = new TCanvas("c" + i_str, "c" + i_str);
        c[i]->SetFillColor(kWhite);

        
        stacks1D[i] = new THStack("stacks1D" + i_str, i_str);

        TH1D *tmpHist;

        // j list MC samples for this channel i
        //std::cout << "debug: number of MC samples (i=" << i << "): " << allMCSamples1D[i]->GetEntries() << std::endl;

        /*
        std::cout << "list of all objects in allMCSamples1D[i=" << i << "]" << std::endl;
        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {
            //std::cout << (((TH1D*)allMCSamples1D[i])->At(j))->GetName() << std::endl;
            tmpHist_draw1D = (TH1D*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw1D->GetName();
            std::cout << tmpName << std::endl;
        }
        */

        // allMCSamples1D[0] contains objects such as: "zr96_rot_k40_2e_P2"
        
        

        /*
        std::cout << "contents of map: MCNameToParamNumberMap" << std::endl;
        for(auto it = MCNameToParamNumberMap.cbegin(); it != MCNameToParamNumberMap.cend(); ++ it)
        {
            std::cout << it->first << " -> " << it->second << std::endl;
        }
        */


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

                    // set branching ratio fraction
                    /*
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
                    
                        //std::cout << "paramNumber=" << paramNumber << " -> tmp_sample_name=" << tmp_sample_name << " ~> tmpName=" << tmpName << std::endl;                    
                        //which_param = paramNumber;
                        which_param = paramNumberToMinuitParamNumberMap.at(paramNumber);
                        found_param = true;
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


            /* if I rebuild the MC in the fitHistograms loop
             * then I don't need to reweight it here
              */

                // no error thrown, which_param is presumably the correct index
                Double_t activity_scale = params.at(which_param); // * activity_scale_branching_ratio;
                tmpHist->Scale(activity_scale);

                if(which_param == 1)
                {
                    std::cout << "error " << __func__ << std::endl;
                    throw "error";
                }

                //TH1D *tmpHist_drawpointer = tmpHist_draw1D;

                
                if(tmpHist->Integral() > 0)
                {
                    stacks1D[i]->Add((TH1D*)tmpHist->Clone());
                    //stacks1D[i]->Add(tmpHist_drawpointer);

                    stacker_helper(tmpHist,
                                 h_2nubb[i],
                                 h_tl208_int[i],
                                 h_bi214_int[i],
                                 h_bi207_int[i],
                                 h_internal[i],
                                 h_neighbours[i],
                                 h_radon[i],
                                 h_external[i],
                                 h_other[i]);


                    if(j == 0)
                    {
                        hAllMC1D[i] = (TH1D*)tmpHist->Clone("Total MC");
                    }
                    else
                    {
                        hAllMC1D[i]->Add((TH1D*)tmpHist);
                    }

                
	            }
                else
                {
                    //std::cout << "not adding to stack, Integral() <= 0: " << tmpHist_draw1D->GetName() << std::endl;
                }
            }
            else
            {
                std::cout << __func__ << " error could not find histogram: tmpName=" << tmpName << std::endl;
            } 

        }


        THStack *stacks1D_major[number1DHists];
        stacks1D_major[i] = new THStack("stacks1D_major" + i_str, i_str);
        stacker_helper_2(stacks1D_major[i],
                         h_2nubb[i],
                         h_tl208_int[i],
                         h_bi214_int[i],
                         h_bi207_int[i],
                         h_internal[i],
                         h_neighbours[i],
                         h_radon[i],
                         h_external[i],
                         h_other[i]);


        double PAD_U_Y_MIN = 0.0;
        double PAD_U_Y_MAX = 500.0;
        double PAD_L_Y_MAX = 3.0; // 1.2
        double PAD_L_Y_MIN = 0.0; // 0.8

        if(i == 0)
        {
            PAD_U_Y_MAX = 350.0;
            PAD_U_Y_MAX = 300.0;
        }
        else if(i == 1)
        {
            PAD_U_Y_MAX = 1000.0;
        }
        else if(i == 2)
        {
            PAD_U_Y_MAX = 450.0;
        }
        else if(i == 3)
        {
            PAD_U_Y_MAX = 1000.0;
        }
        else if(i == 4)
        {
            PAD_U_Y_MAX = 400.0;
        }
        else if(i == 5)
        {
            PAD_U_Y_MAX = 600.0;
        }
        else
        {
            PAD_U_Y_MAX = 350.0;
        }
        //stacks1D_major[i]->Draw("hist");
        //stacks1D_major[i]->SetMaximum(PAD_U_Y_MAX);
        //stacks1D_major[i]->SetMinimum(PAD_U_Y_MIN);
        //stacks1D_major[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
        // TODO moved

        //hAllMC1D[i]->SetMaximum(PAD_U_Y_MAX);
        //hAllMC1D[i]->SetMinimum(PAD_U_Y_MIN);
        hAllMC1D[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
        //data1D[i]->SetMaximum(PAD_U_Y_MAX);
        //data1D[i]->SetMinimum(PAD_U_Y_MIN);
        data1D[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
        if(g_mode_fake_data == true)
        {
            //fakeData1D[i]->SetMaximum(PAD_U_Y_MAX);
            //fakeData1D[i]->SetMinimum(PAD_U_Y_MIN);
            fakeData1D[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
        }


        if(g_mode_fake_data == false)
        {
            hRatio[i] = (TH1D*)data1D[i]->Clone();
        }
        if(g_mode_fake_data == true)
        {
            hRatio[i] = (TH1D*)fakeData1D[i]->Clone();
        }
        hRatio[i]->Sumw2();
        hRatio[i]->Divide(hAllMC1D[i]);
        hRatio[i]->SetTitle("");

        std::cout << "draw-> MC Integral() = " << hAllMC1D[i]->Integral() << std::endl;
        if(g_mode_fake_data == true)
        {
            std::cout << "draw-> fakeData Integral() = " << fakeData1D[i]->Integral() << std::endl;
        }


    // TODO: some other stuff goes here





        if(i == 0)
        {
            hRatio[i]->GetXaxis()->SetTitle("2e Electron Energy [MeV]");

            //PAD_L_Y_MAX = 10.0;
            //PAD_L_Y_MIN = -10.0;
        }
        else if(i == 1)
        {
            hRatio[i]->GetXaxis()->SetTitle("Single Electron Energy [MeV]");

            //PAD_L_Y_MAX = 10.0;
            //PAD_L_Y_MIN = -10.0;
        }
        else if(i == 2)
        {
            hRatio[i]->GetXaxis()->SetTitle("High Energy Electron [MeV]");

            //PAD_L_Y_MAX = 10.0;
            //PAD_L_Y_MIN = -10.0;
        }
        else if(i == 3)
        {
            hRatio[i]->GetXaxis()->SetTitle("Low Energy Electron [MeV]");

            //PAD_L_Y_MAX = 10.0;
            //PAD_L_Y_MIN = -10.0;
        }
        //hRatio[i]->SetMaximum(PAD_L_Y_MAX);
        //hRatio[i]->SetMinimum(PAD_L_Y_MIN);
        hRatio[i]->GetYaxis()->SetRangeUser(PAD_L_Y_MIN, PAD_L_Y_MAX);

        hRatio[i]->GetXaxis()->SetTickSize(0.1);

        //hRatio[i]->GetXaxis()->SetTitle("Electron Energy [MeV]");
        hRatio[i]->GetXaxis()->SetTitleFont(43);
        hRatio[i]->GetXaxis()->SetTitleSize(20);
        hRatio[i]->GetXaxis()->SetTitleOffset(3.0);

        hRatio[i]->GetYaxis()->SetTitle("data / MC");
        hRatio[i]->GetYaxis()->SetTitleFont(43);
        hRatio[i]->GetYaxis()->SetTitleSize(20);
        hRatio[i]->GetYaxis()->SetTitleOffset(1.0);

        hRatio[i]->GetYaxis()->SetLabelFont(43);
        hRatio[i]->GetYaxis()->SetLabelSize(0.0);

        hRatio[i]->GetYaxis()->SetTickSize(0.0);

        hRatio[i]->GetXaxis()->SetLabelFont(43);
        hRatio[i]->GetXaxis()->SetLabelSize(15);

        // in code copying from, canvas alloc. here
        c[i]->cd();
        c[i]->SetBottomMargin(0.0);
        p0[i] = new TPad("pad0", "pad0", 0.0, 0.3, 1.0, 1.0);
        p0[i]->SetBottomMargin(0.0);
        //p0[i]->SetGridx(1);
        //p0[i]->SetGridy(1);
        p0[i]->SetGrid(0, 0);
        p0[i]->SetTicks(2, 2);
        p0[i]->Draw();

        c[i]->cd();
        p1[i] = new TPad("pad1", "pad1", 0.0, 0.0, 1.0, 0.3);
        p1[i]->SetTopMargin(0.0);
        p1[i]->SetBottomMargin(0.4);
        //p1[i]->SetGridx(1);
        //p1[i]->SetGridy(1);
        p1[i]->SetGrid(0, 0);
        p1[i]->SetTicks(2, 2);
        p1[i]->Draw();

        p0[i]->cd();

        // draw regular pad1
        // leave alone for now
        // TODO: axis
        // copy code from other file
        hAllMC1D[i]->SetTitle("");


        stacks1D_major[i]->Draw("hist");
        stacks1D_major[i]->SetMaximum(PAD_U_Y_MAX);
        stacks1D_major[i]->SetMinimum(PAD_U_Y_MIN);
        stacks1D_major[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);

        stacks1D_major[i]->SetTitle("");

        stacks1D_major[i]->GetYaxis()->SetTickSize(0.0);
        //stacks1D_major[i]->GetXaxis()->SetTickSize(0.0);

        stacks1D_major[i]->GetYaxis()->SetTitle("Spectrum");
        stacks1D_major[i]->GetYaxis()->SetTitleSize(20);
        stacks1D_major[i]->GetYaxis()->SetTitleFont(43);
        stacks1D_major[i]->GetYaxis()->SetTitleOffset(1.0);

        stacks1D_major[i]->GetYaxis()->SetLabelFont(43);
        stacks1D_major[i]->GetYaxis()->SetLabelSize(0.0);

        stacks1D_major[i]->GetXaxis()->SetTitle("");
        stacks1D_major[i]->GetXaxis()->SetTitleSize(0);
        stacks1D_major[i]->GetXaxis()->SetTitleFont(43);
        stacks1D_major[i]->GetXaxis()->SetTitleOffset(1.0);

        stacks1D_major[i]->GetXaxis()->SetLabelSize(0);
        stacks1D_major[i]->GetXaxis()->SetLabelFont(63);


        hAllMC1D[i]->SetLineWidth(2);
        hAllMC1D[i]->SetLineColor(kBlack);
        hAllMC1D[i]->SetFillColor(kWhite);
        hAllMC1D[i]->SetFillStyle(0);
        //hAllMC1D[i]->Sumw2();
        //hAllMC1D[i]->Draw("hist sames");
        TString Nmc_str;
        Nmc_str.Form("%i", (int)hAllMC1D[i]->Integral()); // TODO: float?
        hAllMC1D[i]->SetTitle("Total MC (" + Nmc_str + ")");
        hAllMC1D[i]->Draw("hist same");
        hAllMC1D[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
        data1D[i]->SetLineWidth(2);
        data1D[i]->SetMarkerStyle(20);
        data1D[i]->SetMarkerSize(1.0);
        data1D[i]->SetLineColor(kBlack);
        data1D[i]->SetMarkerColor(kBlack);
        data1D[i]->SetFillColor(kBlack);
        if(g_mode_fake_data == true)
        {
            fakeData1D[i]->SetLineWidth(2);
            fakeData1D[i]->SetMarkerStyle(20);
            fakeData1D[i]->SetMarkerSize(1.0);
            fakeData1D[i]->SetLineColor(kBlack);
            fakeData1D[i]->SetMarkerColor(kBlack);
            fakeData1D[i]->SetFillColor(kBlack);
        }
        TString Ndata_str;
        Ndata_str.Form("%i", (int)data1D[i]->Integral()); // TODO: float?
        data1D[i]->SetTitle("Data (" + Ndata_str + ")");
        TString Nfakedata_str;
        if(g_mode_fake_data == true)
        {
            Nfakedata_str.Form("%i", (int)fakeData1D[i]->Integral()); // TODO: float?
            fakeData1D[i]->SetTitle("Fake Data (" + Ndata_str + ")");
        }
        if(g_mode_fake_data == false)
        {
            //data1D[i]->Draw("PEsames");
            data1D[i]->Draw("PEsame");
            data1D[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
        }
        if(g_mode_fake_data == true)
        {
            fakeData1D[i]->Draw("PEsame");
            fakeData1D[i]->GetYaxis()->SetRangeUser(PAD_U_Y_MIN, PAD_U_Y_MAX);
        }

        //double chi2;
        int ndf = -1;
        ndf = 0;
        get_ndf_1D();
        //int ndf;
        //int igood;
        //TString chi2_str;
        TString ndf_str;

        // TODO: should chisquare value include the constraints? because at
        // the moment it does not

        // TODO: chi2 value is different from fit_2e code
        //double prob = data1D[i]->Chi2TestX(hAllMC1D[i], chi2, ndf, igood, "UW");
        // TODO this number produced by Chi2TestX is nonsense
        // TODO: fakedata

        // TODO: check if I can get fcn value from the minuit fit object
        //chi2_str.Form("%4.3f", chi2);
        TString fval_str;
        fval_str.Form("%4.3f", fval);
        ndf_str.Form("%i", ndf);


        
        TGaxis *axis = new TGaxis(0.0, PAD_U_Y_MIN + 0.01, 0.0, PAD_U_Y_MAX, PAD_U_Y_MIN + 0.01, PAD_U_Y_MAX, 510, "");
        axis->SetLabelFont(43);
        axis->SetLabelSize(15);
        axis->Draw();

        TGaxis *axis2 = new TGaxis(5.0, PAD_U_Y_MIN + 0.01, 5.0, PAD_U_Y_MAX, PAD_U_Y_MIN + 0.01, PAD_U_Y_MAX, 510, "+");
        axis2->SetLabelFont(43);
        axis2->SetLabelSize(0);
        axis2->Draw();

        TLegend *leg = new TLegend(0.6, 0.1, 0.85, 0.85);
        if(g_mode_fake_data == false)
        {
            leg->AddEntry(data1D[i], "Data (" + Ndata_str + ")", "PE");
        }
        if(g_mode_fake_data == true)
        {
            leg->AddEntry(fakeData1D[i], "Fake Data (" + Nfakedata_str + ")", "PE");
        }
        leg->AddEntry(hAllMC1D[i], "Total MC (" + Nmc_str + ")", "EP");
        leg->AddEntry(h_2nubb[i], "2#nu#beta#beta", "F");
        leg->AddEntry(h_tl208_int[i], "^{208}Tl Internal", "F");
        leg->AddEntry(h_bi214_int[i], "^{214}Bi Internal", "F");
        leg->AddEntry(h_bi207_int[i], "^{207}Bi Internal", "F");
        leg->AddEntry(h_internal[i], "Internal", "F");
        leg->AddEntry(h_neighbours[i], "Neighbour Foil", "F");
        leg->AddEntry(h_radon[i], "Radon", "F");
        leg->AddEntry(h_external[i], "External", "F");
        //leg->AddEntry((TObject*)nullptr, "#chi^{2}/ndf=" + chi2_str + "/" + ndf_str, "");
        leg->AddEntry((TObject*)nullptr, "#chi^{2}/ndf=" + fval_str + "/" + ndf_str, "");
        //leg->AddEntry(h_other[i], "other", "f");
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        //leg->SetTextFont(62);
        leg->SetTextFont(63);
        //leg->SetTextSize(0.035);
        leg->SetTextSize(15);
        leg->SetShadowColor(kGray + 2);
        leg->SetBorderSize(5);
        leg->Draw("BR");

        //TLatex latexlabel;
        //latexlabel.SetNDC();
        //latexlabel.SetTextFont(62);
        //latexlabel.SetTextSize(0.035);
        //latexlabel.DrawLatex(0.63, 0.23, "#frac{#chi^{2}}{ndf} = #frac{" + chi2_str + "}{" + ndf_str + "}");


        // second pad
        p1[i]->cd();
        hRatio[i]->SetMarkerStyle(20);
        hRatio[i]->SetMarkerSize(1.0);
        hRatio[i]->SetLineWidth(2);
        hRatio[i]->Draw("EP");



        TGaxis *axis3 = new TGaxis(0.0, PAD_L_Y_MIN, 0.0, PAD_L_Y_MAX - 0.0001, PAD_L_Y_MIN, PAD_L_Y_MAX - 0.0001, 503, "");
        axis3->SetLabelFont(43);
        axis3->SetLabelSize(15);
        axis3->Draw();

        TGaxis *axis4 = new TGaxis(5.0, PAD_L_Y_MIN, 5.0, PAD_L_Y_MAX - 0.0001, PAD_L_Y_MIN, PAD_L_Y_MAX - 0.0001, 503, "+");
        axis4->SetLabelFont(43);
        axis4->SetLabelSize(0);
        axis4->Draw();

    


        c[i]->Show();
        //c[i]->Draw();
        c[i]->Update();

    
        if(saveas_filename != "NOSAVE")
        {
            std::string base_name;
            std::string extension;
            filename_split_extension(saveas_filename, base_name, extension);

            //std::string dir = base_name + "_c" + "_" + std::string(i_str);
            std::string dir = saveas_dir;
            std::string name = base_name + "_c" + "_" + std::string(i_str) + extension;
            canvas_saveas_helper(dir, name, c[i]);
        }
    
        delete c[i];
        c[i] = nullptr;




    }




        // channel 2 = high energy
        // channel 3 = low energy
        if(i == 2)
        {
            hHighEnergy_allMC_out = (TH1D*)hAllMC1D[i]->Clone();
            if(g_mode_fake_data == false)
            {
                hHighEnergy_data_out = (TH1D*)data1D[i]->Clone();
            }
            if(g_mode_fake_data)
            {
                hHighEnergy_data_out = (TH1D*)fakeData1D[i]->Clone();
            }
        }
        else if(i == 3)
        {
            hLowEnergy_allMC_out = (TH1D*)hAllMC1D[i]->Clone();
            if(g_mode_fake_data == false)
            {
                hLowEnergy_data_out = (TH1D*)data1D[i]->Clone();
            }
            if(g_mode_fake_data)
            {
                hLowEnergy_data_out = (TH1D*)fakeData1D[i]->Clone();
            }
        }
    }

}
#endif








// TODO: number job id
void draw_covariance_matrix(const double * const CovMatrix, const int number_free_params, const std::string& saveas_filename)
{


    ///////////////////////////////////////////////////////////////////////////
    // draw covarience matrix
    ///////////////////////////////////////////////////////////////////////////

    std::cout << "draw: covariance matrix" << std::endl;

    /*
    std::cout << "contents of free_params:" << std::endl;
    for(int i = 0; i < free_params.size(); ++ i)
    {
        std::cout << free_params.at(i) << std::endl;
    }
    */

    // do not store names elsewhere, create them here
#if 0
    TH2D* hCorrMatrix = new TH2D("hCorrMatrix", "Correlation Matrix",
                                 number_free_params, 0, (double)number_free_params,
                                 number_free_params, 0, (double)number_free_params);

    std::cout << "number_free_params=" << number_free_params << std::endl;
    for(int i = 0; i < number_free_params; i++)
    {
        for(int j = 0; j < number_free_params; j++)
        {
            //std::cout << index_free_params.at(i) << "  " << index_free_params.at(j)<< std::endl;
            //CovMatrix[i][j] = minuit->GetCovarianceMatrixElement(i,j);

            // std::cout << CovMatrix[i][j] << std::endl;
            //Double_t value = CovMatrix[i][j] / (std::sqrt(CovMatrix[i][i]) * std::sqrt(CovMatrix[j][j]));
            //
            Double_t CovMatrix_i_j = CovMatrix[i * number_free_params + j];
            Double_t CovMatrix_i_i = CovMatrix[i * number_free_params + i];
            Double_t CovMatrix_j_j = CovMatrix[j * number_free_params + j];

            Double_t value = CovMatrix_i_j / (std::sqrt(CovMatrix_i_i) * std::sqrt(CovMatrix_j_j));
            //hCorrMatrix->Fill(free_params_names.at(i), free_params_names.at(j), value);
            //hCorrMatrix->Fill(i, j, value);
            
            // need to ignore all HARD constrained (minuit->FixParameter)
            // parameters
            // I don't know if I like this solution, might be better to
            // iterate over the vectors themselves
            // TODO: what order is this vector in and how does that affect
            // output?
            // TODO: change to another incremental index, and have an
            // if(find(fixed params)) statement
            //std::cout << "i=" << i << " j=" << j << std::endl;

            int i_free_param_index = free_params.at(i);
            int j_free_param_index = free_params.at(j);
            //std::cout << "i_free_param_index=" << i_free_param_index << " j_free_param_index=" << j_free_param_index << std::endl;

            // convert minuit/internal parameter
            // ("free parameter" meaning not SOFT or FREE constrain mode)
            // to external parameter number
            //int i_external = minuitParamNumberToParamNumberMap.at(i_free_param_index);
            //int j_external = minuitParamNumberToParamNumberMap.at(j_free_param_index);
            //std::cout << "i_external=" << i_external << " j_external=" << j_external << std::endl;
            // NOTE: using method of accessing vector elements in free_params
            // gives index which is already in external index format
            int i_external = i_free_param_index;
            int j_external = j_free_param_index;

            // convert external parameter number to parameter name
            TString free_param_name_i = paramNameMap[i_external];
            TString free_param_name_j = paramNameMap[j_external];
            //std::cout << "free_param_name_i=" << free_param_name_i << " free_param_name_j=" << free_param_name_j << std::endl;

            // convert parameter name to human readable name
            //TString free_param_human_name_i = paramNameToHumanParamNameMap.at(free_param_name_i);
            //TString free_param_human_name_j = paramNameToHumanParamNameMap.at(free_param_name_j);
            //TString free_param_human_name_i = paramNameToHumanReadableParamNameMap.at(free_param_name_i);
            //TString free_param_human_name_j = paramNameToHumanReadableParamNameMap.at(free_param_name_j);
            //std::cout << "free_param_human_name_i=" << free_param_human_name_i << " free_param_human_name_j=" << free_param_human_name_j << std::endl;

            TString free_param_human_name_i = "NONAME";
            TString free_param_human_name_j = "NONAME";
            // if available, use name from manually assigned names in parameter_names.lst file
            if(paramNumberToHumanReadableParamNameMap.count(i_external) == 1)
            {
                free_param_human_name_i = paramNumberToHumanReadableParamNameMap.at(i_external);
            }
            else
            {
                free_param_human_name_i = paramNameToHumanReadableParamNameMap.at(free_param_name_i);
            }
            if(paramNumberToHumanReadableParamNameMap.count(j_external) == 1)
            {
                free_param_human_name_j = paramNumberToHumanReadableParamNameMap.at(j_external);
            }
            else
            {
                free_param_human_name_j = paramNameToHumanReadableParamNameMap.at(free_param_name_j);
            }

            // fill
            hCorrMatrix->Fill(free_param_human_name_i, free_param_human_name_j, value);

        }
    }
    // TODO: re-enable this check chi2 value is the same
    // remove global variable for chi2?
    // these were commented out for some reason, I did not remove them
    /*
    double chi2, edm, errdef; 
    int nvpar, nparx;
    minuit->GetStats(chi2, edm, errdef, nvpar, nparx);
    int ndf = npfits-nvpar;
    */


    // TODO: margin
    //c->SetLeftMargin(0.15)
    //c->SetTopMargin(0.15);
    
    // create a function to draw the final stacked histograms
    gStyle->SetPaintTextFormat("4.2f");
    Style_t gridstyle = gStyle->GetGridStyle();
    gStyle->SetGridStyle(0);

    hCorrMatrix->SetTitle("Parameter Correlation Matrix");
    hCorrMatrix->LabelsDeflate();
    hCorrMatrix->GetXaxis()->LabelsOption("v");
    hCorrMatrix->GetXaxis()->SetLabelSize(0.035);
    hCorrMatrix->GetYaxis()->SetLabelSize(0.035);
    hCorrMatrix->GetZaxis()->SetRangeUser(-1.0, 1.0);

    //TCanvas *c2 = new TCanvas("c2");
    TCanvas *c2 = new TCanvas("c_cormatrix", "Parameter Correlation Matrix", 200, 10, 1013, 885);
    c2->SetLeftMargin(0.18);
    c2->SetRightMargin(0.14);
    c2->SetTopMargin(0.14);
    c2->SetBottomMargin(0.21);
    c2->SetFillColor(kWhite);
    c2->SetGrid();
    hCorrMatrix->Draw("colz text");
   
    // clean
    delete [] CovMatrix;

    // TODO: save plots

    //if(saveas_filename.size() > 0)
    //{
    //    c2->SaveAs(saveas_filename.c_str());
    //}
    if(saveas_filename.size() > 0)
    {
        std::size_t length = saveas_filename.size();
        if(length >= 2)
        {
            if((saveas_filename[length - 1] == '*') && (saveas_filename[length - 2] == '.'))
            {
                std::string base_name = saveas_filename.substr(0, length - 2);
                std::vector<std::string> extensions;
                extensions.push_back("C");
                extensions.push_back("png");
                extensions.push_back("eps");
                extensions.push_back("pdf");
                for(auto it{extensions.begin()}; it != extensions.end(); ++ it)
                {
                    std::string name = base_name + "." + *it;
                    c2->SaveAs(name.c_str());
                }
            }
        }
        else
        {
            c2->SaveAs(saveas_filename.c_str());
        }
    }
    
    gStyle->SetGridStyle(gridstyle);

#endif  
}








#endif // NEWLOGLIKFITTER_DRAW_H
