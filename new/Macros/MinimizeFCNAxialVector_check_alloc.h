#ifndef MINIMIZEFCNAXIALVECTOR_CHECK_ALLOC_H
#define MINIMIZEFCNAXIALVECTOR_CHECK_ALLOC_H

// not used

#if 0
void
MinimizeFCNAxialVector::check_alloc_V_PHYS_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS
    ///////////////////////////////////////////////////////////////////

    if(V_PHYS_1D_P1_data[0] == nullptr)
    {
        std::cout << "Alloc V_PHYS" << std::endl;

        // 1D
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            std::cout << "ALLOC " << __func__ << std::endl;
        }

        // 2D
        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("V_PHYS_2D_P1_CH%d", ch);
            V_PHYS_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("V_PHYS_2D_P2_CH%d", ch);
            V_PHYS_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */
    }

}
#endif

void
MinimizeFCNAxialVector::check_alloc_V_PHYS_STAT_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS_STAT
    ///////////////////////////////////////////////////////////////////
    //if(V_PHYS_STAT_1D_P1[0] == nullptr)
    
    if(V_PHYS_STAT_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc V_PHYS_STAT" << std::endl;
        //V_PHYS_STAT = new TH2D("V_PHYS_STAT", "V_PHYS_STAT",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);
        std::cout << "Alloc V_PHYS_STAT" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            V_PHYS_STAT_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            V_PHYS_STAT_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
            //std::cout << "ALLOC " << __func__ << std::endl;
        }

        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("V_PHYS_STAT_2D_P1_CH%d", ch);
            //V_PHYS_STAT_2D_P1[ch] = new TH2D(hname, hname,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("V_PHYS_STAT_2D_P2_CH%d", ch);
            //V_PHYS_STAT_2D_P2[ch] = new TH2D(hname, hname,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
            //                  NUM_BINS_XY, 0.0, NUM_BINS_XY);

            V_PHYS_STAT_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY * NUM_BINS_XY, 0.0);
        }
        */
    }

}


void
MinimizeFCNAxialVector::set_V_PHYS_STAT_data() const
{

    ///////////////////////////////////////////////////////////////////
    // V_PHYS_STAT
    ///////////////////////////////////////////////////////////////////

    // done in set_V_PHYS_MATHMORE

}




void MinimizeFCNAxialVector::check_alloc_D() const
{

    ///////////////////////////////////////////////////////////////////
    // V_MATRIX D (data)
    ///////////////////////////////////////////////////////////////////

    //if(D_1D_P1[0] == nullptr)
    if(D_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc D" << std::endl;
        //D = new TH2D("D", "D",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  1, 0.0, 1.0);

        // clear the contents of D
        //for(Int_t bin_i{1}; bin_i <= D->GetNbinsX(); ++ bin_i)
        //{
        //    const Double_t zero = 0.0;
        //    D->SetBinContent(bin_i, 1, zero);
        //}
        std::cout << "Alloc V_D" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            D_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            D_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            //std::cout << "ALLOC " << __func__ << std::endl;
        }

        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("D_2D_P1_CH%d", ch);
            D_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("D_2D_P2_CH%d", ch);
            D_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */
    }

}



void MinimizeFCNAxialVector::check_alloc_M() const
{

    ///////////////////////////////////////////////////////////////////
    // V_MATRIX M (MC)
    ///////////////////////////////////////////////////////////////////

    //if(M_1D_P1[0] == nullptr)
    if(M_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc M" << std::endl;
        //M = new TH2D("M", "M",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  1, 0.0, 1.0);

        // clear the contents of M
        //for(Int_t bin_i{1}; bin_i <= M->GetNbinsX(); ++ bin_i)
        //{
        //    const Double_t zero = 0.0;
        //    M->SetBinContent(bin_i, 1, zero);
        //}

        // 1D
        std::cout << "Alloc V_M" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            M_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            M_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            //std::cout << "ALLOC " << __func__ << std::endl;
        }

        // 2D
        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("M_2D_P1_CH%d", ch);
            M_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("M_2D_P2_CH%d", ch);
            M_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */
    }

}



void MinimizeFCNAxialVector::check_alloc_D_minus_M() const
{

    ///////////////////////////////////////////////////////////////////
    // V_MATRIX D_minus_M (Data - MC)
    ///////////////////////////////////////////////////////////////////

    //if(D_minus_M_1D_P1[0] == nullptr)
    if(D_minus_M_1D_P1_data[0] == nullptr)
    {
        //std::cout << "Alloc D_minus_M" << std::endl;
        //D_minus_M = new TH2D("D_minus_M", "D_minus_M",
        //                  NUM_BINS_XY, 0.0, NUM_BINS_XY,
        //                  1, 0.0, 1.0);

        // clear the contents of M
        //for(Int_t bin_i{1}; bin_i <= D_minus_M->GetNbinsX(); ++ bin_i)
        //{
        //    const Double_t zero = 0.0;
        //    D_minus_M->SetBinContent(bin_i, 1, zero);
        //}

        // 1D
        std::cout << "Alloc D_minus_M" << std::endl;
        for(int ch = 0; ch < number1DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50;
            D_minus_M_1D_P1_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            D_minus_M_1D_P2_data[ch] = new std::vector<double>(NUM_BINS_XY, 0.0);
            //std::cout << "ALLOC " << __func__ << std::endl;
        }


        // 2D
        /*
        for(int ch = 0; ch < number2DHists; ++ ch)
        {
            const Int_t NUM_BINS_XY = 50 * 50;
            TString hname;

            hname.Form("D_minus_M_2D_P1_CH%d", ch);
            D_minus_M_2D_P1[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);

            hname.Form("D_minus_M_2D_P2_CH%d", ch);
            D_minus_M_2D_P2[ch] = new TH2D(hname, hname,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY,
                              NUM_BINS_XY, 0.0, NUM_BINS_XY);
        }
        */

    }

}



// not used
#if 0
void
MinimizeFCNAxialVector::zero_V_PHYS_data() const
{
    // Set to zero
    for(int ch = 0; ch < number1DHists; ++ ch)
    {
        if(channel_enable_1D[ch] == 1)
        {
            /*
            for(Int_t bi = 1; bi <= V_PHYS_1D_P1[ch]->GetNbinsX(); ++ bi)
            {
                const Double_t zero = 0.0;
                V_PHYS_1D_P1[ch]->SetBinContent(bi, 1, zero);
            }

            for(Int_t bi = 1; bi <= V_PHYS_1D_P2[ch]->GetNbinsX(); ++ bi)
            {
                const Double_t zero = 0.0;
                V_PHYS_1D_P2[ch]->SetBinContent(bi, 1, zero);
            }
            */

            const int NUM_BINS_XY = 50;
            for(int i = 0; i < NUM_BINS_XY * NUM_BINS_XY; ++ i)
            {
                const double zero = 0.0;
                #if VECTOR_RANGE_CHECK
                V_PHYS_1D_P1_data[ch]->at(i) = zero;
                V_PHYS_1D_P2_data[ch]->at(i) = zero;
                #else
                V_PHYS_1D_P1_data[ch]->operator[](i) = zero;
                V_PHYS_1D_P2_data[ch]->operator[](i) = zero;
                #endif
            }
        }
    }

    /*
    // Set to zero
    for(int ch = 0; ch < number2DHists; ++ ch)
    {
        if(channel_enable_2D[ch] == 1)
        {
            for(Int_t bj = 1; bj <= V_PHYS_2D_P1[ch]->GetNbinsY(); ++ bj)
            {
                for(Int_t bi = 1; bi <= V_PHYS_2D_P1[ch]->GetNbinsX(); ++ bi)
                {
                    const Double_t zero = 0.0;
                    V_PHYS_2D_P1[ch]->SetBinContent(bi, bj, zero);
                }
            }
            for(Int_t bj = 1; bj <= V_PHYS_2D_P2[ch]->GetNbinsY(); ++ bj)
            {
                for(Int_t bi = 1; bi <= V_PHYS_2D_P2[ch]->GetNbinsX(); ++ bi)
                {
                    const Double_t zero = 0.0;
                    V_PHYS_2D_P2[ch]->SetBinContent(bi, bj, zero);
                }
            }
        }
    }
    */
}
#endif





#endif //MINIMIZEFCNAXIALVECTOR_CHECK_ALLOC_H
