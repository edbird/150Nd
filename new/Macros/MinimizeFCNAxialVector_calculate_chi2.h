#ifndef MINIMIZEFCNAXIALVECTOR_CALCULATE_CHI2_H
#define MINIMIZEFCNAXIALVECTOR_CALCULATE_CHI2_H


void
MinimizeFCNAxialVector::calculate_chi2_P1(double &chi2_P1, int &nch_P1) const
{
    chi2_P1 = 0.0;
    nch_P1 = 0;

    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(channel_enable_1D[channel] == true)
        {
            //std::cout << "CALCULATE CHI2 FOR PHASE 1" << std::endl;

            // calculate (M - D) x V_PHYS x Transpose(M - D)
            //std::cout << "starting matrix calculations" << std::endl;
            double chi2_1D_P1 = 0.0;

            int l_counter = 0;
            //for(Int_t l{1}; l <= V_PHYS_STAT_1D_P1[channel]->GetNbinsY(); ++ l)
            for(Int_t l{0}; l < M_1D_P1_data[channel]->size(); ++ l)
            {
                if(V_ENABLE_BIN_1D_P1[channel]->at(l) == true)
                {
                    // do nothing
                }
                else
                {
                    continue;
                }

                int m_counter = 0;
                //for(Int_t m{1}; m <= V_PHYS_STAT_1D_P1[channel]->GetNbinsX(); ++ m)
                for(Int_t m{0}; m < M_1D_P1_data[channel]->size(); ++ m)
                // transpose, should be applied to (M-D)->GetNbinsY() but this object does not exist
                {
                    if(V_ENABLE_BIN_1D_P1[channel]->at(m) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    //double D_content_1 = D_1D_P1[channel]->GetBinContent(l, 1);
                    //double M_content_1 = M_1D_P1[channel]->GetBinContent(l, 1);
                    #if VECTOR_RANGE_CHECK
                    //double D_content_1 = D_1D_P1_data[channel]->at(l);
                    //double M_content_1 = M_1D_P1_data[channel]->at(l);
                    double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->at(l);
                    #else
                    //double D_content_1 = D_1D_P1_data[channel]->operator[](l);
                    //double M_content_1 = M_1D_P1_data[channel]->operator[](l);
                    double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->operator[](l);
                    #endif
                    //double delta_1 = D_content_1 - M_content_1;
                    double delta_1 = D_minus_M_content_1;
                    //double V_CHEN_content = V_CHEN->GetBinContent(k, l);
                    //double V_PHYS_STAT_content = V_PHYS_STAT_1D_P1[channel]->GetBinContent(l, m);

                    double V_PHYS_content = 0.0;
                    //std::cout << "l=" << l << " m=" << m << " l_counter=" << l_counter << " m_counter=" << m_counter << std::endl;
                    V_PHYS_content = V_PHYS_1D_P1_MATHMORE[channel]->operator[](m_counter).operator[](l_counter);
                    
                    //double D_content_2 = D_1D_P1[channel]->GetBinContent(m, 1);
                    //double M_content_2 = M_1D_P1[channel]->GetBinContent(m, 1); 
                    #if VECTOR_RANGE_CHECK
                    //double D_content_2 = D_1D_P1_data[channel]->at(m);
                    //double M_content_2 = M_1D_P1_data[channel]->at(m);
                    double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(m);
                    #else
                    //double D_content_2 = D_1D_P1_data[channel]->operator[](m);
                    //double M_content_2 = M_1D_P1_data[channel]->operator[](m);
                    double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->operator[](m);
                    #endif
                    //double delta_2 = D_content_2 - M_content_2;
                    double delta_2 = D_minus_M_content_2;
                    //double next = delta_1 * V_CHEN_content * V_PHYS_STAT_content * delta_2;
                    double next = delta_1 * V_PHYS_content * delta_2;

                    //if(next != 0.0)
                    //{
                    //    std::cout << "m=" << m << " l=" << l << " content=" << V_PHYS_content << " " << D_minus_M_content_1 << " " << D_minus_M_content_2 << " -> " << next << std::endl;
                    //}

                    if(std::isnan(next))
                    {
                        //std::cout << "NAN: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "NAN: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    else if(std::isinf(next))
                    {
                        //std::cout << "INF: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "INF: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    chi2_1D_P1 += next;
                    
                    ++ m_counter;
                }
                //std::cout << "after adding strip: " << chi2_1D_P1 << std::endl;

                ++ l_counter;
                ++ nch_P1;
            }

            chi2_P1 += chi2_1D_P1;
            //std::cout << "chi2_1D_P1=" << chi2_1D_P1 << std::endl;
        }

    }

    //std::cout << "nch_P1=" << nch_P1 << std::endl;

}



void
MinimizeFCNAxialVector::calculate_chi2_P2(double &chi2_P2, int &nch_P2) const
{
    chi2_P2 = 0.0;
    nch_P2 = 0;

    for(int channel = 0; channel < number1DHists; ++ channel)
    {
        if(channel_enable_1D[channel] == true)
        {
            //std::cout << "CALCULATE CHI2 FOR PHASE 2" << std::endl;

            double chi2_1D_P2 = 0.0;
            
            int l_counter = 0;
            //for(Int_t l{1}; l <= V_PHYS_STAT_1D_P2[channel]->GetNbinsY(); ++ l)
            for(Int_t l{0}; l < 50; ++ l)
            {
                if(V_ENABLE_BIN_1D_P2[channel]->at(l) == true)
                {
                    // do nothing
                }
                else
                {
                    continue;
                }

                int m_counter = 0;
                //for(Int_t m{1}; m <= V_PHYS_STAT_1D_P2[channel]->GetNbinsX(); ++ m)
                for(Int_t m{0}; m < 50; ++ m)
                // transpose, should be applied to (M-D)->GetNbinsY() but this object does not exist
                {
                    if(V_ENABLE_BIN_1D_P2[channel]->at(m) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    //double D_content_1 = D_1D_P2[channel]->GetBinContent(l, 1);
                    //double M_content_1 = M_1D_P2[channel]->GetBinContent(l, 1);
                    #if VECTOR_RANGE_CHECK
                    //double D_content_1 = D_1D_P2_data[channel]->at(l);
                    //double M_content_1 = M_1D_P2_data[channel]->at(l);
                    double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->at(l);
                    #else
                    //double D_content_1 = D_1D_P2_data[channel]->operator[](l);
                    //double M_content_1 = M_1D_P2_data[channel]->operator[](l);
                    double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->operator[](l);
                    #endif
                    //double delta_1 = D_content_1 - M_content_1;
                    double delta_1 = D_minus_M_content_1;
                    //double V_CHEN_content = V_CHEN->GetBinContent(k, l);
                    //double V_PHYS_STAT_content = V_PHYS_STAT_1D_P2[channel]->GetBinContent(l, m);
                    double V_PHYS_content = 0.0;

                    //std::cout << "l=" << l << " m=" << m << " l_counter=" << l_counter << " m_counter=" << m_counter << std::endl;
                    V_PHYS_content = V_PHYS_1D_P2_MATHMORE[channel]->operator[](m_counter).operator[](l_counter);


                    //double D_content_2 = D_1D_P2[channel]->GetBinContent(m, 1);
                    //double M_content_2 = M_1D_P2[channel]->GetBinContent(m, 1);
                    #if VECTOR_RANGE_CHECK
                    //double D_content_2 = D_1D_P2_data[channel]->at(m);
                    //double M_content_2 = M_1D_P2_data[channel]->at(m);
                    double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->at(m);
                    #else
                    //double D_content_2 = D_1D_P2_data[channel]->operator[](m);
                    //double M_content_2 = M_1D_P2_data[channel]->operator[](m);
                    double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->operator[](m);
                    #endif
                    //double delta_2 = D_content_2 - M_content_2;
                    double delta_2 = D_minus_M_content_2;
                    //double next = delta_1 * V_CHEN_content * V_PHYS_STAT_content * delta_2;
                    double next = delta_1 * V_PHYS_content * delta_2;

                    //if(next != 0.0)
                    //{
                    //    std::cout << "m=" << m << " l=" << l << " content=" << V_PHYS_content << " " << D_minus_M_content_1 << " " << D_minus_M_content_2 << " -> " << next << std::endl;
                    //}

                    if(std::isnan(next))
                    {
                        //std::cout << "NAN: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "NAN: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    else if(std::isinf(next))
                    {
                        //std::cout << "INF: next=" << next << " k=" << k << " l=" << l << " m=" << m << std::endl;
                        std::cout << "INF: next=" << next << " l=" << l << " m=" << m << std::endl;
                        //std::cout << V_CHEN_content << std::endl;
                        std::cout << V_PHYS_content << std::endl;
                        std::cin.get();
                    }
                    chi2_1D_P2 += next;
                    
                    ++ m_counter;
                }
                //std::cout << "after adding strip: " << chi2_1D_P2 << std::endl;

                ++ l_counter;
                ++ nch_P2;
            }

            chi2_P2 += chi2_1D_P2;
            //std::cout << "chi2_1D_P2=" << chi2_1D_P2 << std::endl;
        }
    }

    //std::cout << "nch_P2=" << nch_P2 << std::endl;

}

#endif // MINIMIZEFCNAXIALVECTOR_CALCULATE_CHI2_H
