#ifndef MINIMIZEFCNAXIALVECTOR_SET_V_MATRIX_H
#define MINIMIZEFCNAXIALVECTOR_SET_V_MATRIX_H


void
MinimizeFCNAxialVector::set_V_MATRIX() const
{

    // TODO: speed up this code by only building one matrix object
    // V_PHYS, and building it element by element, using symmetry
    // optimizations if possible. this avoids copying several V_MATRIX
    // objects into one final object
    // 
    // then further optimize by storing only D - M and computing the
    // matrix multiplication in place
    // this might not be possible because I need to invert the matrix
    // unless the inversion can be done in place? (might be slower
    // even if I can invent such an algorithm)


    // A funny thing COULD happen here
    // (but it is not very likely)
    // If the MATHMORE matrix becomes singular due to adding the components
    // of all the V_PHYS matricies
    // But the STAT matrix is non-singular
    // then the algorithm will blow up
    //
    // To fix this: Remove all V_PHYS matricies, and enter values
    // straight into a single V_PHYS matrix, or MATHMORE matrix
    // (need to do former to count number of required bins?)

    for(int channel = 0; channel < number1DHists; ++ channel)
    {

        if(channel_enable_1D[channel] == 1)
        {

            if(recalculate_V_PHYS_xD_Px_MATHMORE == true)
            {


                // set V_PHYS_SYS1
                // NOTE: moved to init, because values do not change

                #if 0
                // TODO: symmetry optimization
                for(Int_t biny{0}; biny < M_1D_P1_data[channel]->size(); ++ biny)
                {
                    //for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                    for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                    {
 
                        // P1
                        {
                            #if VECTOR_RANGE_CHECK
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->at(biny);
                            V_PHYS_SYS1_1D_P1_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                            #else
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P1[channel]->operator[](biny);
                            V_PHYS_SYS1_1D_P1_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                            #endif
                        }

                        // P2
                        {
                            #if VECTOR_RANGE_CHECK
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->at(biny);
                            V_PHYS_SYS1_1D_P2_data[channel]->at(biny * 50 + binx) = coeff1 * coeff2;
                            #else
                            double coeff1 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](binx);
                            double coeff2 = systematic_offset_V_MATRIX_coeff_1D_P2[channel]->operator[](biny);
                            V_PHYS_SYS1_1D_P2_data[channel]->operator[](biny * 50 + binx) = coeff1 * coeff2;
                            #endif
                        }


                    }
                }
                #endif


                // set V_PHYS_STAT
                int counter_P1 = 0;
                int counter_P2 = 0;
                //for(Int_t binx{1}; binx <= M_1D_P1[channel]->GetNbinsX(); ++ binx)
                for(Int_t binx{0}; binx < M_1D_P1_data[channel]->size(); ++ binx)
                {

                    // P1
                    {
                        #if VECTOR_RANGE_CHECK
                        Double_t content_M = M_1D_P1_data[channel]->at(binx);
                        Double_t content_D = D_1D_P1_data[channel]->at(binx);
                        #else
                        Double_t content_M = M_1D_P1_data[channel]->operator[](binx);
                        Double_t content_D = D_1D_P1_data[channel]->operator[](binx);
                        #endif
                        //Double_t content_M = M_1D_P1[channel]->GetBinContent(binx, 1);
                        //Double_t content_D = D_1D_P1[channel]->GetBinContent(binx, 1);
                        //if(content == 0.0) continue;
                        //if(content < 0.0) continue;
                        //if(content <= 0.0)
                        //{
                            //V_PHYS_STAT_1D_P1->SetBinContent(binx, binx, 1.0);
                        //}
                        //else
                        //{

                        Double_t stat = 1.0; // default value
                        if(content_M <= 0.0)
                        {
                            V_ENABLE_BIN_1D_P1[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P1[channel]->push_back(true);

                            Double_t sigma_M = std::sqrt(content_M);
                            Double_t sigma_D = std::sqrt(content_D);
                            stat = sigma_M * sigma_M;

                            ++ counter_P1;

                        }
                            #if VECTOR_RANGE_CHECK
                            V_PHYS_STAT_1D_P1_data[channel]->at(binx * 50 + binx) = stat;
                            #else
                            V_PHYS_STAT_1D_P1_data[channel]->operator[](binx * 50 + binx) = stat;
                            #endif

                        //////Double_t sigma_M = std::sqrt(content_M);
                        //////Double_t sigma_D = std::sqrt(content_D);
                        //Double stat = 1.0 / (sigma * sigma);
                        //Double_t stat = 1.0 / std::abs(content);
                        // TODO; re-enable
                        ////Double_t stat = sigma_M * sigma_M; // + sigma_D * sigma_D;
                        //V_PHYS_STAT_1D_P1->SetBinContent(binx, binx, stat);
                        ////if(stat != 0.0) stat = 1.0 / stat;
                        ////else stat = 1.0;
                        //////Double_t stat = sigma_M * sigma_M;

                        /*
                        if(std::isnan(stat))
                        {
                            std::cout << "binx=" << binx << " sigma_M=" << sigma_M << " content_M=" << content_M << std::endl;
                        }
                        */

                        //if(stat == 0.0)
                        //{
                        //    stat = 1.0;
                        //    V_ENABLE_BIN_1D_P1[channel]->push_back(false);
                        //}
                        //else
                        //{
                        //    V_ENABLE_BIN_1D_P1[channel]->push_back(true);
                        //    ++ counter_P1;
                        //}
                        //V_PHYS_STAT_1D_P1[channel]->SetBinContent(binx, binx, stat);
                        //#if VECTOR_RANGE_CHECK
                        //V_PHYS_STAT_1D_P1_data[channel]->at(binx * 50 + binx) = stat;
                        //#else
                        //V_PHYS_STAT_1D_P1_data[channel]->operator[](binx * 50 + binx) = stat;
                        //#endif
                        //}
                    }

                    // P2
                    {
                        //Double_t content_M = M_1D_P2[channel]->GetBinContent(binx, 1);
                        //Double_t content_D = D_1D_P2[channel]->GetBinContent(binx, 1);
                        #if VECTOR_RANGE_CHECK
                        Double_t content_M = M_1D_P2_data[channel]->at(binx);
                        Double_t content_D = D_1D_P2_data[channel]->at(binx);
                        #else
                        Double_t content_M = M_1D_P2_data[channel]->operator[](binx);
                        Double_t content_D = D_1D_P2_data[channel]->operator[](binx);
                        #endif

                        #if 0
                        //if(content == 0.0) continue;
                        //if(content < 0.0) continue;
                        //if(content <= 0.0)
                        //{
                            //V_PHYS_STAT_1D_P2->SetBinContent(binx, binx, 1.0);
                        //}
                        //else
                        //{
                        Double_t sigma_M = std::sqrt(content_M);
                        Double_t sigma_D = std::sqrt(content_D);
                        //Double stat = 1.0 / (sigma * sigma);
                        //Double_t stat = 1.0 / std::abs(content);
                        // TODO: re-enable
                        ////Double_t stat = sigma_M * sigma_M; // + sigma_D * sigma_D;
                        //V_PHYS_STAT_1D_P2->SetBinContent(binx, binx, stat);
                        ////if(stat != 0.0) stat = 1.0 / stat;
                        ////else stat = 1.0;
                        Double_t stat = sigma_M * sigma_M;
                        // TODO: this might not be the best way to detect signularity
                        // may be more sensible to use if D==0 or M==0, and then or sigmaM==0
                        if(stat == 0.0)
                        {
                            stat = 1.0;
                            V_ENABLE_BIN_1D_P2[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(true);
                            ++ counter_P2;
                        }
                        //V_PHYS_STAT_1D_P2[channel]->SetBinContent(binx, binx, stat);
                        #if VECTOR_RANGE_CHECK
                        V_PHYS_STAT_1D_P2_data[channel]->at(binx * 50 + binx) =  stat;
                        #else
                        V_PHYS_STAT_1D_P2_data[channel]->operator[](binx * 50 + binx) =  stat;
                        #endif
                        //}
                        #endif

                        Double_t stat = 1.0;
                        // TODO: this might not be the best way to detect signularity
                        // may be more sensible to use if D==0 or M==0, and then or sigmaM==0
                        if(content_M <= 0.0)
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(true);
                        
                            Double_t sigma_M = std::sqrt(content_M);
                            Double_t sigma_D = std::sqrt(content_D);
                            stat = sigma_M * sigma_M;

                            ++ counter_P2;

                        }
                            #if VECTOR_RANGE_CHECK
                            V_PHYS_STAT_1D_P2_data[channel]->at(binx * 50 + binx) =  stat;
                            #else
                            V_PHYS_STAT_1D_P2_data[channel]->operator[](binx * 50 + binx) =  stat;
                            #endif
                    }
                }
                //std::cout << "P1: number of enabled bins = " << counter_P1 << std::endl;
                //std::cout << "P2: number of enabled bins = " << counter_P2 << std::endl;
                ///////////////////////////////////////////////////////////////////
                // CERN ROOT MathMore Matrix Lib Objects
                ///////////////////////////////////////////////////////////////////


                //if(V_PHYS_1D_P1_MATHMORE[0] == nullptr)
                //if(recalculate_V_PHYS_xD_Px_MATHMORE == true)
                //{
                //std::cout << "Alloc V_PHYS_MATHMORE channel=" << channel << std::endl;

                //const Int_t NUM_BINS_XY = 50;

                //for(int ch = 0; ch < number1DHists; ++ ch)
                //{


                #if MEASURE_FUNCTION_CALL_TIME
                std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();
                #endif
                {
                    if(V_PHYS_1D_P1_MATHMORE[channel] != nullptr)
                    {
                        if(counter_P1 != V_PHYS_1D_P1_MATHMORE[channel]->GetNrows())
                        {
                            // nrows always == ncols
                            delete V_PHYS_1D_P1_MATHMORE[channel];
                            V_PHYS_1D_P1_MATHMORE[channel] = new TMatrixD(counter_P1, counter_P1); // NUM_BINS_XY, NUM_BINS_XY);
                        }
                    }
                    else
                    {
                    std::cout << "MATHMORE: P1: " << V_PHYS_1D_P1_MATHMORE[0] << std::endl;
                        V_PHYS_1D_P1_MATHMORE[channel] = new TMatrixD(counter_P1, counter_P1); // NUM_BINS_XY, NUM_BINS_XY);
                        std::cout << "ALLOC " << __func__ << std::endl;
                    }

                    if(V_PHYS_1D_P2_MATHMORE[channel] != nullptr)
                    {
                        if(counter_P2 != V_PHYS_1D_P2_MATHMORE[channel]->GetNrows())
                        {
                            delete V_PHYS_1D_P2_MATHMORE[channel];
                            V_PHYS_1D_P2_MATHMORE[channel] = new TMatrixD(counter_P2, counter_P2); //NUM_BINS_XY, NUM_BINS_XY);
                        }
                    }
                    else
                    {
                    std::cout << "MATHMORE: P2: " << V_PHYS_1D_P2_MATHMORE[0] << std::endl;
                        V_PHYS_1D_P2_MATHMORE[channel] = new TMatrixD(counter_P2, counter_P2); //NUM_BINS_XY, NUM_BINS_XY);
                        std::cout << "ALLOC " << __func__ << std::endl;
                    }
                    // TODO: must be able to optimize this
                }
                #if MEASURE_FUNCTION_CALL_TIME
                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                std::cout << "V_PHYS_1D_Px_MATHMORE allocation time: " << 1.0e+06 * runtime_sec.count() << " microsecond" << std::endl;
                #endif

                //std::cout << "size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                //std::cout << "size: " << V_PHYS_1D_P2_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P2_MATHMORE[channel]->GetNcols() << std::endl;
                //std::cin.get();
                //}

                /*
                for(int ch = 0; ch < number1DHists; ++ ch)
                {
                    for(Int_t iy = 0; iy < V_PHYS_1D_P1_MATHMORE[ch]->GetNrows(); ++ iy)
                    {
                        for(Int_t ix = 0; ix < V_PHYS_1D_P1_MATHMORE[ch]->GetNcols(); ++ ix)
                        {
                            V_PHYS_1D_P1_MATHMORE[ch]->operator[](iy).operator[](ix) = 0.0;
                            V_PHYS_1D_P2_MATHMORE[ch]->operator[](iy).operator[](ix) = 0.0;
                        }
                    }
                }
                */
                //}


                // TODO: fix problems here
                // if channel is disabled, ignore it
                // should be done in a big output loop?
                // DONE
                // there is a reference to channel here
                // and yet there is a loop over ch above
                // how can this make sense? is it a bug?
                // FIXED
                // TODO: ALL CHANNELS ARE CURRENTLY ENABLED?
   
                #define DRAWVPHYSMATRIX 1

                #if DRAWVPHYSMATRIX
                TH2D *draw_V_PHYS_STAT_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS1_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS2_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS3_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS4_P1 = nullptr;
                TH2D *draw_V_PHYS_P1 = nullptr;
                TH2D *draw_V_PHYSINV_P1 = nullptr;
                TH2D *draw_V_PHYS_STAT_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS1_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS2_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS3_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS4_P2 = nullptr;
                TH2D *draw_V_PHYS_P2 = nullptr;
                TH2D *draw_V_PHYSINV_P2 = nullptr;
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        draw_V_PHYS_STAT_P1 = new TH2D("draw_V_PHYS_STAT_P1", "draw_V_PHYS_STAT_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS1_P1 = new TH2D("draw_V_PHYS_SYS1_P1", "draw_V_PHYS_SYS1_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS2_P1 = new TH2D("draw_V_PHYS_SYS2_P1", "draw_V_PHYS_SYS2_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS3_P1 = new TH2D("draw_V_PHYS_SYS3_P1", "draw_V_PHYS_SYS3_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS4_P1 = new TH2D("draw_V_PHYS_SYS4_P1", "draw_V_PHYS_SYS4_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_P1 =      new TH2D("draw_V_PHYS_P1", "draw_V_PHYS_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYSINV_P1 =   new TH2D("draw_V_PHYSINV_P1", "draw_V_PHYSINV_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_STAT_P2 = new TH2D("draw_V_PHYS_STAT_P2", "draw_V_PHYS_STAT_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS1_P2 = new TH2D("draw_V_PHYS_SYS1_P2", "draw_V_PHYS_SYS1_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS2_P2 = new TH2D("draw_V_PHYS_SYS2_P2", "draw_V_PHYS_SYS2_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS3_P2 = new TH2D("draw_V_PHYS_SYS3_P2", "draw_V_PHYS_SYS3_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYS4_P2 = new TH2D("draw_V_PHYS_SYS4_P2", "draw_V_PHYS_SYS4_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_P2 =      new TH2D("draw_V_PHYS_P2", "draw_V_PHYS_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYSINV_P2 =   new TH2D("draw_V_PHYSINV_P2", "draw_V_PHYSINV_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        std::cout << "ALLOC " << __func__ << std::endl;
                    }
                }
                #endif


                ///////////////////////////////////////////////////////////
                // Phase 1

                int j_counter = 0;
                //for(Int_t j = 0; j < V_PHYS_STAT_1D_P1[channel]->GetNbinsY(); ++ j)
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P1[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P1[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        //Double_t content = V_PHYS_STAT_1D_P1[channel]->GetBinContent(i + 1, j + 1);
                        #if VECTOR_RANGE_CHECK
                        //Double_t content = V_PHYS_STAT_1D_P1_data[channel]->at(i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P1_data[channel]->at(i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P1_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P1_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P1_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P1_data[channel]->at(i + j * 50);
                            }
                        }
                        Double_t content = c1 + c2 + c3 + c4 + c5;
                        #else
                        //Double_t content = V_PHYS_STAT_1D_P1_data[channel]->operator[](i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P1_data[channel]->operator[](i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                        }
                        Double_t content = c1 + c2 + c3 + c4 + c5;
                        #endif
                        //std::cout << "i=" << i << " j=" << j << " i_counter=" << i_counter << " j_counter=" << j_counter << " content=" << content << std::endl;
                        V_PHYS_1D_P1_MATHMORE[channel]->operator[](j_counter).operator[](i_counter) = content;
                        //std::cout << "j=" << j << " i=" << i << " " << content << std::endl;

                        /*
                        if(channel == 1)
                        {
                            std::cout << "P1" << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c4=" << c4 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                            std::cin.get();
                        }
                        */

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYS_STAT_P1->SetBinContent(i + 1, j + 1, c1);
                                draw_V_PHYS_SYS1_P1->SetBinContent(i + 1, j + 1, c2);
                                draw_V_PHYS_SYS2_P1->SetBinContent(i + 1, j + 1, c3);
                                draw_V_PHYS_SYS3_P1->SetBinContent(i + 1, j + 1, c4);
                                draw_V_PHYS_SYS4_P1->SetBinContent(i + 1, j + 1, c5);
                                draw_V_PHYS_P1->SetBinContent(i + 1, j + 1, content);

                            /*
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                                //std::cout << "channel=" << channel << " i=" << i << " j=" << j << " content=" << content << std::endl;

                                std::cout << "V_ENABLE_SYSALL=" << V_ENABLE_SYSALL << std::endl;
                                std::cout << "V_ENABLE_SYS1=" << V_ENABLE_SYS1 << std::endl;
                                std::cout << "V_ENABLE_SYS2=" << V_ENABLE_SYS2 << std::endl;
                                std::cout << "V_ENABLE_SYS3=" << V_ENABLE_SYS3 << std::endl;
                                std::cout << "V_ENABLE_SYS4=" << V_ENABLE_SYS4 << std::endl;
                                //std::cin.get();
                            */
                            }
                        }
                        #endif


                        ++ i_counter;
                    }

                    ++ j_counter;
                }

                ///////////////////////////////////////////////////////////
                // Phase 2

                j_counter = 0;
                //for(Int_t j = 0; j < V_PHYS_STAT_1D_P2[channel]->GetNbinsY(); ++ j)
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P2[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P2[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        //Double_t content = V_PHYS_STAT_1D_P2[channel]->GetBinContent(i + 1, j + 1);
                        #if VECTOR_RANGE_CHECK
                        //Double_t content = V_PHYS_STAT_1D_P2_data[channel]->at(i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P2_data[channel]->at(i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P2_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P2_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P2_data[channel]->at(i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P2_data[channel]->at(i + j * 50);
                            }
                        }
                        double content = c1 + c2 + c3 + c4 + c5;
                        #else
                        //Double_t content = V_PHYS_STAT_1D_P2_data[channel]->operator[](i + j * 50);
                        double c1 = 0.0;
                        double c2 = 0.0;
                        double c3 = 0.0;
                        double c4 = 0.0;
                        double c5 = 0.0;
                        if(V_ENABLE_STAT == true)
                        {
                            c1 = V_PHYS_STAT_1D_P2_data[channel]->operator[](i + j * 50);
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            if(V_ENABLE_SYS1 == true)
                            {
                                c2 = V_PHYS_SYS1_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS2 == true)
                            {
                                c3 = V_PHYS_SYS2_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS3 == true)
                            {
                                c4 = V_PHYS_SYS3_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                            if(V_ENABLE_SYS4 == true)
                            {
                                c5 = V_PHYS_SYS4_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                        }
                        double content = c1 + c2 + c3 + c4 + c5;
                        #endif
                        //std::cout << "i=" << i << " j=" << j << " i_counter=" << i_counter << " j_counter=" << j_counter << " content=" << content << std::endl;
                        V_PHYS_1D_P2_MATHMORE[channel]->operator[](j_counter).operator[](i_counter) = content;

                        /*
                        if(channel == 1)
                        {
                            std::cout << "P2" << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c4=" << c4 << std::endl;
                            std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                            std::cin.get();
                        }
                        */

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYS_STAT_P2->SetBinContent(i + 1, j + 1, c1);
                                draw_V_PHYS_SYS1_P2->SetBinContent(i + 1, j + 1, c2);
                                draw_V_PHYS_SYS2_P2->SetBinContent(i + 1, j + 1, c3);
                                draw_V_PHYS_SYS3_P2->SetBinContent(i + 1, j + 1, c4);
                                draw_V_PHYS_SYS4_P2->SetBinContent(i + 1, j + 1, c5);
                                draw_V_PHYS_P2->SetBinContent(i + 1, j + 1, content);

                            /*
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c1=" << c1 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c2=" << c2 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c3=" << c3 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c4=" << c4 << std::endl;
                                std::cout << "channel=" << channel << " i=" << i << " j=" << j << " c5=" << c5 << std::endl;
                                //std::cout << "channel=" << channel << " i=" << i << " j=" << j << " content=" << content << std::endl;
                            */
                            }
                        }
                        #endif


                        ++ i_counter;
                    }

                    ++ j_counter;
                }

                //std::cout << "AT END OF SET V_MATRIX_?_MATHMORE" << std::endl;
                //std::cin.get();

                #if DRAWVPHYSMATRIX
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        TCanvas *canvas_V_PHYS_STAT_P1 = new TCanvas("canvas_V_PHYS_STAT_P1", "canvas_V_PHYS_STAT_P1");
                        draw_V_PHYS_STAT_P1->Draw("colz");
                        canvas_V_PHYS_STAT_P1->SaveAs("draw_V_PHYS_STAT_P1.png");

                        TCanvas *canvas_V_PHYS_SYS1_P1 = new TCanvas("canvas_V_PHYS_SYS1_P1", "canvas_V_PHYS_SYS1_P1");
                        draw_V_PHYS_SYS1_P1->Draw("colz");
                        canvas_V_PHYS_SYS1_P1->SaveAs("draw_V_PHYS_SYS1_P1.png");

                        TCanvas *canvas_V_PHYS_SYS2_P1 = new TCanvas("canvas_V_PHYS_SYS2_P1", "canvas_V_PHYS_SYS2_P1");
                        draw_V_PHYS_SYS2_P1->Draw("colz");
                        canvas_V_PHYS_SYS2_P1->SaveAs("draw_V_PHYS_SYS2_P1.png");

                        TCanvas *canvas_V_PHYS_SYS3_P1 = new TCanvas("canvas_V_PHYS_SYS3_P1", "canvas_V_PHYS_SYS3_P1");
                        draw_V_PHYS_SYS3_P1->Draw("colz");
                        canvas_V_PHYS_SYS3_P1->SaveAs("draw_V_PHYS_SYS3_P1.png");

                        TCanvas *canvas_V_PHYS_SYS4_P1 = new TCanvas("canvas_V_PHYS_SYS4_P1", "canvas_V_PHYS_SYS4_P1");
                        draw_V_PHYS_SYS4_P1->Draw("colz");
                        canvas_V_PHYS_SYS4_P1->SaveAs("draw_V_PHYS_SYS4_P1.png");

                        TCanvas *canvas_V_PHYS_P1 = new TCanvas("canvas_V_PHYS_P1", "canvas_V_PHYS_P1");
                        draw_V_PHYS_P1->Draw("colz");
                        canvas_V_PHYS_P1->SaveAs("draw_V_PHYS_P1.png");

                        TCanvas *canvas_V_PHYS_STAT_P2 = new TCanvas("canvas_V_PHYS_STAT_P2", "canvas_V_PHYS_STAT_P2");
                        draw_V_PHYS_STAT_P2->Draw("colz");
                        canvas_V_PHYS_STAT_P2->SaveAs("draw_V_PHYS_STAT_P2.png");

                        TCanvas *canvas_V_PHYS_SYS1_P2 = new TCanvas("canvas_V_PHYS_SYS1_P2", "canvas_V_PHYS_SYS1_P2");
                        draw_V_PHYS_SYS1_P2->Draw("colz");
                        canvas_V_PHYS_SYS1_P2->SaveAs("draw_V_PHYS_SYS1_P2.png");

                        TCanvas *canvas_V_PHYS_SYS2_P2 = new TCanvas("canvas_V_PHYS_SYS2_P2", "canvas_V_PHYS_SYS2_P2");
                        draw_V_PHYS_SYS2_P2->Draw("colz");
                        canvas_V_PHYS_SYS2_P2->SaveAs("draw_V_PHYS_SYS2_P2.png");

                        TCanvas *canvas_V_PHYS_SYS3_P2 = new TCanvas("canvas_V_PHYS_SYS3_P2", "canvas_V_PHYS_SYS3_P2");
                        draw_V_PHYS_SYS3_P2->Draw("colz");
                        canvas_V_PHYS_SYS3_P2->SaveAs("draw_V_PHYS_SYS3_P2.png");

                        TCanvas *canvas_V_PHYS_SYS4_P2 = new TCanvas("canvas_V_PHYS_SYS4_P2", "canvas_V_PHYS_SYS4_P2");
                        draw_V_PHYS_SYS4_P2->Draw("colz");
                        canvas_V_PHYS_SYS4_P2->SaveAs("draw_V_PHYS_SYS4_P2.png");

                        TCanvas *canvas_V_PHYS_P2 = new TCanvas("canvas_V_PHYS_P2", "canvas_V_PHYS_P2");
                        draw_V_PHYS_P2->Draw("colz");
                        canvas_V_PHYS_P2->SaveAs("draw_V_PHYS_P2.png");
                        std::cout << "ALLOC " << __func__ << std::endl;

                        //std::cout << "printed all canvas" << std::endl;
                        //std::cin.get();
                    }
                }
                #endif

                /*
                for(Int_t iy = 0; iy < V_PHYS_1D_P1_MATHMORE[channel]->GetNrows(); ++ iy)
                {
                    for(Int_t ix = 0; ix < V_PHYS_1D_P1_MATHMORE[channel]->GetNcols(); ++ ix)
                    {
                        Double_t vP1 = V_PHYS_1D_P1_MATHMORE[channel]->operator[](iy).operator[](ix);
                        if(vP1 != 0.0)
                        {
                            std::cout << "P1 ix=" << ix << " iy=" << iy << " " << vP1 << std::endl;
                        }
                    }
                }
                for(Int_t iy = 0; iy < V_PHYS_1D_P2_MATHMORE[channel]->GetNrows(); ++ iy)
                {
                    for(Int_t ix = 0; ix < V_PHYS_1D_P2_MATHMORE[channel]->GetNcols(); ++ ix)
                    {
                        Double_t vP2 = V_PHYS_1D_P2_MATHMORE[channel]->operator[](iy).operator[](ix);
                        if(vP2 != 0.0)
                        {
                            std::cout << "P2 ix=" << ix << " iy=" << iy << " " << vP2 << std::endl;
                        }
                    }
                }
                */
                //std::cout << "Matrix size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " x " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                //std::cout << "Matrix size: " << V_PHYS_1D_P2_MATHMORE[channel]->GetNrows() << " x " << V_PHYS_1D_P2_MATHMORE[channel]->GetNcols() << std::endl;
                #define MEASURE_FUNCTION_CALL_TIME_MATRIX_INVERT 0
                #if MEASURE_FUNCTION_CALL_TIME_MATRIX_INVERT
                std::cout << "Start Invert" << std::endl;
                std::chrono::system_clock::time_point start_time = std::chrono::high_resolution_clock::now();
                #endif
                //std::cout << "P1 V_MATRIX size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                if(V_PHYS_1D_P1_MATHMORE[channel]->Determinant() == 0.0)
                {
                    std::cout << "P1 V_MATRIX size: " << V_PHYS_1D_P1_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P1_MATHMORE[channel]->GetNcols() << std::endl;
                    for(int j = 0; j < V_PHYS_1D_P1_MATHMORE[channel]->GetNcols(); ++ j)
                    {
                        for(int i = 0; i < V_PHYS_1D_P1_MATHMORE[channel]->GetNrows(); ++ i)
                        {
                            std::cout << "i=" << i << " j=" << j << " " << V_PHYS_1D_P1_MATHMORE[channel]->operator[](j).operator[](i) << std::endl;
                        }
                    }
                }

                //TMatrixD copy(*V_PHYS_1D_P1_MATHMORE[channel]);
                V_PHYS_1D_P1_MATHMORE[channel]->Invert();
                //TMatrixD multiplymat = copy * *V_PHYS_1D_P1_MATHMORE[channel];
                //std::cout << "Elements after mult" << std::endl;
                //    for(int j = 0; j < V_PHYS_1D_P1_MATHMORE[channel]->GetNcols(); ++ j)
                //    {
                //        for(int i = 0; i < V_PHYS_1D_P1_MATHMORE[channel]->GetNrows(); ++ i)
                //        {
                //            std::cout << "i=" << i << " j=" << j << " " << multiplymat.operator[](j).operator[](i) << std::endl;
                //        }
                //    }
                //    std::cin.get();
                //std::cout << "Next Invert" << std::endl;
                //std::cout << "P2 V_MATRIX size: " << V_PHYS_1D_P2_MATHMORE[channel]->GetNrows() << " " << V_PHYS_1D_P2_MATHMORE[channel]->GetNcols() << std::endl;
                V_PHYS_1D_P2_MATHMORE[channel]->Invert();
                #if MEASURE_FUNCTION_CALL_TIME_MATRIX_INVERT
                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                std::cout << "Done Invert, time=" << 1.0e+06 * runtime_sec.count() << " microsecond" << std::endl;
                #endif




                j_counter = 0;
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P1[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P1[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        double content = V_PHYS_1D_P1_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYSINV_P1->SetBinContent(i + 1, j + 1, content);
                                //std::cout << content << "\t";
                            }
                        }
                        #endif
                        ++ i_counter;
                    }
                    #if DRAWVPHYSMATRIX
                    if(DRAW_V_PHYS_MATRIX == true)
                    {
                        //std::cout << std::endl;
                    }
                    #endif

                    ++ j_counter;
                }


                j_counter = 0;
                for(Int_t j = 0; j < 50; ++ j)
                {
                    if(V_ENABLE_BIN_1D_P2[channel]->at(j) == true)
                    {
                        // do nothing
                    }
                    else
                    {
                        continue;
                    }

                    int i_counter = 0;
                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        if(V_ENABLE_BIN_1D_P2[channel]->at(i) == true)
                        {
                            // do nothing
                        }
                        else
                        {
                            continue;
                        }

                        double content = V_PHYS_1D_P2_MATHMORE[channel]->operator[](j_counter).operator[](i_counter);

                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYSINV_P2->SetBinContent(i + 1, j + 1, content);
                                //std::cout << content << "\t";
                            }
                        }
                        #endif
                        ++ i_counter;
                    }
                        #if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            //std::cout << std::endl;
                        }
                        #endif

                    ++ j_counter;
                }


                #if DRAWVPHYSMATRIX
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        std::cout << "ALLOC " << __func__ << std::endl;
                        TCanvas *canvas_V_PHYSINV_P1 = new TCanvas("canvas_V_PHYSINV_P1", "canvas_V_PHYSINV_P1");
                        draw_V_PHYSINV_P1->Draw("colz");
                        canvas_V_PHYSINV_P1->SaveAs("draw_V_PHYSINV_P1.png");

                        TCanvas *canvas_V_PHYSINV_P2 = new TCanvas("canvas_V_PHYSINV_P2", "canvas_V_PHYSINV_P2");
                        draw_V_PHYSINV_P2->Draw("colz");
                        canvas_V_PHYSINV_P2->SaveAs("draw_V_PHYSINV_P2.png");

                        std::cout << "printed all canvas" << std::endl;
                        //std::cin.get();
                        
                        DRAW_V_PHYS_MATRIX = false;
                    }
                }
                #endif
                
                // disable this to recalculate the V MATRIX each loop
                // if the number of singular elements changes it may crash
                // not sure at the moment
                //recalculate_V_PHYS_xD_Px_MATHMORE = false;
            }
        }
    }
}

#endif // MINIMIZEFCNAXIALVECTOR_SET_V_MATRIX_H
