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

            if(recalculate_V_PHYS_MATHMORE == true)
            {


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

                        Double_t stat = 1.0; // default value
                        if(content_M <= 0.0)
                        {
                            V_ENABLE_BIN_1D_P1[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P1[channel]->push_back(true);

                            //Double_t sigma_M = std::sqrt(content_M);
                            //Double_t sigma_D = std::sqrt(content_D);
                            //stat = sigma_M * sigma_M;
                            stat = content_M;

                            ++ counter_P1;

                        }

                        #if VECTOR_RANGE_CHECK
                        V_PHYS_STAT_1D_P1_data[channel]->at(binx * 50 + binx) = stat;
                        #else
                        V_PHYS_STAT_1D_P1_data[channel]->operator[](binx * 50 + binx) = stat;
                        #endif

                        // TODO; re-enable
                        ////Double_t stat = sigma_M * sigma_M; // + sigma_D * sigma_D;
                        //V_PHYS_STAT_1D_P1->SetBinContent(binx, binx, stat);
                    }

                    // P2
                    {
                        #if VECTOR_RANGE_CHECK
                        Double_t content_M = M_1D_P2_data[channel]->at(binx);
                        Double_t content_D = D_1D_P2_data[channel]->at(binx);
                        #else
                        Double_t content_M = M_1D_P2_data[channel]->operator[](binx);
                        Double_t content_D = D_1D_P2_data[channel]->operator[](binx);
                        #endif

                        // TODO: this might not be the best way to detect signularity
                        // may be more sensible to use if D==0 or M==0, and then or sigmaM==0
                        Double_t stat = 1.0;
                        if(content_M <= 0.0)
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(false);
                        }
                        else
                        {
                            V_ENABLE_BIN_1D_P2[channel]->push_back(true);
                        
                            //Double_t sigma_M = std::sqrt(content_M);
                            //Double_t sigma_D = std::sqrt(content_D);
                            //stat = sigma_M * sigma_M;
                            stat = content_M;

                            ++ counter_P2;

                        }

                        #if VECTOR_RANGE_CHECK
                        V_PHYS_STAT_1D_P2_data[channel]->at(binx * 50 + binx) = stat;
                        #else
                        V_PHYS_STAT_1D_P2_data[channel]->operator[](binx * 50 + binx) = stat;
                        #endif
                    }
                }
                //std::cout << "P1: number of enabled bins = " << counter_P1 << std::endl;
                //std::cout << "P2: number of enabled bins = " << counter_P2 << std::endl;
                ///////////////////////////////////////////////////////////////////
                // CERN ROOT MathMore Matrix Lib Objects
                ///////////////////////////////////////////////////////////////////


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
                        std::cout << "ALLOC " << __func__ << std::endl;
                        V_PHYS_1D_P1_MATHMORE[channel] = new TMatrixD(counter_P1, counter_P1); // NUM_BINS_XY, NUM_BINS_XY);
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
                        std::cout << "ALLOC " << __func__ << std::endl;
                        V_PHYS_1D_P2_MATHMORE[channel] = new TMatrixD(counter_P2, counter_P2); //NUM_BINS_XY, NUM_BINS_XY);
                    }
                    // TODO: must be able to optimize this
                }
                #if MEASURE_FUNCTION_CALL_TIME
                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                std::cout << "V_PHYS_1D_Px_MATHMORE allocation time: " << 1.0e+06 * runtime_sec.count() << " microsecond" << std::endl;
                #endif

   
                #define DRAWVPHYSMATRIX 1
                #if DRAWVPHYSMATRIX
                TH2D *draw_V_PHYS_STAT_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS1_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS2_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS3_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS4_P1 = nullptr;
                TH2D *draw_V_PHYS_SYS5_P1 = nullptr;
                TH2D *draw_V_PHYS_SYSALL_P1 = nullptr;
                /*
                TH2D *draw_V_PHYS_DmMSYSALL_P1 = nullptr;
                TH2D *draw_V_PHYS_SYSALLDmM_P1 = nullptr;
                TH2D *draw_V_PHYS_DmMSYSALLDmM_P1 = nullptr;*/
                TH2D *draw_V_PHYS_P1 = nullptr;
                TH2D *draw_V_PHYSINV_P1 = nullptr;
                TH2D *draw_V_DmMPHYSINVDmM_P1 = nullptr;
                TH2D *draw_V_PHYS_STAT_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS1_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS2_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS3_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS4_P2 = nullptr;
                TH2D *draw_V_PHYS_SYS5_P2 = nullptr;
                TH2D *draw_V_PHYS_SYSALL_P2 = nullptr;
                /*TH2D *draw_V_PHYS_DmMSYSALL_P2 = nullptr;
                TH2D *draw_V_PHYS_SYSALLDmM_P2 = nullptr;
                TH2D *draw_V_PHYS_DmMSYSALLDmM_P2 = nullptr;*/
                TH2D *draw_V_PHYS_P2 = nullptr;
                TH2D *draw_V_PHYSINV_P2 = nullptr;
                TH2D *draw_V_DmMPHYSINVDmM_P2 = nullptr;
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
                        draw_V_PHYS_SYS5_P1 = new TH2D("draw_V_PHYS_SYS5_P1", "draw_V_PHYS_SYS5_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYSALL_P1 = new TH2D("draw_V_PHYS_SYSALL_P1", "draw_V_PHYS_SYSALL_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        /*
                        draw_V_PHYS_DmMSYSALL_P1 = new TH2D("draw_V_PHYS_DmMSYSALL_P1", "draw_V_PHYS_DmMSYSALL_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYSALLDmM_P1 = new TH2D("draw_V_PHYS_SYSALLDmM_P1", "draw_V_PHYS_SYSALLDmM_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_DmMSYSALLDmM_P1 = new TH2D("draw_V_PHYS_DmMSYSALLDmM_P1", "draw_V_PHYS_DmMSYSALLDmM_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        */
                        draw_V_PHYS_P1 =      new TH2D("draw_V_PHYS_P1", "draw_V_PHYS_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYSINV_P1 =   new TH2D("draw_V_PHYSINV_P1", "draw_V_PHYSINV_P1",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_DmMPHYSINVDmM_P1 =   new TH2D("draw_V_DmMPHYSINVDmM_P1", "draw_V_DmMPHYSINVDmM_P1",
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
                        draw_V_PHYS_SYS5_P2 = new TH2D("draw_V_PHYS_SYS5_P2", "draw_V_PHYS_SYS5_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYSALL_P2 = new TH2D("draw_V_PHYS_SYSALL_P2", "draw_V_PHYS_SYSALL_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        /*
                        draw_V_PHYS_DmMSYSALL_P2 = new TH2D("draw_V_PHYS_DmMSYSALL_P2", "draw_V_PHYS_DmMSYSALL_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_SYSALLDmM_P2 = new TH2D("draw_V_PHYS_SYSALLDmM_P2", "draw_V_PHYS_SYSALLDmM_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYS_DmMSYSALLDmM_P2 = new TH2D("draw_V_PHYS_DmMSYSALLDmM_P2", "draw_V_PHYS_DmMSYSALLDmM_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        */
                        draw_V_PHYS_P2 =      new TH2D("draw_V_PHYS_P2", "draw_V_PHYS_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_PHYSINV_P2 =   new TH2D("draw_V_PHYSINV_P2", "draw_V_PHYSINV_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        draw_V_DmMPHYSINVDmM_P2 =   new TH2D("draw_V_DmMPHYSINVDmM_P2", "draw_V_DmMPHYSINVDmM_P2",
                                                             50, 0.0, 5.0,
                                                             50, 0.0, 5.0);
                        std::cout << "ALLOC " << __func__ << std::endl;


                        draw_V_PHYS_STAT_P1->SetContour(1000);
                        draw_V_PHYS_SYS1_P1->SetContour(1000);
                        draw_V_PHYS_SYS2_P1->SetContour(1000);
                        draw_V_PHYS_SYS3_P1->SetContour(1000);
                        draw_V_PHYS_SYS4_P1->SetContour(1000);
                        draw_V_PHYS_SYS5_P1->SetContour(1000);
                        draw_V_PHYS_SYSALL_P1->SetContour(1000);
                        /*draw_V_PHYS_DmMSYSALL_P1->SetContour(1000);
                        draw_V_PHYS_SYSALLDmM_P1->SetContour(1000);
                        draw_V_PHYS_DmMSYSALLDmM_P1->SetContour(1000);*/
                        draw_V_PHYS_P1->SetContour(1000);
                        draw_V_PHYSINV_P1->SetContour(1000);
                        draw_V_DmMPHYSINVDmM_P1->SetContour(1000);
                        draw_V_PHYS_STAT_P2->SetContour(1000);
                        draw_V_PHYS_SYS1_P2->SetContour(1000);
                        draw_V_PHYS_SYS2_P2->SetContour(1000);
                        draw_V_PHYS_SYS3_P2->SetContour(1000);
                        draw_V_PHYS_SYS4_P2->SetContour(1000);
                        draw_V_PHYS_SYS5_P2->SetContour(1000);
                        draw_V_PHYS_SYSALL_P2->SetContour(1000);
                        /*draw_V_PHYS_DmMSYSALL_P2->SetContour(1000);
                        draw_V_PHYS_SYSALLDmM_P2->SetContour(1000);
                        draw_V_PHYS_DmMSYSALLDmM_P2->SetContour(1000);*/
                        draw_V_PHYS_P2->SetContour(1000);
                        draw_V_PHYSINV_P2->SetContour(1000);
                        draw_V_DmMPHYSINVDmM_P2->SetContour(1000);
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

                        double cstat = 0.0;
                        double csys = 0.0;
                        double content = 0.0;
                        #if VECTOR_RANGE_CHECK
                        if(V_ENABLE_STAT == true)
                        {
                            // TODO: move optimization to new loop
                            if(i == j)
                            {
                                cstat = V_PHYS_STAT_1D_P1_data[channel]->at(i + j * 50);
                            }
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            csys = V_PHYS_SYSALL_1D_P1_data[channel]->at(i + j * 50);
                        }
                        #else
                        if(V_ENABLE_STAT == true)
                        {
                            // TODO: move optimization to new loop
                            if(i == j)
                            {
                                cstat = V_PHYS_STAT_1D_P1_data[channel]->operator[](i + j * 50);
                            }
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            csys = V_PHYS_SYSALL_1D_P1_data[channel]->operator[](i + j * 50);
                        }
                        #endif
                        content = cstat + csys;
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

                                double sys1 = 0.0;
                                double sys2 = 0.0;
                                double sys3 = 0.0;
                                double sys4 = 0.0;
                                double sys5 = 0.0;

                                //if(V_ENABLE_SYS1 == true)
                                //{
                                    sys1 = V_PHYS_SYS1_1D_P1_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS2 == true)
                                //{
                                    sys2 = V_PHYS_SYS2_1D_P1_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS3 == true)
                                //{
                                    sys3 = V_PHYS_SYS3_1D_P1_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS4 == true)
                                //{
                                    sys4 = V_PHYS_SYS4_1D_P1_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS5 == true)
                                //{
                                    sys5 = V_PHYS_SYS5_1D_P1_data[channel]->operator[](i + j * 50);
                                //}

                                draw_V_PHYS_STAT_P1->SetBinContent(i + 1, j + 1, cstat);
                                draw_V_PHYS_SYS1_P1->SetBinContent(i + 1, j + 1, sys1);
                                draw_V_PHYS_SYS2_P1->SetBinContent(i + 1, j + 1, sys2);
                                draw_V_PHYS_SYS3_P1->SetBinContent(i + 1, j + 1, sys3);
                                draw_V_PHYS_SYS4_P1->SetBinContent(i + 1, j + 1, sys4);
                                draw_V_PHYS_SYS5_P1->SetBinContent(i + 1, j + 1, sys5);
                                draw_V_PHYS_SYSALL_P1->SetBinContent(i + 1, j + 1, csys);
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

                        double cstat = 0.0;
                        double csys = 0.0;
                        double content = 0.0;
                        #if VECTOR_RANGE_CHECK
                        if(V_ENABLE_STAT == true)
                        {
                            // TODO: move optimization to new loop
                            if(i == j)
                            {
                                cstat = V_PHYS_STAT_1D_P2_data[channel]->at(i + j * 50);
                            }
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            csys = V_PHYS_SYSALL_1D_P2_data[channel]->at(i + j * 50);
                        }
                        #else
                        if(V_ENABLE_STAT == true)
                        {
                            // TODO: move optimization to new loop
                            if(i == j)
                            {
                                cstat = V_PHYS_STAT_1D_P2_data[channel]->operator[](i + j * 50);
                            }
                        }
                        if(V_ENABLE_SYSALL == true)
                        {
                            csys = V_PHYS_SYSALL_1D_P2_data[channel]->operator[](i + j * 50);
                        }
                        #endif
                        content = cstat + csys;
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
                                double sys1 = 0.0;
                                double sys2 = 0.0;
                                double sys3 = 0.0;
                                double sys4 = 0.0;
                                double sys5 = 0.0;

                                // fill figures regardless of whether systematic is enabled
                                //if(V_ENABLE_SYS1 == true)
                                //{
                                    sys1 = V_PHYS_SYS1_1D_P2_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS2 == true)
                                //{
                                    sys2 = V_PHYS_SYS2_1D_P2_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS3 == true)
                                //{
                                    sys3 = V_PHYS_SYS3_1D_P2_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS4 == true)
                                //{
                                    sys4 = V_PHYS_SYS4_1D_P2_data[channel]->operator[](i + j * 50);
                                //}
                                //if(V_ENABLE_SYS5 == true)
                                //{
                                    sys5 = V_PHYS_SYS5_1D_P2_data[channel]->operator[](i + j * 50);
                                //}

                                draw_V_PHYS_STAT_P2->SetBinContent(i + 1, j + 1, cstat);
                                draw_V_PHYS_SYS1_P2->SetBinContent(i + 1, j + 1, sys1);
                                draw_V_PHYS_SYS2_P2->SetBinContent(i + 1, j + 1, sys2);
                                draw_V_PHYS_SYS3_P2->SetBinContent(i + 1, j + 1, sys3);
                                draw_V_PHYS_SYS4_P2->SetBinContent(i + 1, j + 1, sys4);
                                draw_V_PHYS_SYS5_P2->SetBinContent(i + 1, j + 1, sys5);
                                draw_V_PHYS_SYSALL_P2->SetBinContent(i + 1, j + 1, csys);
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
                    /*
                    for(Int_t j = 0; j < 50; ++ j)
                    {
                        for(Int_t i = 0; i < 50; ++ i)
                        {
                            {
                                #if VECTOR_RANGE_CHECK
                                double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->at(i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P1_data[channel]->at(i + j * 50);
                                #else
                                double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->operator[](i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P1_data[channel]->operator[](i + j * 50);
                                #endif
                                double content = D_minus_M_content_1 * V_PHYS_SYSALL_content;
                                draw_V_PHYS_DmMSYSALL_P1->SetBinContent(i + 1, j + 1, content);
                            }

                            {
                                #if VECTOR_RANGE_CHECK
                                double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->at(i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P2_data[channel]->at(i + j * 50);
                                #else
                                double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->operator[](i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P2_data[channel]->operator[](i + j * 50);
                                #endif
                                double content = D_minus_M_content_1 * V_PHYS_SYSALL_content;
                                draw_V_PHYS_DmMSYSALL_P2->SetBinContent(i + 1, j + 1, content);
                            }
                        }
                    }


                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        for(Int_t j = 0; j < 50; ++ j)
                        {
                            {
                                #if VECTOR_RANGE_CHECK
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P1_data[channel]->at(i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(j);
                                #else
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P1_data[channel]->operator[](i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->operator[](j);
                                #endif
                                double content = V_PHYS_SYSALL_content * D_minus_M_content_2;
                                draw_V_PHYS_SYSALLDmM_P1->SetBinContent(i + 1, j + 1, content);
                            }

                            {
                                #if VECTOR_RANGE_CHECK
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P2_data[channel]->at(i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(j);
                                #else
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P2_data[channel]->operator[](i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->operator[](j);
                                #endif
                                double content = V_PHYS_SYSALL_content * D_minus_M_content_2;
                                draw_V_PHYS_SYSALLDmM_P2->SetBinContent(i + 1, j + 1, content);
                            }
                        }
                    }


                    for(Int_t i = 0; i < 50; ++ i)
                    {
                        for(Int_t j = 0; j < 50; ++ j)
                        {
                            {
                                #if VECTOR_RANGE_CHECK
                                double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->at(i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P1_data[channel]->at(i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(j);
                                #else
                                double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->operator[](i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P1_data[channel]->operator[](i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->operator[](j);
                                #endif
                                double content = D_minus_M_content_1 * V_PHYS_SYSALL_content * D_minus_M_content_2;
                                draw_V_PHYS_DmMSYSALLDmM_P1->SetBinContent(i + 1, j + 1, content);
                            }

                            {
                                #if VECTOR_RANGE_CHECK
                                double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->at(i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P2_data[channel]->at(i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->at(j);
                                #else
                                double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->operator[](i);
                                double V_PHYS_SYSALL_content = V_PHYS_SYSALL_1D_P2_data[channel]->operator[](i + j * 50);
                                double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->operator[](j);
                                #endif
                                double content = D_minus_M_content_1 * V_PHYS_SYSALL_content * D_minus_M_content_2;
                                draw_V_PHYS_DmMSYSALLDmM_P2->SetBinContent(i + 1, j + 1, content);
                            }
                        }
                    }
                    */
                }
                #endif


                #if DRAWVPHYSMATRIX
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        TString dir;
                        if(V_ENABLE_SYS1 == true)
                        {
                            dir += "SYS1";
                        }
                        if(V_ENABLE_SYS2 == true)
                        {
                            dir += "SYS2";
                        }
                        if(V_ENABLE_SYS3 == true)
                        {
                            dir += "SYS3";
                        }
                        if(V_ENABLE_SYS4 == true)
                        {
                            dir += "SYS4";
                        }
                        if(V_ENABLE_SYS5 == true)
                        {
                            dir += "SYS5";
                        }

                        TCanvas *canvas_V_PHYS_STAT_P1 = new TCanvas("canvas_V_PHYS_STAT_P1", "canvas_V_PHYS_STAT_P1");
                        draw_V_PHYS_STAT_P1->Draw("colz");
                        canvas_V_PHYS_STAT_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_STAT_P1.png");

                        TCanvas *canvas_V_PHYS_SYS1_P1 = new TCanvas("canvas_V_PHYS_SYS1_P1", "canvas_V_PHYS_SYS1_P1");
                        draw_V_PHYS_SYS1_P1->Draw("colz");
                        canvas_V_PHYS_SYS1_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS1_P1.png");

                        TCanvas *canvas_V_PHYS_SYS2_P1 = new TCanvas("canvas_V_PHYS_SYS2_P1", "canvas_V_PHYS_SYS2_P1");
                        draw_V_PHYS_SYS2_P1->Draw("colz");
                        canvas_V_PHYS_SYS2_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS2_P1.png");

                        TCanvas *canvas_V_PHYS_SYS3_P1 = new TCanvas("canvas_V_PHYS_SYS3_P1", "canvas_V_PHYS_SYS3_P1");
                        draw_V_PHYS_SYS3_P1->Draw("colz");
                        canvas_V_PHYS_SYS3_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS3_P1.png");

                        TCanvas *canvas_V_PHYS_SYS4_P1 = new TCanvas("canvas_V_PHYS_SYS4_P1", "canvas_V_PHYS_SYS4_P1");
                        draw_V_PHYS_SYS4_P1->Draw("colz");
                        canvas_V_PHYS_SYS4_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS4_P1.png");

                        TCanvas *canvas_V_PHYS_SYS5_P1 = new TCanvas("canvas_V_PHYS_SYS5_P1", "canvas_V_PHYS_SYS5_P1");
                        draw_V_PHYS_SYS5_P1->Draw("colz");
                        canvas_V_PHYS_SYS5_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS5_P1.png");

                        TCanvas *canvas_V_PHYS_SYSALL_P1 = new TCanvas("canvas_V_PHYS_SYSALL_P1", "canvas_V_PHYS_SYSALL_P1");
                        draw_V_PHYS_SYSALL_P1->Draw("colz");
                        canvas_V_PHYS_SYSALL_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYSALL_P1.png");
                        /*
                        TCanvas *canvas_V_PHYS_DmMSYSALL_P1 = new TCanvas("canvas_V_PHYS_DmMSYSALL_P1", "canvas_V_PHYS_DmMSYSALL_P1");
                        draw_V_PHYS_DmMSYSALL_P1->Draw("colz");
                        canvas_V_PHYS_DmMSYSALL_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_DmMSYSALL_P1.png");

                        TCanvas *canvas_V_PHYS_SYSALLDmM_P1 = new TCanvas("canvas_V_PHYS_SYSALLDmM_P1", "canvas_V_PHYS_SYSALLDmM_P1");
                        draw_V_PHYS_SYSALLDmM_P1->Draw("colz");
                        canvas_V_PHYS_SYSALLDmM_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYSALLDmM_P1.png");

                        TCanvas *canvas_V_PHYS_DmMSYSALLDmM_P1 = new TCanvas("canvas_V_PHYS_DmMSYSALLDmM_P1", "canvas_V_PHYS_DmMSYSALLDmM_P1");
                        draw_V_PHYS_DmMSYSALLDmM_P1->Draw("colz");
                        canvas_V_PHYS_DmMSYSALLDmM_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_DmMSYSALLDmM_P1.png");
                        */
                        TCanvas *canvas_V_PHYS_P1 = new TCanvas("canvas_V_PHYS_P1", "canvas_V_PHYS_P1");
                        draw_V_PHYS_P1->Draw("colz");
                        canvas_V_PHYS_P1->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_P1.png");

                        TCanvas *canvas_V_PHYS_STAT_P2 = new TCanvas("canvas_V_PHYS_STAT_P2", "canvas_V_PHYS_STAT_P2");
                        draw_V_PHYS_STAT_P2->Draw("colz");
                        canvas_V_PHYS_STAT_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_STAT_P2.png");

                        TCanvas *canvas_V_PHYS_SYS1_P2 = new TCanvas("canvas_V_PHYS_SYS1_P2", "canvas_V_PHYS_SYS1_P2");
                        draw_V_PHYS_SYS1_P2->Draw("colz");
                        canvas_V_PHYS_SYS1_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS1_P2.png");

                        TCanvas *canvas_V_PHYS_SYS2_P2 = new TCanvas("canvas_V_PHYS_SYS2_P2", "canvas_V_PHYS_SYS2_P2");
                        draw_V_PHYS_SYS2_P2->Draw("colz");
                        canvas_V_PHYS_SYS2_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS2_P2.png");

                        TCanvas *canvas_V_PHYS_SYS3_P2 = new TCanvas("canvas_V_PHYS_SYS3_P2", "canvas_V_PHYS_SYS3_P2");
                        draw_V_PHYS_SYS3_P2->Draw("colz");
                        canvas_V_PHYS_SYS3_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS3_P2.png");

                        TCanvas *canvas_V_PHYS_SYS4_P2 = new TCanvas("canvas_V_PHYS_SYS4_P2", "canvas_V_PHYS_SYS4_P2");
                        draw_V_PHYS_SYS4_P2->Draw("colz");
                        canvas_V_PHYS_SYS4_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS4_P2.png");

                        TCanvas *canvas_V_PHYS_SYS5_P2 = new TCanvas("canvas_V_PHYS_SYS5_P2", "canvas_V_PHYS_SYS5_P2");
                        draw_V_PHYS_SYS5_P2->Draw("colz");
                        canvas_V_PHYS_SYS5_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYS5_P2.png");

                        TCanvas *canvas_V_PHYS_SYSALL_P2 = new TCanvas("canvas_V_PHYS_SYSALL_P2", "canvas_V_PHYS_SYSALL_P2");
                        draw_V_PHYS_SYSALL_P2->Draw("colz");
                        canvas_V_PHYS_SYSALL_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYSALL_P2.png");
                        /*
                        TCanvas *canvas_V_PHYS_DmMSYSALL_P2 = new TCanvas("canvas_V_PHYS_DmMSYSALL_P2", "canvas_V_PHYS_DmMSYSALL_P2");
                        draw_V_PHYS_DmMSYSALL_P2->Draw("colz");
                        canvas_V_PHYS_DmMSYSALL_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_DmMSYSALL_P2.png");

                        TCanvas *canvas_V_PHYS_SYSALLDmM_P2 = new TCanvas("canvas_V_PHYS_SYSALLDmM_P2", "canvas_V_PHYS_SYSALLDmM_P2");
                        draw_V_PHYS_SYSALLDmM_P2->Draw("colz");
                        canvas_V_PHYS_SYSALLDmM_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_SYSALLDmM_P2.png");

                        TCanvas *canvas_V_PHYS_DmMSYSALLDmM_P2 = new TCanvas("canvas_V_PHYS_DmMSYSALLDmM_P2", "canvas_V_PHYS_DmMSYSALLDmM_P2");
                        draw_V_PHYS_DmMSYSALLDmM_P2->Draw("colz");
                        canvas_V_PHYS_DmMSYSALLDmM_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_DmMSYSALLDmM_P2.png");
                        */
                        TCanvas *canvas_V_PHYS_P2 = new TCanvas("canvas_V_PHYS_P2", "canvas_V_PHYS_P2");
                        draw_V_PHYS_P2->Draw("colz");
                        canvas_V_PHYS_P2->SaveAs(TString("./") + dir + "/" + "draw_V_PHYS_P2.png");
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

                // invert
                V_PHYS_1D_P1_MATHMORE[channel]->Invert();
                V_PHYS_1D_P2_MATHMORE[channel]->Invert();

                #if MEASURE_FUNCTION_CALL_TIME_MATRIX_INVERT
                std::chrono::system_clock::time_point end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> runtime_sec = end_time - start_time;
                std::cout << "Done Invert, time=" << 1.0e+06 * runtime_sec.count() << " microsecond" << std::endl;
                #endif



                #if DRAWVPHYSMATRIX
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

                        //#if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYSINV_P1->SetBinContent(i + 1, j + 1, content);
                                //std::cout << content << "\t";

                                #if VECTOR_RANGE_CHECK
                                double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->at(i);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->at(j);
                                #else
                                double D_minus_M_content_1 = D_minus_M_1D_P1_data[channel]->operator[](i);
                                double D_minus_M_content_2 = D_minus_M_1D_P1_data[channel]->operator[](j);
                                #endif
                                draw_V_PHYS_DmMSYSALLDmM_P1->SetBinContent(i + 1, j + 1, D_minus_M_content_1 * content * D_minus_M_content_2);
                            }
                        }
                        //#endif
                        ++ i_counter;
                    }
                    //#if DRAWVPHYSMATRIX
                    if(DRAW_V_PHYS_MATRIX == true)
                    {
                        //std::cout << std::endl;
                    }
                    //#endif

                    ++ j_counter;
                }
                #endif

                #if DRAWVPHYSMATRIX
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

                        //#if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            if(channel == 1)
                            {
                                draw_V_PHYSINV_P2->SetBinContent(i + 1, j + 1, content);
                                //std::cout << content << "\t";

                                #if VECTOR_RANGE_CHECK
                                double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->at(i);
                                double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->at(j);
                                #else
                                double D_minus_M_content_1 = D_minus_M_1D_P2_data[channel]->operator[](i);
                                double D_minus_M_content_2 = D_minus_M_1D_P2_data[channel]->operator[](j);
                                #endif
                                draw_V_PHYS_DmMSYSALLDmM_P2->SetBinContent(i + 1, j + 1, D_minus_M_content_1 * content * D_minus_M_content_2);
                            }
                        }
                        //#endif
                        ++ i_counter;
                    }
                        //#if DRAWVPHYSMATRIX
                        if(DRAW_V_PHYS_MATRIX == true)
                        {
                            //std::cout << std::endl;
                        }
                        //#endif

                    ++ j_counter;
                }
                #endif


                #if DRAWVPHYSMATRIX
                if(DRAW_V_PHYS_MATRIX == true)
                {
                    if(channel == 1)
                    {
                        std::cout << "ALLOC " << __func__ << std::endl;
                        TCanvas *canvas_V_PHYSINV_P1 = new TCanvas("canvas_V_PHYSINV_P1", "canvas_V_PHYSINV_P1");
                        draw_V_PHYSINV_P1->Draw("colz");
                        canvas_V_PHYSINV_P1->SaveAs("draw_V_PHYSINV_P1.png");

                        TCanvas *canvas_V_DmMPHYSINVDmM_P1 = new TCanvas("canvas_V_DmMPHYSINVDmM_P1", "canvas_V_DmMPHYSINVDmM_P1");
                        draw_V_DmMPHYSINVDmM_P1->Draw("colz");
                        canvas_V_DmMPHYSINVDmM_P1->SaveAs("draw_V_DmMPHYSINVDmM_P1.png");

                        TCanvas *canvas_V_PHYSINV_P2 = new TCanvas("canvas_V_PHYSINV_P2", "canvas_V_PHYSINV_P2");
                        draw_V_PHYSINV_P2->Draw("colz");
                        canvas_V_PHYSINV_P2->SaveAs("draw_V_PHYSINV_P2.png");

                        TCanvas *canvas_V_DmMPHYSINVDmM_P2 = new TCanvas("canvas_V_DmMPHYSINVDmM_P2", "canvas_V_DmMPHYSINVDmM_P2");
                        draw_V_DmMPHYSINVDmM_P2->Draw("colz");
                        canvas_V_DmMPHYSINVDmM_P2->SaveAs("draw_V_DmMPHYSINVDmM_P2.png");

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
