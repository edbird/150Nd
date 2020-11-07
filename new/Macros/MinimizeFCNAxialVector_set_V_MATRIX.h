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

                        ++ i_counter;
                    }

                    ++ j_counter;
                }

                //std::cout << "AT END OF SET V_MATRIX_?_MATHMORE" << std::endl;
                //std::cin.get();



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


                
                // disable this to recalculate the V MATRIX each loop
                // if the number of singular elements changes it may crash
                // not sure at the moment
                //recalculate_V_PHYS_xD_Px_MATHMORE = false;
            }
        }
    }
}

#endif // MINIMIZEFCNAXIALVECTOR_SET_V_MATRIX_H
