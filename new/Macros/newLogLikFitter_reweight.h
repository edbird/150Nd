#ifndef NEWLOGLIKFITTER_REWEIGHT_H
#define NEWLOGLIKFITTER_REWEIGHT_H



void convert_data_to_histogram_format(const std::vector<std::vector<double>> &data_nEqNull,
                                      const std::vector<std::vector<double>> &data_nEqTwo,
                                      TH2D * &h_nEqNull_return,
                                      TH2D * &h_nEqTwo_return,
                                      const Int_t dimension_xy)
{

    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////

    //const Int_t dimension_xy{3371};
    //const Int_t dimension_xy{1001};
    // don't plot raw data
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    const Double_t keV = 1.0e-03;
    const Double_t E_min = 0.38e-03;
    const Double_t limit_lower  = E_min - 0.5 * keV;
    const Double_t limit_higher = E_min - 0.5 * keV + dimension_xy * keV;
    TH2D *h_nEqNull =
        new TH2D("h_nEqNull", "nEqNull",
                 dimension_xy, limit_lower, limit_higher,
                 dimension_xy, limit_lower, limit_higher);
    TH2D *h_nEqTwo =
        new TH2D("h_nEqTwo", "nEqTwo",
                 dimension_xy, limit_lower, limit_higher,
                 dimension_xy, limit_lower, limit_higher);
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    h_nEqNull->SetStats(0);
    h_nEqTwo->SetStats(0);
    //h_ratio->SetStats(0);

    // required because Nd-150 data is triangular
    std::size_t rd_ix{0};
    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            if(i < dimension_xy - j)
            {
                h_nEqNull->SetBinContent(i, j, data_nEqNull.at(rd_ix)[2]);
                h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(rd_ix)[2]);
            
                ++ rd_ix;
            }
            else
            {
                h_nEqNull->SetBinContent(i, j, 0.0);
                h_nEqTwo->SetBinContent(i, j, 0.0);
            }
            //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            //if(i + j < dimension_xy - 1)
            //{
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            //}
        }
    }
    std::cout << "rd_ix_max=" << rd_ix << std::endl;
    std::cout << "Finished constructing input data histograms" << std::endl;
    
    h_nEqNull_return = h_nEqNull;
    h_nEqTwo_return = h_nEqTwo;

}


// TODO: this is probably equivalent to my previous "manual"
// reading code, but should check
// 2020-10-03: NOTE: This function is not used for anything.
// It gave different results to loading data manually and I did not
// find out why.
void read_data(TTree *tree, std::vector<std::vector<double>>& data)
{
    
    data.clear();

    Float_t electronEnergy1;
    Float_t electronEnergy2;
    Float_t electronWeight;

    tree->SetBranchAddress("electronEnergy1", &electronEnergy1);
    tree->SetBranchAddress("electronEnergy2", &electronEnergy2);
    tree->SetBranchAddress("electronWeight", &electronWeight);

    for(Long64_t event = 0; event < tree->GetEntries(); ++ event)
    {

        tree->GetEntry(event);


        std::cout << std::scientific;

        std::cout << "electronEnergy1=" << electronEnergy1
                  << " electronEnergy2=" << electronEnergy2
                  << " electronWeight=" << electronWeight << std::endl;
        std::cin.get();


        data.emplace_back(std::vector<double>());
        data.back().emplace_back(electronEnergy1);
        data.back().emplace_back(electronEnergy2);
        data.back().emplace_back(electronWeight);

    }

}


std::size_t LD_REWEIGHT_DATA_2(
    TH2D *&h_out,
    const std::string &filename,
    const std::string& h_name,
    const std::string& name//,
    //const double low,
    //const double high
    )
{

    

    //TTree *t_G0 = new TTree();
    //Long64_t count_G0 = t_G0->ReadFile("../../data/150Nd-data/dG150Nd/G0/dG0.dat", "electronEnergy1:electronEnergy2:electronWeight");
    //std::cout << "count_G0=" << count_G0 << std::endl;

    //Long64_t count_G0 = 0;
    std::vector<std::vector<double>> data;
    std::size_t read_count = 0;
    std::ifstream if_Gn(filename);
    if(!if_Gn.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return 0;
    }
    double electronEnergy1, electronEnergy2, electronWeight;
    while(if_Gn.good())
    {
        if_Gn >> electronEnergy1 >> electronEnergy2 >> electronWeight;
        //++ count_G0;
        ++ read_count;

       
        //if(read_count < 10)
        //{
        //    std::cout << "electronEnergy1=" << electronEnergy1
        //          << " electronEnergy2=" << electronEnergy2
        //          << " electronWeight=" << electronWeight << std::endl;
        //}
        //std::cin.get();


        data.emplace_back(std::vector<double>());
        data.back().emplace_back(electronEnergy1);
        data.back().emplace_back(electronEnergy2);
        data.back().emplace_back(electronWeight);
    }
    std::cout << "read_count=" << read_count << std::endl;

    //TTree *t_G2 = new TTree();
    //Long64_t count_G2 = t_G2->ReadFile("../../data/150Nd-data/dG150Nd/G2/dG2.dat", "electronEnergy1:electronEnergy2:electronWeight");
    //std::cout << "count_G2=" << count_G2 << std::endl;

    //output = out;
    //outputG0 = t_G0;
    //outputG2 = t_G2;
    //countG0 = count_G0;
    //countG2 = count_G2;

    // solve quadratic
    //const double A = 0.5;
    //const double B = 0.5;
    //double C = -(double)count;
    //double dimension_xy = (-B + std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
    Long64_t count = read_count;
    double dimension_xy = (int)(-0.5 + std::sqrt(0.25 - 2.0 * (-(double)count)));
    if(dimension_xy != 3371)
    {
        std::cout << "error! dimension_xy != 3371, dimension_xy=" << dimension_xy << std::endl;
    }
    //const double dimension_xy{3371};

    std::cout << "dimension_xy=" << dimension_xy << std::endl;


    ////////////////////////////////////////////////////////////////////////////
    // CONVERT INPUT DATA TO HISTOGRAM FORMAT
    // Note: Added 2018-04-23 (After INTERMEDIATE DATA below)
    // Note: These histograms do NOT have the phase space variable included
    ////////////////////////////////////////////////////////////////////////////
                                      
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "nEqNull", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "nEqTwo", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);


    //const Int_t dimension_xy{3371};
    //const Int_t dimension_xy{1001};
    // don't plot raw data
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_nEqNull = new TH2D("h_nEqNull", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_nEqTwo = new TH2D("h_nEqTwo", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //TH2D *h_ratio = new TH2D("h_ratio", "", dimension_xy, 0.0, 1.0, dimension_xy, 0.0, 1.0);
    //h_nEqNull = new TH2D("h_nEqNull", "nEqNull", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //h_nEqTwo = new TH2D("h_nEqTwo", "nEqTwo", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    //TH2D *h_ = new TH2D(hname.c_str(), name.c_str(), dimension_xy, low, high, dimension_xy, low, high);
    const Double_t keV = 1.0e-03;
    const Double_t E_min = 0.38e-03;
    const Double_t limit_lower  = E_min - 0.5 * keV;
    const Double_t limit_higher = E_min - 0.5 * keV + dimension_xy * keV;
    TH2D *h_ = new TH2D(h_name.c_str(), h_name.c_str(),
                        dimension_xy, limit_lower, limit_higher,
                        dimension_xy, limit_lower, limit_higher);
    //h_nEqTwo = new TH2D("h_nEqTwo", "nEqTwo", dimension_xy, 0.0, bb_Q, dimension_xy, 0.0, bb_Q);
    h_->SetStats(0);
    //h_nEqTwo->SetStats(0);
    //h_ratio->SetStats(0);

    // required because Nd-150 data is triangular
    std::size_t rd_ix{0};
    for(std::size_t i{0}; i < dimension_xy; ++ i)
    {
        for(std::size_t j{0}; j < dimension_xy; ++ j)
        {
            if(i < dimension_xy - j)
            {
                //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(rd_ix)[2]);
                //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(rd_ix)[2]);
                h_->SetBinContent(i, j, data.at(rd_ix)[2]);
                //if(i < 5 && j < 10)
                //{
                //    std::cout << "fill: " << i << " " << j << " " << data.at(rd_ix)[2] << std::endl;
                //}
            
                ++ rd_ix;
            }
            else
            {
                //h_nEqNull->SetBinContent(i, j, 0.0);
                //h_nEqTwo->SetBinContent(i, j, 0.0);
                h_->SetBinContent(i, j, 0.0);
            }
            //h_nEqNull->SetBinContent(i, j, data_nEqNull.at(i * dimension_xy + j)[2]);
            //h_nEqTwo->SetBinContent(i, j, data_nEqTwo.at(i * dimension_xy + j)[2]);
            //if(i < dimension_xy - j)
            //if(i + j < dimension_xy - 1)
            //{
                // TODO: move above lines to inside this if
                //h_ratio->SetBinContent(i, j, ratio.at(i * dimension_xy + j)[2]);
            //}
        }
    }
    std::cout << "rd_ix_max=" << rd_ix << std::endl;
    std::cout << "Finished constructing input data histogram" << std::endl;
    
    //h_nEqNull_return = h_nEqNull;
    //h_nEqTwo_return = h_nEqTwo;
    h_out = h_;

    //TCanvas *c_ = new TCanvas("c_", "c_");
    //h_->Draw("colz");
    //c_->SaveAs("hout.png");
    //std::cout << "histogram done" << std::endl;
    //std::cin.get();

    return read_count;

}


void LD_REWEIGHT_DATA(TTree *&outputG0, TTree *&outputG2, Long64_t &countG0, Long64_t &countG2)
{
    //TH1D* out = nullptr;

    TTree *t_G0 = new TTree();
    Long64_t count_G0 = t_G0->ReadFile("../../data/150Nd-data/dG150Nd/G0/dG0.dat", "electronEnergy1:electronEnergy2:electronWeight");
    std::cout << "count_G0=" << count_G0 << std::endl;

    TTree *t_G2 = new TTree();
    Long64_t count_G2 = t_G2->ReadFile("../../data/150Nd-data/dG150Nd/G2/dG2.dat", "electronEnergy1:electronEnergy2:electronWeight");
    std::cout << "count_G2=" << count_G2 << std::endl;

    //output = out;
    outputG0 = t_G0;
    outputG2 = t_G2;
    countG0 = count_G0;
    countG2 = count_G2;

    return;
}

// function dec
Double_t ReWeight3(
    const Double_t T1,
    const Double_t T2,
    const Double_t epsilon_baseline,
    const Double_t epsilon,
    const TH2D* const h_nEqNull,
    const TH2D* const h_nEqTwo,
    const Double_t psiN0,
    const Double_t psiN2,
    const std::string& debug);


// used to draw hTotalE_ after each call to reweight
// for debugging purposes
static int ctmp_counter = 0;







// electron energies: T1, T2 (in units of Q value)
// theory parameter: epsilon_baseline
// theory parameter: epsilon
// data: G0, G2 (units of Q value)
//Double_t ReWeight(const Double_t T1, const Double_t T2, const Double_t epsilon,
//                  const std::vector<std::vector<double>>& data_nEqNull,
//                  const std::vector<std::vector<double>>& data_nEqTwo)
Double_t ReWeight3(const Double_t T1, const Double_t T2,
                  const Double_t epsilon_baseline,
                  const Double_t epsilon,
                  const TH2D* const h_nEqNull,
                  const TH2D* const h_nEqTwo,
                  const Double_t psiN0, const Double_t psiN2,
                  const std::string& debug = "")
{


    /*
    const double E_min = 0.38e-03; // keV
    const double E_bin_width = 1.0e-03; // 1 keV bin width
    double T1_bin = (T1 - E_min) / E_bin_width;
    double T2_bin = (T2 - E_min) / E_bin_width;
    int T1_bin_ix = (int)std::round(T1_bin);
    int T2_bin_ix = (int)std::round(T2_bin);
    */

    // TODO: NO INTERPOLATION DONE YET

    if(T1 < 0.38e-03 - 0.5e-03)
    {
        std::cout << "problem: lower limit" << std::endl;
        throw "problem";
    }
    if(T1 > 0.38e-03 - 0.5e-03 + 1.0e-03 * 3371.0)
    {
        std::cout << "problem: upper limit" << std::endl;
        throw "problem";
    }

    // find bin corresponding to energies T1, T2
    Int_t bin_x{h_nEqNull->GetXaxis()->FindBin(T1)};
    Int_t bin_y{h_nEqNull->GetYaxis()->FindBin(T2)};
    // TODO: do this conversion manually
    //Int_t bin_x{T1_bin_ix};
    //Int_t bin_y{T2_bin_ix};

    //std::cout << "T1=" << T1 << " T2=" << T2 << std::endl;
    //std::cout << "the bin is " << bin_x << " " << bin_y << std::endl;
    //std::cin.get();

    // bin content corresponding to energies T1, T2
    Double_t h_nEqNull_c{h_nEqNull->GetBinContent(bin_x, bin_y)};
    Double_t h_nEqTwo_c{h_nEqTwo->GetBinContent(bin_x, bin_y)};
  
    h_nEqNull_c = ((TH2*)h_nEqNull)->Interpolate(T1, T2);
    h_nEqTwo_c = ((TH2*)h_nEqTwo)->Interpolate(T1, T2);


    // get the weight for this T1, T2
    // the input data is for epsilon = 0.0
    //Double_t phase_1{1.0 / psiN0};
    //Double_t weight_1{phase_1 * h_nEqNull->GetBinContent(bin_x, bin_y)};
    Double_t phase_1{1.0 / (psiN0 + epsilon_baseline * psiN2)};
    Double_t weight_1{phase_1 * (h_nEqNull_c + epsilon_baseline * h_nEqTwo_c)};

    // get the weight for this T1, T2
    // the input data is for epsilon = some arbitary value
    Double_t phase_2{1.0 / (psiN0 + epsilon * psiN2)};
    Double_t weight_2{phase_2 * (h_nEqNull_c + epsilon * h_nEqTwo_c)};

    if(debug == "true")
    {
        if(std::isnan(weight_2 / weight_1))
        {
            std::cout << std::scientific;
            std::cout << "epsilon=" << epsilon << " epsilon_baseline=" << epsilon_baseline << std::endl;
            std::cout << "T1=" << T1 << " T2=" << T2 << std::endl;
            std::cout << "bin_x=" << bin_x << " bin_y=" << bin_y << std::endl;
            std::cout << "h_nEqNull->GetBinContent(bin_x, bin_y)=" << h_nEqNull->GetBinContent(bin_x, bin_y) << std::endl;
            std::cout << "h_nEqTwo->GetBinContent(bin_x, bin_y)=" << h_nEqTwo->GetBinContent(bin_x, bin_y) << std::endl;
            std::cout << "weight_1=" << weight_1 << " weight_2=" << weight_2 << std::endl;
            std::cout << "phase_1=" << phase_1 << " phase_2=" << phase_2 << std::endl;
            std::cout << "psiN0=" << psiN0 << " psiN2=" << psiN2 << std::endl;
            std::cin.get();
        }
    }


    Double_t weight = weight_2 / weight_1;

/*
    if(std::abs(weight - 1.0) > 1.0e-5)
    {
        std::cout << std::scientific;
        std::cout << "weight=" << weight << " weight_1=" << weight_1 << " weight_2=" << weight_2 << std::endl;
        std::cout << h_nEqNull << " " << h_nEqTwo << std::endl;
        std::cout << "c: " << h_nEqNull_c << ", " << h_nEqTwo_c << std::endl;
        std::cout << "eps: " << epsilon << ", " << epsilon_baseline << std::endl;
        std::cout << "psiN0=" << psiN0 << " psiN2=" << psiN2 << std::endl;
        std::cout << "phase_1=" << phase_1 << " phase_2=" << phase_2 << std::endl;
        std::cout << "T1=" << T1 << " T2=" << T2 << " bin_x=" << bin_x << " bin_y=" << bin_y << std::endl;
        std::cin.get();
    }
*/

    // return re-weighting factor
    return weight;


}

// TODO: move elsewhere
//void reweight_apply(TH1D *output, const TH1D *const input, ...)
//{
//
//    for(Int_t bin_ix{1}; bin_ix <= input->GetNBinsX(); ++ bin_ix)
//    {
//        Double_t bin_energy = ???
//        // >>> cannot reweight in this manner
//        // >>> do not know T1, T2, need to rebuild histogram by reading from TTree
//        // or by storing as 2d histogram
//    }
//
//    ReWeight3(T1, T2, epsilon_baseline, h_nEqNull, h_nEqTwo, psiN0, psiN2, "false");
//    
//
//}

#endif // NEWLOGLIKFITTER_REWEIGHT_H
