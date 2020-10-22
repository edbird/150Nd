#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H



class Systematics
{


    public:


    Systematics()
    {
        systematic_energy_offset = 0.0;
        systematic_energy_offset_last = -std::numeric_limits<double>::infinity();
        systematic_energy_scale = 0.0;
        systematic_energy_scale_last = -std::numeric_limits<double>::infinity();
        systematic_efficiency = 0.0;
        systematic_efficiency_last = -std::numeric_limits<double>::infinity();
        systematic_enrichment = 0.0;
        systematic_enrichment_last = -std::numeric_limits<double>::infinity();
        systematic_energy_offsetsmall = 0.0;
        systematic_energy_offsetsmall_last = -std::numeric_limits<double>::infinity();
        systematic_foil_thickness = 0.0;
        systematic_foil_thickness_last = -std::numeric_limits<double>::infinity();
        systematic_dEdX = 0.0;
        systematic_dEdX_last = -std::numeric_limits<double>::infinity();
        systematic_brem = 0.0;
        systematic_brem_last = -std::numeric_limits<double>::infinity();
        systematic_optical_correction = 0.0;
        systematic_optical_correction_last = -std::numeric_limits<double>::infinity();

        aux_data_is_loaded = false;

        h_systematic_foil_thickness_h = nullptr;
        //h_systematic_foil_thickness_l = nullptr;

        h_systematic_dEdX = nullptr;
        h_systematic_brem = nullptr;

        load_aux_data();
    }


    void reset()
    { 
        systematic_energy_offset = 0.0;
        systematic_energy_scale = 0.0;
        systematic_efficiency = 0.0;
        systematic_enrichment = 0.0;
        systematic_energy_offsetsmall = 0.0;
        systematic_foil_thickness = 0.0;
        systematic_dEdX = 0.0;
        systematic_brem = 0.0;
        systematic_optical_correction = 0.0;
    }


    void load_aux_data()
    {
        if(aux_data_is_loaded == false)
        {
            TFile *f1 = new TFile("systematic_foil_thickness.root");
            h_systematic_foil_thickness_h = (TH2D*)f1->Get("h_systematic_foil_thickness");
            //h_systematic_foil_thickness_l = (TH2D*)f1->Get("h_systematic_foil_thickness_l");
            if(h_systematic_foil_thickness_h == nullptr)
            {
                std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness.root" << std::endl;
            }

            TFile *f2 = new TFile("systematic_dEdX.root");
            h_systematic_dEdX = (TH2D*)f2->Get("h_systematic_dEdX");
            if(h_systematic_dEdX == nullptr)
            {
                std::cout << "Error: could not load systematic data from file " << "systematic_dEdX.root" << std::endl;
            }

            TFile *f3 = new TFile("systematic_brem.root");
            h_systematic_brem = (TH2D*)f3->Get("h_systematic_brem");
            if(h_systematic_brem == nullptr)
            {
                std::cout << "Error: could not load systematic data from file " << "systematic_brem.root" << std::endl;
            }
            
            aux_data_is_loaded = true;
        }
        
    }


    double systematic_energy_offset;
    double systematic_energy_offset_last; // think not used
    double systematic_energy_scale;
    double systematic_energy_scale_last; // think not used
    double systematic_efficiency;
    double systematic_efficiency_last; // think not used
    double systematic_enrichment;
    double systematic_enrichment_last; // think not used
    double systematic_energy_offsetsmall;
    double systematic_energy_offsetsmall_last; // think not used
    double systematic_foil_thickness;
    double systematic_foil_thickness_last;
    double systematic_dEdX;
    double systematic_dEdX_last;
    double systematic_brem;
    double systematic_brem_last;
    double systematic_optical_correction;
    double systematic_optical_correction_last;

    bool aux_data_is_loaded;

    TH2D *h_systematic_foil_thickness_h;
    //TH2D *h_systematic_foil_thickness_l;
    TH2D *h_systematic_dEdX;
    TH2D *h_systematic_brem;


}gSystematics;


#endif // SYSTEMATICS_H
