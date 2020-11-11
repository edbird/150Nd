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
        systematic_foil_thickness_virtual = 0.0;
        systematic_foil_thickness_virtual_last = -std::numeric_limits<double>::infinity();
        systematic_dEdX_virtual = 0.0;
        systematic_dEdX_virtual_last = -std::numeric_limits<double>::infinity();
        systematic_brem_virtual = 0.0;
        systematic_brem_virtual_last = -std::numeric_limits<double>::infinity();
        systematic_foil_thickness_nominal = 0.0;
        systematic_foil_thickness_nominal_last = -std::numeric_limits<double>::infinity();
        systematic_dEdX_nominal = 0.0;
        systematic_dEdX_nominal_last = -std::numeric_limits<double>::infinity();
        systematic_brem_nominal = 0.0;
        systematic_brem_nominal_last = -std::numeric_limits<double>::infinity();
        systematic_optical_correction = 0.0;
        systematic_optical_correction_last = -std::numeric_limits<double>::infinity();

        aux_data_is_loaded = false;

        h_systematic_foil_thickness_virtual_h = nullptr;
        h_systematic_foil_thickness_virtual_l = nullptr;
        h_systematic_dEdX_virtual_h = nullptr;
        h_systematic_dEdX_virtual_l = nullptr;
        h_systematic_brem_virtual_h = nullptr;
        h_systematic_brem_virtual_l = nullptr;
        h_systematic_foil_thickness_nominal_h = nullptr;
        h_systematic_foil_thickness_nominal_l = nullptr;
        h_systematic_dEdX_nominal_h = nullptr;
        h_systematic_dEdX_nominal_l = nullptr;
        h_systematic_brem_nominal_h = nullptr;
        h_systematic_brem_nominal_l = nullptr;

        load_aux_data();
    }



    void reset()
    { 
        systematic_energy_offset = 0.0;
        systematic_energy_scale = 0.0;
        systematic_efficiency = 0.0;
        systematic_enrichment = 0.0;
        systematic_energy_offsetsmall = 0.0;
        systematic_foil_thickness_virtual = 0.0;
        systematic_dEdX_virtual = 0.0;
        systematic_brem_virtual = 0.0;
        systematic_foil_thickness_nominal = 0.0;
        systematic_dEdX_nominal = 0.0;
        systematic_brem_nominal = 0.0;
        systematic_optical_correction = 0.0;
    }


    void load_aux_data()
    {
        if(aux_data_is_loaded == false)
        {
            {
                TFile *f1 = new TFile("systematic_foil_thickness_virtual.root");
                h_systematic_foil_thickness_virtual_h = (TH2D*)f1->Get("h_systematic_foil_thickness_h");
                h_systematic_foil_thickness_virtual_l = (TH2D*)f1->Get("h_systematic_foil_thickness_l");
                if(h_systematic_foil_thickness_virtual_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_virtual.root" << std::endl;
                }
                if(h_systematic_foil_thickness_virtual_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_virtual.root" << std::endl;
                }
            }

            {
                TFile *f2 = new TFile("systematic_dEdX_virtual.root");
                h_systematic_dEdX_virtual_h = (TH2D*)f2->Get("h_systematic_dEdX_h");
                h_systematic_dEdX_virtual_l = (TH2D*)f2->Get("h_systematic_dEdX_l");
                if(h_systematic_dEdX_virtual_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_virtual.root" << std::endl;
                }
                if(h_systematic_dEdX_virtual_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_virtual.root" << std::endl;
                }
            }

            {
                TFile *f3 = new TFile("systematic_brem_virtual.root");
                h_systematic_brem_virtual_h = (TH2D*)f3->Get("h_systematic_brem_h");
                h_systematic_brem_virtual_l = (TH2D*)f3->Get("h_systematic_brem_l");
                if(h_systematic_brem_virtual_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_virtual.root" << std::endl;
                }
                if(h_systematic_brem_virtual_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_virtual.root" << std::endl;
                }
            }

            {
                TFile *f4 = new TFile("systematic_foil_thickness_nominal.root");
                h_systematic_foil_thickness_nominal_h = (TH2D*)f4->Get("h_systematic_foil_thickness_h");
                h_systematic_foil_thickness_nominal_l = (TH2D*)f4->Get("h_systematic_foil_thickness_l");
                if(h_systematic_foil_thickness_nominal_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_nominal.root" << std::endl;
                }
                if(h_systematic_foil_thickness_nominal_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_foil_thickness_nominal.root" << std::endl;
                }
            }

            {
                TFile *f5 = new TFile("systematic_dEdX_nominal.root");
                h_systematic_dEdX_nominal_h = (TH2D*)f5->Get("h_systematic_dEdX_h");
                h_systematic_dEdX_nominal_l = (TH2D*)f5->Get("h_systematic_dEdX_l");
                if(h_systematic_dEdX_nominal_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_nominal.root" << std::endl;
                }
                if(h_systematic_dEdX_nominal_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_dEdX_nominal.root" << std::endl;
                }
            }

            {
                TFile *f6 = new TFile("systematic_brem_nominal.root");
                h_systematic_brem_nominal_h = (TH2D*)f6->Get("h_systematic_brem_h");
                h_systematic_brem_nominal_l = (TH2D*)f6->Get("h_systematic_brem_l");
                if(h_systematic_brem_nominal_h == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_nominal.root" << std::endl;
                }
                if(h_systematic_brem_nominal_l == nullptr)
                {
                    std::cout << "Error: could not load systematic data from file " << "systematic_brem_nominal.root" << std::endl;
                }
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
    double systematic_foil_thickness_virtual;
    double systematic_foil_thickness_virtual_last;
    double systematic_dEdX_virtual;
    double systematic_dEdX_virtual_last;
    double systematic_brem_virtual;
    double systematic_brem_virtual_last;
    double systematic_foil_thickness_nominal;
    double systematic_foil_thickness_nominal_last;
    double systematic_dEdX_nominal;
    double systematic_dEdX_nominal_last;
    double systematic_brem_nominal;
    double systematic_brem_nominal_last;
    double systematic_optical_correction;
    double systematic_optical_correction_last;

    bool aux_data_is_loaded;

    TH2D *h_systematic_foil_thickness_virtual_h;
    TH2D *h_systematic_foil_thickness_virtual_l;
    TH2D *h_systematic_dEdX_virtual_h;
    TH2D *h_systematic_dEdX_virtual_l;
    TH2D *h_systematic_brem_virtual_h;
    TH2D *h_systematic_brem_virtual_l;
    TH2D *h_systematic_foil_thickness_nominal_h;
    TH2D *h_systematic_foil_thickness_nominal_l;
    TH2D *h_systematic_dEdX_nominal_h;
    TH2D *h_systematic_dEdX_nominal_l;
    TH2D *h_systematic_brem_nominal_h;
    TH2D *h_systematic_brem_nominal_l;


}gSystematics;



void gSystematics_print()
{
        std::cout << "SYSTEMATICS: CONSTANT OFFSET: " << gSystematics.systematic_energy_offset << " MeV" << std::endl;
        std::cout << "SYSTEMATICS: CONSTANT SCALE: " << gSystematics.systematic_energy_scale << " MeV" << std::endl;
        std::cout << "SYSTEMATICS: EFFICIENCY: " << gSystematics.systematic_efficiency << "" << std::endl;
        std::cout << "SYSTEMATICS: ENRICHMENT: " << gSystematics.systematic_enrichment << "" << std::endl;
        std::cout << "SYSTEMATICS: OFFSETSMALL: " << gSystematics.systematic_energy_offsetsmall << " MeV" << std::endl;
        std::cout << "SYSTEMATICS: FOIL THICKNESS (V): " << gSystematics.systematic_foil_thickness_virtual << " " << std::endl;
        std::cout << "SYSTEMATICS: ENERGY LOSS (V): " << gSystematics.systematic_dEdX_virtual << " " << std::endl;
        std::cout << "SYSTEMATICS: BREMSSTRAHLUNG (V): " << gSystematics.systematic_brem_virtual << " " << std::endl;
        std::cout << "SYSTEMATICS: FOIL THICKNESS (N): " << gSystematics.systematic_foil_thickness_nominal << " " << std::endl;
        std::cout << "SYSTEMATICS: ENERGY LOSS (N): " << gSystematics.systematic_dEdX_nominal << " " << std::endl;
        std::cout << "SYSTEMATICS: BREMSSTRAHLUNG (N): " << gSystematics.systematic_brem_nominal << " " << std::endl;
}


#endif // SYSTEMATICS_H
