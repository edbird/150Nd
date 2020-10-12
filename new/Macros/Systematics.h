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
    }


    void reset()
    { 
        systematic_energy_offset = 0.0;
        systematic_energy_scale = 0.0;
        systematic_efficiency = 0.0;
        systematic_enrichment = 0.0;
        systematic_energy_offsetsmall = 0.0;
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


}gSystematics;


#endif // SYSTEMATICS_H
