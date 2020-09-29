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
    }

    double systematic_energy_offset;
    double systematic_energy_offset_last; // think not used
    double systematic_energy_scale;
    double systematic_energy_scale_last; // think not used


}gSystematics;


#endif // SYSTEMATICS_H
