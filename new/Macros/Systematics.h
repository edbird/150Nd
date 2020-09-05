#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H


class Systematics
{


    public:


    Systematics()
    {
        systematic_energy_offset = 0.0;
        systematic_energy_offset_last = -std::numeric_limits<double>::infinity();
        systematic_energy_multiply = 0.0;
    }

    double systematic_energy_offset;
    double systematic_energy_offset_last;
    double systematic_energy_multiply;


}gSystematics;


#endif // SYSTEMATICS_H
