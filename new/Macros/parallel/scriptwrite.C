

#include <string>
#include <iostream>
#include <fstream>

void scriptwrite()
{

    const int NUM_CORES = 12;
    const int RES = 301;

    int start[NUM_CORES] = 
    {
        0,  0,  0,  0,
        0,  0,  0,  0,
        0,  0,  0,  0
    };
    int stop[NUM_CORES] = 
    {
        0,  0,  0,  0,
        0,  0,  0,  0,
        0,  0,  0,  0
    };

    for(int i = 0; i < NUM_CORES; ++ i)
    {
        int si = (i * RES) / NUM_CORES;
        int sf = ((i + 1) * RES) / NUM_CORES;
        start[i] = si;
        stop[i] = sf;
        std::cout << "si=" << si << std::endl;
        std::cout << "sf=" << sf << std::endl;
        //std::cout << std::endl;
    }

    // NOTE: use script index starting from 1 so that JID0 is reserved
    // for concatinated output
    for(int i = 0; i < NUM_CORES; ++ i)
    {
        std::string name = std::string("script") + std::to_string(i + 1) + ".txt";
        std::ofstream ofs(name);
        ofs << "NUMBER=" << std::to_string(i + 1) << std::endl;
        //#ofs << "OUTPUT_NAME=batch" << std::to_string(i) << ".txt" << std::endl;
        ofs << "OUTPUT_NAME=mps_resultsmatrix" << std::endl;
        //ofs << "START_INDEX=" << std::to_string(i) << std::endl;
        //ofs << "STOP_INDEX=" << std::to_string(i + 1) << std::endl;
        ofs << "START_INDEX=" << std::to_string(start[i]) << std::endl;
        ofs << "STOP_INDEX=" << std::to_string(stop[i]) << std::endl;
        ofs << "RUNNING=false" << std::endl;
        ofs.close();
    }

    /*
    for(int i = 0; i < 12; ++ i)
    {
        std::string name = std::string("script") + std::to_string(i) + ".txt";
        std::ofstream ofs(name);
        ofs << "NUMBER=" << std::to_string(i) << std::endl;
        //#ofs << "OUTPUT_NAME=batch" << std::to_string(i) << ".txt" << std::endl;
        ofs << "OUTPUT_NAME=mps_resultsmatrix" << std::endl;
        //ofs << "START_INDEX=" << std::to_string(i) << std::endl;
        //ofs << "STOP_INDEX=" << std::to_string(i + 1) << std::endl;
        ofs << "START_INDEX=" << std::to_string(start[i]) << std::endl;
        ofs << "STOP_INDEX=" << std::to_string(stop[i]) << std::endl;
        ofs << "RUNNING=false" << std::endl;
        ofs.close();
    }
    */
}
