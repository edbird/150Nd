

#include <string>
#include <iostream>
#include <fstream>

void scriptwrite()
{

    /*
    int start[12] = 
    {
          0,  25,  50,  75,
        100, 125, 150, 175,
        200, 225, 250, 275,
    };
    int stop[12] = 
    {
         25,  50,  75, 100,
        125, 150, 175, 200,
        225, 250, 275, 301
    };
    */
    int start[12] = 
    {
          51,  71,  91,  111,
        131, 141, 150, 175,
        200, 225, 250, 275,
    };
    int stop[12] = 
    {
         71,  91,  111, 131,
        141, 150, 175, 200,
        225, 250, 275, 301
    };

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

}
