

#include <string>
#include <iostream>
#include <fstream>

void scriptwrite()
{

    for(int i = 0; i < 301; ++ i)
    {
        std::string name = std::string("script") + std::to_string(i) + ".txt";
        std::ofstream ofs(name);
        ofs << "NUMBER=" << std::to_string(i) << std::endl;
        ofs << "OUTPUT_NAME=batch" << std::to_string(i) << ".txt" << std::endl;
        ofs << "START_INDEX=" << std::to_string(i) << std::endl;
        ofs << "STOP_INDEX=" << std::to_string(i + 1) << std::endl;
        ofs << "RUNNING=false" << std::endl;
        ofs.close();
    }

}
