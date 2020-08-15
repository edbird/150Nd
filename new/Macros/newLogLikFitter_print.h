#ifndef NEWLOGLIKFITTER_PRINT_H
#define NEWLOGLIKFITTER_PRINT_H


///////////////////////////////////////////////////////////////////////////////
// add timestamp to ofstream object
///////////////////////////////////////////////////////////////////////////////

void timestamp(std::ostream& os)
{
    time_t rawtime;
    time(&rawtime);
    struct tm *info;
    info = localtime(&rawtime);
    char timestringbuf[100];
    size_t ncount = strftime(timestringbuf, sizeof(timestringbuf), "%F %T", info);
    std::string datetimestamp(timestringbuf);
    os << datetimestamp << std::endl;
}



///////////////////////////////////////////////////////////////////////////////
// print map to std::cout
///////////////////////////////////////////////////////////////////////////////

template<typename T, typename U>
void print_map(std::map<T, U> &map, const std::string& name)
{
    std::cout << "contents of map " << name << std::endl;
    for(auto it = map.cbegin(); it != map.cend(); ++ it)
    {
        std::cout << it->first << " -> " << it->second << std::endl;
    }
    
    /*
    std::cout << "contents of map paramNameToNumberMap:" << std::endl;
    for(auto it = paramNameToNumberMap.cbegin(); it != paramNameToNumberMap.cend(); ++ it)
    {
        std::cout << it->first << " -> " << it->second << std::endl;
    }
    */
}


///////////////////////////////////////////////////////////////////////////////
// print parameters map summary
///////////////////////////////////////////////////////////////////////////////

// TODO: change name


#if 0
void print_paramNameMap()
{

    TString sample_names[numberParams];
    // read samples names out again just to make sure its correct (useful when changing the number of params)
    for(int i = 0; i < numberParams; i++)
    {
        // tmpStr doesn't appear to do anything ?
        // TODO
        TString tmpStr = "";

        TString i_str;
        i_str.Form("%i", i);

        sample_names[i] = i_str;
        
        std::cout << "param " << i << ": sample_names[" << i << "]=" << sample_names[i] << " -> ";
        std::cout << paramNameMap[i] << " : ";
        bool first = true;
        for(auto it = MCNameToParamNameMap.cbegin(); it != MCNameToParamNameMap.cend(); ++ it)
        {
            if(it->second == paramNameMap[i])
            {
                // match for MC name
                if(!first) std::cout << ", ";
                std::cout << it->first;
                first = false;
            }
        }
        std::cout << std::endl;

        /*
        for(int j = 0; j < paramNameMap[i].size(); j++)
        { 
            std::cout << paramNameMap[i].at(j);
            tmpStr += paramNameMap[i].at(j);
            
            if(j + 1 < paramNameMap[i].size())
            {
                std::cout << ", ";
                tmpStr += ",";
            }
        }
        */
        
        std::cout << std::endl;
    }

}
#endif


#endif //NEWLOGLIKFITTER_PRINT_H
