#ifndef NEWLOGLIKFITTER_LOGPOISSON_H
#define NEWLOGLIKFITTER_LOGPOISSON_H


///////////////////////////////////////////////////////////////////////////////
// function declarations
///////////////////////////////////////////////////////////////////////////////



/*
double logpoisson(const double nData, const double nMC)
{
    double x = nData;
    double par = nMC;
    if(par < 0.0)
    {
        par = 0.0;
    }
    double ret = 0.0;
    if(x > 0.0)
    {
        ret = x * std::log(par) - par - TMath::LnGamma(x + 1.0);
    }
    else
    {
        ret = -par;
    }
    return ret;
}

*/


//double logpoisson_sterling(const double nData, const double nMC)
double logpoisson(const double nData, const double nMC)
{
    //std::cout << __func__ << "(nData=" << nData << ", nMC=" << nMC << ")" << std::endl;

    if(nData < 0.0) throw "error";

    double mnu = nMC;
    double dnu = nData;
    //std::cout << __func__ << "(dnu=" << dnu << ", mnu=" << mnu << ")" << std::endl;

    if(mnu < 0.0)
    {
        std::cout << "setting mnu=0" << std::endl;
        mnu = 0.0;
    }
    double ret = 0.0;
    if(dnu > 0.0)
    {   
        //std::cout << "dnu > 0.0" << std::endl;
        ret = -1.0 * (mnu - dnu + dnu * std::log(dnu / mnu));
        //std::cout << "ret=" << ret << std::endl;
    }
    else
    {
        //std::cout << "dnu ! > 0.0" << std::endl;
        //ret = -1.0 * (mnu - dnu);
        ret = -mnu;
        //std::cout << "ret=" << ret << std::endl;
    }
    //std::cout << "ll: nData=" << dnu << " nMC=" << mnu << " ll=" << ret << std::endl;
    return ret;
}


/*
double logpoisson(const double nData, const double nMC)
{
    double mnu = nMC;
    double dnu = nData;

    if(mnu < 0.0)
    {
        mnu = 0.0;
    }
    double ret = 0.0;
    if(dnu > 0.0)
    {   
        ret = TMath::Log(TMath::Poisson(nData, nMC));
    }
    else
    {
        //ret = -1.0 * (mnu - dnu);
        ret = TMath::Log(TMath::Poisson(nData, nMC));
    }
    //std::cout << "ll: nData=" << dnu << " nMC=" << mnu << " ll=" << ret << std::endl;
    return ret;
}
*/




#endif // NEWLOGLIKFITTER_LOGPOISSON_H
