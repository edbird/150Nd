#ifndef INPUTCOLORDEF_H
#define INPUTCOLORDEF_H


#include "InputNumberDef.h"


Color_t RadonBkgColor = 5;
Color_t ExternalBkgColor = 2;
Color_t InternalBkgColor = 6; //for grouping
Color_t NeighbourColor = 3;
Color_t Nd150Color = 4;
Color_t DataColor = kBlack;



// 9
Color_t Rn222BkgColors[nRn222Bkgs] =
{
kYellow,    // Bi 214 air
kYellow,    // Pb 214 air (same color)
kGreen,     // Bi 214 SFoil
kGreen,     // Pb 214 SFoil (same color)
kGreen+1,   // Bi 210 SFoil
kCyan,      // Bi 214 SWire
kCyan,      // Pb 214 SWire (same color)
kCyan+1,    // Bi 210 SWire
kGray,      // Bi 210 SScin
kBlue+3,    // Bi 214 Mylar
kBlue+3     // Pb 214 Mylar (same color)
};



// 11
Color_t ExternalBkgColors[nExternalBkgs] =
{
kYellow,    // Bi 214 Fe Shield
kYellow,    // Pb 214 Fe Shield
kYellow+1,  // Tl 208 Fe Shield

kOrange-6,  // Bi 214 PMT
kGreen+2,   // Ac 228 PMT
kGreen+2,   // Tl 208 PMT
kGreen+4,   // K 40 PMT

kTeal+5,   // K 40 Scint IN
kTeal+5,   // K 40 Scint OUT
kTeal+5,   // K 40 Scint PET

kCyan+4,    // Pa 234m SScin
kOrange+9//,  // Co 60 Cu Tower
//kYellow     // Neutrons
};



// 0 
Color_t Rn220BkgColors[nRn220Bkgs] =
{
kBlack      // Tl 208 Air
};


// 12
Color_t InternalBkgColors[nInternalBkgs] =
{
kAzure-2,   // Bi 214
kAzure-2,   // Pb 214 (same color)
kViolet,      // Ac 228
kViolet,    // Tl 208

kViolet,      // Bi 212 (same color as Ac 228)
kViolet+1,  // Bi 207
kPink+7,    // Eu 152
kMagenta+2, // Eu 154

kViolet-1,  // K 40
kBlue,    // Pa 234m
}; 


Color_t NeighbourColors[nNeighbours] =
{
kRed+2,     // 100 Mo
kRed+2,     // 100 Mo
kRed+2,     // 100 Mo
kRed+2,     // 100 Mo
kOrange+7,  // 48 Ca
kOrange+7,  // 48 Ca (will be disabled, same colror)
//kOrange+7,  // 48 Ca
kOrange+7,  // 48 Ca
kOrange+7,  // 48 Ca
kOrange,    // 96 Zr
kOrange,    // 96 Zr
kOrange,    // 96 Zr
kOrange     // 96 Zr
};

Color_t Nd150Colors[1] =
{
kAzure + 1
};




#endif
