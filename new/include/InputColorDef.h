#ifndef INPUTCOLORDEF_H
#define INPUTCOLORDEF_H


#include "InputNumberDef.h"


//Color_t Nd150Color = GetRGB(184, 67, 34); //kAzure + 1;
//Color_t tl208InternalBkgColor = GetRGB(196, 119, 55); //kGreen;
//Color_t bi214InternalBkgColor = GetRGB(209, 186, 88); //kBlue;
//Color_t bi207InternalBkgColor = GetRGB(218, 202, 127); //kViolet;
//Color_t InternalBkgColor = TColor::GetRGB(210, 231, 210); //kOrange; //for grouping
//Color_t NeighbourColor = TColor::GetRGB(153, 204, 223); //kCyan;
//Color_t RadonBkgColor = TColor::GetRGB(102, 149, 223); //kYellow;
//Color_t ExternalBkgColor = TColor::GetRGB(46, 70, 227); //kRed;
Color_t DataColor = kBlack;

Int_t ci_c0 = TColor::GetFreeColorIndex();
TColor *c_c0 = new TColor(ci_c0, 184.0 / 255.0, 67.0 / 255.0, 34.0 / 255.0);

Int_t ci_c1 = TColor::GetFreeColorIndex(); 
TColor *c_c1 = new TColor(ci_c1, 196.0 / 255.0, 119.0 / 255.0, 55.0 / 255.0);
//TColor *c_c1 = new TColor(ci_c1, 190.0 / 255.0, 92.0 / 255.0, 44.0 / 255.0);

Int_t ci_c2 = TColor::GetFreeColorIndex();
TColor *c_c2 = new TColor(ci_c2 , 209.0 / 255.0, 168.0 / 255.0, 88.0 / 255.0);
//TColor *c_c2 = new TColor(ci_c2 , 196.0 / 255.0, 119.0 / 255.0, 55.0 / 255.0);

Int_t ci_c3 = TColor::GetFreeColorIndex();
//TColor *c_c3 = new TColor(ci_c3 , 220.0 / 255.0, 223.0 / 255.0, 171.0 / 255.0);
TColor *c_c3 = new TColor(ci_c3 , 221.0 / 255.0, 215.0 / 255.0, 150.0 / 255.0);

// TODO: change above color it's grim

Int_t ci_c4 = TColor::GetFreeColorIndex();
//TColor *c_c4 = new TColor(ci_c4 , 210.0 / 255.0, 231.0 / 255.0, 210.0 / 255.0);
//TColor *c_c4 = new TColor(ci_c4 , 190.0 / 255.0, 223.0 / 255.0, 215.0 / 255.0);
TColor *c_c4 = new TColor(ci_c4 , 210.0 / 255.0, 231.0 / 255.0, 210.0 / 255.0);

Int_t ci_c5 = TColor::GetFreeColorIndex();
//TColor *c_c5 = new TColor(ci_c5, 153.0 / 255.0, 204.0 / 255.0, 223.0 / 255.0);
//TColor *c_c5 = new TColor(ci_c5, 123.0 / 255.0, 171.0 / 255.0, 228.0 / 255.0);
TColor *c_c5 = new TColor(ci_c5, 138.0 / 255.0, 187.0 / 255.0, 226.0 / 255.0);

Int_t ci_c6 = TColor::GetFreeColorIndex();
//TColor *c_c6 = new TColor(ci_c6, 102.0 / 255.0, 149.0 / 255.0, 223.0 / 255.0);
TColor *c_c6 = new TColor(ci_c6, 82.0 / 255.0, 129.0 / 255.0, 217.0 / 255.0);

Int_t ci_c7 = TColor::GetFreeColorIndex();
TColor *c_c7 = new TColor(ci_c7, 46.0 / 255.0, 70.0 / 255.0, 227.0 / 255.0);

Color_t Nd150Color = ci_c7;
Color_t tl208InternalBkgColor = ci_c6;
Color_t bi214InternalBkgColor = ci_c5;
Color_t bi207InternalBkgColor = ci_c4;
Color_t InternalBkgColor = ci_c3;
Color_t NeighbourColor = ci_c2;
Color_t RadonBkgColor = ci_c1;
Color_t ExternalBkgColor = ci_c0;




// 9
/*
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
*/

Color_t Rn222BkgColors[nRn222Bkgs] =
{
kYellow,    // Bi 214 air
kYellow,    // Pb 214 air (same color)
kGreen,     // Bi 214 SFoil
kGreen,     // Bi 214 SFoil
kGreen,     // Pb 214 SFoil (same color)
kGreen,     // Pb 214 SFoil (same color)
kGreen+1,   // Bi 210 SFoil
kCyan,      // Bi 214 SWire
kCyan,      // Bi 214 SWire
kCyan,      // Bi 214 SWire
kCyan,      // Pb 214 SWire (same color)
kCyan,      // Pb 214 SWire (same color)
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

Color_t Nd150Colors[nNd150Samples] =
{
kAzure + 1
};

Color_t Nd150SystematicColors[nNd150Systematics] =
{
kBlack,
kBlack,
kBlack,
kBlack,
kBlack,
kBlack,
kBlack
};




#endif
