#ifndef NEWLOGLIKFITTER_STACKER_HELPER_H
#define NEWLOGLIKFITTER_STACKER_HELPER_H

void stacker_helper(TH1F* &tmpHist_draw1D,
                    TH1F* &h_2nubb,
                    TH1F* &h_tl208_int,
                    TH1F* &h_bi214_int,
                    TH1F* &h_bi207_int,
                    TH1F* &h_internal,
                    TH1F* &h_neighbours,
                    TH1F* &h_radon,
                    TH1F* &h_external,
                    TH1F* &h_other)
{

    TString hname = tmpHist_draw1D->GetName();


    // 150 Nd
    if(hname.Contains("150"))
    {
        if(h_2nubb == nullptr)
        {
            h_2nubb = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_2nubb->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // tl208 internal
    else if(hname.Contains("tl208_int"))// ||
            //hname.Contains("ac228_int") ||
            //hname.Contains("bi212_int"))
    {
        if(h_tl208_int == nullptr)
        {
            h_tl208_int = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_tl208_int->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // bi 214 internal
    else if(hname.Contains("bi214_int") ||
            hname.Contains("pb214_int")
            )
    {
        //stacks1D_bi214_int->Add((TH1F*)tmpHist_draw1D->Clone());
        if(h_bi214_int == nullptr)
        {
            h_bi214_int = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_bi214_int->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // bi 207 internal
    else if(hname.Contains("bi207_int"))
    {
        //stacks1D_bi207_int->Add((TH1F*)tmpHist_draw1D->Clone());
        if(h_bi207_int == nullptr)
        {
            h_bi207_int = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_bi207_int->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // internals (other)
    else if(
        hname.Contains("int") ||
        hname.Contains("mylar")
        )
    {
        //stacks1D_internal->Add((TH1F*)tmpHist_draw1D->Clone());
        if(h_internal == nullptr)
        {
            h_internal = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_internal->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // neighbours
    else if(
        hname.Contains("mo100") ||
        hname.Contains("zr96") ||
        hname.Contains("ca48")
        )
    {
        //stacks1D_neighbours->Add((TH1F*)tmpHist_draw1D->Clone());
        if(h_neighbours == nullptr)
        {
            h_neighbours = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_neighbours->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // radon
    else if(
        hname.Contains("swire") ||
        hname.Contains("sfoil")
        )
    {
        //stacks1D_radon->Add((TH1F*)tmpHist_draw1D->Clone());
        if(h_radon == nullptr)
        {
            h_radon = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_radon->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // external
    else if(
        hname.Contains("feShield") ||
        hname.Contains("pmt") ||
        hname.Contains("cuTower") ||
        hname.Contains("sscin") ||
        hname.Contains("scint") ||
        hname.Contains("air") // TODO: air now included with externals
        )
    {
        //stacks1D_external->Add((TH1F*)tmpHist_draw1D->Clone());
        if(h_external == nullptr)
        {
            h_external = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_external->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }
    // everything else
    else
    {
        std::cout << "adding " << tmpHist_draw1D->GetName() << " to others TODO: FIX" << std::endl;
        //stacks1D_other->Add((TH1F*)tmpHist_draw1D->Clone());
        if(h_other == nullptr)
        {
            h_other = (TH1F*)tmpHist_draw1D->Clone();
        }
        else
        {
            h_other->Add((TH1F*)tmpHist_draw1D->Clone());
        }
    }



}



void stacker_helper_2(THStack* &stacks1D_major,
                      TH1F* &h_2nubb,
                      TH1F* &h_tl208_int,
                      TH1F* &h_bi214_int,
                      TH1F* &h_bi207_int,
                      TH1F* &h_internal,
                      TH1F* &h_neighbours,
                      TH1F* &h_radon,
                      TH1F* &h_external,
                      TH1F* &h_other)
{

    if(h_external != nullptr)
    {
        //h_external->SetLineWidth(1);
        //h_external->SetLineColor(kBlack);
        h_external->SetFillColor(ExternalBkgColor);
        h_external->SetLineStyle(0);
        h_external->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_external);
    }
    if(h_radon != nullptr)
    {
        //h_radon->SetLineWidth(1);
        //h_radon->SetLineColor(kBlack);
        h_radon->SetFillColor(RadonBkgColor);
        h_radon->SetLineStyle(0);
        h_radon->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_radon);
    }
    if(h_neighbours != nullptr)
    {
        h_neighbours->SetFillColor(NeighbourColor);
        h_neighbours->SetLineStyle(0);
        h_neighbours->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_neighbours);
    }
    if(h_internal != nullptr)
    {
        //h_internal->SetLineWidth(1);
        //h_internal->SetLineColor(kBlack);
        h_internal->SetFillColor(InternalBkgColor);
        h_internal->SetLineStyle(0);
        h_internal->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_internal);
    }
    if(h_bi207_int != nullptr)
    {
        h_bi207_int->SetFillColor(bi207InternalBkgColor);
        h_bi207_int->SetLineStyle(0);
        h_bi207_int->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_bi207_int);
    }
    if(h_bi214_int != nullptr)
    {
        h_bi214_int->SetFillColor(bi214InternalBkgColor);
        h_bi214_int->SetLineStyle(0);
        h_bi214_int->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_bi214_int);
    }
    if(h_tl208_int != nullptr)
    {
        h_tl208_int->SetFillColor(tl208InternalBkgColor);
        h_tl208_int->SetLineStyle(0);
        h_tl208_int->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_tl208_int);
    }
    if(h_2nubb != nullptr)
    {
        h_2nubb->SetFillColor(Nd150Color);
        h_2nubb->SetLineStyle(0);
        h_2nubb->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_2nubb);
    }
    if(h_other != nullptr)
    {
        h_other->SetLineStyle(0);
        h_other->SetLineWidth(0);
        stacks1D_major->Add((TH1F*)h_other);
        std::cout << "h_other is non zero" << std::endl;
    }

}

#endif //NEWLOGLIKFITTER_STACKER_HELPER_H
