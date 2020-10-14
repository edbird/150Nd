#ifndef NEWLOGLIKFITTER_STACKER_HELPER_H
#define NEWLOGLIKFITTER_STACKER_HELPER_H


void new_stacker_helper
(
    TH1D *tmpHist1D,
    TH1D *h_stack_external,
    TH1D *h_stack_radon,
    TH1D *h_stack_neighbour,
    TH1D *h_stack_internal,
    TH1D *h_stack_bi207_int,
    TH1D *h_stack_bi214_int,
    TH1D *h_stack_tl208_int,
    TH1D *h_stack_2nubb,
    TH1D *h_stack_total_MC,
    const double scale_factor
)
{

    //std::cout << __func__ << std::endl;

    TString hname = tmpHist1D->GetName();


    // 150 Nd
    if(hname.Contains("150"))
    {
        h_stack_2nubb->Add(tmpHist1D, scale_factor);
        h_stack_total_MC->Add(tmpHist1D, scale_factor);
    }
    // tl208 internal
    else if(hname.Contains("tl208_int"))// ||
            //hname.Contains("ac228_int") ||
            //hname.Contains("bi212_int"))
    {
        h_stack_tl208_int->Add(tmpHist1D, scale_factor);
        h_stack_2nubb->Add(tmpHist1D, scale_factor);
        h_stack_total_MC->Add(tmpHist1D, scale_factor);
    }
    // bi 214 internal
    else if(hname.Contains("bi214_int") ||
            hname.Contains("pb214_int")
            )
    {
        h_stack_bi214_int->Add(tmpHist1D, scale_factor);
        h_stack_tl208_int->Add(tmpHist1D, scale_factor);
        h_stack_2nubb->Add(tmpHist1D, scale_factor);
        h_stack_total_MC->Add(tmpHist1D, scale_factor);
    }
    // bi 207 internal
    else if(hname.Contains("bi207_int"))
    {
        h_stack_bi207_int->Add(tmpHist1D, scale_factor);    // 1
        h_stack_bi214_int->Add(tmpHist1D, scale_factor);
        h_stack_tl208_int->Add(tmpHist1D, scale_factor);    // 3
        h_stack_2nubb->Add(tmpHist1D, scale_factor);
        h_stack_total_MC->Add(tmpHist1D, scale_factor);     // 5
    }
    // internals (other)
    else if(
        hname.Contains("int") ||
        hname.Contains("mylar")
        )
    {
        h_stack_internal->Add(tmpHist1D, scale_factor);     // 1
        h_stack_bi207_int->Add(tmpHist1D, scale_factor);
        h_stack_bi214_int->Add(tmpHist1D, scale_factor);    // 3
        h_stack_tl208_int->Add(tmpHist1D, scale_factor);
        h_stack_2nubb->Add(tmpHist1D, scale_factor);        // 5
        h_stack_total_MC->Add(tmpHist1D, scale_factor);
    }
    // neighbours
    else if(
        hname.Contains("mo100") ||
        hname.Contains("zr96") ||
        hname.Contains("ca48")
        )
    {
        h_stack_neighbour->Add(tmpHist1D, scale_factor);    // 1
        h_stack_internal->Add(tmpHist1D, scale_factor);
        h_stack_bi207_int->Add(tmpHist1D, scale_factor);    // 3
        h_stack_bi214_int->Add(tmpHist1D, scale_factor);
        h_stack_tl208_int->Add(tmpHist1D, scale_factor);    // 5
        h_stack_2nubb->Add(tmpHist1D, scale_factor);
        h_stack_total_MC->Add(tmpHist1D, scale_factor);     // 7
    }
    // radon
    else if(
        hname.Contains("swire") ||
        hname.Contains("sfoil")
        )
    {
        h_stack_radon->Add(tmpHist1D, scale_factor);        // 1
        h_stack_neighbour->Add(tmpHist1D, scale_factor);
        h_stack_internal->Add(tmpHist1D, scale_factor);     // 3
        h_stack_bi207_int->Add(tmpHist1D, scale_factor);
        h_stack_bi214_int->Add(tmpHist1D, scale_factor);    // 5
        h_stack_tl208_int->Add(tmpHist1D, scale_factor);
        h_stack_2nubb->Add(tmpHist1D, scale_factor);        // 7
        h_stack_total_MC->Add(tmpHist1D, scale_factor);
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
        h_stack_external->Add(tmpHist1D, scale_factor);     // 1
        h_stack_radon->Add(tmpHist1D, scale_factor);
        h_stack_neighbour->Add(tmpHist1D, scale_factor);    // 3
        h_stack_internal->Add(tmpHist1D, scale_factor);
        h_stack_bi207_int->Add(tmpHist1D, scale_factor);    // 5
        h_stack_bi214_int->Add(tmpHist1D, scale_factor);
        h_stack_tl208_int->Add(tmpHist1D, scale_factor);    // 7
        h_stack_2nubb->Add(tmpHist1D, scale_factor);
        h_stack_total_MC->Add(tmpHist1D, scale_factor);     // 9
    }
    // everything else
    else
    {
        std::cout << "adding " << tmpHist1D->GetName() << " to others TODO: FIX" << std::endl;
    }


    //std::cout << __func__ << " exit" <<std::endl;

}










void stacker_helper(TH1D* &tmpHist_draw1D,
                    TH1D* &h_2nubb,
                    TH1D* &h_tl208_int,
                    TH1D* &h_bi214_int,
                    TH1D* &h_bi207_int,
                    TH1D* &h_internal,
                    TH1D* &h_neighbours,
                    TH1D* &h_radon,
                    TH1D* &h_external,
                    TH1D* &h_other)
{

    TString hname = tmpHist_draw1D->GetName();


    // 150 Nd
    if(hname.Contains("150"))
    {
        if(h_2nubb == nullptr)
        {
            h_2nubb = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_2nubb->Add((TH1D*)tmpHist_draw1D->Clone());
            h_2nubb->Add((TH1D*)tmpHist_draw1D);
        }
    }
    // tl208 internal
    else if(hname.Contains("tl208_int"))// ||
            //hname.Contains("ac228_int") ||
            //hname.Contains("bi212_int"))
    {
        if(h_tl208_int == nullptr)
        {
            h_tl208_int = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_tl208_int->Add((TH1D*)tmpHist_draw1D->Clone());
            h_tl208_int->Add((TH1D*)tmpHist_draw1D);
        }
    }
    // bi 214 internal
    else if(hname.Contains("bi214_int") ||
            hname.Contains("pb214_int")
            )
    {
        //stacks1D_bi214_int->Add((TH1D*)tmpHist_draw1D->Clone());
        if(h_bi214_int == nullptr)
        {
            h_bi214_int = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_bi214_int->Add((TH1D*)tmpHist_draw1D->Clone());
            h_bi214_int->Add((TH1D*)tmpHist_draw1D);
        }
    }
    // bi 207 internal
    else if(hname.Contains("bi207_int"))
    {
        //stacks1D_bi207_int->Add((TH1D*)tmpHist_draw1D->Clone());
        if(h_bi207_int == nullptr)
        {
            h_bi207_int = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_bi207_int->Add((TH1D*)tmpHist_draw1D->Clone());
            h_bi207_int->Add((TH1D*)tmpHist_draw1D);
        }
    }
    // internals (other)
    else if(
        hname.Contains("int") ||
        hname.Contains("mylar")
        )
    {
        //stacks1D_internal->Add((TH1D*)tmpHist_draw1D->Clone());
        if(h_internal == nullptr)
        {
            h_internal = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_internal->Add((TH1D*)tmpHist_draw1D->Clone());
            h_internal->Add((TH1D*)tmpHist_draw1D);
        }
    }
    // neighbours
    else if(
        hname.Contains("mo100") ||
        hname.Contains("zr96") ||
        hname.Contains("ca48")
        )
    {
        //stacks1D_neighbours->Add((TH1D*)tmpHist_draw1D->Clone());
        if(h_neighbours == nullptr)
        {
            h_neighbours = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_neighbours->Add((TH1D*)tmpHist_draw1D->Clone());
            h_neighbours->Add((TH1D*)tmpHist_draw1D);
        }
    }
    // radon
    else if(
        hname.Contains("swire") ||
        hname.Contains("sfoil")
        )
    {
        //stacks1D_radon->Add((TH1D*)tmpHist_draw1D->Clone());
        if(h_radon == nullptr)
        {
            h_radon = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_radon->Add((TH1D*)tmpHist_draw1D->Clone());
            h_radon->Add((TH1D*)tmpHist_draw1D);
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
        //stacks1D_external->Add((TH1D*)tmpHist_draw1D->Clone());
        if(h_external == nullptr)
        {
            h_external = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_external->Add((TH1D*)tmpHist_draw1D->Clone());
            h_external->Add((TH1D*)tmpHist_draw1D);
        }
    }
    // everything else
    else
    {
        std::cout << "adding " << tmpHist_draw1D->GetName() << " to others TODO: FIX" << std::endl;
        //stacks1D_other->Add((TH1D*)tmpHist_draw1D->Clone());
        if(h_other == nullptr)
        {
            h_other = (TH1D*)tmpHist_draw1D->Clone();
        }
        else
        {
            //h_other->Add((TH1D*)tmpHist_draw1D->Clone());
            h_other->Add((TH1D*)tmpHist_draw1D);
        }
    }



}



void stacker_helper_2(THStack* &stacks1D_major,
                      TH1D* &h_2nubb,
                      TH1D* &h_tl208_int,
                      TH1D* &h_bi214_int,
                      TH1D* &h_bi207_int,
                      TH1D* &h_internal,
                      TH1D* &h_neighbours,
                      TH1D* &h_radon,
                      TH1D* &h_external,
                      TH1D* &h_other)
{

    if(h_external != nullptr)
    {
        //h_external->SetLineWidth(1);
        //h_external->SetLineColor(kBlack);
        h_external->SetFillColor(ExternalBkgColor);
        h_external->SetLineStyle(0);
        h_external->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_external);
    }
    if(h_radon != nullptr)
    {
        //h_radon->SetLineWidth(1);
        //h_radon->SetLineColor(kBlack);
        h_radon->SetFillColor(RadonBkgColor);
        h_radon->SetLineStyle(0);
        h_radon->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_radon);
    }
    if(h_neighbours != nullptr)
    {
        h_neighbours->SetFillColor(NeighbourColor);
        h_neighbours->SetLineStyle(0);
        h_neighbours->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_neighbours);
    }
    if(h_internal != nullptr)
    {
        //h_internal->SetLineWidth(1);
        //h_internal->SetLineColor(kBlack);
        h_internal->SetFillColor(InternalBkgColor);
        h_internal->SetLineStyle(0);
        h_internal->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_internal);
    }
    if(h_bi207_int != nullptr)
    {
        h_bi207_int->SetFillColor(bi207InternalBkgColor);
        h_bi207_int->SetLineStyle(0);
        h_bi207_int->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_bi207_int);
    }
    if(h_bi214_int != nullptr)
    {
        h_bi214_int->SetFillColor(bi214InternalBkgColor);
        h_bi214_int->SetLineStyle(0);
        h_bi214_int->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_bi214_int);
    }
    if(h_tl208_int != nullptr)
    {
        h_tl208_int->SetFillColor(tl208InternalBkgColor);
        h_tl208_int->SetLineStyle(0);
        h_tl208_int->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_tl208_int);
    }
    if(h_2nubb != nullptr)
    {
        h_2nubb->SetFillColor(Nd150Color);
        h_2nubb->SetLineStyle(0);
        h_2nubb->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_2nubb);
    }
    if(h_other != nullptr)
    {
        h_other->SetLineStyle(0);
        h_other->SetLineWidth(0);
        stacks1D_major->Add((TH1D*)h_other);
        std::cout << "h_other is non zero" << std::endl;
    }

}

#endif //NEWLOGLIKFITTER_STACKER_HELPER_H
