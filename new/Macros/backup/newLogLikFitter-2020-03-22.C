#include "TH1F.h"
#include "TArray.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"

#include "THStack.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <string>
#include "TLegend.h"
#include <vector>
#include <TMinuit.h>

#include "newLogLikFitter.h"

//-----------------------------------------
//    Functions
//----------------------------------------
void loadFiles();

void  book1DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);
void  book2DHistograms(Int_t channel_counter, TString theChannel,TString thePhase, TString theHistogram);
void logLikelihood(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */); 
Double_t getNumberMC1D(Int_t channel, Int_t binx, Double_t *p);
Double_t getNumberMC2D(Int_t channel, Int_t binx, Int_t biny, Double_t *p);

void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err);

void loadFiles() {

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    allDataSamples1D = new TObjArray();
    allDataSamples2D = new TObjArray();

    // First, create a root file to hold all of the histograms
    TFile *myFile = TFile::Open("Nd150_loglikResults.root","RECREATE");
    myFile->Close();

    // Read in all the histograms we want
    //book2DHistograms(0,"1e1g_","P1","hEeVEg_");
    //book2DHistograms(0,"1e1g_","P2","hEeVEg_");
    //book1DHistograms(0,"1e1g_","P1","hClusterEnergy_");
    //book1DHistograms(1,"1e1g_","P2","hClusterEnergy_");
    
    //book1DHistograms(0,"1e_","P1","hEe_");
    //book1DHistograms(0,"1e_","P2","hEe_");
    //book1DHistograms(2,"1e1g_","P1","hTotalE_");
    //book1DHistograms(3,"1e1g_","P2","hTotalE_");
    //book1DHistograms(4,"1e2g_","P1","hClusterEnergyMax_");
    //book1DHistograms(5,"1e2g_","P2","hClusterEnergyMax_");
    //book1DHistograms(2,"1e1gExt_","P1","hClusterEnergy_");
    //book1DHistograms(3,"1e1gExt_","P2","hClusterEnergy_");
    //book1DHistograms(6,"2e_","P1","hTotalE_");
    //book1DHistograms(7,"2e_","P2","hTotalE_");
    //book1DHistograms(4,"OCE_","P1","hTotalE_");
    //book1DHistograms(5,"OCE_","P2","hTotalE_");
    book1DHistograms(0, "2e_", "P1", "hTotalE_");
  
    // Array to hold activity adjustment parameters
    Double_t AdjustActs[numberParams];
    Double_t AdjustActs_Err[numberParams];
    for(Int_t ix{0}; ix < numberParams; ++ ix)
    {
        AdjustActs[ix] = 1.0;
        AdjustActs_Err[ix] = 1.0;
    }
  
    fitBackgrounds(AdjustActs, AdjustActs_Err);

    std::cout << "The following adjustments should be made:" << std::endl;
    for(int i = 0; i < numberParams; i++)
    {
        std::cout << i << " :\t" << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;  
        std::ofstream myFile("fit_results.txt", std::ios::out | std::ios::app);
        myFile << AdjustActs[i] << " +- " << AdjustActs_Err[i] << std::endl;

        myFile.close();
    }
  
    //Now let's draw the results...

    THStack *stacks1D[number1DHists];
    TCanvas *c;
    TH1F *data1D[number1DHists];

    // each channel 1D hists
    std::cout << "debug: " << allDataSamples1D->GetEntries() << std::endl;
    for(int i = 0; i < allDataSamples1D->GetEntries(); i++)
    {
        data1D[i] = (TH1F*)allDataSamples1D->At(i)->Clone();
        TString i_str;
        i_str.Form("%i", i);
        c = new TCanvas("c" + i_str);
        c->SetFillColor(kWhite);
        stacks1D[i] = new THStack("stacks1D" + i_str, i_str);

        TH1F *tmpHist_draw1D;

        // j list MC samples for this channel i
        for(int j = 0; j < allMCSamples1D[i]->GetEntries(); j++)
        {
            TString j_str;
            j_str.Form("%i", j);
            tmpHist_draw1D = (TH1F*)allMCSamples1D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw1D->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            int which_param;
            bool foundParam = false;

            // k searching through array of params for the right one
            for(int k = 0; (k < numberParams) && !foundParam; k++)
            {

	            // match up the isotope with parm
                for(int n = 0; (n < paramNameMap[k].size()) && !foundParam; n++)
                {
	                //std::cout <<"is it...  " <<paramNameMap[k].at(j) << std::endl;

	                if(tmpName.Contains(paramNameMap[k].at(n)))
                    {
	                    foundParam = true;
                        which_param = k;
	                }
	            }
            }

            if(foundParam)
            {
                tmpHist_draw1D->Scale(AdjustActs[which_param]);
                if(tmpHist_draw1D->Integral() > 0)
                {
                    stacks1D[i]->Add(tmpHist_draw1D);
	            }
                std::cout << "found histogram: tmpName=" << tmpName << std::endl;
            }
            else
            {
                std::cout << "error could not find histograms" << std::endl;
                std::cout << "tmpName=" << tmpName << std::endl;
            } 

        }


        stacks1D[i]->Draw("hist");
        data1D[i]->SetLineWidth(2);
        data1D[i]->SetMarkerStyle(20);
        data1D[i]->SetMarkerSize(0.5);
        data1D[i]->Draw("PEsames");

        c->SaveAs("finalHisto1D_" + i_str + ".C");
        c->SaveAs("finalHisto1D_" + i_str + ".eps");
        c->SaveAs("finalHisto1D_" + i_str + ".png");
        std::cout << "saving to file (canvas)" << std::endl;
    
    }

    // each channel 1D hists
    THStack *stacks2D[number2DHists];
    TH2F *data2D[number2DHists];
    for(int i = 0; i < allDataSamples2D->GetEntries(); i++)
    {
        data2D[i] = (TH2F*)allDataSamples2D->At(i)->Clone();
        TString i_str;
        i_str.Form("%i",i);
        c = new TCanvas("c_"+i_str);
        c->SetFillColor(kWhite);
        stacks2D[i] = new THStack("stacks2D"+i_str,i_str);

        // j list MC samples for this channel i
        TH2F *tmpHist_draw2D;
        for(int j = 0; j < allMCSamples2D[i]->GetEntries(); j++)
        {
            TString j_str;
            j_str.Form("%i",j);
            tmpHist_draw2D = (TH2F*)allMCSamples2D[i]->At(j)->Clone();
            TString tmpName = tmpHist_draw2D->GetName();

            //std::cout << "looking for " << tmpName << std::endl;
            
            // k searching through array of params for the right one
            int which_param;
            bool foundParam = false;
            for(int k = 0; (k < numberParams) && !foundParam; k++)
            {
            
                // match up the isotope with parm
                for(int n = 0; (n < paramNameMap[k].size()) && !foundParam; n++)
                {
	                //std::cout <<"is it...  " <<paramNameMap[k].at(j) << std::endl;

	                if(tmpName.Contains(paramNameMap[k].at(n)))
                    {
	                    foundParam = true;
                        which_param = k;
	                }
	            }
            }

            if(foundParam)
            {
                tmpHist_draw2D->Scale(AdjustActs[which_param]);
     
                if(tmpHist_draw2D->Integral() > 0)
                {
                    stacks2D[i]->Add(tmpHist_draw2D);
                }
            }
            else
            {
                std::cout << "error could not find histograms" << std::endl;
            }

        }

        stacks2D[i]->Draw("");
        data2D[i]->SetLineWidth(2);
        data2D[i]->SetMarkerStyle(20);
        data2D[i]->SetMarkerSize(0.5);
        data2D[i]->Draw("PEsames");

        c->SaveAs("finalHisto2D_"+i_str+".C");
        c->SaveAs("finalHisto2D_"+i_str+".eps");
        c->SaveAs("finalHisto2D_"+i_str+".png");
        std::cout << "saved 2D histogram to file" << std::endl;

    }

}

void book1DHistograms(Int_t channel_counter, TString theChannel, TString thePhase, TString theHistogram) {

    std::cout << "booking 1D hists for " << theChannel << " " << thePhase << std::endl;
    allMCSamples1D[channel_counter] = new TObjArray();

    TFile *aFile = TFile::Open("/home/ebirdsall/NEMO3/Nd150_analysis/MeasureStuff/new/Macros/Nd150_"+theChannel+thePhase+".root");
    gDirectory->cd("singleHistos");
    //gDirectory->ls();
    TH1F *tmpHist  = new TH1F("tmpHist_" + theChannel + thePhase, "" , 1, 0, 1);
    for(int i = 0; i < nExternalBkgs_oneParam; i++)
    {
        if(gDirectory->GetListOfKeys()->Contains(theHistogram + ExternalBkgOneParamFiles[i] + "_fit"))
        {
            //check if the histograms exists
            tmpHist = (TH1F*)gDirectory->Get(theHistogram+ExternalBkgOneParamFiles[i]+"_fit")->Clone(ExternalBkgOneParamFiles[i]+"_"+theChannel+thePhase);
            allMCSamples1D[channel_counter]->Add(tmpHist);

            //std::cout << tmpHist->GetName() << std::endl;

        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + ExternalBkgOneParamFiles[i] + "_fit" << std::endl;
        }
    }

    for(int i = 0; i < nInternalBkgs; i++)
    {
        if(gDirectory->GetListOfKeys()->Contains(theHistogram+InternalBkgFiles[i]+"_fit"))
        {
            //check if the histograms exists
            tmpHist = (TH1F*)gDirectory->Get(theHistogram + InternalBkgFiles[i] + "_fit")->Clone(InternalBkgFiles[i] + "_" + theChannel + thePhase);
            allMCSamples1D[channel_counter]->Add(tmpHist);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + InternalBkgFiles[i] + "_fit" << std::endl;
        }
    }

    for(int i = 0; i < nRn222Bkgs; i++)
    {
        if(gDirectory->GetListOfKeys()->Contains(theHistogram + Rn222BkgFiles[i] + "_fit"))
        {
            //check if the histograms exists
            tmpHist = (TH1F*)gDirectory->Get(theHistogram + Rn222BkgFiles[i] + "_fit")->Clone(Rn222BkgFiles[i] + "_" + theChannel + thePhase);
            allMCSamples1D[channel_counter]->Add(tmpHist);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + Rn222BkgFiles[i] + "_fit" << std::endl;
        }
    }

    for(int i = 0; i < nRn220Bkgs; i++)
    {
        if(gDirectory->GetListOfKeys()->Contains(theHistogram + Rn220BkgFiles[i] + "_fit"))
        {
            //check if the histograms exists
            tmpHist = (TH1F*)gDirectory->Get(theHistogram + Rn220BkgFiles[i] + "_fit")->Clone(Rn220BkgFiles[i] + "_" + theChannel + thePhase);
            allMCSamples1D[channel_counter]->Add(tmpHist);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + Rn220BkgFiles[i] + "_fit" << std::endl;
        }
    }

    for(int i = 0; i < nNd150Samples; i++)
    {
        if(gDirectory->GetListOfKeys()->Contains(theHistogram + Nd150Files[i] + "_fit"))
        {
            //check if the histograms exists
            tmpHist = (TH1F*)gDirectory->Get(theHistogram + Nd150Files[i] + "_fit")->Clone(Nd150Files[i] + "_" + theChannel + thePhase);
            allMCSamples1D[channel_counter]->Add(tmpHist);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + Nd150Files[i] + "_fit" << std::endl;
        }
    }

    for(int i = 0; i < nNeighbours; i++)
    {
        if(gDirectory->GetListOfKeys()->Contains(theHistogram + NeighbourFiles[i] + "_fit"))
        {
            //check if the histograms exists
            tmpHist = (TH1F*)gDirectory->Get(theHistogram + NeighbourFiles[i] + "_fit")->Clone(NeighbourFiles[i] + "_" + theChannel + thePhase);
            allMCSamples1D[channel_counter]->Add(tmpHist);
        }
        else
        {
            std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + NeighbourFiles[i] + "_fit" << std::endl;
        }
    }

    if(gDirectory->GetListOfKeys()->Contains(theHistogram + "Data"))
    {
        tmpHist = (TH1F*)gDirectory->Get(theHistogram + "Data")->Clone("Data_" + theChannel + thePhase);
        allDataSamples1D->Add((TH1F*)tmpHist);
    }
    else
    {
        std::cout << "gDirectory->GetListOfKeys() does not contain " << theHistogram + "Data" << std::endl;
    }

    // std::cout << tmpHist->GetName() << std::endl;
    // tmpHist->Delete();
    // aFile->Close();
    // aFile->Delete();
}

void book2DHistograms(Int_t channel_counter, TString theChannel, TString thePhase, TString theHistogram) {

  std::cout << "booking 2D hists for "<< theChannel<<" " <<thePhase<<std::endl;
  allMCSamples2D[channel_counter] = new TObjArray();

  TFile *aFile = TFile::Open("/home/blotsd/NEMO3/Nd150_analysis/MeasureStuff/Macros/Nd150_"+theChannel+thePhase+".root");
  gDirectory->cd("singleHistos");
  // gDirectory->ls();
  TH2F *tmpHist2  = new TH2F("tmpHist2_"+theChannel+thePhase,"",1,0,1,1,0,1);
  for ( int i = 0; i < nExternalBkgs_oneParam; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+ExternalBkgOneParamFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+ExternalBkgOneParamFiles[i]+"_fit")->Clone(ExternalBkgOneParamFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);

      //std::cout << tmpHist->GetName() << std::endl;

    }
  }
  /*
  for ( int i = 0; i < nExternalBkgs_twoParam; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+ExternalBkgTwoParamFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+ExternalBkgTwoParamFiles[i]+"_fit")->Clone(ExternalBkgTwoParamFiles[i]+"_"+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);

      //std::cout << tmpHist->GetName() << std::endl;

    }
  }*/

  for ( int i = 0; i < nInternalBkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+InternalBkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+InternalBkgFiles[i]+"_fit")->Clone(InternalBkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nRn222Bkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Rn222BkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Rn222BkgFiles[i]+"_fit")->Clone(Rn222BkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nRn220Bkgs; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Rn220BkgFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Rn220BkgFiles[i]+"_fit")->Clone(Rn220BkgFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nNd150Samples; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+Nd150Files[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+Nd150Files[i]+"_fit")->Clone(Nd150Files[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  for ( int i = 0; i < nNeighbours; i++ ) {
    if ( gDirectory->GetListOfKeys()->Contains(theHistogram+NeighbourFiles[i]+"_fit") ) {//check if the histograms exists
      tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+NeighbourFiles[i]+"_fit")->Clone(NeighbourFiles[i]+"_"+theChannel+thePhase);
      allMCSamples2D[channel_counter]->Add(tmpHist2);
    }
  }

  if ( gDirectory->GetListOfKeys()->Contains(theHistogram+"Data") ) {
    tmpHist2 = (TH2F*)gDirectory->Get(theHistogram+"Data")->Clone("Data_"+theChannel+thePhase);
    allDataSamples2D->Add((TH2F*)tmpHist2);
  }

  std::cout << tmpHist2->GetName() << std::endl;
  // tmpHist2->Delete();
  //aFile->Close();
  // aFile->Delete();
}


void fitBackgrounds(double *AdjustActs, double *AdjustActs_Err)
{

    //TVirtualFitter::SetDefaultFitter("Minuit2");
    TMinuit * minuit = new TMinuit(numberParams);

    std::cout << "Fit created" << std::endl;

    for(int i = 0; i < numberParams; i++)
    {
        // paramNameMap[i] = std::vector<TString>();
        //paramActMap[i] = std::vector<double>();
        //paramActErrMap[i] = std::vector<double>();

        paramNameMap[i].clear();
        paramActMap[i].clear();
        paramActErrMap[i].clear();
    }

    std::ifstream paramFile;
    paramFile.open("parameter_names.lst");
    while(!paramFile.eof())
    {
        TString name;
        int dummy;
        int paramNumber;
        double activity;
        double activity_error;
        TString fixed_str;

        std::stringstream ss;
        std::string s;
        std::getline(paramFile, s);
        
        if(s.size() > 0)
        {

            ss << s;

            ss >> dummy >> paramNumber >> name >> activity >> activity_error >> fixed_str;
            std::cout <<  paramNumber << "\t" << name << "\t" << activity << "\t" << activity_error << " fixed_str=" << fixed_str << std::endl;

        
            if(fixed_str.CompareTo("fixed") == 0)
            {
                //std::cout << "detected fixed parameter" << std::endl;
                fixed_params.push_back(paramNumber);
            }
            else if(fixed_str.CompareTo("notfixed") == 0 || fixed_str.CompareTo("free") == 0)
            {
                //std::cout << "detected free parameter" << std::endl;
            }
            else
            {
                std::cout << "unknown fixed/free parameter specification" << std::endl;
                std::cout << "fixed_str=" << fixed_str << std::endl;
            }
            // TODO

            paramNameMap[paramNumber].push_back(name);
            paramActMap[paramNumber].push_back(activity);
            paramActErrMap[paramNumber].push_back(activity_error);
        }

    }

  
    TString sample_names[numberParams];
    // read samples names out again just to make sure its correct (useful when changing the number of params)
    for(int i = 0; i < numberParams; i++)
    {
        TString tmpStr = "";
        std::cout << "param " << i << " : ";

        TString i_str;
        i_str.Form("%i",i);

        sample_names[i] = i_str;
        for(int j = 0; j < paramNameMap[i].size(); j++)
        { 
            std::cout << paramNameMap[i].at(j) << ", ";
            tmpStr += paramNameMap[i].at(j)+",";
        }
        
        std::cout << std::endl;
    }

    // Set the parameters for the fit, and give them an arbitrary 10% error to start.
    for(int i = 0; i < numberParams; i++)
    {
        // AdjustActs[i] = 1.;
        //AdjustActs_Err[i] = 1.;
        TString i_str;
        i_str.Form("%i",i);
        minuit->DefineParameter(i, "_" + i_str + "_", 1.0, 0.1, 0.0, 1000.0);

        if(std::find(fixed_params.begin(), fixed_params.end(), i) != fixed_params.end())
        {
            minuit->FixParameter(i);
        }

        //if ((i != 17 ) && (i < 32) || (i > 38)  ) {
        //if ( (i > 32) ) {
        //  minuit->FixParameter(i);
        //  fixed_params.push_back(i);
        //}

        //if ( i == 24 ) {
        //  minuit->FixParameter(i);
        //  fixed_params.push_back(i);
        //}
    }
  

    // Fix the tracker radon activities
    //minuit->FixParameter(34);
    //minuit->FixParameter(35);
    //minuit->FixParameter(23);
    // minuit->FixParameter(24);
    /*
    minuit->FixParameter(1);
    minuit->FixParameter(3);
    minuit->FixParameter(7);
    minuit->FixParameter(11);
    minuit->FixParameter(21);
    minuit->FixParameter(22);
    minuit->FixParameter(29);
    minuit->FixParameter(30);
    //minuit->FixParameter(31);
    minuit->FixParameter(34);
    minuit->FixParameter(35);
    //minuit->FixParameter(48);
    minuit->FixParameter(49);
    minuit->FixParameter(50);
    minuit->FixParameter(51);
    minuit->FixParameter(52);
    minuit->FixParameter(53);
    minuit->FixParameter(54);

    fixed_params.push_back(1);
    fixed_params.push_back(3);
    fixed_params.push_back(7);
    fixed_params.push_back(11);
    fixed_params.push_back(21);
    fixed_params.push_back(22);
    fixed_params.push_back(29);
    fixed_params.push_back(30);
    //fixed_params.push_back(31);
    fixed_params.push_back(34);
    fixed_params.push_back(35);
    fixed_params.push_back(49);
    fixed_params.push_back(50);
    fixed_params.push_back(51);
    fixed_params.push_back(52);
    fixed_params.push_back(53);
    fixed_params.push_back(54);
    */
    //fixed_params.push_back(48);

    // Sort out fixed vs free params (this can be improved, its too hackish now)
    std::vector<TString> free_params_names;
    for(int i = 0; i < numberParams; i++)
    {
        TString i_str;
        i_str.Form("%i",i);

        bool fixed = false;
        for(int j = 0; j < fixed_params.size(); j++)
        {
            if(i == fixed_params.at(j))
            {
                fixed = true;
            }
        }

        if(fixed)
        {
            continue;
        }
        else
        {
            free_params.push_back(i);
            index_free_params.push_back(i_str);
            free_params_names.push_back(paramNameMap[i].front());
        }
    }


    // Array to pass arguments to fitter
    //double arglist[100];
    //arglist[0] = 0;

    // Set print level
    //  minuit->ExecuteCommand("SET PRINT",arglist,2);
    std::cout << "lets do the fit" << std::endl;
    // Do minimisation
    // arglist[0] = 50000;  // number of function calls
    // arglist[1] = 0.1;  // tolerance

    minuit->SetErrorDef(0.5);
    minuit->SetMaxIterations(50000);
    //give it the function to use
    minuit->SetFCN(logLikelihood);
    minuit->mnsimp();
    minuit->Migrad();
    //minuit->ExecuteCommand("SIMPLEX",arglist,2);

    // Then get results
    for(int i = 0; i < numberParams; i++)
    {
        minuit->GetParameter(i, AdjustActs[i], AdjustActs_Err[i]);
        //AdjustActs_Err[i] = minuit->GetParError(i);
    }


    Double_t CovMatrix[free_params.size()][free_params.size()];
    minuit->mnemat(&CovMatrix[0][0],free_params.size());

  
    /*
    TH2D* hCorrMatrix = new TH2D("hCorrMatrix","Correlation Matrix", free_params.size(), 0, free_params.size(), free_params.size(), 0,free_params.size());

    for(int i = 0; i < free_params.size(); i++)
    {
        for(int j = 0; j < free_params.size(); j++)
        {
          
            //std::cout << index_free_params.at(i) << "  " << index_free_params.at(j)<< std::endl;
            //CovMatrix[i][j] = minuit->GetCovarianceMatrixElement(i,j);

            // std::cout << CovMatrix[i][j] << std::endl;
            Double_t value = CovMatrix[i][j]/(sqrt(CovMatrix[i][i])*sqrt(CovMatrix[j][j]));
            hCorrMatrix->Fill(free_params_names.at(i), free_params_names.at(j), value);

        }
    }
  
    */
    // these were commented out for some reason, I did not remove them
    /*
    double chi2, edm, errdef; 
    int nvpar, nparx;
    minuit->GetStats(chi2, edm, errdef, nvpar, nparx);
    int ndf = npfits-nvpar;
    */
    
    // create a function to draw the final stacked histograms
    /*
    hCorrMatrix->LabelsDeflate();
    hCorrMatrix->GetZaxis()->SetRangeUser(-1,1);
    hCorrMatrix->SetTitle("Correlation Matrix");
    gStyle->SetGridStyle(0);
    hCorrMatrix->GetXaxis()->LabelsOption("v");

    TCanvas *c = new TCanvas("c");
    c->SetFillColor(kWhite);
    c->SetGrid();
    hCorrMatrix->Draw("colz");
    */

}


void logLikelihood(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  ) {

  //std::cout << "entering fit" << std::endl;

  double loglik = 0; 
  //double tmp;
  //   std::cout << "getting 1D histograms" << std::endl;

  TH1F *tmpData1D;
  // std::cout << allDataSamples1D->GetEntries()  << std::endl;
  for ( int i = 0; i < allDataSamples1D->GetEntries() ; i++ ) { // there are i samples for each channel
    TString i_str;
    i_str.Form("%i",i);
    //  std::cout << i << std::endl;

    tmpData1D = (TH1F*)allDataSamples1D->At(i)->Clone("tmpData1D"+i_str+"_");

    // std::cout << tmpData1D->Integral() << std::endl;

    int nBinsX = tmpData1D->GetNbinsX();
    for ( int ix = 1; ix <= nBinsX; ++ix ) {
      Int_t nData = (Int_t)tmpData1D->GetBinContent(ix);
      double nMC = getNumberMC1D(i,ix,p);

      if ( nMC > 0 && TMath::Poisson(nData,nMC) > 0 ) {
	loglik += TMath::Log(TMath::Poisson(nData,nMC));
      } else {
	loglik-=10.;
      }
        
    }//~bins
  }//channels
  //  std::cout << "getting 2D histograms" << std::endl;

  TH2F *tmpData2D;
  // std::cout << allDataSamples->GetEntries()  << std::endl;
  for ( int i = 0; i < allDataSamples2D->GetEntries() ; i++ ) { // there are i samples for each channel
    TString i_str;
    i_str.Form("%i",i);

    tmpData2D = (TH2F*)allDataSamples2D->At(i)->Clone("tmpData2D"+i_str+"_");

    int nBinsX = tmpData2D->GetNbinsX();
    int nBinsY = tmpData2D->GetNbinsY();

    for ( int ix = 1; ix <= nBinsX; ++ix ) {
      for ( int iy = 1; iy <= nBinsY; ++iy ) {
        Int_t nData = (Int_t)tmpData2D->GetBinContent(ix,iy);
	double nMC = getNumberMC2D(i,ix,iy,p);

	if ( nMC > 0 && TMath::Poisson(nData,nMC) > 0 ) {
	  loglik += TMath::Log(TMath::Poisson(nData,nMC));
	} else {
	  loglik-=10.;
	}
      }//~ybins
    }//~xbins
  }//channels

 
  // add constraints to improve likelihood
   //will eventually add gaussian constraint
    for ( int i = 0; i < 14; i++ ) {
    double constraint = 0.;
    //constraint = pow( (p[i] - 1) / (0.01), 2);
    // std::cout << i << "   " << p[i];
    bool fixed = false;
    for ( int j = 0; j < fixed_params.size(); j++ ) {
      if ( i == fixed_params.at(j) ) fixed = true;
    }

    //if ( (i <2) || (i == 7 ) || (i==10) || (i == 14) || (i==15) || (i==16) || (i == 18) || (i == 19) || (i == 21) || 9i==22) ) {
      // constraint = pow( (p[i] - 1.) / (paramActErrMap[i].front() / paramActMap[i].front()), 2);
      constraint = TMath::Gaus(p[i],1.,(paramActErrMap[i].front() / paramActMap[i].front()),true);
      

      // std::cout << "   " << constraint << std::endl;
      //  loglik += log(constraint);  
	// }
  }
  

  
  fval = -2*loglik; 
    
  //tmpData->Delete();

}


Double_t getNumberMC1D(Int_t channel, Int_t binx, Double_t *p) {

  // std::cout <<"getting number of1D MC... "  <<channel << std::endl;

  double nMC = 0;

  // (1) grab a hist from the sample list of this channel
  // (2) figure out which parameter it corresponds to
  TH1F *tmpHist;
  int which_param;
  //  std::cout <<"getting number of MC... "  <<channel << std::endl;

  //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

  for ( int k = 0; k < allMCSamples1D[channel]->GetEntries(); k++ ) {
    tmpHist = (TH1F*)allMCSamples1D[channel]->At(k);

    //if ( tmpHist->Integral() == 0 ) continue;
   
    TString tmpName = tmpHist->GetName();

    //      std::cout << "looking for " << tmpName << std::endl;
    bool foundParam = false;

    for ( int i = 0; (i < numberParams) && !foundParam; i++ ) {

      // std::vector<int>::iterator it = find(fixed_params.begin(), fixed_params.end(), i);
      //if ( it != fixed_params.end() ) continue;      

      for ( int j = 0; (j < paramNameMap[i].size()) && !foundParam; j++ ) {
	//std::cout <<"is it...  " <<paramNameMap[i].at(j) << std::endl;

	if ( tmpName.Contains(paramNameMap[i].at(j)) ) {
	  foundParam = true;
          which_param = i;
	}

      }//~j searching through array of params for the right one
    }//~i list of isotopes with same parameter

    if ( foundParam ) {nMC += p[which_param]*tmpHist->GetBinContent(binx);}
    //  else {std::cout << "error could not find histogram: " << tmpName << std::endl;}	 
  }//each histogram

  //tmpHist->Delete();

  return nMC;
}

Double_t getNumberMC2D(Int_t channel, Int_t binx, Int_t biny, Double_t *p) {

  // std::cout <<"getting number of 2D MC... "  <<channel << std::endl;

  double nMC = 0;

  // (1) grab a hist from the sample list of this channel
  // (2) figure out which parameter it corresponds to
  TH2F *tmpHist2;
  int which_param;
  //  std::cout <<"getting number of MC... "  <<channel << std::endl;

  //std::cout << allMCSamples[channel]->GetEntries() << std::endl;

  for ( int k = 0; k < allMCSamples2D[channel]->GetEntries(); k++ ) {
    tmpHist2 = (TH2F*)allMCSamples2D[channel]->At(k);
    TString tmpName = tmpHist2->GetName();

    //  std::cout << "looking for " << tmpName << std::endl;
    bool foundParam = false;

    for ( int i = 0; (i < numberParams) && !foundParam; i++ ) {

      std::vector<int>::iterator it = find(fixed_params.begin(), fixed_params.end(), i);
      if ( it != fixed_params.end() ) continue;    

      for ( int j = 0; (j < paramNameMap[i].size()) && !foundParam; j++ ) {
	//	   std::cout <<"is it...  " <<paramNameMap[i].at(j) << std::endl;

	if ( tmpName.Contains(paramNameMap[i].at(j)) ) {
	  foundParam = true;
          which_param = i;
	}

      }//~j searching through array of params for the right one
    }//~i list of isotopes with same parameter

    if ( foundParam ) {nMC += p[which_param]*tmpHist2->GetBinContent(binx,biny);}
    // else {std::cout << "error could not find histogram: " << tmpName << std::endl;}	 
  }//each histogram

  return nMC;
}
