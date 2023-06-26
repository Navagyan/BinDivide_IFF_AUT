// This code should divide the NTuple tree into bins with equal number of entries!!!
#include <cmath>
#include "TFile.h"
#include "TH1.h"
// void Bindivide(const char *infile="./Ntuple_V4_P123.root"){
// void Bindivide(const char *infile = "./Ntuple_RawTree_P123_iff_dataV4_noTOF.root")
void Bindivide_5pTBin_9MinvBin(const char *infile = "/Users/nghimire/Research/Run17_AUT/IFF_Analysis/Ntuple_V1/Asym_Code_P20ic.SL22b/Asym_Code_Kaons/Ntuple_RawTree_P123_TPC_Kaons.root")
{
    TFile *f = new TFile(infile);

    TTree *ntuple1 = (TTree *)f->Get("ntuple1");
    TTree *ntuple2 = (TTree *)f->Get("ntuple2");
    TTree *ntuple3 = (TTree *)f->Get("ntuple3");
    TTree *ntuple4 = (TTree *)f->Get("ntuple4");
    TTree *ntuple5 = (TTree *)f->Get("ntuple5");

    float cone, Minv, pT_pair, eta_pair, fitPts_min_pair;

    ntuple1->SetBranchAddress("cone", &cone);
    ntuple2->SetBranchAddress("Minv", &Minv);
    ntuple2->SetBranchAddress("pT_pair", &pT_pair);
    ntuple2->SetBranchAddress("eta_pair", &eta_pair);
    ntuple5->SetBranchAddress("fitPts_min_pair", &fitPts_min_pair);

    Int_t nentries = (Int_t)ntuple1->GetEntries();

    cout << nentries << "nentries before cut" << endl;

    ntuple1->AddFriend("ntuple2");
    ntuple1->AddFriend("ntuple4");
    ntuple1->AddFriend("ntuple5");

    //************************************************* Time Saving Block Starts **********************************************************************//
    // This is the step which takes more time so run the code with this block first only which makes the root file contating histogram which can be used later on
    // TFile *h5pTBinRoot = new TFile("h5pTBin.root", "RECREATE");
    // TH1D *h_5pTBin = new TH1D("h_5pTBin", "h_5pTBin", 1000000, 0, 25);
    // for (int i = 0; i < nentries; i++)
    //// for (int i = 0; i < 100000; i++)
    //{
    //    ntuple1->GetEntry(i);
    //    if (Minv > 4 && cone > 0.7 && fitPts_min_pair < 15 && pT_pair > 25)
    //        continue;
    //    h_5pTBin->Fill(pT_pair);
    //}
    // h5pTBinRoot->Write();
    // h5pTBinRoot->Close();
    //************************************************* Time Saving Block Ends**********************************************************************//

    //  Once you created root file from above time saving block then read the root file and get the histogram
    TFile *h5pTBin_infile = new TFile("h5pTBin.root");
    TH1D *h_5pTBin = (TH1D *)h5pTBin_infile->Get("h_5pTBin")->Clone();

    TCanvas *c = new TCanvas("5pTBin_histo", "5pTBin_histo", 700, 900);
    c->cd()->SetLogy();
    h_5pTBin->Draw();
    c->SaveAs("./5pTBin_histo.pdf");

    TAxis *axis_h = h_5pTBin->GetXaxis();
    double pT15th = (0.2) * h_5pTBin->GetEntries();
    double pT25th = (0.4) * h_5pTBin->GetEntries();
    double pT35th = (0.6) * h_5pTBin->GetEntries();
    double pT45th = (0.8) * h_5pTBin->GetEntries();
    double pT55th = (1) * h_5pTBin->GetEntries();
    cout << pT15th << "\t" << pT25th << "\t" << pT35th << "\t" << pT45th << "\t" << pT55th << endl;

    double xx = 0;
    vector<double> diff_h_pT_1;
    vector<double> diff_h_pT_2;
    vector<double> diff_h_pT_3;
    vector<double> diff_h_pT_4;

    vector<double> diff_h_pT_1_index;
    vector<double> diff_h_pT_2_index;
    vector<double> diff_h_pT_3_index;
    vector<double> diff_h_pT_4_index;

    for (int i = 0; i < h_5pTBin->GetNbinsX(); i++)
    {
        // cout << h1->GetNbinsX() << "Number if h1 bins"<<endl;
        double y = (double)h_5pTBin->GetBinContent(i);
        if (y <= 0)
            continue;
        double dx = axis_h->GetBinWidth(i);
        // cout << xx << "+"<< dx << "*"<< y << "="<<xx+dx*y<< " \t in bin \t "<< i<<endl;
        xx = xx + y;
        // cout << abs((xx-yy)/yy) << "\t diff  \t "<< endl;

        if (((abs(((xx - pT15th) / pT15th)) > 0.0001) && (abs(((xx - pT15th) / pT15th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " 1/5th data"
                 << "on bin"
                 << "from " << pT15th << i << "with pT value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - pT15th) / pT15th) << "\t diff  \t " << endl;
            diff_h_pT_1.push_back(abs((xx - pT15th) / pT15th));
            diff_h_pT_1_index.push_back(i);
        }
        if (((abs(((xx - pT25th) / pT25th)) > 0.0001) && (abs(((xx - pT25th) / pT25th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " 2/5th data"
                 << "on bin"
                 << "from " << pT25th << i << "with pT value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - pT25th) / pT25th) << "\t diff  \t " << endl;
            diff_h_pT_2.push_back(abs((xx - pT25th) / pT25th));
            diff_h_pT_2_index.push_back(i);
        }

        if (((abs(((xx - pT35th) / pT35th)) > 0.0001) && (abs(((xx - pT35th) / pT35th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " 3/5th data"
                 << "on bin"
                 << "from " << pT35th << i << "with pT value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - pT35th) / pT35th) << "\t diff  \t " << endl;
            diff_h_pT_3.push_back(abs((xx - pT35th) / pT35th));
            diff_h_pT_3_index.push_back(i);
        }
        if (((abs(((xx - pT45th) / pT45th)) > 0.0001) && (abs(((xx - pT45th) / pT45th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " 4/5th data"
                 << "on bin"
                 << "from " << pT45th << i << "with pT value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - pT45th) / pT45th) << "\t diff  \t " << endl;
            diff_h_pT_4.push_back(abs((xx - pT45th) / pT45th));
            diff_h_pT_4_index.push_back(i);
        }
    }

    double pTBin_0 = 2.5;
    double pTBin_1 = axis_h->GetBinUpEdge(diff_h_pT_1_index.at(min_element(diff_h_pT_1.begin(), diff_h_pT_1.end()) - diff_h_pT_1.begin()));
    double pTBin_2 = axis_h->GetBinUpEdge(diff_h_pT_2_index.at(min_element(diff_h_pT_2.begin(), diff_h_pT_2.end()) - diff_h_pT_2.begin()));
    double pTBin_3 = axis_h->GetBinUpEdge(diff_h_pT_3_index.at(min_element(diff_h_pT_3.begin(), diff_h_pT_3.end()) - diff_h_pT_3.begin()));
    double pTBin_4 = axis_h->GetBinUpEdge(diff_h_pT_4_index.at(min_element(diff_h_pT_4.begin(), diff_h_pT_4.end()) - diff_h_pT_4.begin()));
    double pTBin_5 = 25;

    double PTBIN[6] = {pTBin_0, pTBin_1, pTBin_2, pTBin_3, pTBin_4, pTBin_5};

    cout << "==========Final Result========" << endl;
    cout << "==========Final Result========" << endl;
    cout << "==========Final Result========" << endl;
    const int PTBIN_size = (sizeof(PTBIN) / sizeof(PTBIN[0]));

    cout << PTBIN_size << "pT Bin size\t" << endl;
    for (int i = 0; i < (sizeof(PTBIN) / sizeof(PTBIN[0])); i++)
    {
        cout << PTBIN[i] << "\t pT value for boundary\t" << i << endl;
    }

    vector<double> diff_Minv[6][10];
    vector<double> diff_Minv_index[6][10];
    vector<double> Minv_bin[10];
    double Minv_j9th[10];
    TH1D *h1_Minv_bin[10];

    //************************ Second Time Saving Block  Starts **************************//
    // TFile *fROOT = new TFile("hist_5pTBin_Minv.root", "RECREATE");
    // TCanvas *canv_Minv = new TCanvas("h1_Minv", "h1_Minv", 900, 700);
    // for (int i = 0; i < (PTBIN_size - 1); i++)
    //{
    //    h1_Minv_bin[i] = new TH1D(Form("h1_Minv_bin%i", i), "", 8000, 0.5, 4);
    //}
    // for (int j = 0; j < (PTBIN_size - 1); j++)
    //{
    //    for (int i = 0; i < nentries; i++)
    //    {
    //        ntuple1->GetEntry(i);
    //        if (pT_pair < 0.5 && Minv > 4 && cone <= 0.7 && fitPts_min_pair < 15 && Minv < 0.5 && Minv > 4)
    //            continue;
    //        if (pT_pair < PTBIN[j] || pT_pair >= PTBIN[j + 1])
    //            continue;
    //        h1_Minv_bin[j]->Fill(Minv);
    //    }
    //    canv_Minv->cd()->SetLogy();
    //    h1_Minv_bin[j]->Draw();
    //    canv_Minv->SaveAs(Form("h1_Minv_bin%i.png", j));
    //}
    // fROOT->Write();
    //************************ Second Time Saving Block  Ends **************************//

    TFile *hist_5pTBin_Minv_infile = new TFile("./hist_5pTBin_Minv.root");
    for (int i = 0; i < (PTBIN_size - 1); i++)
    {
        h1_Minv_bin[i] = (TH1D *)hist_5pTBin_Minv_infile->Get(Form("h1_Minv_bin%i", i))->Clone();
    }

    for (int i = 0; i < (PTBIN_size - 1); i++)
    {
        // divide each Minv bin into 8 pT bin(the 1st and last bin boundaries are fixed, so we need only 8 bin values
        for (int j = 0; j < 8; j++)
        {
            Minv_j9th[j] = ((j + 1) / 9.0) * h1_Minv_bin[i]->GetEntries();
            double xx_Minv = 0;
            for (int k = 0; k < h1_Minv_bin[i]->GetNbinsX(); k++)
            {
                double y_Minv = h1_Minv_bin[i]->GetBinContent(k);
                if (y_Minv <= 0)
                    continue;
                xx_Minv = xx_Minv + y_Minv;
                for (int m = 0; m < 8; m++)
                {
                    if ((abs(xx_Minv - Minv_j9th[m]) / Minv_j9th[m]) > 0.0001 && (abs(xx_Minv - Minv_j9th[m]) / Minv_j9th[m]) < 0.01)
                    {
                        diff_Minv[i][m].push_back(abs(xx_Minv - Minv_j9th[m]) / Minv_j9th[m]);
                        diff_Minv_index[i][m].push_back(k);
                    }
                } // m-loop ended
            }     // loop over bins of individual Minvbin-loop ended;
        }         // 8Minv-loop ended
    }             // 5pTBin-loop ended

    for (int pTbin = 0; pTbin < (PTBIN_size - 1); pTbin++)
    {
        for (int mbin = 0; mbin < 8; mbin++)
        {
            Minv_bin[pTbin].push_back(h1_Minv_bin[pTbin]->GetXaxis()->GetBinUpEdge(diff_Minv_index[pTbin][mbin].at(min_element(diff_Minv[pTbin][mbin].begin(), diff_Minv[pTbin][mbin].end()) - diff_Minv[pTbin][mbin].begin())));

        } // mbin ended
    }     // pTbin ended

    for (int pTbin = 0; pTbin < (PTBIN_size - 1); pTbin++)
    {
        Minv_bin[pTbin].insert(Minv_bin[pTbin].begin(), 0.5);
        Minv_bin[pTbin].push_back(4.0);
        cout << "Minv Bin Boundary for \t " << pTbin << "\t pT bin" << endl;
        for (int mbin = 0; mbin < Minv_bin[pTbin].size(); mbin++)
        {
            if (mbin == 0)
            {
                cout << "{" << Minv_bin[pTbin].at(mbin) << ",";
            }
            if (mbin > 0 && mbin < (Minv_bin[pTbin].size() - 1))
            {
                cout << Minv_bin[pTbin].at(mbin) << ",";
            }
            if (mbin == (Minv_bin[pTbin].size() - 1))
            {
                cout << Minv_bin[pTbin].at(mbin) << "}" << endl;
            }
        } // pT bin ended;
    }     // mbin ended;

} // Main Finction Ends
