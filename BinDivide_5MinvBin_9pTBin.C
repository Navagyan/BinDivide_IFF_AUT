// This code should divide the total data set into 5 Minv Bins and Each Minv Bin is further divided into 9 pT Bins

#include <cmath>
void BinDivide_5MinvBin_9pTBin(const char *infile = "/Users/nghimire/Research/Run17_AUT/IFF_Analysis/Ntuple_V1/Asym_Code_P20ic.SL22b/Asym_Code_Kaons/Ntuple_RawTree_P123_TPC_Kaons.root")
{
    TFile *f = new TFile(infile);

    TTree *ntuple1 = (TTree *)f->Get("ntuple1");
    TTree *ntuple2 = (TTree *)f->Get("ntuple2");
    TTree *ntuple3 = (TTree *)f->Get("ntuple3");
    TTree *ntuple4 = (TTree *)f->Get("ntuple4");
    TTree *ntuple5 = (TTree *)f->Get("ntuple5");
    // TTree *ntuple6 = (TTree*)f->Get("ntuple6");

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
    //************************************************* Time Saving Block Starts**********************************************************************//
    // This is the step which takes more time so run the code with this block first only which makes the root file contating histogram which can be used later on
    // TFile *hMinvRoot = new TFile("hMinv.root", "RECREATE");
    // TH1D *h_Minv = new TH1D("h_Minv", "h_Minv", 100000, 0.5, 4);
    // for (int i = 0; i < nentries; i++)
    //// for (int i = 0; i < 100000; i++)
    //{
    //    ntuple1->GetEntry(i);
    //    if (Minv > 4 && cone > 0.7 && fitPts_min_pair < 15 && pT_pair > 25)
    //        continue;
    //    h_Minv->Fill(Minv);
    //}
    // hMinvRoot->Write();
    // hMinvRoot->Close();
    //************************************************* Time Saving Block Ends**********************************************************************//

    //  Once you created root file from above time saving block then read the root file and get the histogram
    TFile *hMinv_infile = new TFile("hMinv.root");
    TH1D *h_Minv = (TH1D *)hMinv_infile->Get("h_Minv")->Clone();

    cout << h_Minv->GetEntries() << "\tTotal Entries\t" << endl;
    TCanvas *c = new TCanvas("Minv_histo", "Minv_histo", 700, 900);
    c->cd()->SetLogy();
    h_Minv->Draw();
    c->SaveAs("./h_Minv.pdf");
    TAxis *axis_h = h_Minv->GetXaxis();
    double Minv15th = h_Minv->GetEntries() / 5;
    double Minv25th = (0.4) * h_Minv->GetEntries();
    double Minv35th = (0.6) * h_Minv->GetEntries();
    double Minv45th = (0.8) * h_Minv->GetEntries();
    double Minv55th = (1) * h_Minv->GetEntries();

    cout << Minv15th << "\t" << Minv25th << "\t" << Minv35th << "\t" << Minv45th << endl;

    double xx = 0;
    vector<double> diff_h_Minv_1;
    vector<double> diff_h_Minv_2;
    vector<double> diff_h_Minv_3;
    vector<double> diff_h_Minv_4;

    vector<double> diff_h_Minv_1_index;
    vector<double> diff_h_Minv_2_index;
    vector<double> diff_h_Minv_3_index;
    vector<double> diff_h_Minv_4_index;

    for (int i = 0; i < h_Minv->GetNbinsX(); i++)
    {
        double y = (double)h_Minv->GetBinContent(i);
        if (y <= 0)
            continue;
        double dx = axis_h->GetBinWidth(i);
        xx = xx + y;
        // cout << xx << "\tSum of h1_Minv\t" << Minv15th << "\t Minv15th\t" << (xx - Minv15th) / Minv15th << endl;
        // cout << xx << "\tSum of h1_Minv\t" << Minv25th << "\t Minv15th\t" << (xx - Minv25th) / Minv25th << endl;
        // cout << xx << "\tSum of h1_Minv\t" << Minv35th << "\t Minv15th\t" << (xx - Minv35th) / Minv35th << endl;
        // cout << xx << "\tSum of h1_Minv\t" << Minv45th << "\t Minv15th\t" << (xx - Minv45th) / Minv45th << endl;
        if (((abs(((xx - Minv15th) / Minv15th)) > 0.0001) && (abs(((xx - Minv15th) / Minv15th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " 1/5th data"
                 << "on bin"
                 << "from " << Minv15th << i << "with Minv value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - Minv15th) / Minv15th) << "\t diff  \t " << endl;
            diff_h_Minv_1.push_back(abs((xx - Minv15th) / Minv15th));
            diff_h_Minv_1_index.push_back(i);
        }

        if (((abs(((xx - Minv25th) / Minv25th)) > 0.0001) && (abs(((xx - Minv25th) / Minv25th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " 2/5th data"
                 << "on bin"
                 << "from " << Minv25th << i << "with Minv value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - Minv25th) / Minv25th) << "\t diff  \t " << endl;
            diff_h_Minv_2.push_back(abs((xx - Minv25th) / Minv25th));
            diff_h_Minv_2_index.push_back(i);
        }
        if (((abs(((xx - Minv35th) / Minv35th)) > 0.0001) && (abs(((xx - Minv35th) / Minv35th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " \t 3/5 th data"
                 << "\t on bin\t " << i << "\t from " << Minv35th << "\t with pT value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - Minv35th) / Minv35th) << "\t diff  \t " << endl;
            diff_h_Minv_3.push_back(abs((xx - Minv35th) / Minv35th));
            diff_h_Minv_3_index.push_back(i);
        }
        if (((abs(((xx - Minv45th) / Minv45th)) > 0.0001) && (abs(((xx - Minv45th) / Minv45th)) < 0.01)))
        {
            cout << std::setprecision(7) << xx << " \t 4/5 th data"
                 << "\t on bin\t" << i << "\t from \t" << Minv45th << "\t with pT value" << axis_h->GetBinUpEdge(i) << endl;
            cout << abs((xx - Minv45th) / Minv45th) << "\t diff  \t " << endl;
            diff_h_Minv_4.push_back(abs((xx - Minv45th) / Minv45th));
            diff_h_Minv_4_index.push_back(i);
        }
    }

    double Minv_0 = 0.5;
    double Minv_1 = axis_h->GetBinUpEdge(diff_h_Minv_1_index.at(min_element(diff_h_Minv_1.begin(), diff_h_Minv_1.end()) - diff_h_Minv_1.begin()));
    double Minv_2 = axis_h->GetBinUpEdge(diff_h_Minv_2_index.at(min_element(diff_h_Minv_2.begin(), diff_h_Minv_2.end()) - diff_h_Minv_2.begin()));
    double Minv_3 = axis_h->GetBinUpEdge(diff_h_Minv_3_index.at(min_element(diff_h_Minv_3.begin(), diff_h_Minv_3.end()) - diff_h_Minv_3.begin()));
    double Minv_4 = axis_h->GetBinUpEdge(diff_h_Minv_4_index.at(min_element(diff_h_Minv_4.begin(), diff_h_Minv_4.end()) - diff_h_Minv_4.begin()));
    double Minv_5 = 4.0;

    double MINV[6] = {Minv_0, Minv_1, Minv_2, Minv_3, Minv_4, Minv_5};

    cout << "==========Final Result========" << endl;
    cout << "==========Final Result========" << endl;
    cout << "==========Final Result========" << endl;

    cout << "Dividing Entire data set in 9 Minv bin 1st" << endl;

    const int MINBIN_size=(sizeof(MINV) / sizeof(MINV[0]));

    for (int i = 0; i < MINBIN_size; i++)
    {
        cout << MINV[i] << "\t Minv value for boundary\t" << i << endl;
    }

    vector<double> diff_pT[6][10];
    vector<double> diff_pT_index[6][10];
    vector<double> pT_bin[10];
    double pT_j9th[10];
    TH1D *h1_pT_bin[10];

    //************************ Second Time Saving Block  Starts **************************//
   //  TFile *fROOT = new TFile("hist_5MinvBin_pT.root", "RECREATE");
   //  TCanvas *canv_pT = new TCanvas("h1_pT", "h1_pT", 900, 700);
   //  for (int i = 0; i < (MINBIN_size-1); i++)
   // {
   //     h1_pT_bin[i] = new TH1D(Form("h1_pT_bin%i", i), "", 8000, 0.5, 25);
   // }
   //  for (int j = 0; j < (MINBIN_size-1); j++)
   // {
   //     for (int i = 0; i < nentries; i++)
   //     {
   //         ntuple1->GetEntry(i);
   //         if (pT_pair < 0.5 && Minv > 4 && cone <= 0.7 && fitPts_min_pair < 15 && Minv < 0.2 && Minv > 4)
   //             continue;
   //         if (Minv < MINV[j] || Minv >= MINV[j + 1])
   //             continue;
   //         h1_pT_bin[j]->Fill(pT_pair);
   //     }
   //     canv_pT->cd()->SetLogy();
   //     h1_pT_bin[j]->Draw();
   //     canv_pT->SaveAs(Form("h1_pT_bin%i.png", j));
   // }
   //  fROOT->Write();
    //************************ Second Time Saving Block  Ends **************************//

    //  Once you created root file from above 2nd time saving block then read the root file and get the histograms

    TFile *hist_pT_infile = new TFile("./hist_5MinvBin_pT.root");
    for (int i = 0; i < (MINBIN_size-1); i++)
    {
        h1_pT_bin[i] = (TH1D *)hist_pT_infile->Get(Form("h1_pT_bin%i", i))->Clone();
    }

    for (int i = 0; i < (MINBIN_size-1); i++)
    {
        // divide each Minv bin into 8 pT bin(the 1st and last bin boundaries are fixed, so we need only 8 bin values
        for (int j = 0; j < 8; j++)
        {
            pT_j9th[j] = ((j + 1) / 9.0) * h1_pT_bin[i]->GetEntries();
            double xx_pT = 0;
            for (int k = 0; k < h1_pT_bin[i]->GetNbinsX(); k++)
            {
                double y_pT = h1_pT_bin[i]->GetBinContent(k);
                if (y_pT <= 0)
                    continue;
                xx_pT = xx_pT + y_pT;
                for (int m = 0; m < 8; m++)
                {
                    if ((abs(xx_pT - pT_j9th[m]) / pT_j9th[m]) > 0.0001 && (abs(xx_pT - pT_j9th[m]) / pT_j9th[m]) < 0.01)
                    {
                        diff_pT[i][m].push_back(abs(xx_pT - pT_j9th[m]) / pT_j9th[m]);
                        diff_pT_index[i][m].push_back(k);
                    }
                } // m-loop ended
            }     // loop over bins of individual pTbin-loop ended;
        }         // 8pT-loop ended
    }             // Minv-loop ended

    for (int mbin = 0; mbin < (MINBIN_size-1); mbin++)
    {
        for (int pTbin = 0; pTbin < 8; pTbin++)
        {
            pT_bin[mbin].push_back(h1_pT_bin[mbin]->GetXaxis()->GetBinUpEdge(diff_pT_index[mbin][pTbin].at(min_element(diff_pT[mbin][pTbin].begin(), diff_pT[mbin][pTbin].end()) - diff_pT[mbin][pTbin].begin())));

        } // pTbin ended
    }     // mbin ended

    for (int mbin = 0; mbin < (MINBIN_size-1); mbin++)
    {
        pT_bin[mbin].insert(pT_bin[mbin].begin(), 0.5);
        pT_bin[mbin].push_back(25.0);
        cout << "pT Bin Boundary for \t " << mbin << "\t Minv bin" << endl;
        for (int pTbin = 0; pTbin < pT_bin[mbin].size(); pTbin++)
        {
            if (pTbin == 0)
            {
                cout << "{" << pT_bin[mbin].at(pTbin) << ",";
            }
            if (pTbin > 0 && pTbin < (pT_bin[mbin].size() - 1))
            {
                cout << pT_bin[mbin].at(pTbin) << ",";
            }
            if (pTbin == (pT_bin[mbin].size() - 1))
            {
                cout << pT_bin[mbin].at(pTbin) << "}" << endl;
            }
        } // pT bin ended;
    }     // mbin ended;

} // Main function eneded;
