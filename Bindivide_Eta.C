// This code should divide the NTuple tree into bins with equal number of entries!!!
// Need to improve this code look Bindivide.C and Bindivide_pT.C for looping concept

#include <cmath>
#include "TFile.h"
#include "TH1.h"
void Bindivide_Eta(const char *infile = "Ntuple_V1_P123.root")
{
    TFile *f = new TFile(infile);

    TTree *ntuple1 = (TTree *)f->Get("ntuple1");
    TTree *ntuple2 = (TTree *)f->Get("ntuple2");
    TTree *ntuple3 = (TTree *)f->Get("ntuple3");
    TTree *ntuple4 = (TTree *)f->Get("ntuple4");
    TTree *ntuple5 = (TTree *)f->Get("ntuple5");
    TTree *ntuple6 = (TTree *)f->Get("ntuple6");

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

    TH1D *h = new TH1D("eta_pair", "eta_pair", 1000, -1.2, 1.2);

    for (int i = 0; i < nentries; i++)
    // for (int i=0;i<1000;i++)
    {
        ntuple1->GetEntry(i);
        if (Minv > 4 && cone > 0.7 && fitPts_min_pair < 15 && pT_pair > 25)
            continue;
        h->Fill(eta_pair);

        // if (Minv<4 && cone<0.7 && fitPts_min_pair>15 && pT_pair>2.60 && pT_pair<25){h2_Minv->Fill(Minv);}

        // cout << pT_pair << "pT_pair Value"<< endl;
    }
    cout << h->GetEntries() << "Entries" << endl;

    h->Draw();
    gPad->SetLogy();
    gPad->Update();

    double a17th = h->GetEntries() / 7;
    double a27th = (2.0 / 7.0) * h->GetEntries();
    double a37th = (3.0 / 7.0) * h->GetEntries();
    double a47th = (4.0 / 7.0) * h->GetEntries();
    double a57th = (5.0 / 7.0) * h->GetEntries();
    double a67th = (6.0 / 7.0) * h->GetEntries();

    cout << h->GetEntries() << "\t Total Entries" << endl;
    cout << a17th << "\t a17th" << endl;
    cout << a27th << "\t a27th" << endl;
    cout << a37th << "\t a37th" << endl;
    cout << a47th << "\t a47th" << endl;
    cout << a57th << "\t a57th" << endl;
    cout << a67th << "\t a67th" << endl;

    vector<double> diff_h_1;
    vector<double> diff_h_2;
    vector<double> diff_h_3;
    vector<double> diff_h_4;
    vector<double> diff_h_5;
    vector<double> diff_h_6;
    vector<double> diff_h_1_Index;
    vector<double> diff_h_2_Index;
    vector<double> diff_h_3_Index;
    vector<double> diff_h_4_Index;
    vector<double> diff_h_5_Index;
    vector<double> diff_h_6_Index;

    TAxis *axis_h = h->GetXaxis();

    double kk = 0;
    for (int i = 0; i < h->GetNbinsX(); i++)
    {
        double y = (double)h->GetBinContent(i);
        if (y <= 0)
            continue;

        kk = kk + y;
        cout << kk << "\t kk Value" << endl;
        if (((abs(((kk - a17th) / a17th)) > 0.0001) && (abs(((kk - a17th) / a17th)) < 0.01)))
        {
            cout << std::setprecision(7) << kk << " \t 1/7 th data"
                 << "\t on bin\t " << i << "\t from " << a17th << "\t with eta value" << axis_h->GetBinUpEdge(i) << endl;
            diff_h_1.push_back(abs(((kk - a17th) / a17th)));
            diff_h_1_Index.push_back(i);
        }
        if (((abs(((kk - a27th) / a27th)) > 0.0001) && (abs(((kk - a27th) / a27th)) < 0.01)))
        {
            cout << std::setprecision(7) << kk << " \t 2/7 th data"
                 << "\t on bin\t " << i << "\t from " << a27th << "\t with eta value" << axis_h->GetBinUpEdge(i) << endl;
            diff_h_2.push_back(abs(((kk - a27th) / a27th)));
            diff_h_2_Index.push_back(i);
        }

        if (((abs(((kk - a37th) / a37th)) > 0.0001) && (abs(((kk - a37th) / a37th)) < 0.01)))
        {
            cout << std::setprecision(7) << kk << " \t 3/7 th data"
                 << "\t on bin\t " << i << "\t from " << a37th << "\t with eta value" << axis_h->GetBinUpEdge(i) << endl;
            diff_h_3.push_back(abs(((kk - a37th) / a37th)));
            diff_h_3_Index.push_back(i);
        }

        if (((abs(((kk - a47th) / a47th)) > 0.0001) && (abs(((kk - a47th) / a47th)) < 0.01)))
        {
            cout << std::setprecision(7) << kk << " \t 4/7 th data"
                 << "\t on bin\t " << i << "\t from " << a47th << "\t with eta value" << axis_h->GetBinUpEdge(i) << endl;
            diff_h_4.push_back(abs(((kk - a47th) / a47th)));
            diff_h_4_Index.push_back(i);
        }

        if (((abs(((kk - a57th) / a57th)) > 0.0001) && (abs(((kk - a57th) / a57th)) < 0.01)))
        {
            cout << std::setprecision(7) << kk << " \t 5/7 th data"
                 << "\t on bin\t " << i << "\t from " << a57th << "\t with eta value" << axis_h->GetBinUpEdge(i) << endl;
            diff_h_5.push_back(abs(((kk - a57th) / a57th)));
            diff_h_5_Index.push_back(i);
        }
        if (((abs(((kk - a67th) / a67th)) > 0.0001) && (abs(((kk - a67th) / a67th)) < 0.01)))
        {
            cout << std::setprecision(7) << kk << " \t 6/7 th data"
                 << "\t on bin\t " << i << "\t from " << a67th << "\t with eta value" << axis_h->GetBinUpEdge(i) << endl;
            diff_h_6.push_back(abs(((kk - a67th) / a67th)));
            diff_h_6_Index.push_back(i);
        }
    }

    cout << *min_element(diff_h_2.begin(), diff_h_2.end()) << "\t minimun value in diff_h_2" << endl;
    cout << min_element(diff_h_2.begin(), diff_h_2.end()) - diff_h_2.begin() << "\t minimun value index in diff_h_2" << endl;
    cout << diff_h_2_Index.at(min_element(diff_h_2.begin(), diff_h_2.end()) - diff_h_2.begin()) << "\t bin number where diff_h_2 is minimum" << endl;
    cout << axis_h->GetBinUpEdge(diff_h_2_Index.at(min_element(diff_h_2.begin(), diff_h_2.end()) - diff_h_2.begin())) << " 2nd pT limit" << endl;

    double eta_pair_0 = -1.20;
    double eta_pair_1 = axis_h->GetBinUpEdge(diff_h_1_Index.at(min_element(diff_h_1.begin(), diff_h_1.end()) - diff_h_1.begin()));
    double eta_pair_2 = axis_h->GetBinUpEdge(diff_h_2_Index.at(min_element(diff_h_2.begin(), diff_h_2.end()) - diff_h_2.begin()));
    double eta_pair_3 = axis_h->GetBinUpEdge(diff_h_3_Index.at(min_element(diff_h_3.begin(), diff_h_3.end()) - diff_h_3.begin()));
    double eta_pair_4 = axis_h->GetBinUpEdge(diff_h_4_Index.at(min_element(diff_h_4.begin(), diff_h_4.end()) - diff_h_4.begin()));
    double eta_pair_5 = axis_h->GetBinUpEdge(diff_h_5_Index.at(min_element(diff_h_5.begin(), diff_h_5.end()) - diff_h_5.begin()));
    double eta_pair_6 = axis_h->GetBinUpEdge(diff_h_6_Index.at(min_element(diff_h_6.begin(), diff_h_6.end()) - diff_h_6.begin()));
    double eta_pair_7 = 1.20;

    double eta_range[] = {eta_pair_0, eta_pair_1, eta_pair_2, eta_pair_3, eta_pair_4, eta_pair_5, eta_pair_6, eta_pair_7};

    cout << "final Result" << endl;
    cout << "final Result" << endl;
    cout << "final Result" << endl;

    for (int i = 0; i < sizeof(eta_range) / sizeof(eta_range[0]); i++)
    {
        cout << eta_range[i];
    }

} // Main Finction Ends
