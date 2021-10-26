void E_TOF_CURVE(){
//========= Macro generated from object: E_TOF_CURVE_MG1/Graph
//========= by ROOT version6.24/00
   
   Double_t E_TOF_CURVE_MG1_fx2[15] = {
   0.807351,
   0.993092,
   1.25313,
   1.59984,
   1.95894,
   2.40472,
   3.04862,
   3.68014,
   4.47263,
   5.36418,
   6.87487,
   7.79119,
   9.10375,
   18.9356,
   20.4958};
   Double_t E_TOF_CURVE_MG1_fy2[15] = {
   300.19,
   313.193,
   323.397,
   330.982,
   338.115,
   344.256,
   351.028,
   355.362,
   357.53,
   358.884,
   360.419,
   361.142,
   361.593,
   362.586,
   362.857};
   auto* graph_MG1 = new TGraph(15,E_TOF_CURVE_MG1_fx2,E_TOF_CURVE_MG1_fy2);
   graph_MG1->SetName("E_TOF_CURVE_MG1");
   auto* f = new TFile("E_TOF_CURVE_MG1.root", "recreate");
   graph_MG1->Write();

   Double_t E_TOF_CURVE_MG3_fx1[15] = {
   0.757821,
   1.09215,
   1.42649,
   1.79797,
   2.29327,
   2.90003,
   3.4201,
   4.17544,
   5.32703,
   6.78819,
   6.78819,
   8.12552,
   9.77242,
   16.6572,
   19.8519};
   Double_t E_TOF_CURVE_MG3_fy1[15] = {
   301.544,
   317.527,
   329.537,
   336.58,
   344.256,
   350.667,
   354.189,
   357.168,
   359.245,
   360.329,
   360.329,
   360.69,
   361.232,
   362.586,
   363.038};
   TGraph *graph_MG3 = new TGraph(15,E_TOF_CURVE_MG3_fx1,E_TOF_CURVE_MG3_fy1);
   graph_MG3->SetName("E_TOF_CURVE_MG3");
   f = new TFile("E_TOF_CURVE_MG3.root", "recreate");
   graph_MG3->Write();

   Double_t E_TOF_CURVE_MG4_fx3[18] = {
   0.931178,
   1.11692,
   1.11692,
   1.27789,
   1.51316,
   1.51316,
   1.81035,
   2.15706,
   2.66475,
   3.14768,
   3.14768,
   3.74205,
   4.54692,
   5.33942,
   6.30527,
   6.30527,
   9.8591,
   20.4958};
   Double_t E_TOF_CURVE_MG4_fy3[18] = {
   300.142,
   325.555,
   325.555,
   332.65,
   340.39,
   340.39,
   348.646,
   355.095,
   362.448,
   367.092,
   367.092,
   370.833,
   372.897,
   374.445,
   375.219,
   375.219,
   376.122,
   378.057};
   auto* graph_MG4 = new TGraph(18,E_TOF_CURVE_MG4_fx3,E_TOF_CURVE_MG4_fy3);
   graph_MG4->SetName("E_TOF_CURVE_MG4");
   f = new TFile("E_TOF_CURVE_MG4.root", "recreate");
   graph_MG4->Write();


   Double_t E_TOF_CURVE_MG5_fx4[14] = {
   0.398722,
   0.782586,
   1.17883,
   1.58746,
   1.58746,
   2.08277,
   2.60284,
   3.22198,
   3.72967,
   4.49739,
   5.62422,
   6.9244,
   8.36079,
   19.4804};
   Double_t E_TOF_CURVE_MG5_fy4[14] = {
   281.042,
   306.017,
   320.671,
   329.752,
   329.752,
   337.699,
   343.065,
   348.328,
   350.702,
   352.456,
   354.107,
   355.243,
   355.862,
   356.997};
   auto* graph_MG5 = new TGraph(14,E_TOF_CURVE_MG5_fx4,E_TOF_CURVE_MG5_fy4);
   graph_MG5->SetName("E_TOF_CURVE_MG5");
   f = new TFile("E_TOF_CURVE_MG5.root", "recreate");
   graph_MG5->Write();


   Double_t E_TOF_CURVE_MG7_fx5[13] = {
   1.30266,
   1.93418,
   2.44187,
   2.78858,
   3.17245,
   3.61822,
   4.23736,
   5.15368,
   6.55292,
   7.82834,
   8.497,
   18.1184,
   20.2234};
   Double_t E_TOF_CURVE_MG7_fy5[13] = {
   280.217,
   304.985,
   320.671,
   329.236,
   336.254,
   341.207,
   344.303,
   347.193,
   349.257,
   349.979,
   350.186,
   352.147,
   352.559};
   auto* graph_MG7 = new TGraph(13,E_TOF_CURVE_MG7_fx5,E_TOF_CURVE_MG7_fy5);
   graph_MG7->SetName("E_TOF_CURVE_MG7");
   f = new TFile("E_TOF_CURVE_MG7.root", "recreate");
   graph_MG7->Write();


   Double_t E_TOF_CURVE_MG11_fx6[17] = {
   0.633994,
   0.695907,
   0.844499,
   1.00547,
   1.29028,
   1.53793,
   1.53793,
   2.21898,
   3.22198,
   3.22198,
   4.37357,
   5.5623,
   5.5623,
   6.88725,
   8.38556,
   12.7443,
   19.3071};
   Double_t E_TOF_CURVE_MG11_fy6[17] = {
   282.384,
   308.493,
   325.831,
   333.055,
   338.834,
   342.652,
   342.652,
   347.915,
   351.734,
   351.734,
   353.901,
   355.759,
   355.759,
   356.687,
   357.41,
   358.029,
   358.648};
   auto* graph_MG11 = new TGraph(17,E_TOF_CURVE_MG11_fx6,E_TOF_CURVE_MG11_fy6);
   graph_MG11->SetName("E_TOF_CURVE_MG11");
   f = new TFile("E_TOF_CURVE_MG11.root", "recreate");
   graph_MG11->Write();

}
