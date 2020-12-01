void DE_E_cuts(){

   auto* cut1 = new TCutG("DE_E_m4_z2_MM",33);
   cut1->SetVarX("Must2Data.CsI_E+Must2Data.SI_E");
   cut1->SetVarY("Must2Data.MugastData.SI_E");
   cut1->SetTitle("Graph");
   cut1->SetFillStyle(1000);
   cut1->SetPoint(0,25.4355,23.7311);
   cut1->SetPoint(1,25.9009,21.9508);
   cut1->SetPoint(2,27.879,18.8447);
   cut1->SetPoint(3,30.2061,16.9508);
   cut1->SetPoint(4,33.5805,14.9432);
   cut1->SetPoint(5,40.9109,12.6326);
   cut1->SetPoint(6,49.8703,10.8144);
   cut1->SetPoint(7,58.5971,9.45076);
   cut1->SetPoint(8,76.3996,7.82197);
   cut1->SetPoint(9,91.7586,6.83712);
   cut1->SetPoint(10,111.074,5.92803);
   cut1->SetPoint(11,125.618,5.35985);
   cut1->SetPoint(12,132.251,5.01894);
   cut1->SetPoint(13,134.345,5.70076);
   cut1->SetPoint(14,132.832,5.96591);
   cut1->SetPoint(15,125.037,6.38258);
   cut1->SetPoint(16,112.703,6.98864);
   cut1->SetPoint(17,99.6709,7.70833);
   cut1->SetPoint(18,99.6709,7.70833);
   cut1->SetPoint(19,86.639,8.76894);
   cut1->SetPoint(20,86.639,8.76894);
   cut1->SetPoint(21,73.9561,9.86742);
   cut1->SetPoint(22,60.4588,11.3068);
   cut1->SetPoint(23,51.9648,12.4811);
   cut1->SetPoint(24,46.0306,13.8068);
   cut1->SetPoint(25,40.7945,15.3598);
   cut1->SetPoint(26,35.3258,18.0871);
   cut1->SetPoint(27,33.4641,20.0568);
   cut1->SetPoint(28,29.9734,23.2008);
   cut1->SetPoint(29,28.6935,24.7538);
   cut1->SetPoint(30,27.2972,24.9432);
   cut1->SetPoint(31,26.7154,24.5644);
   cut1->SetPoint(32,25.4355,23.7311);
   
   auto* cut2 = new TCutG("DE_E_m3_z2_MM",36);
   cut2->SetVarX("Must2Data.CsI_E+Must2Data.SI_E");
   cut2->SetVarY("Must2Data.MugastData.SI_E");
   cut2->SetTitle("Graph");
   cut2->SetFillStyle(1000);
   cut2->SetPoint(0,21.363,21.1553);
   cut2->SetPoint(1,21.363,19.2235);
   cut2->SetPoint(2,22.992,16.6098);
   cut2->SetPoint(3,27.2972,13.3144);
   cut2->SetPoint(4,31.3697,11.572);
   cut2->SetPoint(5,37.5366,9.94318);
   cut2->SetPoint(6,45.0997,8.54167);
   cut2->SetPoint(7,52.1975,7.5947);
   cut2->SetPoint(8,62.7859,6.57197);
   cut2->SetPoint(9,73.4907,5.92803);
   cut2->SetPoint(10,86.5226,5.32197);
   cut2->SetPoint(11,98.8564,4.94318);
   cut2->SetPoint(12,114.099,4.33712);
   cut2->SetPoint(13,130.389,3.73106);
   cut2->SetPoint(14,132.251,4.10985);
   cut2->SetPoint(15,131.087,4.48864);
   cut2->SetPoint(16,129.225,4.79167);
   cut2->SetPoint(17,114.564,5.32197);
   cut2->SetPoint(18,99.5545,5.92803);
   cut2->SetPoint(19,87.9189,6.45833);
   cut2->SetPoint(20,80.3557,6.95076);
   cut2->SetPoint(21,67.5565,7.85985);
   cut2->SetPoint(22,56.7354,9.07197);
   cut2->SetPoint(23,48.3577,10.322);
   cut2->SetPoint(24,41.26,11.8371);
   cut2->SetPoint(25,41.26,11.8371);
   cut2->SetPoint(26,36.8384,12.9735);
   cut2->SetPoint(27,33.115,14.2235);
   cut2->SetPoint(28,30.5552,15.4735);
   cut2->SetPoint(29,27.879,17.1023);
   cut2->SetPoint(30,26.3664,18.8068);
   cut2->SetPoint(31,25.2028,20.8144);
   cut2->SetPoint(32,24.7374,22.5189);
   cut2->SetPoint(33,23.8065,22.5568);
   cut2->SetPoint(34,22.7593,22.1023);
   cut2->SetPoint(35,21.363,21.1553);
   

   auto* cut3 = new TCutG("DE_E_m3_z1_MM",30);
   cut3->SetVarX("Must2Data.CsI_E+Must2Data.SI_E");
   cut3->SetVarY("Must2Data.MugastData.SI_E");
   cut3->SetTitle("Graph");
   cut3->SetFillStyle(1000);
   cut3->SetPoint(0,8.68019,7.36742);
   cut3->SetPoint(1,11.0073,5.92803);
   cut3->SetPoint(2,13.5672,4.82955);
   cut3->SetPoint(3,13.5672,4.82955);
   cut3->SetPoint(4,17.2906,3.80682);
   cut3->SetPoint(5,23.4574,2.93561);
   cut3->SetPoint(6,29.2753,2.48106);
   cut3->SetPoint(7,35.0931,2.14015);
   cut3->SetPoint(8,42.7726,1.83712);
   cut3->SetPoint(9,50.6848,1.53409);
   cut3->SetPoint(10,57.6662,1.42045);
   cut3->SetPoint(11,66.9747,1.19318);
   cut3->SetPoint(12,75.1197,0.92803);
   cut3->SetPoint(13,86.8717,1.00379);
   cut3->SetPoint(14,84.6609,1.23106);
   cut3->SetPoint(15,84.6609,1.23106);
   cut3->SetPoint(16,64.6476,1.72348);
   cut3->SetPoint(17,53.4774,2.02651);
   cut3->SetPoint(18,47.3105,2.17803);
   cut3->SetPoint(19,39.3983,2.51894);
   cut3->SetPoint(20,34.2786,2.74621);
   cut3->SetPoint(21,27.2972,3.23864);
   cut3->SetPoint(22,28.1117,3.20076);
   cut3->SetPoint(23,22.5266,3.8447);
   cut3->SetPoint(24,18.2214,4.60227);
   cut3->SetPoint(25,16.2434,5.17045);
   cut3->SetPoint(26,14.0326,5.89015);
   cut3->SetPoint(27,11.8218,6.98864);
   cut3->SetPoint(28,9.84375,8.35227);
   cut3->SetPoint(29,8.68019,7.36742);
   
   auto* cut4 = new TCutG("DE_E_m2_z1_MM",27);
   cut4->SetVarX("Must2Data.CsI_E+Must2Data.SI_E");
   cut4->SetVarY("Must2Data.MugastData.SI_E");
   cut4->SetTitle("Graph");
   cut4->SetFillStyle(1000);
   cut4->SetPoint(0,7.40027,6.45833);
   cut4->SetPoint(1,8.68019,5.17045);
   cut4->SetPoint(2,11.4727,3.92045);
   cut4->SetPoint(3,16.5924,2.82197);
   cut4->SetPoint(4,23.2247,1.98864);
   cut4->SetPoint(5,31.3697,1.60985);
   cut4->SetPoint(6,40.7945,1.26894);
   cut4->SetPoint(7,49.5213,1.04167);
   cut4->SetPoint(8,59.6443,0.92803);
   cut4->SetPoint(9,63.9495,0.852273);
   cut4->SetPoint(10,66.393,1.07955);
   cut4->SetPoint(11,64.0658,1.19318);
   cut4->SetPoint(12,58.4807,1.3447);
   cut4->SetPoint(13,49.8703,1.53409);
   cut4->SetPoint(14,43.3544,1.68561);
   cut4->SetPoint(15,35.6749,1.95076);
   cut4->SetPoint(16,30.9043,2.21591);
   cut4->SetPoint(17,26.5991,2.51894);
   cut4->SetPoint(18,19.6177,3.23864);
   cut4->SetPoint(19,16.127,3.80682);
   cut4->SetPoint(20,13.5672,4.48864);
   cut4->SetPoint(21,11.5891,5.24621);
   cut4->SetPoint(22,9.96011,6.23106);
   cut4->SetPoint(23,8.9129,6.875);
   cut4->SetPoint(24,8.0984,6.79924);
   cut4->SetPoint(25,8.0984,6.79924);
   cut4->SetPoint(26,7.40027,6.45833);
   
   auto *cut5 = new TCutG("DE_E_m1_z1_MM",21);
   cut5->SetVarX("Must2Data.CsI_E+Must2Data.SI_E");
   cut5->SetVarY("Must2Data.MugastData.SI_E");
   cut5->SetTitle("Graph");
   cut5->SetFillStyle(1000);
   cut5->SetPoint(0,5.18949,4.56439);
   cut5->SetPoint(1,7.28391,3.04924);
   cut5->SetPoint(2,12.0545,1.91288);
   cut5->SetPoint(3,18.3378,1.26894);
   cut5->SetPoint(4,25.5519,0.852273);
   cut5->SetPoint(5,32.4169,0.738636);
   cut5->SetPoint(6,41.8418,0.662879);
   cut5->SetPoint(7,52.4302,0.549242);
   cut5->SetPoint(8,50.8012,0.814394);
   cut5->SetPoint(9,47.5432,0.92803);
   cut5->SetPoint(10,40.4455,1.11742);
   cut5->SetPoint(11,29.2753,1.45833);
   cut5->SetPoint(12,21.7121,1.79924);
   cut5->SetPoint(13,16.2434,2.32955);
   cut5->SetPoint(14,12.2872,3.01136);
   cut5->SetPoint(15,10.0765,3.54167);
   cut5->SetPoint(16,8.79654,4.18561);
   cut5->SetPoint(17,7.16755,5.01894);
   cut5->SetPoint(18,5.77128,4.98106);
   cut5->SetPoint(19,5.07314,4.64015);
   cut5->SetPoint(20,5.18949,4.56439);

   TFile ff("MUGAST.root", "update");
   cut1->Write();
   cut2->Write();
   cut3->Write();
   cut4->Write();
   cut5->Write();
   ff->Write();
   ff->Close();

}
