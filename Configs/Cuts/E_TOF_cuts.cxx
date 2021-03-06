void E_TOF_cuts(){
   auto *cut1 = new TCutG("E_TOF_m2_z1_MG1",32);
   cut1->SetVarX("E vs TOF of MG1");
   cut1->SetVarY("");
   cut1->SetTitle("Graph");
   cut1->SetFillStyle(1000);
   cut1->SetPoint(0,1.88794,332.837);
   cut1->SetPoint(1,2.31269,339.676);
   cut1->SetPoint(2,2.85205,345.687);
   cut1->SetPoint(3,3.16218,348.692);
   cut1->SetPoint(4,3.54647,350.868);
   cut1->SetPoint(5,4.01167,352.526);
   cut1->SetPoint(6,4.6422,354.495);
   cut1->SetPoint(7,5.94289,356.567);
   cut1->SetPoint(8,7.00214,357.811);
   cut1->SetPoint(9,8.40414,359.884);
   cut1->SetPoint(10,9.37543,361.127);
   cut1->SetPoint(11,10.8627,361.127);
   cut1->SetPoint(12,11.3298,360.816);
   cut1->SetPoint(13,11.4404,359.158);
   cut1->SetPoint(14,8.91035,357.086);
   cut1->SetPoint(15,8.19546,356.049);
   cut1->SetPoint(16,7.47027,355.22);
   cut1->SetPoint(17,6.83657,354.184);
   cut1->SetPoint(18,6.10843,353.044);
   cut1->SetPoint(19,5.23872,351.697);
   cut1->SetPoint(20,4.47687,349.832);
   cut1->SetPoint(21,3.64086,346.308);
   cut1->SetPoint(22,3.05431,341.749);
   cut1->SetPoint(23,2.4745,334.599);
   cut1->SetPoint(24,2.03626,327.345);
   cut1->SetPoint(25,1.6722,320.92);
   cut1->SetPoint(26,1.47668,318.018);
   cut1->SetPoint(27,1.42274,317.707);
   cut1->SetPoint(28,1.3351,319.366);
   cut1->SetPoint(29,1.66545,327.863);
   cut1->SetPoint(30,1.82052,331.49);
   cut1->SetPoint(31,1.88794,332.837);


   auto *cut2 = new TCutG("E_TOF_m1_z1_MG1",37);
   cut2->SetVarX("E vs TOF of MG1");
   cut2->SetVarY("");
   cut2->SetTitle("Graph");
   cut2->SetFillStyle(1000);
   cut2->SetPoint(0,3.89449,352.422);
   cut2->SetPoint(1,4.42412,354.391);
   cut2->SetPoint(2,5.10951,355.738);
   cut2->SetPoint(3,6.40242,357.5);
   cut2->SetPoint(4,7.55513,358.951);
   cut2->SetPoint(5,8.04581,359.78);
   cut2->SetPoint(6,8.17822,360.816);
   cut2->SetPoint(7,8.03803,361.749);
   cut2->SetPoint(8,7.00214,361.438);
   cut2->SetPoint(9,5.51452,361.334);
   cut2->SetPoint(10,4.502,361.127);
   cut2->SetPoint(11,3.91007,361.127);
   cut2->SetPoint(12,3.34929,360.816);
   cut2->SetPoint(13,2.80409,360.091);
   cut2->SetPoint(14,2.62495,359.365);
   cut2->SetPoint(15,2.75736,358.536);
   cut2->SetPoint(16,3.08448,358.433);
   cut2->SetPoint(17,3.7543,358.744);
   cut2->SetPoint(18,3.32592,357.396);
   cut2->SetPoint(19,2.88976,354.702);
   cut2->SetPoint(20,2.65611,352.526);
   cut2->SetPoint(21,2.05638,344.858);
   cut2->SetPoint(22,1.86167,341.334);
   cut2->SetPoint(23,1.51118,335.013);
   cut2->SetPoint(24,1.20743,326.619);
   cut2->SetPoint(25,1.02829,319.987);
   cut2->SetPoint(26,1.02829,316.568);
   cut2->SetPoint(27,1.15291,316.671);
   cut2->SetPoint(28,1.25416,317.915);
   cut2->SetPoint(29,1.41772,322.474);
   cut2->SetPoint(30,1.83051,332.837);
   cut2->SetPoint(31,2.26668,339.573);
   cut2->SetPoint(32,2.95986,347.137);
   cut2->SetPoint(33,3.19352,349.314);
   cut2->SetPoint(34,3.42718,350.661);
   cut2->SetPoint(35,3.62189,351.593);
   cut2->SetPoint(36,3.89449,352.422);


   auto *cut3 = new TCutG("E_TOF_m4_z2_MG1",28);
   cut3->SetVarX("E vs TOF of MG1");
   cut3->SetVarY("");
   cut3->SetTitle("Graph");
   cut3->SetFillStyle(1000);
   cut3->SetPoint(0,1.32458,314.599);
   cut3->SetPoint(1,1.70561,319.884);
   cut3->SetPoint(2,2.2833,325.998);
   cut3->SetPoint(3,3.31578,333.044);
   cut3->SetPoint(4,4.68012,339.884);
   cut3->SetPoint(5,5.94613,344.858);
   cut3->SetPoint(6,7.76526,350.35);
   cut3->SetPoint(7,9.74417,356.257);
   cut3->SetPoint(8,11.5018,358.536);
   cut3->SetPoint(9,12.5097,359.78);
   cut3->SetPoint(10,14.7959,360.091);
   cut3->SetPoint(11,15.0786,358.64);
   cut3->SetPoint(12,14.9926,356.049);
   cut3->SetPoint(13,13.788,355.531);
   cut3->SetPoint(14,11.76,353.355);
   cut3->SetPoint(15,10.9119,351.593);
   cut3->SetPoint(16,9.2771,347.966);
   cut3->SetPoint(17,7.5686,343.303);
   cut3->SetPoint(18,5.46677,335.22);
   cut3->SetPoint(19,3.91806,327.863);
   cut3->SetPoint(20,2.77496,319.366);
   cut3->SetPoint(21,1.76707,309.832);
   cut3->SetPoint(22,1.49666,308.485);
   cut3->SetPoint(23,1.34916,308.381);
   cut3->SetPoint(24,1.18937,309.417);
   cut3->SetPoint(25,1.17708,312.112);
   cut3->SetPoint(26,1.25083,313.77);
   cut3->SetPoint(27,1.32458,314.599);


   auto *cut4 = new TCutG("E_TOF_m2_z1_MG3",42);
   cut4->SetVarX("E vs TOF of MG3");
   cut4->SetVarY("");
   cut4->SetTitle("Graph");
   cut4->SetFillStyle(1000);
   cut4->SetPoint(0,1.20081,317.569);
   cut4->SetPoint(1,1.35846,323.419);
   cut4->SetPoint(2,1.6552,330.413);
   cut4->SetPoint(3,2.04468,336.899);
   cut4->SetPoint(4,2.44343,341.477);
   cut4->SetPoint(5,2.87928,346.182);
   cut4->SetPoint(6,3.1853,349.235);
   cut4->SetPoint(7,3.51914,351.778);
   cut4->SetPoint(8,4.06626,354.194);
   cut4->SetPoint(9,4.59484,355.084);
   cut4->SetPoint(10,5.77255,356.865);
   cut4->SetPoint(11,6.8297,358.009);
   cut4->SetPoint(12,8.0445,359.027);
   cut4->SetPoint(13,9.02747,360.68);
   cut4->SetPoint(14,9.37059,361.061);
   cut4->SetPoint(15,10.233,361.443);
   cut4->SetPoint(16,10.2979,361.443);
   cut4->SetPoint(17,11.0954,361.316);
   cut4->SetPoint(18,11.0583,358.645);
   cut4->SetPoint(19,10.1032,358.136);
   cut4->SetPoint(20,9.73224,357.628);
   cut4->SetPoint(21,9.10166,357.119);
   cut4->SetPoint(22,8.3227,356.102);
   cut4->SetPoint(23,7.39537,354.703);
   cut4->SetPoint(24,6.84825,354.067);
   cut4->SetPoint(25,6.14348,353.05);
   cut4->SetPoint(26,5.49435,352.414);
   cut4->SetPoint(27,4.91013,351.142);
   cut4->SetPoint(28,4.27955,349.107);
   cut4->SetPoint(29,3.89007,347.581);
   cut4->SetPoint(30,3.50986,344.911);
   cut4->SetPoint(31,3.11111,342.113);
   cut4->SetPoint(32,2.68454,337.026);
   cut4->SetPoint(33,2.24869,330.541);
   cut4->SetPoint(34,1.75721,323.038);
   cut4->SetPoint(35,1.47901,317.696);
   cut4->SetPoint(36,1.26572,313.627);
   cut4->SetPoint(37,1.19154,311.719);
   cut4->SetPoint(38,1.16372,310.321);
   cut4->SetPoint(39,1.14517,313.373);
   cut4->SetPoint(40,1.15444,315.026);
   cut4->SetPoint(41,1.20081,317.569);

   auto *cut5 = new TCutG("E_TOF_m1_z1_MG3",35);
   cut5->SetVarX("E vs TOF of MG3");
   cut5->SetVarY("");
   cut5->SetTitle("Graph");
   cut5->SetFillStyle(1000);
   cut5->SetPoint(0,1.07724,313.285);
   cut5->SetPoint(1,1.15846,317.295);
   cut5->SetPoint(2,1.27216,321.59);
   cut5->SetPoint(3,1.38587,324.74);
   cut5->SetPoint(4,1.60515,329.608);
   cut5->SetPoint(5,1.65388,331.135);
   cut5->SetPoint(6,2.02749,337.245);
   cut5->SetPoint(7,2.4417,342.208);
   cut5->SetPoint(8,2.99398,348.031);
   cut5->SetPoint(9,3.49753,352.231);
   cut5->SetPoint(10,4.13915,354.713);
   cut5->SetPoint(11,5.03255,356.24);
   cut5->SetPoint(12,6.33204,358.054);
   cut5->SetPoint(13,7.52594,359.104);
   cut5->SetPoint(14,8.35437,359.772);
   cut5->SetPoint(15,8.46807,360.726);
   cut5->SetPoint(16,8.42746,362.063);
   cut5->SetPoint(17,7.84269,362.063);
   cut5->SetPoint(18,6.41326,361.776);
   cut5->SetPoint(19,5.14626,361.585);
   cut5->SetPoint(20,3.7087,361.299);
   cut5->SetPoint(21,2.49043,360.249);
   cut5->SetPoint(22,2.39297,359.104);
   cut5->SetPoint(23,2.53916,357.958);
   cut5->SetPoint(24,3.60312,358.245);
   cut5->SetPoint(25,2.96149,355.285);
   cut5->SetPoint(26,2.47419,351.085);
   cut5->SetPoint(27,1.9219,344.69);
   cut5->SetPoint(28,1.54018,337.817);
   cut5->SetPoint(29,1.30465,331.995);
   cut5->SetPoint(30,1.14221,325.695);
   cut5->SetPoint(31,1.01226,319.967);
   cut5->SetPoint(32,0.987897,316.531);
   cut5->SetPoint(33,0.979776,313.954);
   cut5->SetPoint(34,1.07724,313.285);


   auto *cut6 = new TCutG("E_TOF_m4_z2_MG3",26);
   cut6->SetVarX("E vs TOF of MG3");
   cut6->SetVarY("");
   cut6->SetTitle("Graph");
   cut6->SetFillStyle(1000);
   cut6->SetPoint(0,1.27366,311.338);
   cut6->SetPoint(1,1.58026,319.604);
   cut6->SetPoint(2,1.9812,325.962);
   cut6->SetPoint(3,2.94818,333.72);
   cut6->SetPoint(4,4.15101,340.714);
   cut6->SetPoint(5,6.55666,349.743);
   cut6->SetPoint(6,9.03307,354.703);
   cut6->SetPoint(7,11.368,358.645);
   cut6->SetPoint(8,13.0897,360.171);
   cut6->SetPoint(9,16.9576,361.316);
   cut6->SetPoint(10,24.1038,362.715);
   cut6->SetPoint(11,24.1038,360.553);
   cut6->SetPoint(12,22.2406,358.772);
   cut6->SetPoint(13,17.6415,358.136);
   cut6->SetPoint(14,14.4576,356.865);
   cut6->SetPoint(15,11.7453,354.703);
   cut6->SetPoint(16,10.2831,351.778);
   cut6->SetPoint(17,8.25477,347.073);
   cut6->SetPoint(18,6.06138,339.697);
   cut6->SetPoint(19,4.50478,332.067);
   cut6->SetPoint(20,3.39629,325.708);
   cut6->SetPoint(21,2.73592,319.604);
   cut6->SetPoint(22,1.65102,311.719);
   cut6->SetPoint(23,1.43875,306.251);
   cut6->SetPoint(24,1.08498,306.633);
   cut6->SetPoint(25,1.27366,311.338);

   auto *cut7 = new TCutG("E_TOF_m2_z1_MG4",40);
   cut7->SetVarX("E vs TOF of MG4");
   cut7->SetVarY("");
   cut7->SetTitle("Graph");
   cut7->SetFillStyle(1000);
   cut7->SetPoint(0,2.2393,350.715);
   cut7->SetPoint(1,2.39118,353.988);
   cut7->SetPoint(2,2.80631,359.562);
   cut7->SetPoint(3,3.25182,363.897);
   cut7->SetPoint(4,3.62645,366.904);
   cut7->SetPoint(5,4.0112,368.851);
   cut7->SetPoint(6,4.79084,371.239);
   cut7->SetPoint(7,5.69198,372.655);
   cut7->SetPoint(8,6.56275,373.362);
   cut7->SetPoint(9,7.40314,374.07);
   cut7->SetPoint(10,7.53476,375.22);
   cut7->SetPoint(11,8.34478,375.84);
   cut7->SetPoint(12,10.3597,375.928);
   cut7->SetPoint(13,10.7444,375.928);
   cut7->SetPoint(14,11.2811,375.928);
   cut7->SetPoint(15,12.1215,374.778);
   cut7->SetPoint(16,11.4532,373.893);
   cut7->SetPoint(17,11.19,373.716);
   cut7->SetPoint(18,10.7849,373.716);
   cut7->SetPoint(19,10.2281,373.539);
   cut7->SetPoint(20,8.98266,373.097);
   cut7->SetPoint(21,7.50439,371.859);
   cut7->SetPoint(22,6.19824,370.443);
   cut7->SetPoint(23,5.28698,368.585);
   cut7->SetPoint(24,4.46684,366.639);
   cut7->SetPoint(25,4.04158,364.25);
   cut7->SetPoint(26,3.42394,358.942);
   cut7->SetPoint(27,2.88731,353.104);
   cut7->SetPoint(28,2.5633,347.53);
   cut7->SetPoint(29,2.20892,341.869);
   cut7->SetPoint(30,1.92542,336.384);
   cut7->SetPoint(31,1.61154,329.837);
   cut7->SetPoint(32,1.37866,327.449);
   cut7->SetPoint(33,1.16603,326.741);
   cut7->SetPoint(34,1.13565,328.51);
   cut7->SetPoint(35,1.32803,332.314);
   cut7->SetPoint(36,1.66216,339.215);
   cut7->SetPoint(37,1.94567,345.23);
   cut7->SetPoint(38,2.05705,347.442);
   cut7->SetPoint(39,2.2393,350.715);

   auto *cut8 = new TCutG("E_TOF_m1_z1_MG4",41);
   cut8->SetVarX("E vs TOF of MG4");
   cut8->SetVarY("");
   cut8->SetTitle("Graph");
   cut8->SetFillStyle(1000);
   cut8->SetPoint(0,1.20653,330.81);
   cut8->SetPoint(1,1.35841,334.526);
   cut8->SetPoint(2,1.64191,339.745);
   cut8->SetPoint(3,1.95579,346.203);
   cut8->SetPoint(4,2.28992,353.546);
   cut8->SetPoint(5,2.58355,357.969);
   cut8->SetPoint(6,2.90756,361.685);
   cut8->SetPoint(7,3.33282,365.577);
   cut8->SetPoint(8,3.7682,368.32);
   cut8->SetPoint(9,4.30483,370.532);
   cut8->SetPoint(10,5.1756,372.124);
   cut8->SetPoint(11,5.89449,373.362);
   cut8->SetPoint(12,6.69438,373.893);
   cut8->SetPoint(13,7.30189,374.247);
   cut8->SetPoint(14,7.39301,375.309);
   cut8->SetPoint(15,7.26139,376.105);
   cut8->SetPoint(16,6.68425,376.282);
   cut8->SetPoint(17,5.46923,375.84);
   cut8->SetPoint(18,4.44659,376.016);
   cut8->SetPoint(19,3.6872,375.486);
   cut8->SetPoint(20,2.79618,374.866);
   cut8->SetPoint(21,2.24942,374.336);
   cut8->SetPoint(22,2.13805,373.982);
   cut8->SetPoint(23,2.10767,372.478);
   cut8->SetPoint(24,2.42155,371.77);
   cut8->SetPoint(25,2.90756,372.212);
   cut8->SetPoint(26,3.61632,373.097);
   cut8->SetPoint(27,3.85933,373.097);
   cut8->SetPoint(28,3.62645,372.566);
   cut8->SetPoint(29,3.15056,370.797);
   cut8->SetPoint(30,2.78606,367.524);
   cut8->SetPoint(31,2.5633,364.781);
   cut8->SetPoint(32,2.20892,359.916);
   cut8->SetPoint(33,2.01655,356.289);
   cut8->SetPoint(34,1.70266,350.892);
   cut8->SetPoint(35,1.44954,344.699);
   cut8->SetPoint(36,1.23691,338.684);
   cut8->SetPoint(37,1.12553,335.587);
   cut8->SetPoint(38,1.08503,332.491);
   cut8->SetPoint(39,1.09515,330.987);
   cut8->SetPoint(40,1.20653,330.81);


   auto *cut9 = new TCutG("E_TOF_m2_z1_MG5",40);
   cut9->SetVarX("E vs TOF of MG5");
   cut9->SetVarY("");
   cut9->SetTitle("Graph");
   cut9->SetFillStyle(1000);
   cut9->SetPoint(0,2.76744,343.02);
   cut9->SetPoint(1,2.65178,342.064);
   cut9->SetPoint(2,2.33634,338.556);
   cut9->SetPoint(3,2.04194,333.772);
   cut9->SetPoint(4,1.77907,329.307);
   cut9->SetPoint(5,1.50569,324.524);
   cut9->SetPoint(6,1.37952,321.016);
   cut9->SetPoint(7,1.19025,315.913);
   cut9->SetPoint(8,1.09562,312.618);
   cut9->SetPoint(9,1.07459,310.492);
   cut9->SetPoint(10,1.25334,312.618);
   cut9->SetPoint(11,1.92628,322.929);
   cut9->SetPoint(12,2.33634,330.264);
   cut9->SetPoint(13,2.70436,335.579);
   cut9->SetPoint(14,3.25112,341.213);
   cut9->SetPoint(15,3.70324,343.977);
   cut9->SetPoint(16,4.39721,346.954);
   cut9->SetPoint(17,5.0386,348.442);
   cut9->SetPoint(18,6.1216,350.249);
   cut9->SetPoint(19,7.93011,351.631);
   cut9->SetPoint(20,9.40216,352.906);
   cut9->SetPoint(21,10.1277,352.694);
   cut9->SetPoint(22,10.6639,352.906);
   cut9->SetPoint(23,11.1581,353.013);
   cut9->SetPoint(24,11.9257,353.863);
   cut9->SetPoint(25,11.9257,355.351);
   cut9->SetPoint(26,11.2002,355.883);
   cut9->SetPoint(27,10.2118,355.564);
   cut9->SetPoint(28,9.89634,355.777);
   cut9->SetPoint(29,9.24444,355.564);
   cut9->SetPoint(30,8.56099,355.139);
   cut9->SetPoint(31,7.90908,354.501);
   cut9->SetPoint(32,7.17306,354.076);
   cut9->SetPoint(33,6.37395,353.438);
   cut9->SetPoint(34,5.36455,352.588);
   cut9->SetPoint(35,4.28155,350.568);
   cut9->SetPoint(36,3.6717,348.867);
   cut9->SetPoint(37,3.06185,345.572);
   cut9->SetPoint(38,2.85156,343.658);
   cut9->SetPoint(39,2.76744,343.02);

   auto *cut10 = new TCutG("E_TOF_m1_z1_MG5",35);
   cut10->SetVarX("E vs TOF of MG5");
   cut10->SetVarY("");
   cut10->SetTitle("Graph");
   cut10->SetFillStyle(1000);
   cut10->SetPoint(0,1.0115,310.705);
   cut10->SetPoint(1,1.25334,319.209);
   cut10->SetPoint(2,1.57929,326.756);
   cut10->SetPoint(3,1.99988,333.666);
   cut10->SetPoint(4,2.36789,339.406);
   cut10->SetPoint(5,2.70436,342.914);
   cut10->SetPoint(6,3.0934,346.422);
   cut10->SetPoint(7,3.5981,349.08);
   cut10->SetPoint(8,4.19743,350.674);
   cut10->SetPoint(9,5.28043,352.8);
   cut10->SetPoint(10,6.06903,353.651);
   cut10->SetPoint(11,7.21512,354.288);
   cut10->SetPoint(12,8.30864,355.458);
   cut10->SetPoint(13,8.40327,356.308);
   cut10->SetPoint(14,8.26658,356.733);
   cut10->SetPoint(15,7.07843,356.521);
   cut10->SetPoint(16,5.44867,355.989);
   cut10->SetPoint(17,4.69162,355.989);
   cut10->SetPoint(18,4.0292,356.202);
   cut10->SetPoint(19,3.20906,355.458);
   cut10->SetPoint(20,2.65178,354.82);
   cut10->SetPoint(21,2.68333,353.013);
   cut10->SetPoint(22,3.29317,353.651);
   cut10->SetPoint(23,3.95559,353.863);
   cut10->SetPoint(24,3.34575,351.525);
   cut10->SetPoint(25,2.87259,349.08);
   cut10->SetPoint(26,2.49406,345.784);
   cut10->SetPoint(27,2.09451,341.745);
   cut10->SetPoint(28,1.54775,332.178);
   cut10->SetPoint(29,1.19025,325.374);
   cut10->SetPoint(30,0.969446,318.465);
   cut10->SetPoint(31,0.811727,313.15);
   cut10->SetPoint(32,0.790698,307.941);
   cut10->SetPoint(33,0.916873,307.941);
   cut10->SetPoint(34,1.0115,310.705);


   auto *cut11 = new TCutG("E_TOF_m4_z2_MG5",24);
   cut11->SetVarX("E vs TOF of MG5");
   cut11->SetVarY("");
   cut11->SetTitle("Graph");
   cut11->SetFillStyle(1000);
   cut11->SetPoint(0,1.27051,310.811);
   cut11->SetPoint(1,1.6356,318.039);
   cut11->SetPoint(2,2.07371,324.099);
   cut11->SetPoint(3,3.33936,330.158);
   cut11->SetPoint(4,4.75104,336.217);
   cut11->SetPoint(5,6.40612,342.808);
   cut11->SetPoint(6,8.01252,347.698);
   cut11->SetPoint(7,10.6168,351.631);
   cut11->SetPoint(8,13.9513,354.501);
   cut11->SetPoint(9,18.6488,356.308);
   cut11->SetPoint(10,25.6342,356.733);
   cut11->SetPoint(11,25.9506,354.607);
   cut11->SetPoint(12,23.6627,353.651);
   cut11->SetPoint(13,19.9145,353.651);
   cut11->SetPoint(14,16.1419,352.162);
   cut11->SetPoint(15,12.6127,349.717);
   cut11->SetPoint(16,10.1057,346.209);
   cut11->SetPoint(17,7.86648,341.426);
   cut11->SetPoint(18,5.55424,333.453);
   cut11->SetPoint(19,4.1669,326.225);
   cut11->SetPoint(20,2.77955,319.528);
   cut11->SetPoint(21,2.02503,313.894);
   cut11->SetPoint(22,1.53825,309.535);
   cut11->SetPoint(23,1.27051,310.811);


   auto *cut12 = new TCutG("E_TOF_m2_z1_MG7",46);
   cut12->SetVarX("E vs TOF of MG7");
   cut12->SetVarY("");
   cut12->SetTitle("Graph");
   cut12->SetFillStyle(1000);
   cut12->SetPoint(0,2.08348,295.234);
   cut12->SetPoint(1,2.01584,295.972);
   cut12->SetPoint(2,2.02551,298.608);
   cut12->SetPoint(3,2.24775,304.405);
   cut12->SetPoint(4,2.39269,309.043);
   cut12->SetPoint(5,2.55696,312.943);
   cut12->SetPoint(6,2.7792,319.9);
   cut12->SetPoint(7,3.00144,325.909);
   cut12->SetPoint(8,3.30099,332.55);
   cut12->SetPoint(9,3.61986,337.715);
   cut12->SetPoint(10,4.01603,340.772);
   cut12->SetPoint(11,4.4122,342.564);
   cut12->SetPoint(12,5.06927,344.672);
   cut12->SetPoint(13,6.02588,346.253);
   cut12->SetPoint(14,6.60564,346.885);
   cut12->SetPoint(15,7.57191,347.834);
   cut12->SetPoint(16,8.15168,348.677);
   cut12->SetPoint(17,8.5865,349.626);
   cut12->SetPoint(18,9.21458,350.153);
   cut12->SetPoint(19,10.3725,350.469);
   cut12->SetPoint(20,11.2828,350.891);
   cut12->SetPoint(21,12.9321,351.102);
   cut12->SetPoint(22,13.157,349.942);
   cut12->SetPoint(23,13.0606,348.994);
   cut12->SetPoint(24,12.5787,348.256);
   cut12->SetPoint(25,11.8076,347.94);
   cut12->SetPoint(26,11.0579,347.412);
   cut12->SetPoint(27,10.6081,347.202);
   cut12->SetPoint(28,9.95861,346.991);
   cut12->SetPoint(29,9.24356,346.148);
   cut12->SetPoint(30,8.44156,344.988);
   cut12->SetPoint(31,7.6782,344.145);
   cut12->SetPoint(32,6.85687,343.407);
   cut12->SetPoint(33,5.98722,342.037);
   cut12->SetPoint(34,5.15623,340.455);
   cut12->SetPoint(35,4.72141,339.296);
   cut12->SetPoint(36,4.28658,336.661);
   cut12->SetPoint(37,3.89041,332.339);
   cut12->SetPoint(38,3.59087,326.752);
   cut12->SetPoint(39,3.37829,321.06);
   cut12->SetPoint(40,3.14638,314.63);
   cut12->SetPoint(41,2.75021,304.089);
   cut12->SetPoint(42,2.34438,296.71);
   cut12->SetPoint(43,2.14146,295.234);
   cut12->SetPoint(44,2.1318,295.234);
   cut12->SetPoint(45,2.08348,295.234);

   auto *cut13 = new TCutG("E_TOF_m1_z1_MG7",46);
   cut13->SetVarX("E vs TOF of MG7");
   cut13->SetVarY("");
   cut13->SetTitle("Graph");
   cut13->SetFillStyle(1000);
   cut13->SetPoint(0,1.83721,293.653);
   cut13->SetPoint(1,1.97643,298.397);
   cut13->SetPoint(2,2.24417,305.775);
   cut13->SetPoint(3,2.47977,312.943);
   cut13->SetPoint(4,2.7368,320.322);
   cut13->SetPoint(5,2.9724,326.752);
   cut13->SetPoint(6,3.22942,332.128);
   cut13->SetPoint(7,3.56141,337.609);
   cut13->SetPoint(8,3.86127,340.35);
   cut13->SetPoint(9,4.3539,342.88);
   cut13->SetPoint(10,4.96433,344.777);
   cut13->SetPoint(11,5.51051,345.937);
   cut13->SetPoint(12,6.14236,346.991);
   cut13->SetPoint(13,6.84917,347.518);
   cut13->SetPoint(14,7.57741,348.045);
   cut13->SetPoint(15,8.01649,348.677);
   cut13->SetPoint(16,8.21997,349.415);
   cut13->SetPoint(17,8.23068,350.786);
   cut13->SetPoint(18,7.93082,351.207);
   cut13->SetPoint(19,7.1062,350.891);
   cut13->SetPoint(20,5.88534,351.102);
   cut13->SetPoint(21,4.56809,350.469);
   cut13->SetPoint(22,3.50787,349.837);
   cut13->SetPoint(23,2.92956,349.415);
   cut13->SetPoint(24,2.71538,348.994);
   cut13->SetPoint(25,2.59757,348.256);
   cut13->SetPoint(26,2.58687,347.202);
   cut13->SetPoint(27,2.60828,346.148);
   cut13->SetPoint(28,2.92956,346.042);
   cut13->SetPoint(29,3.78631,346.464);
   cut13->SetPoint(30,4.31107,347.307);
   cut13->SetPoint(31,3.93624,346.042);
   cut13->SetPoint(32,3.64709,345.093);
   cut13->SetPoint(33,3.30439,342.142);
   cut13->SetPoint(34,2.96169,337.715);
   cut13->SetPoint(35,2.59757,330.125);
   cut13->SetPoint(36,2.27629,323.274);
   cut13->SetPoint(37,1.96572,315.368);
   cut13->SetPoint(38,1.76225,306.83);
   cut13->SetPoint(39,1.62303,301.875);
   cut13->SetPoint(40,1.55877,295.656);
   cut13->SetPoint(41,1.63374,291.545);
   cut13->SetPoint(42,1.73012,291.229);
   cut13->SetPoint(43,1.79438,292.178);
   cut13->SetPoint(44,1.79438,292.178);
   cut13->SetPoint(45,1.83721,293.653);

   auto *cut14 = new TCutG("E_TOF_m4_z2_MG7",28);
   cut14->SetVarX("E vs TOF of MG7");
   cut14->SetVarY("");
   cut14->SetTitle("Graph");
   cut14->SetFillStyle(1000);
   cut14->SetPoint(0,2.0911,285.806);
   cut14->SetPoint(1,2.51879,292.867);
   cut14->SetPoint(2,4.32253,308.576);
   cut14->SetPoint(3,6.14487,321.258);
   cut14->SetPoint(4,7.46513,330.481);
   cut14->SetPoint(5,9.77095,340.425);
   cut14->SetPoint(6,12.114,346.334);
   cut14->SetPoint(7,14.6615,349.504);
   cut14->SetPoint(8,16.3537,350.945);
   cut14->SetPoint(9,18.5107,352.242);
   cut14->SetPoint(10,20.519,352.386);
   cut14->SetPoint(11,21.3558,350.225);
   cut14->SetPoint(12,21.5046,347.342);
   cut14->SetPoint(13,20.9467,344.748);
   cut14->SetPoint(14,19.0872,343.739);
   cut14->SetPoint(15,16.019,342.154);
   cut14->SetPoint(16,13.9921,341.29);
   cut14->SetPoint(17,11.9838,337.398);
   cut14->SetPoint(18,9.26888,326.59);
   cut14->SetPoint(19,7.40935,314.34);
   cut14->SetPoint(20,6.07049,304.973);
   cut14->SetPoint(21,4.82461,292.867);
   cut14->SetPoint(22,3.57872,279.465);
   cut14->SetPoint(23,2.76053,275.141);
   cut14->SetPoint(24,2.03532,272.835);
   cut14->SetPoint(25,1.7006,276.726);
   cut14->SetPoint(26,1.83077,280.762);
   cut14->SetPoint(27,2.0911,285.806);

   auto *cut15 = new TCutG("E_TOF_m4_z2_MG4",28);
   cut15->SetVarX("E vs TOF of MG4");
   cut15->SetVarY("");
   cut15->SetTitle("Graph");
   cut15->SetFillStyle(1000);
   cut15->SetPoint(0,1.33817,322.586);
   cut15->SetPoint(1,1.85738,330.612);
   cut15->SetPoint(2,2.91942,339.709);
   cut15->SetPoint(3,4.17026,347.2);
   cut15->SetPoint(4,5.94032,356.831);
   cut15->SetPoint(5,7.09676,361.647);
   cut15->SetPoint(6,9.12643,368.603);
   cut15->SetPoint(7,10.8257,372.348);
   cut15->SetPoint(8,12.5486,374.489);
   cut15->SetPoint(9,15.4987,376.094);
   cut15->SetPoint(10,19.5108,377.164);
   cut15->SetPoint(11,22.5081,376.986);
   cut15->SetPoint(12,22.6497,375.916);
   cut15->SetPoint(13,22.4609,375.024);
   cut15->SetPoint(14,17.8587,373.062);
   cut15->SetPoint(15,13.4218,370.743);
   cut15->SetPoint(16,11.3213,367.711);
   cut15->SetPoint(17,8.56001,360.398);
   cut15->SetPoint(18,6.36514,349.697);
   cut15->SetPoint(19,4.92549,340.957);
   cut15->SetPoint(20,3.72185,333.466);
   cut15->SetPoint(21,2.3294,322.051);
   cut15->SetPoint(22,1.52697,312.42);
   cut15->SetPoint(23,1.22016,305.642);
   cut15->SetPoint(24,1.00775,309.031);
   cut15->SetPoint(25,1.07856,313.312);
   cut15->SetPoint(26,1.24376,319.732);
   cut15->SetPoint(27,1.33817,322.586);

   auto *cut16 = new TCutG("E_TOF_m2_z1_MG11",47);
   cut16->SetVarX("E vs TOF of MG11");
   cut16->SetVarY("");
   cut16->SetTitle("Graph");
   cut16->SetFillStyle(1000);
   cut16->SetPoint(0,0.739561,322.175);
   cut16->SetPoint(1,0.779307,324.173);
   cut16->SetPoint(2,0.842903,326.805);
   cut16->SetPoint(3,0.898548,328.167);
   cut16->SetPoint(4,0.978042,330.437);
   cut16->SetPoint(5,1.06549,332.616);
   cut16->SetPoint(6,1.15293,334.431);
   cut16->SetPoint(7,1.32781,336.701);
   cut16->SetPoint(8,1.54245,339.697);
   cut16->SetPoint(9,1.70144,341.876);
   cut16->SetPoint(10,1.82068,342.965);
   cut16->SetPoint(11,2.26584,345.416);
   cut16->SetPoint(12,2.55997,347.141);
   cut16->SetPoint(13,2.98924,348.594);
   cut16->SetPoint(14,3.4185,349.774);
   cut16->SetPoint(15,4.126,350.954);
   cut16->SetPoint(16,4.81759,352.316);
   cut16->SetPoint(17,5.82717,353.496);
   cut16->SetPoint(18,7.25011,354.404);
   cut16->SetPoint(19,8.52996,354.858);
   cut16->SetPoint(20,8.76939,356.583);
   cut16->SetPoint(21,9.1977,356.946);
   cut16->SetPoint(22,9.75416,357.218);
   cut16->SetPoint(23,10.838,357.309);
   cut16->SetPoint(24,11.107,356.855);
   cut16->SetPoint(25,11.0881,355.221);
   cut16->SetPoint(26,10.0136,355.039);
   cut16->SetPoint(27,9.18413,354.222);
   cut16->SetPoint(28,8.25968,353.405);
   cut16->SetPoint(29,7.9099,353.405);
   cut16->SetPoint(30,7.06727,352.679);
   cut16->SetPoint(31,5.85896,351.499);
   cut16->SetPoint(32,4.63476,349.683);
   cut16->SetPoint(33,3.76033,347.504);
   cut16->SetPoint(34,3.13233,346.052);
   cut16->SetPoint(35,2.58382,343.51);
   cut16->SetPoint(36,2.1625,341.512);
   cut16->SetPoint(37,1.78093,338.517);
   cut16->SetPoint(38,1.5027,334.431);
   cut16->SetPoint(39,1.36756,331.889);
   cut16->SetPoint(40,1.20062,328.803);
   cut16->SetPoint(41,1.10523,325.716);
   cut16->SetPoint(42,0.985991,322.629);
   cut16->SetPoint(43,0.803156,320.087);
   cut16->SetPoint(44,0.723662,321.086);
   cut16->SetPoint(45,0.707763,321.449);
   cut16->SetPoint(46,0.739561,322.175);

   auto* cut17 = new TCutG("E_TOF_m1_z1_MG11",49);
   cut17->SetVarX("E vs TOF of MG11");
   cut17->SetVarY("");
   cut17->SetTitle("Graph");
   cut17->SetFillStyle(1000);
   cut17->SetPoint(0,0.880085,327.804);
   cut17->SetPoint(1,0.946065,330.255);
   cut17->SetPoint(2,1.08745,333.523);
   cut17->SetPoint(3,1.21941,335.793);
   cut17->SetPoint(4,1.37022,337.609);
   cut17->SetPoint(5,1.53046,339.969);
   cut17->SetPoint(6,1.74725,342.602);
   cut17->SetPoint(7,1.91691,343.691);
   cut17->SetPoint(8,2.32222,346.052);
   cut17->SetPoint(9,2.76522,348.14);
   cut17->SetPoint(10,3.61354,350.591);
   cut17->SetPoint(11,4.09425,351.226);
   cut17->SetPoint(12,4.8483,352.679);
   cut17->SetPoint(13,5.5081,353.224);
   cut17->SetPoint(14,5.96996,353.859);
   cut17->SetPoint(15,6.69574,354.404);
   cut17->SetPoint(16,7.06334,354.495);
   cut17->SetPoint(17,7.53463,354.676);
   cut17->SetPoint(18,8.13787,355.039);
   cut17->SetPoint(19,8.43949,355.221);
   cut17->SetPoint(20,8.43949,356.31);
   cut17->SetPoint(21,8.26041,356.764);
   cut17->SetPoint(22,7.99649,357.672);
   cut17->SetPoint(23,7.62888,357.854);
   cut17->SetPoint(24,7.16703,357.491);
   cut17->SetPoint(25,6.88425,357.4);
   cut17->SetPoint(26,6.29044,357.127);
   cut17->SetPoint(27,5.5081,357.037);
   cut17->SetPoint(28,4.87658,356.401);
   cut17->SetPoint(29,4.16023,356.129);
   cut17->SetPoint(30,3.84918,356.22);
   cut17->SetPoint(31,3.26479,355.584);
   cut17->SetPoint(32,2.86891,355.039);
   cut17->SetPoint(33,2.60499,354.495);
   cut17->SetPoint(34,2.19968,353.405);
   cut17->SetPoint(35,2.08657,352.225);
   cut17->SetPoint(36,2.1337,351.68);
   cut17->SetPoint(37,2.66154,352.588);
   cut17->SetPoint(38,2.99144,352.679);
   cut17->SetPoint(39,2.56729,351.317);
   cut17->SetPoint(40,2.22796,349.502);
   cut17->SetPoint(41,1.85093,347.323);
   cut17->SetPoint(42,1.52103,343.964);
   cut17->SetPoint(43,1.31367,341.059);
   cut17->SetPoint(44,1.02147,336.247);
   cut17->SetPoint(45,0.889511,332.343);
   cut17->SetPoint(46,0.814105,328.53);
   cut17->SetPoint(47,0.785828,326.351);
   cut17->SetPoint(48,0.880085,327.804);

   auto* cut18 = new TCutG("E_TOF_m4_z2_MG11",31);
   cut18->SetVarX("E vs TOF of MG11");
   cut18->SetVarY("");
   cut18->SetTitle("Graph");
   cut18->SetFillStyle(1000);
   cut18->SetPoint(0,1.30396,329.347);
   cut18->SetPoint(1,1.51986,334.159);
   cut18->SetPoint(2,1.79745,337.7);
   cut18->SetPoint(3,1.93625,339.243);
   cut18->SetPoint(4,2.59939,342.965);
   cut18->SetPoint(5,3.27795,345.598);
   cut18->SetPoint(6,4.28036,347.777);
   cut18->SetPoint(7,5.35989,350.137);
   cut18->SetPoint(8,6.42399,351.317);
   cut18->SetPoint(9,7.65774,352.861);
   cut18->SetPoint(10,9.77053,354.313);
   cut18->SetPoint(11,11.2202,355.13);
   cut18->SetPoint(12,12.6081,356.583);
   cut18->SetPoint(13,13.4101,357.218);
   cut18->SetPoint(14,18.2063,358.852);
   cut18->SetPoint(15,18.5764,358.035);
   cut18->SetPoint(16,18.6072,357.037);
   cut18->SetPoint(17,18.2525,355.312);
   cut18->SetPoint(18,12.7161,354.041);
   cut18->SetPoint(19,10.4028,352.407);
   cut18->SetPoint(20,7.99702,350.046);
   cut18->SetPoint(21,6.03845,348.049);
   cut18->SetPoint(22,4.75844,345.689);
   cut18->SetPoint(23,4.06446,344.327);
   cut18->SetPoint(24,3.2471,342.148);
   cut18->SetPoint(25,2.49143,338.789);
   cut18->SetPoint(26,1.99794,335.611);
   cut18->SetPoint(27,1.61239,331.617);
   cut18->SetPoint(28,1.38106,328.621);
   cut18->SetPoint(29,1.38106,328.621);
   cut18->SetPoint(30,1.30396,329.347);

   auto* cut19 = new TCutG("E_TOF_m2_z1_MG1_2",23);
   cut19->SetVarX("MugastData.SI_E");
   cut19->SetVarY("MugastData.T");
   cut19->SetTitle("Graph");
   cut19->SetFillStyle(1000);
   cut19->SetPoint(0,1.36758,319.929);
   cut19->SetPoint(1,1.47375,323.268);
   cut19->SetPoint(2,1.69347,328.674);
   cut19->SetPoint(3,1.89402,333.509);
   cut19->SetPoint(4,2.14471,337.498);
   cut19->SetPoint(5,2.4883,341.851);
   cut19->SetPoint(6,2.913,346.634);
   cut19->SetPoint(7,3.16369,348.739);
   cut19->SetPoint(8,3.68423,351.637);
   cut19->SetPoint(9,4.0352,352.793);
   cut19->SetPoint(10,4.2269,353.417);
   cut19->SetPoint(11,4.41713,353.105);
   cut19->SetPoint(12,4.41713,351);
   cut19->SetPoint(13,4.25639,350.168);
   cut19->SetPoint(14,3.6341,345.997);
   cut19->SetPoint(15,3.12387,342.436);
   cut19->SetPoint(16,2.81272,338.551);
   cut19->SetPoint(17,2.43079,333.924);
   cut19->SetPoint(18,2.14471,329.506);
   cut19->SetPoint(19,1.78343,321.943);
   cut19->SetPoint(20,1.61384,319.539);
   cut19->SetPoint(21,1.38233,319.318);
   cut19->SetPoint(22,1.36758,319.929);

   auto* cut20 = new TCutG("E_TOF_m2_z1_MG3_2",22);
   cut20->SetVarX("MugastData.SI_E");
   cut20->SetVarY("MugastData.T");
   cut20->SetTitle("Graph");
   cut20->SetFillStyle(1000);
   cut20->SetPoint(0,1.65466,322.42);
   cut20->SetPoint(1,2.01533,328.723);
   cut20->SetPoint(2,2.28754,332.126);
   cut20->SetPoint(3,2.6346,336.16);
   cut20->SetPoint(4,3.12457,341.58);
   cut20->SetPoint(5,3.61454,346.37);
   cut20->SetPoint(6,3.64857,347.378);
   cut20->SetPoint(7,3.77106,349.773);
   cut20->SetPoint(8,3.83911,351.79);
   cut20->SetPoint(9,3.67579,352.672);
   cut20->SetPoint(10,3.36275,351.16);
   cut20->SetPoint(11,2.96805,347.756);
   cut20->SetPoint(12,2.53252,343.218);
   cut20->SetPoint(13,2.16504,338.555);
   cut20->SetPoint(14,1.8384,333.513);
   cut20->SetPoint(15,1.54561,328.14);
   cut20->SetPoint(16,1.28394,321.72);
   cut20->SetPoint(17,1.10424,314.143);
   cut20->SetPoint(18,1.09663,310.193);
   cut20->SetPoint(19,1.15788,307.798);
   cut20->SetPoint(20,1.23954,310.95);
   cut20->SetPoint(21,1.65466,322.42);

   auto* cut21 = new TCutG("E_TOF_m2_z1_MG4_2",29);
   cut21->SetVarX("MugastData.SI_E");
   cut21->SetVarY("MugastData.T");
   cut21->SetTitle("Graph");
   cut21->SetFillStyle(1000);
   cut21->SetPoint(0,1.4437,328.393);
   cut21->SetPoint(1,1.54577,330.588);
   cut21->SetPoint(2,1.54577,330.588);
   cut21->SetPoint(3,1.73632,333.708);
   cut21->SetPoint(4,1.73632,333.708);
   cut21->SetPoint(5,2.09699,341.218);
   cut21->SetPoint(6,2.30795,345.147);
   cut21->SetPoint(7,2.54613,348.845);
   cut21->SetPoint(8,2.54613,348.845);
   cut21->SetPoint(9,2.91361,353.697);
   cut21->SetPoint(10,2.91361,353.697);
   cut21->SetPoint(11,3.26067,357.511);
   cut21->SetPoint(12,3.62815,360.515);
   cut21->SetPoint(13,3.88675,364.097);
   cut21->SetPoint(14,4.02966,365.945);
   cut21->SetPoint(15,3.94119,366.754);
   cut21->SetPoint(16,3.94119,366.754);
   cut21->SetPoint(17,3.38317,362.017);
   cut21->SetPoint(18,2.92423,357.514);
   cut21->SetPoint(19,2.6544,354.36);
   cut21->SetPoint(20,2.4642,351.994);
   cut21->SetPoint(21,2.19907,348.729);
   cut21->SetPoint(22,1.93367,344.8);
   cut21->SetPoint(23,1.93367,344.8);
   cut21->SetPoint(24,1.52536,337.405);
   cut21->SetPoint(25,1.37564,333.939);
   cut21->SetPoint(26,1.21913,330.011);
   cut21->SetPoint(27,1.21913,330.011);
   cut21->SetPoint(28,1.4437,328.393);

   auto* cut22 = new TCutG("E_TOF_m2_z1_MG5_2",29);
   cut22->SetVarX("MugastData.SI_E");
   cut22->SetVarY("MugastData.T");
   cut22->SetTitle("Graph");
   cut22->SetFillStyle(1000);
   cut22->SetPoint(0,1.28524,314.357);
   cut22->SetPoint(1,1.61619,318.979);
   cut22->SetPoint(2,1.61619,318.979);
   cut22->SetPoint(3,1.91562,323.286);
   cut22->SetPoint(4,2.15989,327.382);
   cut22->SetPoint(5,2.51447,332.529);
   cut22->SetPoint(6,2.51447,332.529);
   cut22->SetPoint(7,2.86117,337.887);
   cut22->SetPoint(8,2.86117,337.887);
   cut22->SetPoint(9,3.34183,341.668);
   cut22->SetPoint(10,3.34183,341.668);
   cut22->SetPoint(11,3.87765,344.084);
   cut22->SetPoint(12,4.34255,346.5);
   cut22->SetPoint(13,4.34255,346.5);
   cut22->SetPoint(14,4.31891,348.496);
   cut22->SetPoint(15,4.31891,348.496);
   cut22->SetPoint(16,3.86977,347.55);
   cut22->SetPoint(17,3.46791,345.765);
   cut22->SetPoint(18,2.97937,342.718);
   cut22->SetPoint(19,2.82178,341.773);
   cut22->SetPoint(20,2.41203,337.361);
   cut22->SetPoint(21,2.41203,337.361);
   cut22->SetPoint(22,1.99441,332.95);
   cut22->SetPoint(23,1.99441,332.95);
   cut22->SetPoint(24,1.79742,330.218);
   cut22->SetPoint(25,1.58467,325.807);
   cut22->SetPoint(26,1.41132,322.235);
   cut22->SetPoint(27,1.31676,318.454);
   cut22->SetPoint(28,1.28524,314.357);

    TFile ff("MUGAST.root", "recreate");
    cut1->Write(); 
    cut2->Write(); 
    cut3->Write(); 
    cut4->Write(); 
    cut5->Write(); 
    cut6->Write(); 
    cut7->Write(); 
    cut8->Write(); 
    cut9->Write(); 
    cut10->Write(); 
    cut11->Write(); 
    cut12->Write(); 
    cut13->Write(); 
    cut14->Write(); 
    cut15->Write(); 
    cut16->Write(); 
    cut17->Write(); 
    cut18->Write(); 
    cut19->Write(); 
    cut20->Write(); 
    cut21->Write(); 
    cut22->Write(); 

    ff.Write();
    ff.Close();
}
