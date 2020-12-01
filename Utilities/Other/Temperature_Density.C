#include <map>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TFrame.h>
#include <TPad.h>
#include <TGraph2D.h>
#include <TGraphDelaunay.h>
#include <TMultiGraph.h>


void Temperature_Density(){
    
    std::map<double,TGraph*> graph;

    graph.emplace(3.33, new TGraph());
    auto it = graph.find(3.33);

    it->second->SetPoint(1	,0.298	,880.61);
    it->second->SetPoint(2	,0.328	,784.44);
    it->second->SetPoint(3	,0.384	,646.17);
    it->second->SetPoint(4	,0.445	,564.87);
    it->second->SetPoint(5	,0.503	,475.11);
    it->second->SetPoint(6	,0.558	,421.95);
    it->second->SetPoint(7	,0.618	,365.03);
    it->second->SetPoint(8	,0.673	,323.88);
    it->second->SetPoint(9	,0.726	,290.32);
    it->second->SetPoint(10	,0.772	,263.65);
    it->second->SetPoint(11	,0.807	,247.60);
    it->second->SetPoint(12	,0.851	,224.28);
    it->second->SetPoint(13	,0.899	,205.31);
    it->second->SetPoint(14	,0.949	,183.44);
    it->second->SetPoint(15	,0.992	,168.10);
    it->second->SetPoint(16	,1.027	,153.61);
    it->second->SetPoint(17	,1.056	,143.15);
    it->second->SetPoint(18	,1.085	,130.66);
    it->second->SetPoint(19	,1.118	,116.44);
    it->second->SetPoint(20	,1.144	,99.95);
    it->second->SetPoint(21	,1.154	,87.92);
    it->second->SetPoint(22	,1.156	,81.64);
    it->second->SetPoint(23	,1.157	,77.27);
    it->second->SetPoint(24	,1.159	,72.44);
    it->second->SetPoint(25	,1.161	,68.20);
    it->second->SetPoint(26	,1.165	,63.56);
    it->second->SetPoint(27	,1.172	,59.94);
    it->second->SetPoint(28	,1.199	,55.45);
    it->second->SetPoint(29	,1.227	,53.38);
    it->second->SetPoint(30	,1.272	,51.16);
    it->second->SetPoint(31	,1.333	,49.41);
    it->second->SetPoint(32	,1.408	,47.92);
    it->second->SetPoint(33	,1.495	,46.72);
    it->second->SetPoint(34	,1.654	,45.05);
    it->second->SetPoint(35	,1.755	,44.22);
    it->second->SetPoint(36	,1.943	,43.82);
    it->second->SetPoint(37	,2.198	,41.73);
    it->second->SetPoint(38	,2.590	,40.58);
    it->second->SetPoint(39	,2.992	,39.32);
    it->second->SetPoint(40	,3.418	,38.22);
    it->second->SetPoint(41	,3.800	,37.57);
    it->second->SetPoint(42	,4.069	,37.07);
    it->second->SetPoint(43	,5.79	,34.83);
    it->second->RemovePoint(0);


    it = graph.emplace(3.5, new TGraph()).first;

	it->second->SetPoint(1	,0.277	,1000.36);
	it->second->SetPoint(2	,0.322	,868.69);
	it->second->SetPoint(3	,0.411	,658.44);
	it->second->SetPoint(4	,0.488	,533.41);
	it->second->SetPoint(5	,0.559	,447.16);
	it->second->SetPoint(6	,0.629	,393.09);
	it->second->SetPoint(7	,0.693	,339.18);
	it->second->SetPoint(8	,0.748	,304.50);
	it->second->SetPoint(9	,0.808	,276.33);
	it->second->SetPoint(10	,0.869	,246.79);
	it->second->SetPoint(11	,0.931	,219.90);
	it->second->SetPoint(12	,0.984	,200.68);
	it->second->SetPoint(13	,1.036	,183.18);
	it->second->SetPoint(14	,1.093	,166.54);
	it->second->SetPoint(15	,1.147	,146.64);
	it->second->SetPoint(16	,1.174	,138.13);
	it->second->SetPoint(17	,1.206	,127.94);
	it->second->SetPoint(18	,1.243	,117.90);
	it->second->SetPoint(19	,1.274	,107.88);
	it->second->SetPoint(20	,1.306	,98.02);
	it->second->SetPoint(21	,1.331	,87.54);
	it->second->SetPoint(22	,1.349	,78.99);
	it->second->SetPoint(23	,1.373	,73.76);
	it->second->SetPoint(24	,1.381	,68.34);
	it->second->SetPoint(25	,1.400	,64.38);
	it->second->SetPoint(26	,1.422	,59.92);
	it->second->SetPoint(27	,1.458	,57.00);
	it->second->SetPoint(28	,1.487	,55.17);
	it->second->SetPoint(29	,1.649	,49.85);
	it->second->SetPoint(30	,1.816	,47.09);
	it->second->SetPoint(31	,2.216	,43.70);
	it->second->SetPoint(32	,2.622	,41.69);
	it->second->SetPoint(33	,3.004	,40.43);
	it->second->SetPoint(34	,3.404	,39.41);
	it->second->SetPoint(35	,3.810	,38.34);
	it->second->SetPoint(36	,4.183	,37.71);
	it->second->SetPoint(37	,4.427	,37.24);
	it->second->SetPoint(38	,4.672	,36.88);
	it->second->SetPoint(39	,4.931	,36.52);
	it->second->SetPoint(40	,5.70	,35.64);
	it->second->SetPoint(41	,7.05	,34.11);
	it->second->SetPoint(42	,8.65	,32.92);
	it->second->SetPoint(43	,10.20	,32.01);
	it->second->SetPoint(44	,11.70	,31.12);
	it->second->SetPoint(45	,13.30	,30.55);
	it->second->SetPoint(46	,14.85	,29.96);
	it->second->SetPoint(47	,16.45	,29.40);
	it->second->SetPoint(48	,20.20	,28.35);
	it->second->SetPoint(49	,23.50	,27.66);
	it->second->SetPoint(50	,26.55	,27.02);
	it->second->SetPoint(51	,30.00	,26.49);
	it->second->SetPoint(52	,33.95	,26.01);
	it->second->SetPoint(53	,36.25	,25.55);
	it->second->SetPoint(54	,39.65	,25.15);
	it->second->SetPoint(55	,43.10	,24.76);
	it->second->SetPoint(56	,46.20	,24.45);
	it->second->SetPoint(57	,49.20	,24.16);
	it->second->SetPoint(58	,52.30	,23.89);
	it->second->SetPoint(59	,55.85	,23.61);
	it->second->SetPoint(60	,59.10	,23.34);
	it->second->SetPoint(61	,62.30	,23.07);
	it->second->SetPoint(62	,65.55	,22.88);
	it->second->SetPoint(63	,68.70	,22.67);
	it->second->SetPoint(64	,71.70	,22.49);
	it->second->SetPoint(65	,74.90	,22.29);
	it->second->SetPoint(66	,78.30	,22.10);
	it->second->SetPoint(67	,81.35	,21.94);
	it->second->SetPoint(68	,84.75	,21.77);
	it->second->SetPoint(69	,87.95	,21.61);
	it->second->SetPoint(70	,91.20	,21.46);
	it->second->SetPoint(71	,94.50	,21.32);
	it->second->SetPoint(72	,96.80	,21.20);
    it->second->RemovePoint(0);

    it = graph.emplace(3.75, new TGraph()).first;

	it->second->SetPoint(1	,0.485	,580.75);
	it->second->SetPoint(2	,0.535	,527.44);
	it->second->SetPoint(3	,0.645	,413.51);
	it->second->SetPoint(4	,0.776	,326.02);
	it->second->SetPoint(5	,0.888	,272.12);
	it->second->SetPoint(6	,0.984	,234.56);
	it->second->SetPoint(7	,1.050	,214.24);
	it->second->SetPoint(8	,1.105	,196.79);
	it->second->SetPoint(9	,1.168	,177.27);
	it->second->SetPoint(10	,1.258	,156.30);
	it->second->SetPoint(11	,1.457	,112.12);
	it->second->SetPoint(12	,1.339	,137.87);
	it->second->SetPoint(13	,1.399	,123.42);
	it->second->SetPoint(14	,1.512	,101.98);
	it->second->SetPoint(15	,1.541	,96.53);
	it->second->SetPoint(16	,1.580	,88.59);
	it->second->SetPoint(17	,1.599	,84.51);
	it->second->SetPoint(18	,1.614	,82.30);
	it->second->SetPoint(19	,1.628	,79.31);
	it->second->SetPoint(20	,1.645	,77.14);
	it->second->SetPoint(21	,1.667	,74.07);
	it->second->SetPoint(22	,1.679	,71.66);
	it->second->SetPoint(23	,1.715	,67.30);
	it->second->SetPoint(24	,1.802	,60.24);
	it->second->SetPoint(25	,1.880	,55.60);
	it->second->SetPoint(26	,1.977	,52.72);
	it->second->SetPoint(27	,2.190	,48.81);
	it->second->SetPoint(28	,2.606	,44.70);
	it->second->SetPoint(29	,3.027	,42.43);
	it->second->SetPoint(30	,3.400	,41.00);
	it->second->SetPoint(31	,3.806	,39.88);
	it->second->SetPoint(32	,4.193	,38.84);
	it->second->SetPoint(33	,4.617	,38.04);
	it->second->SetPoint(34	,5.03	,37.25);
	it->second->SetPoint(35	,5.46	,36.67);
	it->second->SetPoint(36	,6.10	,35.81);
	it->second->SetPoint(37	,7.20	,34.69);
	it->second->SetPoint(38	,8.00	,33.94);
	it->second->SetPoint(39	,8.95	,33.20);
	it->second->SetPoint(40	,10.20	,32.41);
	it->second->SetPoint(41	,11.70	,31.55);
	it->second->SetPoint(42	,13.20	,30.88);
	it->second->SetPoint(43	,14.85	,30.21);
	it->second->SetPoint(44	,16.30	,29.66);
	it->second->SetPoint(45	,17.50	,29.32);
	it->second->SetPoint(46	,20.15	,28.58);
	it->second->SetPoint(47	,23.50	,27.81);
	it->second->SetPoint(48	,26.80	,27.15);
	it->second->SetPoint(49	,30.10	,26.59);
	it->second->SetPoint(50	,33.30	,26.08);
	it->second->SetPoint(51	,36.60	,25.62);
	it->second->SetPoint(52	,39.85	,25.23);
	it->second->SetPoint(53	,43.20	,24.86);
	it->second->SetPoint(54	,46.30	,24.53);
	it->second->SetPoint(55	,49.45	,24.22);
	it->second->SetPoint(56	,52.70	,23.94);
	it->second->SetPoint(57	,55.95	,23.65);
	it->second->SetPoint(58	,59.15	,23.39);
	it->second->SetPoint(59	,62.50	,23.12);
	it->second->SetPoint(60	,65.70	,22.91);
	it->second->SetPoint(61	,68.90	,22.70);
	it->second->SetPoint(62	,72.00	,22.50);
	it->second->SetPoint(63	,74.95	,22.31);
	it->second->SetPoint(64	,78.55	,22.13);
	it->second->SetPoint(65	,81.65	,21.95);
	it->second->SetPoint(66	,84.75	,21.80);
	it->second->SetPoint(67	,88.15	,21.64);
	it->second->SetPoint(68	,91.15	,21.49);
	it->second->SetPoint(69	,92.40	,21.45);
    it->second->RemovePoint(0);

    it = graph.emplace(4.21, new TGraph()).first;

	it->second->SetPoint(1	,0.398	,796.46);
	it->second->SetPoint(2	,0.506	,610.23);
	it->second->SetPoint(3	,0.636	,475.11);
	it->second->SetPoint(4	,0.783	,388.64);
	it->second->SetPoint(5	,0.984	,296.30);
	it->second->SetPoint(6	,1.080	,263.65);
	it->second->SetPoint(7	,1.136	,238.14);
	it->second->SetPoint(8	,1.206	,219.53);
	it->second->SetPoint(9	,1.300	,197.38);
	it->second->SetPoint(10	,1.378	,180.87);
	it->second->SetPoint(11	,1.468	,163.27);
	it->second->SetPoint(12	,1.562	,147.44);
	it->second->SetPoint(13	,1.741	,123.42);
	it->second->SetPoint(14	,1.741	,125.15);
	it->second->SetPoint(15	,1.562	,148.50);
	it->second->SetPoint(16	,1.649	,135.63);
	it->second->SetPoint(17	,1.836	,112.12);
	it->second->SetPoint(18	,1.949	,98.98);
	it->second->SetPoint(19	,2.046	,87.54);
	it->second->SetPoint(20	,2.152	,79.58);
	it->second->SetPoint(21	,2.246	,72.96);
	it->second->SetPoint(22	,2.351	,67.10);
	it->second->SetPoint(23	,2.496	,61.00);
	it->second->SetPoint(24	,2.620	,57.64);
	it->second->SetPoint(25	,2.820	,52.88);
	it->second->SetPoint(26	,3.014	,50.35);
	it->second->SetPoint(27	,3.411	,46.49);
	it->second->SetPoint(28	,3.817	,43.80);
	it->second->SetPoint(29	,4.214	,42.16);
	it->second->SetPoint(30	,4.621	,40.88);
	it->second->SetPoint(31	,5.04	,39.92);
	it->second->SetPoint(32	,5.28	,39.32);
	it->second->SetPoint(33	,5.55	,38.67);
	it->second->SetPoint(34	,5.79	,38.22);
	it->second->SetPoint(35	,6.40	,37.20);
	it->second->SetPoint(36	,7.30	,35.98);
	it->second->SetPoint(37	,8.30	,34.89);
	it->second->SetPoint(38	,9.25	,34.06);
	it->second->SetPoint(39	,10.45	,33.22);
	it->second->SetPoint(40	,11.75	,32.33);
	it->second->SetPoint(41	,13.30	,31.50);
	it->second->SetPoint(42	,14.85	,30.80);
	it->second->SetPoint(43	,16.30	,30.19);
	it->second->SetPoint(44	,18.70	,29.41);
	it->second->SetPoint(45	,21.10	,28.77);
	it->second->SetPoint(46	,23.80	,28.08);
	it->second->SetPoint(47	,26.80	,27.47);
	it->second->SetPoint(48	,30.00	,26.86);
	it->second->SetPoint(49	,33.20	,26.33);
	it->second->SetPoint(50	,36.50	,25.86);
	it->second->SetPoint(51	,39.85	,25.44);
	it->second->SetPoint(52	,43.10	,25.04);
	it->second->SetPoint(53	,46.30	,24.71);
	it->second->SetPoint(54	,49.40	,24.40);
	it->second->SetPoint(55	,52.70	,24.09);
	it->second->SetPoint(56	,55.75	,23.82);
	it->second->SetPoint(57	,59.05	,23.54);
	it->second->SetPoint(58	,62.20	,23.27);
	it->second->SetPoint(59	,65.55	,23.03);
	it->second->SetPoint(60	,68.55	,22.83);
	it->second->SetPoint(61	,71.95	,22.63);
	it->second->SetPoint(62	,75.15	,22.43);
	it->second->SetPoint(63	,78.35	,22.26);
	it->second->SetPoint(64	,81.60	,22.06);
	it->second->SetPoint(65	,84.70	,21.91);
	it->second->SetPoint(66	,88.00	,21.73);
	it->second->SetPoint(67	,91.10	,21.59);
	it->second->SetPoint(68	,94.25	,21.44);
	it->second->SetPoint(69	,95.15	,21.40);
    it->second->RemovePoint(0);


    it = graph.emplace(4.52, new TGraph()).first;

	it->second->SetPoint(1	,0.412	,898.98);
	it->second->SetPoint(2	,0.664	,507.74);
	it->second->SetPoint(3	,1.048	,297.15);
	it->second->SetPoint(4	,1.195	,252.19);
	it->second->SetPoint(5	,1.63	,166.59);
	it->second->SetPoint(6	,1.99	,121.33);
	it->second->SetPoint(7	,2.43	,84.54);
	it->second->SetPoint(8	,2.91	,62.95);
	it->second->SetPoint(9	,3.44	,52.70);
	it->second->SetPoint(10	,3.93	,47.68);
	it->second->SetPoint(11	,4.97	,42.35);
	it->second->SetPoint(12	,5.53	,40.74);
	it->second->SetPoint(13	,6.89	,37.81);
	it->second->SetPoint(14	,7.73	,36.61);
	it->second->SetPoint(15	,9.29	,34.88);
	it->second->SetPoint(16	,10.39	,33.93);
	it->second->SetPoint(17	,11.96	,32.83);
	it->second->SetPoint(18	,13.39	,31.97);
	it->second->SetPoint(19	,13.70	,31.80);
	it->second->SetPoint(20	,15.12	,31.15);
	it->second->SetPoint(21	,17.67	,30.10);
	it->second->SetPoint(22	,20.79	,29.06);
	it->second->SetPoint(23	,23.50	,28.30);
	it->second->SetPoint(24	,26.34	,27.72);
	it->second->SetPoint(25	,29.70	,26.98);
	it->second->SetPoint(26	,33.21	,26.37);
	it->second->SetPoint(27	,33.35	,26.42);
	it->second->SetPoint(28	,36.60	,25.87);
	it->second->SetPoint(29	,36.70	,25.90);
	it->second->SetPoint(30	,39.10	,25.60);
	it->second->SetPoint(31	,39.50	,25.48);
	it->second->SetPoint(32	,43.45	,25.03);
	it->second->SetPoint(33	,46.30	,24.74);
	it->second->SetPoint(34	,49.85	,24.40);
	it->second->SetPoint(35	,52.95	,24.13);
	it->second->SetPoint(36	,55.60	,23.79);
	it->second->SetPoint(37	,56.05	,23.84);
	it->second->SetPoint(38	,59.25	,23.59);
	it->second->SetPoint(39	,62.50	,23.34);
	it->second->SetPoint(40	,64.90	,23.05);
	it->second->SetPoint(41	,65.85	,23.10);
	it->second->SetPoint(42	,68.85	,22.87);
	it->second->SetPoint(43	,72.00	,22.68);
	it->second->SetPoint(44	,74.30	,22.42);
	it->second->SetPoint(45	,75.20	,22.46);
	it->second->SetPoint(46	,78.45	,22.28);
	it->second->SetPoint(47	,81.60	,22.10);
	it->second->SetPoint(48	,84.40	,21.94);
	it->second->SetPoint(49	,87.40	,21.81);
	it->second->SetPoint(50	,91.55	,21.60);
	it->second->SetPoint(51	,97.05	,21.32);
    it->second->RemovePoint(0);

    it = graph.emplace(4.97, new TGraph()).first;

	it->second->SetPoint(1	,0.475	,821.84);
	it->second->SetPoint(2	,0.703	,544.09);
	it->second->SetPoint(3	,1.096	,325.70);
	it->second->SetPoint(4	,1.325	,262.85);
	it->second->SetPoint(5	,1.64	,192.46);
	it->second->SetPoint(6	,1.94	,158.08);
	it->second->SetPoint(7	,2.41	,114.78);
	it->second->SetPoint(8	,2.96	,83.83);
	it->second->SetPoint(9	,3.44	,67.02);
	it->second->SetPoint(10	,3.94	,57.40);
	it->second->SetPoint(11	,4.68	,49.30);
	it->second->SetPoint(12	,5.34	,45.33);
	it->second->SetPoint(13	,6.33	,41.71);
	it->second->SetPoint(14	,7.35	,39.31);
	it->second->SetPoint(15	,8.29	,37.67);
	it->second->SetPoint(16	,9.23	,36.46);
	it->second->SetPoint(17	,11.09	,34.61);
	it->second->SetPoint(18	,12.78	,33.29);
	it->second->SetPoint(19	,14.35	,32.28);
	it->second->SetPoint(20	,16.13	,31.42);
	it->second->SetPoint(21	,18.56	,30.40);
	it->second->SetPoint(22	,21.74	,29.27);
	it->second->SetPoint(23	,25.44	,28.25);
	it->second->SetPoint(24	,28.39	,27.57);
	it->second->SetPoint(25	,31.06	,27.04);
	it->second->SetPoint(26	,34.40	,26.46);
	it->second->SetPoint(27	,37.52	,25.96);
	it->second->SetPoint(28	,43.31	,25.21);
	it->second->SetPoint(29	,47.50	,24.76);
	it->second->SetPoint(30	,57.60	,23.82);
	it->second->SetPoint(31	,62.15	,23.47);
	it->second->SetPoint(32	,67.10	,23.06);
	it->second->SetPoint(33	,71.90	,22.75);
	it->second->SetPoint(34	,76.61	,22.44);
	it->second->SetPoint(35	,81.40	,22.12);
	it->second->SetPoint(36	,86.10	,21.87);
	it->second->SetPoint(37	,91.28	,21.58);
    it->second->RemovePoint(0);

    it = graph.emplace(5.42, new TGraph()).first;

	it->second->SetPoint(1	,0.487	,867.44);
	it->second->SetPoint(2	,0.633	,656.86);
	it->second->SetPoint(3	,0.802	,504.85);
	it->second->SetPoint(4	,0.999	,401.35);
	it->second->SetPoint(5	,1.104	,363.01);
	it->second->SetPoint(6	,1.23	,317.91);
	it->second->SetPoint(7	,1.47	,255.46);
	it->second->SetPoint(8	,1.71	,214.78);
	it->second->SetPoint(9	,1.97	,185.52);
	it->second->SetPoint(10	,2.44	,135.85);
	it->second->SetPoint(11	,2.92	,106.66);
	it->second->SetPoint(12	,3.40	,84.95);
	it->second->SetPoint(13	,3.89	,71.33);
	it->second->SetPoint(14	,4.38	,61.69);
	it->second->SetPoint(15	,4.86	,55.83);
	it->second->SetPoint(16	,5.36	,51.28);
	it->second->SetPoint(17	, 5.80	, 48.36);
	it->second->SetPoint(18	,6.55	,44.90);
	it->second->SetPoint(19	,7.31	,42.35);
	it->second->SetPoint(20	,8.29	,40.04);
	it->second->SetPoint(21	,9.47	,37.94);
	it->second->SetPoint(22	,10.71	,36.28);
	it->second->SetPoint(23	,12.83	,34.40);
	it->second->SetPoint(24	,14.77	,33.00);
	it->second->SetPoint(25	,17.14	,31.71);
	it->second->SetPoint(26	,19.84	,30.60);
	it->second->SetPoint(27	,20.28	,30.39);
	it->second->SetPoint(28	,23.48	,29.29);
	it->second->SetPoint(29	,26.84	,28.36);
	it->second->SetPoint(30	,30.30	,27.58);
	it->second->SetPoint(31	,33.80	,26.87);
	it->second->SetPoint(32	,39.90	,25.92);
	it->second->SetPoint(33	,44.15	,25.35);
	it->second->SetPoint(34	,48.20	,24.92);
	it->second->SetPoint(35	,53.05	,34.42);
	it->second->SetPoint(36	,57.80	,23.99);
	it->second->SetPoint(37	,62.80	,23.56);
	it->second->SetPoint(38	,67.25	,23.22);
	it->second->SetPoint(39	,71.90	,22.89);
	it->second->SetPoint(40	,76.75	,22.53);
	it->second->SetPoint(41	,81.60	,22.27);
	it->second->SetPoint(42	,86.35	,21.99);
	it->second->SetPoint(43	,92.15	,21.68);
	it->second->SetPoint(44	,96.70	,21.43);
    it->second->RemovePoint(0);

    it = graph.emplace(5.93, new TGraph()).first;

	it->second->SetPoint(1	,0.341	,1416.4);
	it->second->SetPoint(2	,0.692	,656.86);
	it->second->SetPoint(3	,0.989	,450.70);
	it->second->SetPoint(4	,1.126	,392.84);
	it->second->SetPoint(5	,1.17	,377.31);
	it->second->SetPoint(6	,1.56	,272.41);
	it->second->SetPoint(7	,1.93	,213.80);
	it->second->SetPoint(8	,2.42	,162.15);
	it->second->SetPoint(9	,2.79	,134.57);
	it->second->SetPoint(10	,3.80	,88.02);
	it->second->SetPoint(11	,4.84	,66.18);
	it->second->SetPoint(12	,5.78	,55.54);
	it->second->SetPoint(13	,6.80	,48.68);
	it->second->SetPoint(14	, 7.67	, 44.99);
	it->second->SetPoint(15	,7.76	,44.59);
	it->second->SetPoint(16	,8.75	,41.80);
	it->second->SetPoint(17	,9.71	,39.95);
	it->second->SetPoint(18	,10.58	,38.48);
	it->second->SetPoint(19	,11.65	,37.11);
	it->second->SetPoint(20	,13.79	,35.01);
	it->second->SetPoint(21	,15.20	,33.92);
	it->second->SetPoint(22	,17.57	,32.47);
	it->second->SetPoint(23	,19.89	,31.35);
	it->second->SetPoint(24	,21.67	,30.60);
	it->second->SetPoint(25	,24.01	,29.76);
	it->second->SetPoint(26	,26.60	,29.02);
	it->second->SetPoint(27	,29.86	,28.15);
	it->second->SetPoint(28	,33.30	,27.42);
	it->second->SetPoint(29	,36.45	,26.86);
	it->second->SetPoint(30	,39.18	,26.36);
	it->second->SetPoint(31	,42.90	,25.81);
	it->second->SetPoint(32	,42.95	,25.84);
	it->second->SetPoint(33	,46.25	,25.43);
	it->second->SetPoint(34	,46.25	,25.36);
	it->second->SetPoint(35	,49.50	,25.03);
	it->second->SetPoint(36	,49.50	,25.07);
	it->second->SetPoint(37	,52.40	,24.78);
	it->second->SetPoint(38	,52.60	,24.73);
	it->second->SetPoint(39	,55.85	,24.42);
	it->second->SetPoint(40	,59.25	,24.04);
	it->second->SetPoint(41	,62.65	,23.80);
	it->second->SetPoint(42	,65.60	,23.55);
	it->second->SetPoint(43	,68.60	,23.32);
	it->second->SetPoint(44	,71.90	,23.08);
	it->second->SetPoint(45	,75.30	,22.89);
	it->second->SetPoint(46	,78.25	,22.67);
	it->second->SetPoint(47	,81.60	,22.50);
	it->second->SetPoint(48	,84.45	,22.28);
	it->second->SetPoint(49	,87.70	,22.12);
	it->second->SetPoint(50	,91.30	,21.93);
	it->second->SetPoint(51	,94.60	,21.76);
	it->second->SetPoint(52	,97.75	,21.60);
    it->second->RemovePoint(0);


    it = graph.emplace(6.94, new TGraph()).first;

	it->second->SetPoint(1	,0.255	,2572.0);
	it->second->SetPoint(2	,0.329	,1839.6);
	it->second->SetPoint(3	, 0.623	, 892.86);
	it->second->SetPoint(4	,0.873	,622.06);
	it->second->SetPoint(5	,1.172	,451.85);
	it->second->SetPoint(6	,1.43	,365.56);
	it->second->SetPoint(7	,1.97	,257.60);
	it->second->SetPoint(8	,2.44	,202.81);
	it->second->SetPoint(9	,2.88	,167.18);
	it->second->SetPoint(10	,3.88	,115.81);
	it->second->SetPoint(11	,4.85	,89.51);
	it->second->SetPoint(12	,5.80	,72.67);
	it->second->SetPoint(13	,6.79	,61.86);
	it->second->SetPoint(14	,7.75	,54.76);
	it->second->SetPoint(15	,8.66	,50.07);
	it->second->SetPoint(16	,9.72	,46.11);
	it->second->SetPoint(17	,12.03	,40.71);
	it->second->SetPoint(18	,14.49	,37.42);
	it->second->SetPoint(19	,20.10	,33.04);
	it->second->SetPoint(20	,23.41	,31.44);
	it->second->SetPoint(21	,26.66	,30.27);
	it->second->SetPoint(22	,30.05	,29.29);
	it->second->SetPoint(23	,33.25	,28.40);
	it->second->SetPoint(24	,36.40	,27.74);
	it->second->SetPoint(25	,40.05	,27.03);
	it->second->SetPoint(26	,43.15	,26.50);
	it->second->SetPoint(27	,43.15	,26.41);
	it->second->SetPoint(28	,46.40	,25.98);
	it->second->SetPoint(29	,46.40	,26.02);
	it->second->SetPoint(30	,48.45	,25.75);
	it->second->SetPoint(31	,49.55	,25.49);
	it->second->SetPoint(32	,52.85	,25.08);
	it->second->SetPoint(33	,55.90	,24.90);
	it->second->SetPoint(34	,59.30	,24.49);
	it->second->SetPoint(35	,62.50	,24.25);
	it->second->SetPoint(36	,65.70	,23.99);
	it->second->SetPoint(37	,69.15	,23.64);
	it->second->SetPoint(38	,71.70	,23.48);
	it->second->SetPoint(39	,75.10	,23.25);
	it->second->SetPoint(40	,78.50	,22.98);
	it->second->SetPoint(41	,81.60	,22.73);
	it->second->SetPoint(42	,84.90	,22.58);
	it->second->SetPoint(43	,88.10	,22.39);
	it->second->SetPoint(44	,91.25	,22.21);
	it->second->SetPoint(45	, 94.30	, 22.04);
	it->second->SetPoint(46	,97.70	, 21.84);
    it->second->RemovePoint(0);

    it = graph.emplace(8.99, new TGraph()).first;

	it->second->SetPoint(1	,0.683	,1036.0);
	it->second->SetPoint(2	,1.083	,648.67);
	it->second->SetPoint(3	,1.635	,423.02);
	it->second->SetPoint(4	,2.30	,301.76);
	it->second->SetPoint(5	,3.11	,217.13);
	it->second->SetPoint(6	,3.91	,169.05);
	it->second->SetPoint(7	,4.88	,132.95);
	it->second->SetPoint(8	,5.85	,109.35);
	it->second->SetPoint(9	,6.84	,92.08);
	it->second->SetPoint(10	,7.76	,80.46);
	it->second->SetPoint(11	,8.88	,70.51);
	it->second->SetPoint(12	,9.81	,64.05);
	it->second->SetPoint(13	,11.22	,56.85);
	it->second->SetPoint(14	,13.95	,47.97);
	it->second->SetPoint(15	,17.04	,42.04);
	it->second->SetPoint(16	,20.79	,37.84);
	it->second->SetPoint(17	,25.46	,34.43);
	it->second->SetPoint(18	,25.39	,34.46);
	it->second->SetPoint(19	,28.92	,32.59);
	it->second->SetPoint(20	,32.30	,31.29);
	it->second->SetPoint(21	,40.15	,28.95);
	it->second->SetPoint(22	,40.15	,28.99);
	it->second->SetPoint(23	,44.70	,28.00);
	it->second->SetPoint(24	,44.70	,27.96);
	it->second->SetPoint(25	,49.60	,27.03);
	it->second->SetPoint(26	,50.10	,27.00);
	it->second->SetPoint(27	,52.40	,26.63);
	it->second->SetPoint(28	,56.00	,26.14);
	it->second->SetPoint(29	,59.15	,25.68);
	it->second->SetPoint(30	,62.45	,25.36);
	it->second->SetPoint(31	,65.65	,24.99);
	it->second->SetPoint(32	,68.60	,24.74);
	it->second->SetPoint(33	,71.95	,24.38);
	it->second->SetPoint(34	,75.40	,24.09);
	it->second->SetPoint(35	,78.45	,23.87);
	it->second->SetPoint(36	,81.35	,23.66);
	it->second->SetPoint(37	,84.50	,23.39);
	it->second->SetPoint(38	,86.70	,23.24);
	it->second->SetPoint(39	,88.90 	,23.09 );
    it->second->RemovePoint(0);

    it = graph.emplace(13., new TGraph()).first;

	it->second->SetPoint(1	,0.658	,1656.8);
	it->second->SetPoint(2	,0.822	,1251.9);
	it->second->SetPoint(3	,1.074	,973.53);
	it->second->SetPoint(4	,1.16	,911.75);
	it->second->SetPoint(5	,2.06	,500.91);
	it->second->SetPoint(6	,3.44	,299.75);
	it->second->SetPoint(7	,4.24	,242.85);
	it->second->SetPoint(8	,5.05	,204.66);
	it->second->SetPoint(9	,5.86	,177.22);
	it->second->SetPoint(10	,7.31	,140.69);
	it->second->SetPoint(11	,8.73	,118.17);
	it->second->SetPoint(12	,9.95	,103.14);
	it->second->SetPoint(13	,11.68	,87.91);
	it->second->SetPoint(14	,13.62	,76.25);
	it->second->SetPoint(15	,15.70	,67.22);
	it->second->SetPoint(16	,17.94	,59.83);
	it->second->SetPoint(17	,20.27	,54.13);
	it->second->SetPoint(18	,23.79	,47.95);
	it->second->SetPoint(19	,26.86	,43.97);
	it->second->SetPoint(20	,33.00	,38.73);
	it->second->SetPoint(21	,36.75	,36.52);
	it->second->SetPoint(22	,43.72	,33.40);
	it->second->SetPoint(23	,47.80	,32.02);
	it->second->SetPoint(24	,51.90	,30.84);
	it->second->SetPoint(25	,56.70	,29.69);
	it->second->SetPoint(26	,56.70	,29.80);
	it->second->SetPoint(27	,61.85	,28.72);
	it->second->SetPoint(28	,66.10	,27.96);
	it->second->SetPoint(29	,66.31	,27.88);
	it->second->SetPoint(30	,69.10	,27.45);
	it->second->SetPoint(31	,72.25	,27.01);
	it->second->SetPoint(32	,75.55	,26.50);
	it->second->SetPoint(33	,78.45	,26.15);
	it->second->SetPoint(34	,81.75	,25.78);
	it->second->SetPoint(35	,85.00	,25.46);
	it->second->SetPoint(36	,88.00	,25.14);
	it->second->SetPoint(37	,90.50	,24.92);
	it->second->SetPoint(38	,93.50	,24.68);
	it->second->SetPoint(39	,97.05	,24.41);
    it->second->RemovePoint(0);
    


    std::map<double,TGraph*> graph2;

    it = graph2.emplace(0.987, new TGraph()).first;

	it->second->SetPoint(1	,4.28	,291.13);
	it->second->SetPoint(2	,4.41	,306.95);
	it->second->SetPoint(3	,4.48	,313.09);
	it->second->SetPoint(4	,4.63	,327.97);
	it->second->SetPoint(5	,4.89	,353.07);
	it->second->SetPoint(6	,5.18	,384.33);
	it->second->SetPoint(7	,5.43	,405.05);
	it->second->SetPoint(8	,5.71	,430.38);
	it->second->SetPoint(9	,5.94	,452.25);
	it->second->SetPoint(10	,6.30	,482.64);
	it->second->SetPoint(11	,6.96	,540.68);
	it->second->SetPoint(12	,7.83	,614.58);
	it->second->SetPoint(13	,8.82	,710.04);
	it->second->SetPoint(14	,9.87	,800.26);
	it->second->SetPoint(15	,10.92	,894.18);
	it->second->SetPoint(16	,11.83	,973.53);
	it->second->SetPoint(17	,13.01	,1074.8);
	it->second->SetPoint(18	,14	    ,1191.1);
	it->second->SetPoint(19	,15.83	,1254.5);
    it->second->RemovePoint(0);

    it = graph2.emplace(1.926, new TGraph()).first;

	it->second->SetPoint(1	,4.30	,110.60);
	it->second->SetPoint(2	,4.40	,118.52);
	it->second->SetPoint(3	,4.62	,135.31);
	it->second->SetPoint(4	,4.87	,152.70);
	it->second->SetPoint(5	,4.97	,158.85);
	it->second->SetPoint(6	,5.18	,171.85);
	it->second->SetPoint(7	,5.42	,185.68);
	it->second->SetPoint(8	,5.70	,200.84);
	it->second->SetPoint(9	,5.94	,212.84);
	it->second->SetPoint(10	,6.30	,231.82);
	it->second->SetPoint(11	,6.97	,264.86);
	it->second->SetPoint(12	,7.84	,305.67);
	it->second->SetPoint(13	,8.84	,352.12);
	it->second->SetPoint(14	,9.90	,403.82);
	it->second->SetPoint(15	,10.95	,455.05);
	it->second->SetPoint(16	,11.88	,484.18);
	it->second->SetPoint(17	,13.00	,546.95);
	it->second->SetPoint(18	,14.55	,607.40);
	it->second->SetPoint(19	,15.26	,639.73);
    it->second->RemovePoint(0);

    it = graph2.emplace(2.87, new TGraph()).first;

	it->second->SetPoint(1	,4.32	,55.45);
	it->second->SetPoint(2	,4.83	,80.27);
	it->second->SetPoint(3	,5.38	,106.04);
	it->second->SetPoint(4	,5.82	,123.56);
	it->second->SetPoint(5	,6.35	,143.90);
	it->second->SetPoint(6	,6.89	,164.30);
	it->second->SetPoint(7	,7.42	,183.77);
	it->second->SetPoint(8	,7.91	,200.21);
	it->second->SetPoint(9	,9.03	,237.91);
	it->second->SetPoint(10	,9.96	,270.26);
	it->second->SetPoint(11	,11.04	,309.24);
	it->second->SetPoint(12	,11.98	,336.19);
	it->second->SetPoint(13	,12.84	,362.40);
	it->second->SetPoint(14	,13.52	,382.77);
	it->second->SetPoint(15	,15.36	,439.41);
	it->second->SetPoint(16	,17.97	,518.83);
	it->second->SetPoint(17	,20.38	,591.22);
    it->second->RemovePoint(0);

    it = graph2.emplace(5.55, new TGraph()).first;

	it->second->SetPoint(1	,4.23	,38.78);
	it->second->SetPoint(2	,4.43	,39.95);
	it->second->SetPoint(3	,4.50	,40.37);
	it->second->SetPoint(4	,4.64	,41.38);
	it->second->SetPoint(5	,4.79	,42.87);
	it->second->SetPoint(6	,4.88	,43.75);
	it->second->SetPoint(7	,5.00	,44.71);
	it->second->SetPoint(8	,5.22	,47.04);
	it->second->SetPoint(9	,5.46	,50.11);
	it->second->SetPoint(10	,5.74	,54.13);
	it->second->SetPoint(11	,5.98	,58.16);
	it->second->SetPoint(12	,6.35	,64.71);
	it->second->SetPoint(13	,7.02	,77.84);
	it->second->SetPoint(14	,7.89	,95.46);
	it->second->SetPoint(15	,8.89	,113.77);
	it->second->SetPoint(16	,9.96	,133.05);
	it->second->SetPoint(17	,11.04	,151.72);
	it->second->SetPoint(18	,11.98	,167.68);
	it->second->SetPoint(19	,13.16	,189.08);
	it->second->SetPoint(20	,13.86	,199.36);
	it->second->SetPoint(21	,15.52	,226.94);
	it->second->SetPoint(22	,17.18	,253.510);
    it->second->RemovePoint(0);


    it = graph2.emplace(9.97, new TGraph()).first;

	it->second->SetPoint(1	,4.29	,33.65);
	it->second->SetPoint(2	,5.45	,37.29);
	it->second->SetPoint(3	,5.73	,38.50);
	it->second->SetPoint(4	,6.43	,42.04);
	it->second->SetPoint(5	,6.98	,45.59);
	it->second->SetPoint(6	,7.54	,49.77);
	it->second->SetPoint(7	,8.03	,53.98);
	it->second->SetPoint(8	,9.20	,65.37);
	it->second->SetPoint(9	,10.03	,73.81);
	it->second->SetPoint(10	,11.26	,86.34);
	it->second->SetPoint(11	,12.26	,95.92);
	it->second->SetPoint(12	,13.42	,106.79);
	it->second->SetPoint(13	,13.95	,111.66);
	it->second->SetPoint(14	,15.97	,130.23);
    it->second->RemovePoint(0);

    it = graph2.emplace(15.91, new TGraph()).first;

	it->second->SetPoint(1	,4.30	,30.54);
	it->second->SetPoint(2	,4.42	,30.69);
	it->second->SetPoint(3	,4.49	,30.78);
	it->second->SetPoint(4	,4.62	,30.97);
	it->second->SetPoint(5	,4.82	,31.22);
	it->second->SetPoint(6	,4.89	,31.36);
	it->second->SetPoint(7	,4.98	,31.52);
	it->second->SetPoint(8	,5.21	,31.88);
	it->second->SetPoint(9	,5.46	,32.31);
	it->second->SetPoint(10	,5.74	,32.84);
	it->second->SetPoint(11	,5.98	,33.30);
	it->second->SetPoint(12	,6.34	,34.28);
	it->second->SetPoint(13	,7.16	,36.55);
	it->second->SetPoint(14	,7.92	,39.30);
	it->second->SetPoint(15	,8.93	,43.61);
	it->second->SetPoint(16	,9.98	,48.92);
	it->second->SetPoint(17	,11.08	,55.24);
	it->second->SetPoint(18	,12.01	,60.70);
	it->second->SetPoint(19	,13.20	,67.40);
	it->second->SetPoint(20	,14.28	,73.73);
	it->second->SetPoint(21	,15.58	,81.46);
	it->second->SetPoint(22	,17.24	,91.30);
    it->second->RemovePoint(0);

    const double gmol = 3.016029;
    bool use_rho = true;

    auto*c1 = new TCanvas();
    c1->SetLogy();
    c1->SetLogx();
    Color_t col = 1;

    auto* gr = new TMultiGraph("gr", "Volume at different temperatures;P[atm];V[cm^3/mol]");
    auto* gr2 = new TMultiGraph("gr2", "Volume at different pressures;T[K];V[cm^3/mol]");

    if (use_rho){
        gr->SetTitle("Density at different temperatures;P[atm];rho[g/cm^3]");
        gr2->SetTitle("Density at different pressures;T[K];rho[g/cm^3]");
        for(auto& it2: graph){
            for (int i=0; i<it2.second->GetN(); ++i){
                it2.second->SetPoint(i, it2.second->GetPointX(i), gmol/it2.second->GetPointY(i));
            }   
        }
        for(auto& it2: graph2){
            for (int i=0; i<it2.second->GetN(); ++i){
                it2.second->SetPoint(i, it2.second->GetPointX(i), gmol/it2.second->GetPointY(i));
            }   
        }
    }

    auto* gr_tot = new TGraph2D();
    auto* gr_2d = new TGraph2D();
    for(auto& it2: graph){
        col+=2;
        it2.second->SetMarkerColor(col);
        it2.second->SetTitle(Form("T = %f K", it2.first));
        it2.second->GetYaxis()->SetTitleOffset(1.5);
        
        gr->Add(it2.second);
        for (int i=0; i<it2.second->GetN(); ++i){
            gr_2d->SetPoint(gr_2d->GetN(),it2.second->GetPointX(i),it2.second->GetPointY(i),it2.first);
            gr_tot->SetPoint(gr_tot->GetN(),it2.second->GetPointX(i),it2.second->GetPointY(i),it2.first);
        }
    }

    gr->Draw("ap*");
    c1->BuildLegend();

    auto*c2 = new TCanvas();
    c2->SetLogy();
    c2->SetLogx();
    col = 1;
    auto* gr2_2d = new TGraph2D();
    for(auto& it2: graph2){
        col+=2;
        it2.second->SetMarkerColor(col);
        it2.second->SetTitle(Form("P = %f atm", it2.first));
        it2.second->GetYaxis()->SetTitleOffset(1.5);
        gr2->Add(it2.second);
        for (int i=0; i<it2.second->GetN(); ++i){
            gr2_2d->SetPoint(gr2_2d->GetN(),it2.second->GetPointX(i),it2.second->GetPointY(i),it2.first);
            gr_tot->SetPoint(gr_tot->GetN(),it2.first, it2.second->GetPointY(i),it2.second->GetPointX(i));
        }
    }
    gr2->Draw("ap*");
    c2->BuildLegend();

    auto* c3 = new TCanvas();
    c3->SetLogy();
    c3->SetLogx();
    c3->SetLogz();

    gr_2d->Draw();

    auto* c4 = new TCanvas();
    c4->SetLogy();
    c4->SetLogx();
    c4->SetLogz();
    gr2_2d->Draw();
    
    auto* c5 = new TCanvas();
    c5->SetLogy();
    c5->SetLogx();
    c5->SetLogz();
    gr_tot->Draw();
    //gr_tot->Draw();
    
    gStyle->SetPalette(kLightTemperature);
    auto* c6 = new TCanvas();
    c6->Divide(1,2, 0.01, 0.01);
    auto* pad = c6->cd(1);
    pad->SetLogy();
    pad->SetLogx();
    gr->Draw("apl* PMC PLC");
    auto * legend = c1->BuildLegend(0.6, 0.6, pad->GetFrame()->GetX2(), pad->GetFrame()->GetY1());
    legend->Draw();
    pad = c6->cd(2);
    pad->SetLogy();
    pad->SetLogx();
    gr2->Draw("apl* PMC PLC");
    legend = pad->BuildLegend(0.6, 0, 1,0.4);
    legend->Draw();


//    int nsteps =1000;
//    double startx = 0;
//    double starty = 1;
//    double endx = 100;
//    double endy = 2000;    
//    double stepx = (endx-startx)/nsteps;
//    double stepy = (endy-starty)/nsteps;
//    auto* gr_final = new TGraph2D();
//    TGraphDelaunay interp(gr_tot);
//    for(int i=0; i<nsteps; ++i){
//        for(int j=0; j<nsteps; ++j){
//            double res = interp.ComputeZ(startx+i*stepx, starty+j*stepy);
//            if (res >0)
//
//                gr_final->SetPoint(gr_final->GetN(), startx+i*stepx, starty+j*stepy, res);
//        }
//    }
//    auto* c6 = new TCanvas();
//    c6->SetLogy();
//    c6->SetLogx();
//    c6->SetLogz();
//    gr_final->Draw("surf"); 
}   
