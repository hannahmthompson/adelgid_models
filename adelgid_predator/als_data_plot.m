%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. tsugae, L. nigrinus, S. tsugae model
% Plot of data (over model time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data
ll_data_time = [21	23	24	73	127	128	129	130	131	132];
ll_data = [3.666666667	1.75	1.666666667	4	3	4.75	7.8	4.6	3	1.5];
la_data_time = [0 1	3	3.5	11	18	19	20	21	23	49	51	52	55	56	57	58	62	72	73	74	75	103	104	105	106	107	108	109	111	112	113	114	115	117	118	119	120	121	123	124	125	127	129	130	131	132];
la_data = [3 1.25	2	1.5	1	1	2.5	1	1.333333333	1.333333333	1	1.4	1	1	1	1	1	1	2.6	2.5	2	1	3.428571429	4.714285714	5.111111111	3.1	3.571428571	4.333333333	4	3.5	3.25	1	1	3.333333333	1	1.833333333	1.714285714	5	2	1	1	2	2	2	1	1	1];
sl_data_time = [78	79	80	82	84	134	135	136	138	139	140	144];
sl_data = [1	1	1	1	1	2	4.333333333	4	1.5	1.5	2.25	1];
sa_data_time = [19	21	32	33	35	57	80	83	84	85	88	89	90	94	97	104	106	125	127	130	133	135	136	138	139	140	141	142	143	144];
sa_data = [1	1	1	1	3	1	1	2	2	1	1.5	1.5	1	1	1	1	1	1	2	1	1.25	1	1	1.75	2.6	1.5	1.666666667	1	1	1];
a_data_time = [14	50];
a_data = [0.519444444	4.4125];

% plot
figure
hold on
ax = gca;
ax.FontSize = 16;
% adelgid data in black
plot(a_data_time ./ 52, a_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0]);
% LN larva data in orange, open diamonds
plot(ll_data_time ./ 52, ll_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'LineWidth', 2, 'MarkerFaceColor', [1, 1, 1]);
% LN adult data in orange, closed diamonds
plot(la_data_time ./ 52, la_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [0.8, 0.4, 0]);
% ST larva data in blue, open diamonds
plot(sl_data_time ./ 52, sl_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'LineWidth', 2, 'MarkerFaceColor', [1, 1, 1]);
% ST adult data in blue, closed diamonds
plot(sa_data_time ./ 52, sa_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [0, 0.45, 0.7]); 
ylim([0, 8.3])
xlabel('Time (years)', 'FontSize', 16)
ylabel('Density', 'FontSize', 16)
% with default plot size and default legend placement, legend may overlap data
legend('Adelgid data, A', 'L. nigirinus larvae data, Ll', 'L. nigrinus adult data, La', 'S. tsugae larvae data, Sl', 'S. tsugae adult data, Sa')
