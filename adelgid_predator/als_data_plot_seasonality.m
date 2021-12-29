%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. tsugae, L. nigrinus, S. tsugae model
% Plot of data to show seasonality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data- with each year of data for each class saved seperately
ll_1_time = [12	14	15];
ll_1 = [3.666666667	1.75	1.666666667];
ll_2_time = 12;
ll_2 = 4;
ll_3_time = [14	15	16	17	18	19];
ll_3 = [3	4.75	7.8	4.6	3	1.5];
% two set of LN adults in data for the beginning and end of calendar year
% (lab used for the end of the calendar year)
lab_0_time = [43	44	46	46.5];
lab_0 = [3	1.25	2	1.5];
la_1_time = [2	9	10	11	12	14];
la_1 = [1	1	2.5	1	1.333333333	1.333333333];
lab_1_time = [40	42	43	46	47	48	49];
lab_1 = [1	1.4	1	1	1	1	1];
la_2_time = [1 11	12	13	14];
la_2 = [1	2.6	2.5	2	1];
lab_2_time = [42	43	44	45	46	47	48	50	51	52];
lab_2 = [3.428571429	4.714285714	5.111111111	3.1	3.571428571	4.333333333	4	3.5	3.25	1];
la_3_time = [1	2	4	5	6	7	8	10	11	12	14	16	17	18	19];
la_3 = [1	3.333333333	1	1.833333333	1.714285714	5	2	1	1	2	2	2	1	1	1];
sl_2_time = [17	18	19	21	23];
sl_2 = [1	1	1	1	1];
sl_3_time = [21	22	23	25	26	27	31];
sl_3 = [2	4.333333333	4	1.5	1.5	2.25	1];
sa_1_time = [10	12	23	24	26	48];
sa_1 = [1	1	1	1	3	1];
sa_2_time = [19	22	23	24	27	28	29	33	36	43	45];
sa_2 = [1	2	2	1	1.5	1.5	1	1	1	1	1];
sa_3_time = [12	14	17	20	22	23	25	26	27	28	29	30	31];
sa_3 = [1	2	1	1.25	1	1	2.25	2.6	1.5	1.666666667	1	1	1];

% plot
figure
subplot(4, 1, 1)
hold on
ax = gca;
ax.FontSize = 13;
% LN larva data in shades of orange, open diamonds
plot(ll_1_time, ll_1, 'd-', 'Color', [0.8, 0.4, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [1, 1, 1]);
plot(ll_2_time, ll_2, 'd-', 'Color', [0.6, 0.2, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.6, 0.2, 0], 'MarkerFaceColor', [1, 1, 1]);
plot(ll_3_time, ll_3, 'd-', 'Color', [0.4, 0, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.4, 0, 0], 'MarkerFaceColor', [1, 1, 1]);
xlabel('Time (weeks in calendar year)', 'FontSize', 13)
ylabel('L. nigrinus larva density', 'FontSize', 13)
xlim([0, 52])
ylim([0, 9])
legend('2011', '2012', '2013')

subplot(4, 1, 2)
hold on
ax = gca;
ax.FontSize = 13;
% LN adult data in shades of orange, closed diamonds
% save plot objects to use in legend, so each year only appears in legend once
p(1) = plot(lab_0_time, lab_0, 'd-', 'Color', [1, 0.6, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [1, 0.6, 0], 'MarkerFaceColor', [1, 0.6, 0]);
p(2) = plot(la_1_time, la_1, 'd-', 'Color', [0.8, 0.4, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [0.8, 0.4, 0]);
p(3) = plot(lab_1_time, lab_1, 'd-', 'Color', [0.8, 0.4, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [0.8, 0.4, 0]);
p(4) = plot(la_2_time, la_2, 'd-', 'Color', [0.6, 0.2, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.6, 0.2, 0], 'MarkerFaceColor', [0.6, 0.2, 0]);
p(5) = plot(lab_2_time, lab_2, 'd-', 'Color', [0.6, 0.2, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.6, 0.2, 0], 'MarkerFaceColor', [0.6, 0.2, 0]);
p(6) = plot(la_3_time, la_3, 'd-', 'Color', [0.4, 0, 0], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0.4, 0, 0], 'MarkerFaceColor', [0.4, 0, 0]);
xlabel('Time (weeks in calendar year)', 'FontSize', 13)
ylabel('L. nigrinus adult density', 'FontSize', 13)
xlim([0, 52])
ylim([0, 6])
legend([p(1), p(2), p(4), p(6)], '2010', '2011', '2012', '2013')

subplot(4, 1, 3)
hold on
ax = gca;
ax.FontSize = 13;
% ST larva data in shades of blue, open diamonds
plot(sl_2_time, sl_2, 'd-', 'Color', [0, 0.45, 0.7], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [1, 1, 1]);
plot(sl_3_time, sl_3, 'd-', 'Color', [0, 0.25, 0.5], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.25, 0.5], 'MarkerFaceColor', [1, 1, 1]);
xlabel('Time (weeks in calendar year)', 'FontSize', 13)
ylabel('S. tsugae larva density', 'FontSize', 13)
xlim([0, 52])
ylim([0, 5])
legend('2012', '2013')

subplot( 4, 1, 4 )
hold on
ax = gca;
ax.FontSize = 13;
% ST adult data in shades of blue, closed diamonds
plot(sa_1_time, sa_1, 'd-', 'Color', [0, 0.65, 0.9], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.65, 0.9], 'MarkerFaceColor', [0, 0.65, 0.9]);%LN larvae data in orange
plot(sa_2_time, sa_2, 'd-', 'Color', [0, 0.45, 0.7], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [0, 0.45, 0.7]);%LN larvae data in orange
plot(sa_3_time, sa_3, 'd-', 'Color', [0, 0.25, 0.5], 'LineWidth', 3, 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.25, 0.5], 'MarkerFaceColor', [0, 0.25, 0.5]);%LN larvae data in orange
xlabel('Time (weeks in calendar year)', 'FontSize', 13)
ylabel('S. tsugae adult density', 'FontSize', 13)
xlim([0, 52])
ylim([0, 4])
legend('2011', '2012', '2013')
