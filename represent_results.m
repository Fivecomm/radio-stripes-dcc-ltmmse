% This sript generates figures from the paper.
%
% This is version 1.0 (Last edited: 2025-04-29)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.
%%
close all; clear; clc;

fprintf('\t\t\t\t LP-MMSE (DCC) \t LTMMSE (All) \t Proposed\n')

scenarios = {'scenario_05', 'scenario_07a', 'scenario_07', ...
    'scenario_07b', 'scenario_05b'};

for sc = 1:length(scenarios)

    load(['results/', scenarios{sc}, '.mat'], 'SE_LPMMSE_DCC', ...
        'SE_LTMMSE_all', 'SE_LTMMSE_DCC', 'R_LTMMSE_all', ...
        'time_LPMMSE_DCC', 'time_LTMMSE_all', 'time_LTMMSE_DCC');
    
    [K, nbrOfSetups] = size(SE_LPMMSE_DCC); 
    
    figure_handle = figure;
    set(figure_handle, 'Position', [100, 100, 430, 400]);
    hold on; box on; grid on;
    set(gca,'LineWidth', 1.5,'fontsize',13, 'FontName', 'Times New Roman');
    plot(sort(R_LTMMSE_all(:)),linspace(0,1,K*nbrOfSetups),'k-', ...
        'LineWidth',2);
    plot(sort(SE_LTMMSE_all(:)),linspace(0,1,K*nbrOfSetups),'b--', ...
        'LineWidth',2);
    plot(sort(SE_LTMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r--', ...
        'LineWidth',2);
    plot(sort(SE_LPMMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'g--', ...
        'LineWidth',2);
    xlabel('SE per UE [bit/s/Hz]', 'FontName', 'Times New Roman');
    ylabel('CDF', 'FontName', 'Times New Roman');
    hLegend = legend({'Baseline', 'LTMMSE (All)', 'Proposed', ...
        'LP-MMSE (DCC)'}, 'Location', 'northoutside', 'FontName', ...
        'Times New Roman');
    set(hLegend, 'Orientation', 'horizontal', 'NumColumns', 2, 'Box', ...
        'off');
    % Customize the x-axis ticks to include 4 and 8
    ax = gca;
    ax.XTick = [0, 2, 4, 6, 8, 10, 12, 14];
    % Turn on the grid
    xlim([0 14]);
    ylim([-0.05 1.05]);
    
    % Print the messages in the command window
    avg_time_LPMMSE_DCC = mean(time_LPMMSE_DCC);
    avg_time_LTMMSE_all = mean(time_LTMMSE_all);
    avg_time_LTMMSE_DCC = mean(time_LTMMSE_DCC);
    time_saving_LTMMSE_DCC = ...
        (avg_time_LPMMSE_DCC - avg_time_LTMMSE_all) / avg_time_LPMMSE_DCC;
    time_saving_TMMSE_DCC_from_LPMMSE = ...
        (avg_time_LTMMSE_all - avg_time_LTMMSE_DCC) / avg_time_LTMMSE_all;
    time_saving_TMMSE_DCC_from_all = ...
        (avg_time_LPMMSE_DCC - avg_time_LTMMSE_DCC) / avg_time_LPMMSE_DCC;
    fprintf([scenarios{sc}, ...
        ' \t %.3f \t\t \t %.3f (%.2f) \t %.3f (%.2f) (%.2f) \n'], ...
        avg_time_LPMMSE_DCC, avg_time_LTMMSE_all, time_saving_LTMMSE_DCC, ...
        avg_time_LTMMSE_DCC, time_saving_TMMSE_DCC_from_LPMMSE, ...
        time_saving_TMMSE_DCC_from_all);
    
    % Save the figure as an EPS file
    saveas(gcf, ['figs_paper/', scenarios{sc}, '.eps'], 'epsc');    

end

