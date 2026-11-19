clear; clc; close all;

%% 
data_dir = fullfile(pwd, 'data');
out_dir  = fullfile(pwd, 'figures_final');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultAxesLineWidth', 1.0);
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultLegendFontSize', 9);
style = algo_style();   
%% ===== Load data=====
mat_path = fullfile(data_dir, 'engineering_data.mat');
if exist(mat_path, 'file')
    S = load(mat_path);
    fprintf('Load successfully, variable list：\n');
    disp(fieldnames(S));
else
    error('Cannot find %s. Please check the .mat file name in the data directory', mat_path);
end


%% ==========
key_probs = [1, 4, 5, 6];  
figure('Position', [100 100 900 600]);

for k = 1:length(key_probs)
    subplot(2, 2, k);
    hold on;
    p_idx = key_probs(k);
    for a = 1:length(style)
        
    end
    title(prob_names{p_idx}, 'FontSize', 11);
    xlabel('Evaluation frequency（FEs）');
    ylabel('Optimal objective value');
    set(gca, 'YScale', 'log'); 
    grid on; box on;
end

subplot(2, 2, 4);
legend({style.name}, 'Location', 'best', 'NumColumns', 2);

save_final(gcf, fullfile(out_dir, 'eng_conv_key'));

%% ==========
figure('Position', [100 100 1000 700]);
for p = 1:length(prob_names)
    subplot(2, 3, p);
    hold on;
    title(prob_names{p}, 'FontSize', 11);
    ylabel('objective function value');
    set(gca, 'XTickLabelRotation', 45);
    grid on; box on;
end
save_final(gcf, fullfile(out_dir, 'eng_box_all'));

%% ==========
figure('Position', [100 100 900 400]);
for k = 1:length(key_probs)
    subplot(1, length(key_probs), k);
    hold on;
    p_idx = key_probs(k);
    % boxplot(squeeze(box_data(p_idx, :, :)), ...);
    title(prob_names{p_idx}, 'FontSize', 11);
    ylabel('objective function value');
    set(gca, 'XTickLabelRotation', 45);
    grid on; box on;
end
save_final(gcf, fullfile(out_dir, 'eng_box_key'));

%% ==========
figure('Position', [100 100 700 400]);
hold on;
b = bar(rank_data, 'FaceColor', 'flat');
for a = 1:length(style)
    b.CData(a, :) = style(a).color;
end
set(gca, 'XTickLabel', {style.name}, 'XTickLabelRotation', 45);
ylabel('Average Friedman rank');
ylim([0 max(rank_data)*1.15]);
grid on; box on;
for a = 1:length(style)
    text(a, rank_data(a)+0.15, sprintf('%.2f', rank_data(a)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end
save_final(gcf, fullfile(out_dir, 'eng_rank_bar'));

%% ==========
for p = 1:length(prob_names)
    figure('Position', [100 100 700 450]);
    hold on;
    for a = 1:length(style)
        % y = squeeze(conv_data(p, a, :));
        % plot(...);
    end
    title(prob_names{p}, 'FontSize', 12);
    xlabel('Evaluation frequency（FEs）');
    ylabel('Optimal objective value');
    set(gca, 'YScale', 'log');
    legend({style.name}, 'Location', 'best');
    grid on; box on;
    save_final(gcf, fullfile(out_dir, sprintf('eng_conv_all_%s', prob_names{p})));
    close(gcf);
end

disp(' ');
disp('===All images have been generated ===');
disp(['Output directory：', out_dir]);
disp('.fig The file can be opened using openfig() and then the legend, coordinates, and colors can be directly edited');
disp('The .pdf file is a 300dpi vector image that can be directly inserted into the paper');

%% =========
function save_final(fig_handle, base_path)
    savefig(fig_handle, [base_path '.fig'])
    print(fig_handle, [base_path '.pdf'], '-dpdf', '-r300', '-bestfit');
    fprintf('Save：%s.fig and %s.pdf\n', basename(base_path), basename(base_path));
end

function n = basename(p)
    [~, n, ~] = fileparts(p);
end