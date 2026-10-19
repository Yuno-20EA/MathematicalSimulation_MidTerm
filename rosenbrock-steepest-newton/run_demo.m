%% RUN_DEMO  Steepest descent vs Newton + backtracking line search tren ham Rosenbrock.
%  Chay:  >> run_demo         (MATLAB hoac Octave; khong can toolbox)
%  Ket qua: bang tren console, 3 hinh, va (neu opts.saveFigures) file trong results/.
clear; clc; close all;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
opts = default_options();
outDir = fullfile(root, opts.outDir);
if ~exist(fullfile(outDir, 'images'), 'dir'), mkdir(fullfile(outDir, 'images')); end

logFile = fullfile(outDir, 'console_output.txt');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);

nStart = size(opts.starts, 2);
R = cell(nStart, 2);
for i = 1:nStart
    R{i,1} = line_search_descent(opts.starts(:,i), 'steepest', opts);
    R{i,2} = line_search_descent(opts.starts(:,i), 'newton',   opts);
    print_tables(opts.starts(:,i), R{i,1}, R{i,2}, opts);
end
diary off;

figs = {plot_surface3d(opts.starts, R, opts), ...
        plot_zoom(opts.starts, R, opts), ...
        plot_convergence(opts.starts, R, opts)};
names = {'fig1_3d_surface', 'fig2_zoom', 'fig3_convergence'};

if opts.saveFigures
    for k = 1:numel(figs)
        try
            print(figs{k}, '-dpng', '-r150', fullfile(outDir, 'images', [names{k} '.png']));
        catch err
            warning('Khong luu duoc %s: %s', names{k}, err.message);
        end
    end
end
