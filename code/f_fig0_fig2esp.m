% Set the folder path containing the .fig files
path_figure = 'D:\SetoLab\Phenology\figure\Revision'; % <-- change this to your actual path

% Get list of all .fig files in the folder
fig_files = dir(fullfile(path_figure, '*.fig'));

% Loop through each .fig file
for k =13:14% 1:length(fig_files)
    fig_name = fig_files(k).name;
    
    % Open the .fig file
    fig_path = fullfile(path_figure, fig_name);
    openfig(fig_path, 'invisible');  % 'invisible' to avoid UI pop-ups

    ax = gca;                       % Get current axes
    if k~=13 && k~= 14
        ax.YAxisLocation = 'right';    % Move Y-axis to right
    end
    yline(ax.YLim, 'Color', 'k', 'LineWidth', 0.5);
    xline(ax.XLim(1), 'Color', 'k', 'LineWidth', 0.5);
    ax.Box = 'off';                % Turn off top and right axis lines
    ax.TickDir = 'out';            % Ticks pointing outward   

    % Get base name without extension
    [~, base_name, ~] = fileparts(fig_name);
    
    % Construct new filename
    eps_filename = sprintf('%s/ESP/%s.eps', path_figure, base_name);
    
    % Export current figure to .eps
    exportgraphics(gcf, eps_filename, 'BackgroundColor', 'none', 'ContentType', 'vector');
    
    % Close the figure
    close(gcf);
    
    fprintf('Saved: %s\n', eps_filename);
end
