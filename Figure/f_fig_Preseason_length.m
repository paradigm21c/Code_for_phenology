% % Start a parallel pool (if not already started)
% if isempty(gcp('nocreate'))
%     parpool('Processes',4); % Adjust the number of workers as needed
% end

product='SDC_full';
phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask'};

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask',...
    'Bare_Soil','Building','Freshwater_Aquatic_Vegetation','Freshwater_Wetland',...
    'Other_Paved_Surface','Other_Water','Road_Railroad','Saltwater_Aquatic_Vegetation',...
    'Tidal_Wetland','All_layer','ECM_tree'};

newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;

path_home = 'D:/SetoLab/Phenology';
path_mask = sprintf('%s/mask',path_home);
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));
mkn = length(files_mask);

vin = {'EVI2','NDVI','NIRv'};

path_cal = sprintf('%s/data_cal/Parks_polygon/ParkDS/',path_home);
path_figure = sprintf('%s/figure/Preseason_%s_DS/all_years/',path_home,product);
% vi_range=importdata(sprintf('%s/vi_range_threshold_%s.mat',path_cal,product));

ps_stack_full = importdata(sprintf('%s/PS_%s_DS_mv_full.mat',path_cal,product));

shrub =[2 6];
tree = [1 3 4 5];
dayn=1;
winsize_all = [7 10 12 15];

for winsize = 10 %winsize_all
    ps_stack = importdata(sprintf('%s/PS_%s_DS_mv%d.mat',path_cal,product,winsize));
    
    for vi =[3]%1 2] %1:length(vin)
        for dayn = [1,2] %,4]
            for phn=[1] % 4 5 7]

                ps_max = ps_stack{dayn,vi,phn,1} ;
                A = ps_max(:,phn);
                An = repmat("2000-2011",length(A),1);

                ps_max = ps_stack{dayn,vi,phn,23-winsize} ;
                B = ps_max(:,phn);
                Bn = repmat("2011-2022",length(B),1);

                ps_max_full = ps_stack_full{dayn,vi,phn} ;
                C = ps_max_full(:,phn);
                Cn = repmat("2000-2022",length(C),1);

                table_ps = table([A;B;C],[An;Bn;Cn],'VariableNames',{'PS','MV'});
                table_ps.MV = categorical(table_ps.MV ,{'2000-2011','2011-2022','2000-2022'});

                figure('Position',[600,200,500,300]);
                hold on

                h = boxchart(table_ps.PS,'GroupByColor',table_ps.MV);
                h(1).BoxFaceColor = newmap(1,:);
                h(1).MarkerColor = newmap(1,:);
                h(1).MarkerStyle = '.';
                h(1).BoxWidth = 0.75;

                h(2).BoxFaceColor = newmap(end,:);
                h(2).MarkerColor = newmap(end,:);
                h(2).MarkerStyle = '.';
                h(2).BoxWidth = 0.75;

                h(3).BoxFaceColor = newmap(end/2,:);
                h(3).MarkerColor = newmap(end/2,:);
                h(3).MarkerStyle = '.';
                h(3).BoxWidth = 0.75;

                ax = gca;
                ax.FontName = fig_font;
                % ax.FontWeight = 'bold';
                ax.FontSize = fig_font_size-2;
                ax.Box = 'on';
                ax.XTickLabel = 'Data period';
                legend;
                ylabel('Preseason length (Day)')
    
                savefig(gcf,sprintf('%s/Boxplot_Win10_%s_X_Year_Y_PS_%s_%s.fig',path_figure,daymet_names{dayn}, product,vin{vi}))
                saveas(gcf,sprintf('%s/Boxplot_Win10_%s_X_Year_Y_PS_%s_%s.png',path_figure,daymet_names{dayn},product,vin{vi}))

            end
        end
    end
end
