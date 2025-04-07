product='SDC_full';
% product='MODIS';

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask'};

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask',...
    'Bare_Soil','Building','Freshwater_Aquatic_Vegetation','Freshwater_Wetland',...
    'Other_Paved_Surface','Other_Water','Road_Railroad','Saltwater_Aquatic_Vegetation',...
    'Tidal_Wetland','All_layer'};

phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};

parks_names = {'Botanical_garden_and_zoo' , 'Central_Park', 'Dyker_beach_park',...
    'Forest_Park', 'Highland_Park', 'Kissena_Park', 'Prospect_Park',...
    'Van_corlandt_Park', 'combined_mask'};

daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};

path_mask = sprintf('D:/SetoLab/Phenology/mask/');
% files_mask = dir(sprintf('%s/parks_arc_buf_%s/*.tif',path_mask,product));
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));

% park_area = importdata(sprintf('park_area_%s.mat',product));
% park_name = importdata(sprintf('park_name_%s.mat',product));

park_name_SDC = importdata(sprintf('park_name_SDC.mat'));
park_name_MODIS = importdata(sprintf('park_name_MODIS.mat'));

path_cal = 'D:/SetoLab/Phenology/data_cal/Parks_polygon/ParkDS/HUMID';
% 
% for pn = 1:length(park_name_MODIS)
%     commask_list(pn) = find(park_name_SDC==park_name_MODIS(pn));
% end
% 
% small_list = 1:length(common_list);
% 
% for cm = 1:length(commask_list)
%     small_list(small_list == commask_list(cm))=[]; 
% end
% mkn = length(commask_list);
mkn = length(files_mask);

vin = {'EVI2','NDVI','NIRv'};
etn={'Tree'};

path_figure = sprintf('D:/SetoLab/Phenology/figure/Preseason_%s_DS/all_years/HUMID',product);
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
fig_y_axis_loc = 'right';
fig_font_weight = "normal"; % "bold"

park_eco_ratio = importdata("park_eco_ratio_DS.mat");
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);

% vi_range=importdata('vi_range_30_180_SDC.mat');
vi_range=importdata(sprintf('%s/vi_range_threshold_%s.mat',path_cal,product));
file_dates=importdata(sprintf('%s/file_dates_%s.mat',path_cal,product));

%
shrub =[2 6];
% tree = [4 5];
tree = [1 3 4 5];

vi=1; etype = tree; phn=1; dayn=1; 

winsize = [8 10 12];
mark = {'v','^','<'};

color_interval = [1 256/2 256];

for dayn =[1] %,2,4]
    for vi =1:length(vin)
        for phn=1 %[1 4 5 7] %:length(phenology_names)
            figure('Position',[600,200,500,200]);
            hold on
            for mv =1:length(winsize)
                pos_stack=importdata(sprintf('%s/Neg_%s_DS_mv%d_2018.mat',path_cal ,product,winsize(mv)));
                percent_save_pos = pos_stack{dayn,vi,phn,end};       

                X=2000+floor(winsize(mv)/2):2018-ceil((winsize(mv)/2));
                Y=percent_save_pos*100;
         
                h=plot(X,Y);
                h.Marker = mark{mv};
                h.MarkerFaceColor = newmap(color_interval(mv),:);
                h.MarkerEdgeColor = 'none';
                h.LineStyle = "none";

                mln= fitlm(X,Y);
                

                significance_value_tau = 0.05;
                significance_value_ac = 0.05;
                gpu_shift_critical_size = 550;
                [tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(X, Y, significance_value_tau, significance_value_ac, gpu_shift_critical_size);
                
                % Estimate Theil-Sen parameters.
                % plims = [min(X), max(X)]';
                [est_ts, ~] = TheilSen(X',Y);
                Yts = X.*est_ts(2) + est_ts(1);

                % h=plot(X,mln.Fitted);
                h=plot(X,Yts);
                h.LineWidth=2;
                h.Color = newmap(color_interval(mv),:);

                
                ax = gca;
                ax.FontName = fig_font;
                ax.FontWeight = fig_font_weight;
                ax.FontSize = fig_font_size;
                ax.Box = 'on';
                ax.YAxisLocation=fig_y_axis_loc;
    
                clear rho_max pval_max
            end
            
            l=legend('mv=8','','mv=10','','mv=12','');
            l.Location='northeast';
            % l.Location='southwest';
            
            xlim([2000 2018])
            % ylim([0 5])
            hold off
            
            savefig(gcf,sprintf('%s/Plot_X_year_Y_%s_%s_tree_SigPer_neg_Fusion_%s.fig',path_figure, phenology_names{phn},daymet_names{dayn},vin{vi}))
            saveas(gcf,sprintf('%s/Plot_X_year_Y_%s_%s_tree_SigPer_neg_Fusion_%s.png',path_figure, phenology_names{phn},daymet_names{dayn},vin{vi}))
            fprintf('Plot_X_year_Y_%s_%s_tree_SigPer_neg_Fusion_%s.png\n', phenology_names{phn},daymet_names{dayn},vin{vi})
        
            % close all
        end
    end
end