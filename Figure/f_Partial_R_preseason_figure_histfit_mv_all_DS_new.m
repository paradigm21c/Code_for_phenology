% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
    parpool('Processes',28); % Adjust the number of workers as needed
end


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

path_home = '/home/jk2954/gibbs_pi/Product/';
path_mask = sprintf('%s/mask',path_home);
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));

path_cal = sprintf('%s/Parks_polygon/ParkDS/',path_home);

mkn = length(files_mask);

vin = {'EVI2','NDVI','NIRv'};

path_figure = sprintf('%s/Parks_polygon/Preseason_%s_DS',path_home,product);

newmap = importdata(sprintf('%s/Parks_polygon/batlow.mat',path_home));
fig_font = 'Arial';
fig_font_size = 14;


vi_range=importdata(sprintf('%s/vi_range_threshold_%s.mat',path_cal,product));
file_dates=importdata(sprintf('%s/file_dates_%s.mat',path_cal,product));

park_eco_ratio = importdata(sprintf('%s/Parks_polygon/park_eco_ratio_DS.mat',path_home));
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);


%
shrub =[2 6];
tree = [1 3 4 5];

winsize_all = [7 10 12 15];

% winsize =15;
% color_interval = floor(length(newmap)./(23-winsize));
% percent_save_neg= nan(size(1:23-winsize,2),1);
% percent_save_pos= nan(size(1:23-winsize,2),1);

% vi=3; 
etype = tree;
% phn=1; dayn=1; 

num_bin = 40;
for winsize = winsize_all
    color_interval = floor(length(newmap)./(23-winsize));
    percent_save_neg= nan(size(1:23-winsize,2),1);
    percent_save_pos= nan(size(1:23-winsize,2),1);

for dayn =[1,2,4]
    for vi =1:length(vin)
        for phn=[1 4 5 7]
            parfor ps = 0:120
                vi_daymet = vi_range{ps+1}; %(yr,ec,mk,phn)
                [rho_full, pval_full] = f_cal_partial_preseason_weight(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, 1, 23, park_veg_ratio);                        
                rho_ps_full(:,phn,ps+1)=rho_full; % (mk,phn)
                pval_ps_full(:,phn,ps+1)=pval_full; % (mk,phn)
            end

            parfor mk=1:mkn
                if phn == 1 || phn == 4 
                    [~, ind_all] = min(squeeze(rho_ps_full(mk,phn,:)),[],"all");
                elseif phn ==7 || phn == 5
                    [~, ind_all] = max(squeeze(rho_ps_full(mk,phn,:)),[],"all");
                end

                rho_max_full(mk,phn) = rho_ps_full(mk,phn,ind_all);
                pval_max_full(mk,phn) = pval_ps_full(mk,phn,ind_all);
                ps_max_full(mk,phn) = ind_all;
            end


            figure;
            hold on
            for i = [1,23-winsize]
                % Tmax / Tmin
                parfor ps = 0:120
                    % for phn=1:length(phenology_names) 
                        vi_daymet = vi_range{ps+1}; %(yr,ec,mk,phn)
                        % [rho, pval] = f_cal_partial_preseason(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize-1);
                        [rho, pval] = f_cal_partial_preseason_weight(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize-1, park_veg_ratio);
                        % [rho, pval] = f_cal_partial_preseason_common(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, 1, 23,commask_list);
                        % [rho, pval] = f_cal_partial_preseason(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, 1, 23);
                        rho_ps(:,phn,ps+1)=rho; % (mk,phn)
                        pval_ps(:,phn,ps+1)=pval; % (mk,phn)                   
                    % end
                end
            
                temp = vi_range{1}; 
                temp = temp{vi,1};
                ph_date(:,etype,:,phn) = temp(1:23,etype,:,phn);
            
                parfor mk=1:mkn
                    % for phn=7%:length(phenology_names)
                        if phn == 1 || phn == 4 
                            [~, ind] = min(squeeze(rho_ps(mk,phn,:)),[],"all");
                        elseif phn ==7 || phn == 5
                            [~, ind] = max(squeeze(rho_ps(mk,phn,:)),[],"all");   
                        end
            
                        rho_max(mk,phn) = rho_ps(mk,phn,ind);
                        pval_max(mk,phn) = pval_ps(mk,phn,ind);
                        ps_max(mk,phn) = ind;
                    % end
                end

                hf = histfit(rho_max_full(:,phn),num_bin);
                hf(1).FaceColor = newmap(1+color_interval*(i-1),:);
                hf(1).EdgeColor = 'none';
                hf(2).Color = newmap(1+color_interval*(i-1),:);
                hf(1).FaceAlpha = 0.2;

            
                percent_pos=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) > 0)./sum(~isnan(rho_max(:,phn)));
                percent_neg=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) < 0)./sum(~isnan(rho_max(:,phn)));
            
                percent_save_pos(i) = percent_pos;
                percent_save_neg(i) = percent_neg;
            
                % annotation(gcf,'textbox',...
                % [0.655357142857142 0.534238094284421 0.166785717836448 0.0638095247631981],...
                % 'String',{sprintf('%.3g (%.3g)',round(percent_pos,3),round(percent_neg,3))});
            
                ylabel('Frequency');
                xlabel(sprintf('R_{(%s - %s)}',phenology_names{phn}, daymet_names{dayn}))
            
                hy=xline(0);
                hy.LineWidth =2;
                hy.Color = 'k';
            
                ax = gca;
                ax.FontName = fig_font;
                ax.FontWeight = 'bold';
                ax.FontSize = fig_font_size;
                ax.Box = 'on';
            
                xlim([-1 1])
                ylim([0 35])
            
                % rho_stack{dayn,vi,phn,i} = rho_max;
                % pval_stack{dayn,vi,phn,i} = pval_max;
                % ps_stack{dayn,vi,phn,i} = ps_max;
            end
            
            hf = histfit(rho_max_full(:,phn));
            hf(1).FaceColor = 'none'; %newmap(color_interval(w_s),:);
            hf(1).EdgeColor = 'none';
            % hf(1).FaceColor = newmap(1,:);
            hf(2).Color = newmap(end/2,:) ;%[1 0 0];

            % rho_stack_full{dayn,vi,phn} = rho_max_full;
            % pval_stack_full{dayn,vi,phn} = pval_max_full;
            % ps_stack_full{dayn,vi,phn} = ps_max_full;
        
            savefig(gcf,sprintf('%s/Histfit_X_%s_%s_tree_Fusion_%s_mv%d.fig',path_figure, phenology_names{phn},daymet_names{dayn},vin{vi},winsize))
            saveas(gcf,sprintf('%s/Histfit_X_%s_%s_tree_Fusion_%s_mv%d.png',path_figure, phenology_names{phn},daymet_names{dayn},vin{vi},winsize))

            close all
        end
    end
end

% save(sprintf('%s/Rho_%s_DS_mv%d.mat',path_figure,product,winsize),'rho_stack');
% save(sprintf('%s/Pval_%s_DS_mv%d.mat',path_figure,product,winsize),'pval_stack');
% save(sprintf('%s/PS_%s_DS_mv%d.mat',path_figure,product,winsize),'ps_stack');
end

% save(sprintf('%s/Rho_%s_DS_mv_full.mat',path_figure,product),'rho_stack_full');
% save(sprintf('%s/Pval_%s_DS_mv_full.mat',path_figure,product),'pval_stack_full');
% save(sprintf('%s/PS_%s_mv_DS_full.mat',path_figure,product),'ps_stack_full');
