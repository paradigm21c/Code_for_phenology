% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
    parpool('Processes',28); % Adjust the number of workers as needed
end


product='SDC_full';
phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask',...
    'Bare_Soil','Building','Freshwater_Aquatic_Vegetation','Freshwater_Wetland',...
    'Other_Paved_Surface','Other_Water','Road_Railroad','Saltwater_Aquatic_Vegetation',...
    'Tidal_Wetland','All_layer','ECM_tree'};

path_home = '/home/jk2954/gibbs_pi/Product/';
path_mask = sprintf('%s/mask',path_home);
files_mask = dir(sprintf('%s/parks_arc_%s/*.tif',path_mask,product));
path_cal = sprintf('%s/Parks_polygon/ParkOG/',path_home);

mkn = length(files_mask);

vin = {'EVI2','NDVI','NIRv'};


path_figure = sprintf('%s/Parks_polygon/Preseason_%s',path_home,product);

newmap = importdata(sprintf('%s/Parks_polygon/batlow.mat',path_home));
fig_font = 'Arial';
fig_font_size = 14;

vi_range=importdata(sprintf('%s/vi_range_threshold_%s.mat',path_cal,product));
file_dates=importdata(sprintf('%s/file_dates_%s.mat',path_cal,product));

park_eco_ratio = importdata(sprintf('%s/Parks_polygon/park_eco_ratio.mat',path_home));
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);

%
shrub =[2 6];
% tree = [4 5];
tree = [1 3 4 5];

vi=1; etype = tree; phn=1; dayn=1; 

winsize = [7 10 12 15];
mark = {'v','^','<','>'};

% winsize = [5 7 10 12 15];
% mark = {'o','v','^','<','>'};

% winsize =15;
color_interval = floor(length(newmap)./(length(winsize)));

for dayn =[1,2,4]
    for vi =1:length(vin)
        for phn=[1 4 5 7] %:length(phenology_names)
            figure;
            hold on
            for mv =1:length(winsize)
            
                percent_save_neg= nan(size(1:23-winsize(mv)+1,2),1);
                percent_save_pos= nan(size(1:23-winsize(mv)+1,2),1);
                for i = 1:23-winsize(mv)+1
                    % Tmax / Tmin
                    parfor ps = 0:120
                        % for phn=1:length(phenology_names) 
                            vi_daymet = vi_range{ps+1}; %(yr,ec,mk,phn)
                            % [rho, pval] = f_cal_partial_preseason(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize(mv)-1);
                            [rho, pval] = f_cal_partial_preseason_weight(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize(mv)-1, park_veg_ratio);
                            % [rho, pval] = f_cal_partial_preseason_common(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, i, i+winsize(mv)-1,commask_list);
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
                
                    percent_pos=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) > 0)./sum(~isnan(rho_max(:,phn)));
                    percent_neg=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) < 0)./sum(~isnan(rho_max(:,phn)));
                
                    percent_save_pos(i) = percent_pos;
                    percent_save_neg(i) = percent_neg;
                
                    rho_stack(:,:,i) = rho_max(:,phn);
                end
                
                
                h=plot(2000+floor(winsize(mv)/2):2022-(winsize(mv)/2)+1,percent_save_neg*100);
                h.Marker = mark{mv};
                h.MarkerFaceColor = newmap(1+color_interval*(mv-1),:);
                h.MarkerEdgeColor = 'none';
                h.LineStyle = "none";
                mln= fitlm(2000+floor(winsize(mv)/2):2022-(winsize(mv)/2)+1,percent_save_neg*100);
                h=plot(2000+floor(winsize(mv)/2):2022-(winsize(mv)/2)+1,mln.Fitted);
                h.LineWidth=2;
                h.Color = newmap(1+color_interval*(mv-1),:);
                
                ax = gca;
                ax.FontName = fig_font;
                ax.FontWeight = 'bold';
                ax.FontSize = fig_font_size;
                ax.Box = 'on';
                
                
                clear rho_max pval_max
            end
            
            l=legend('mv=7','','mv=10','','mv=12','','mv=15','');
            l.Location='northeast';
            % l.Location='southwest';
            
            xlim([2000 2022])
            ylim([0 50])
            hold off
            
            savefig(gcf,sprintf('%s/Plot_X_year_Y_%s_%s_tree_SigPer_neg_Fusion_%s.fig',path_figure ,phenology_names{phn},daymet_names{dayn},vin{vi}))
            saveas(gcf,sprintf('%s/Plot_X_year_Y_%s_%s_tree_SigPer_neg_Fusion_%s.png',path_figure ,phenology_names{phn},daymet_names{dayn},vin{vi}))
            fprintf('Plot_X_year_Y_%s_%s_tree_SigPer_pos_Fusion_%s.png \n', phenology_names{phn},daymet_names{dayn},vin{vi})
            
            close all
        end
    end
end