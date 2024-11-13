% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
    parpool('Processes',30); % Adjust the number of workers as needed
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

path_home = '/home/jk2954/gibbs_pi/Product';
path_mask = sprintf('%s/mask',path_home);
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));
mkn = length(files_mask);

vin = {'EVI2','NDVI','NIRv'};

path_cal = sprintf('%s/Parks_polygon/',path_home);
% path_figure = sprintf('%s/figure/Preseason_%s/all_years/',path_home,product);
vi_range=importdata(sprintf('%s/ParkDS/vi_range_threshold_%s.mat',path_cal,product));

ps_stack_full = importdata(sprintf('%s/Preseason_%s_DS/PS_%s_DS_mv_full.mat',path_cal,product,product));

park_eco_ratio = importdata(sprintf('%s/park_eco_ratio_DS.mat',path_cal));
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);


%
shrub =[2 6];
tree = [1 3 4 5];

winsize_all = [10 7 12 15];

% vi=3; 
etype = tree;
% phn=1; dayn=1; 
for vi =[3 ] %1:length(vin) 1 2
    for winsize = [10 12 15 7]%winsize_all 7 12 
        ps_stack = importdata(sprintf('%s/Preseason_%s_DS/PS_%s_DS_mv%d.mat',path_cal,product,product,winsize));
        percent_save_neg= nan(size(1:23-winsize,2),1);
        percent_save_pos= nan(size(1:23-winsize,2),1);
        for dayn = [1,2,4]
            for phn=[1 4 5 7]
                ps_max_full = ps_stack_full{dayn,vi,phn} ;
                pr_table_full = [];

                if winsize ==10
                    parfor mk=1:mkn
                    % for mk=1:mkn
                        ind_all=ps_max_full(mk,phn);
                        vi_daymet = vi_range{ind_all}; %(yr,ec,mk,phn)
                        [pr_table_temp] = f_cal_partial_preseason_fitlm(vi_daymet, daymet_names, vi, phn, dayn, mk, etype, 1, 23, park_veg_ratio);                        
                        pr_table_full= [pr_table_full;pr_table_temp];
                    end
                fitlm_stack_full = pr_table_full;
                fprintf('Winsize full: %s daymet %s %s Done \n',...
                        daymet_names{dayn},vin{vi},phenology_names{phn})
                end
    
                for i = 1:23-winsize %[1,23-winsize]      
                    ps_max = ps_stack{dayn,vi,phn,i} ;
                    pr_table = [];
    
                    parfor mk=1:mkn
                    % for mk=1:mkn
                        ind_all=ps_max(mk,phn);
                        vi_daymet = vi_range{ind_all}; %(yr,ec,mk,phn)
                        [pr_table_temp] = f_cal_partial_preseason_fitlm(vi_daymet, daymet_names, vi, phn, dayn, mk, etype,  i, i+winsize-1, park_veg_ratio);                        
                        pr_table= [pr_table;pr_table_temp];
                    end
    
                    fitlm_stack{i} = pr_table;
                    fprintf('Winsize %d (%d th): %s daymet %s %s Done \n',...
                        winsize,i,daymet_names{dayn},vin{vi},phenology_names{phn})
                end
               save(sprintf('%s/ParkDS/Fitlm_%s_%s_%s_%s_DS_mv%d.mat',...
                   path_cal,vin{vi},daymet_names{dayn},phenology_names{phn},product,winsize), ...
                   'fitlm_stack','-v7.3');
               if winsize ==10
                   save(sprintf('%s/ParkDS/Fitlm_%s_%s_%s_%s_DS_mv_full.mat', ...
                       path_cal,vin{vi},daymet_names{dayn},phenology_names{phn},product), ...
                       'fitlm_stack_full','-v7.3');
               end  
            end
        end       
    end   
end
    