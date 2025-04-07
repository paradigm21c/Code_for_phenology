% 
% Start a parallel pool (if not already started)
% if isempty(gcp('nocreate'))
%     parpool('Processes',30); % Adjust the number of workers as needed
% end

product='SDC_full';
% product='MODIS';

file_path = '/gpfs/gibbs/pi/seto/jk2954/Product/Parks_polygon/ParkDS/HUMID/';

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask',...
    'Bare_Soil','Building','Freshwater_Aquatic_Vegetation','Freshwater_Wetland',...
    'Other_Paved_Surface','Other_Water','Road_Railroad','Saltwater_Aquatic_Vegetation',...
    'Tidal_Wetland','All_layer','ECM_tree'};

phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};

daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation'};

files_mask = dir(sprintf('/home/jk2954/gibbs_pi/Product/mask/parks_arc_%s_DS/*.tif',product));
mkn = length(files_mask);

% vin = {'EVI2'} ;% ,'NIRv'};
vin = {'EVI2','NDVI','NIRv'};

etn={'Shrub','Tree'};
% etn={'Tree'};

% path_mask = sprintf('D:/SetoLab/Phenology/mask/');
% files_mask = dir(sprintf('%s/Parks/*_%s.tif',path_mask,product));

if strcmp(product, 'SDC_full')
    evi2_all=importdata(sprintf('%s/evi2_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    ndvi_all=importdata(sprintf('%s/ndvi_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    nirv_all=importdata(sprintf('%s/nirv_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    
    tmax_all=importdata(sprintf('%s/tmax_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    tmin_all=importdata(sprintf('%s/tmin_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    prcp_all=importdata(sprintf('%s/prcp_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    drad_all=importdata(sprintf('%s/drad_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    % dayl_all=importdata(sprintf('%s/dayl_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
    vp_all=importdata(sprintf('%s/q2_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)

elseif strcmp(product,'MODIS')
    evi2_all=importdata(sprintf('evi2_all_arc_%s.mat',product));  %(fn,ec,mk)
    ndvi_all=importdata(sprintf('ndvi_all_arc_%s.mat',product));  %(fn,ec,mk)
    nirv_all=importdata(sprintf('nirv_all_arc_%s.mat',product));  %(fn,ec,mk)
    
    tmax_all=importdata(sprintf('tmax_all_arc_%s.mat',product));  %(fn,ec,mk)
    tmin_all=importdata(sprintf('tmin_all_arc_%s.mat',product));  %(fn,ec,mk)
    prcp_all=importdata(sprintf('prcp_all_arc_%s.mat',product));  %(fn,ec,mk)
    drad_all=importdata(sprintf('drad_all_arc_%s.mat',product));  %(fn,ec,mk)
    % dayl_all=importdata(sprintf('dayl_all_arc_%s.mat',product));  %(fn,ec,mk)
    vp_all=importdata(sprintf('%s/q2_merged_arc_%s.mat',file_path,product));  %(fn,ec,mk)
end
file_dates=importdata(sprintf('%s/file_dates_%s.mat',file_path,product));

%
years = unique(file_dates.Year);
% for vi =1:length(vin)   
%     for mk=1:mkn %length(files_mask)
%         for yr = 1:length(years)
%             for ec = 1:length(ecotype)
%                 if vi == 1
%                     vi_year = evi2_all(file_dates.Year == years(yr),ec,mk);
%                     amp_threshold = 0.10;
%                 elseif vi ==2
%                     vi_year = ndvi_all(file_dates.Year == years(yr),ec,mk);
%                     amp_threshold = 0.10;
%                 elseif vi ==3
%                     vi_year = nirv_all(file_dates.Year == years(yr),ec,mk);
%                     amp_threshold = 0.07;
%                 end
% 
%                 date_year = file_dates(file_dates.Year == years(yr));
% 
%                 tmax_year = tmax_all(file_dates.Year == years(yr),ec,mk);
%                 tmin_year = tmin_all(file_dates.Year == years(yr),ec,mk);
%                 prcp_year = prcp_all(file_dates.Year == years(yr),ec,mk);
%                 drad_year = drad_all(file_dates.Year == years(yr),ec,mk);
%                 vp_year = vp_all(file_dates.Year == years(yr),ec,mk);
% 
%                 p = prctile(vi_year,[10 90],"all");
%                 vi_year(vi_year<=p(1))=p(1);
% 
%                 vi_smooth = sgolayfilt(vi_year,2,9);
%                 vi_smooth = sgolayfilt(vi_smooth,6,9);
%                 % vi_smooth = smooth(vi_year, 12, 'sgolay', 3);
% 
%                 vi_max = max(vi_smooth,[],1);
%                 ind_max = find(vi_smooth==vi_max);
% 
%                 vi_min_sos = min(vi_smooth(1:ind_max),[],1);
%                 amp_sos = vi_max - vi_min_sos;
%                 val_15_sos = vi_min_sos+(amp_sos.*0.15);
%                 val_50_sos = vi_min_sos+(amp_sos.*0.50);
%                 val_90_sos = vi_min_sos+(amp_sos.*0.90);
% 
%                 vi_min_eos = min(vi_smooth(ind_max:end),[],1);
%                 amp_eos = vi_max - vi_min_eos;
%                 val_15_eos = vi_min_eos+(amp_eos.*0.15);
%                 val_50_eos = vi_min_eos+(amp_eos.*0.50);
%                 val_90_eos = vi_min_eos+(amp_eos.*0.90);
% 
%                 for phn =1:length(phenology_names) %{Greenup, MidGreenUp, Maturity, Peak, Senescence, MidGreenDown, Dormancy}
%                     if amp_sos >= amp_threshold
%                         if phn ==1
%                             ind_sos15 = find(vi_smooth>=val_15_sos);
%                             temp = date_year(ind_sos15(1));
%                         elseif phn==2
%                             ind_sos50 = find(vi_smooth>=val_50_sos);
%                             temp = date_year(ind_sos50(1));
%                         elseif phn==3
%                             ind_sos90 = find(vi_smooth>=val_90_sos);
%                             temp = date_year(ind_sos90(1));
%                         elseif phn==4
%                             temp =  date_year(ind_max);
%                         elseif phn==5
%                             ind_eos90 = find(vi_smooth>=val_90_eos);
%                             temp = date_year(ind_eos90(end));
%                         elseif phn==6
%                             ind_eos50 = find(vi_smooth>=val_50_eos);
%                             temp = date_year(ind_eos50(end));
%                         elseif phn==7
%                             ind_eos15 = find(vi_smooth>=val_15_eos);
%                             temp = date_year(ind_eos15(end));
%                         end
%                         vi_date(yr,ec,mk,phn) = temp;
%                         Tmax_date(yr,ec,mk,phn) = nan;
%                         Tmin_date(yr,ec,mk,phn) = nan;
%                         Prcp_date(yr,ec,mk,phn) = nan;
%                         Drad_date(yr,ec,mk,phn) = nan;
%                         % Dayl_date(yr,ec,mk,phn) = nan;
%                         Vp_date(yr,ec,mk,phn) = nan;
% 
%                     else
%                         vi_date(yr,ec,mk,phn) = NaT;
%                         Tmax_date(yr,ec,mk,phn) = nan;
%                         Tmin_date(yr,ec,mk,phn) = nan;
%                         Prcp_date(yr,ec,mk,phn) = nan;
%                         Drad_date(yr,ec,mk,phn) = nan;
%                         % Dayl_date(yr,ec,mk,phn) = nan;
%                         Vp_date(yr,ec,mk,phn) = nan;
% 
%                     end
% 
%                 end
%             end
%         end
%     end
%     vi_daymet{vi,1} = vi_date;
%     vi_daymet{vi,2} = Tmax_date;
%     vi_daymet{vi,3} = Tmin_date;
%     vi_daymet{vi,4} = Prcp_date;
%     vi_daymet{vi,5} = Vp_date;
%     % vi_daymet{vi,5} = Drad_date;
%     % vi_daymet{vi,6} = Vp_date;
% 
% end
% 
% save(sprintf('%s/vi_HUMID_%s.mat',file_path,product),'vi_daymet','-v7.3');

vi_daymet=importdata(sprintf('%s/vi_HUMID_%s.mat',file_path,product));  

parfor ps = 0:150 % Preseason days
    vi_daymet_pre = cell(length(vin), 6); 
    for vi =1:length(vin)
        vi_date=vi_daymet{vi,1};
        Tmax_date = nan(size(vi_date));
        Tmin_date = nan(size(vi_date));
        Prcp_date = nan(size(vi_date));
        Drad_date = nan(size(vi_date));
        Vp_date = nan(size(vi_date));
        % Dayl_date = nan(size(vi_date));

        for mk=1:mkn %length(files_mask)
            for yr = 1:length(years)
                for ec = 1:length(ecotype)
                    for phn =[1,4,5,7] % 1:length(phenology_names) %{Greenup, MidGreenUp, Maturity, Peak, Senescence, MidGreenDown, Dormancy}         
                        if vi_date(yr,ec,mk,phn) ~= NaT
                            temp_date = vi_date(yr,ec,mk,phn);
                            peak_date = vi_date(yr,ec,mk,4);
                            ind_daymet = find(file_dates == temp_date);

                            temp_date.Format='DDD';
                            temp_date = str2double(string(temp_date));

                            peak_date.Format='DDD';
                            peak_date = str2double(string(peak_date));

                            if phn==1 & ind_daymet>ps & temp_date>30 & temp_date < 180
                                Tmax_date(yr,ec,mk,phn) = mean(tmax_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Tmin_date(yr,ec,mk,phn) = mean(tmin_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');      
                                Prcp_date(yr,ec,mk,phn) = sum(prcp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Drad_date(yr,ec,mk,phn) = mean(drad_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Vp_date(yr,ec,mk,phn) = mean(vp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                % Dayl_date(yr,ec,mk,phn) = mean(dayl_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                            end

                            % if phn==2 & temp_date>100 & temp_date < 240
                            %     Tmax_date(yr,ec,mk,phn) = mean(tmax_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                            %     Tmin_date(yr,ec,mk,phn) = mean(tmin_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');      
                            %     Prcp_date(yr,ec,mk,phn) = sum(prcp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                            %     Drad_date(yr,ec,mk,phn) = mean(drad_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                            % end

                            if phn==4 & temp_date>150 & temp_date < 240
                                Tmax_date(yr,ec,mk,phn) = mean(tmax_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Tmin_date(yr,ec,mk,phn) = mean(tmin_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');      
                                Prcp_date(yr,ec,mk,phn) = sum(prcp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Drad_date(yr,ec,mk,phn) = mean(drad_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Vp_date(yr,ec,mk,phn) = mean(vp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                % Dayl_date(yr,ec,mk,phn) = mean(dayl_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                            end

                            if phn==5 & temp_date>180 & temp_date < 300
                                Tmax_date(yr,ec,mk,phn) = mean(tmax_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Tmin_date(yr,ec,mk,phn) = mean(tmin_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');      
                                Prcp_date(yr,ec,mk,phn) = sum(prcp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Drad_date(yr,ec,mk,phn) = mean(drad_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Vp_date(yr,ec,mk,phn) = mean(vp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                % Dayl_date(yr,ec,mk,phn) = mean(dayl_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                            end

                            if phn==7 & ind_daymet>ps & (ps< temp_date-peak_date) & temp_date>240 & temp_date < 350
                                Tmax_date(yr,ec,mk,phn) = mean(tmax_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Tmin_date(yr,ec,mk,phn) = mean(tmin_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');      
                                Prcp_date(yr,ec,mk,phn) = sum(prcp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Drad_date(yr,ec,mk,phn) = mean(drad_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                Vp_date(yr,ec,mk,phn) = mean(vp_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                                % Dayl_date(yr,ec,mk,phn) = mean(dayl_all(ind_daymet-ps:ind_daymet,ec,mk),'omitnan');
                            end
                        end    
                    end
                end
            end
        end
        vi_daymet_pre{vi,1} = vi_daymet{vi,1};
        vi_daymet_pre{vi,2} = Tmax_date;
        vi_daymet_pre{vi,3} = Tmin_date;
        vi_daymet_pre{vi,4} = Prcp_date;
        vi_daymet_pre{vi,5} = Drad_date;
        vi_daymet_pre{vi,6} = Vp_date;
        % vi_daymet_pre{vi,6} = Dayl_date;
    end
    vi_range_threshold{ps+1} = vi_daymet_pre;
    fprintf('ps %d done\n',ps+1)
end

save(sprintf('%s/vi_range_threshold_%s.mat',file_path,product),'vi_range_threshold','-v7.3');