
% % Start a parallel pool (if not already started)
% if isempty(gcp('nocreate'))
%    parpool('Processes',30); % Adjust the number of workers as needed
% end

% % Pixel-level pheonolgy
% product ='SDC_full';
% VIs = {'EVI2','NDVI','NIRv'};
% maps_names = {'max_value', 'min_value', 'max_date', 'min_date', 'val_15_value', 'val_50_value', 'val_90_value'...
%     'greenup_date','dormancy_date','midgreenup_date','midgreendown_date','maturity_date','senescence_date',...
%     'val_sos','sos_date','val_eos','eos_date'};
% 
% path_maps=  'D:\SetoLab\Phenology\map\SDC_full\';
% path_park_mask = 'D:\SetoLab\Phenology\mask\parks_arc_buf_SDC_full\';
% path_ecm_mask = 'D:\SetoLab\Phenology\mask\SDC_full\';
% 
% files_park = dir([path_park_mask,'*.tif']);
% files_ecm = dir([path_ecm_mask,'*.tif']);
% 
% 
% years = 2000:2022;
% 
% maps_vi = cell(1, length(VIs));
% 
% for vi = 1:length(VIs)
%     maps = nan(length(files_park),length(years),length(files_ecm),length(maps_names));
%     for mp = 1:length(maps_names)
%         for yr = 1:length(years)
%             A_map = readgeoraster(sprintf('%s/%s_%s_%s_%d.tif',path_maps,VIs{vi},product,maps_names{mp},years(yr)));
% 
%             parfor pk =1:length(files_park)
%                 A_park = readgeoraster(sprintf('%s/%s',files_park(pk).folder,files_park(pk).name));
%                 A_park(A_park==0)=nan; A_park(~isnan(A_park))=1;
% 
%                 for ec = 1:length(files_ecm)
%                     A_ecm = readgeoraster(sprintf('%s/%s',files_ecm(ec).folder,files_ecm(ec).name));
%                     maps(pk,yr,ec,mp) = mean(A_map.*A_park,'all','omitnan');
%                 end  
% 
%             end
%         end
%     end
%     maps_vi{vi} = maps;
% end



%% Par

% Pixel-level phenology
product = 'SDC_full';
VIs = {'EVI2', 'NDVI', 'NIRv'};
maps_names = {'max_value', 'min_value', 'max_date', 'min_date', 'val_15_value', 'val_50_value', 'val_90_value', ...
    'greenup_date', 'dormancy_date', 'midgreenup_date', 'midgreendown_date', 'maturity_date', 'senescence_date', ...
    'val_sos', 'sos_date', 'val_eos', 'eos_date'};

path_maps = 'D:\SetoLab\Phenology\map\SDC_full\';
path_park_mask = 'D:\SetoLab\Phenology\mask\parks_arc_buf_SDC_full\';
path_ecm_mask = 'D:\SetoLab\Phenology\mask\SDC_full\';

files_park = dir([path_park_mask, '*.tif']);
files_ecm = dir([path_ecm_mask, '*.tif']);

years = 2000:2022;

maps_vi = cell(1, length(VIs));

for vi = 1:length(VIs)
    maps = nan(length(files_park), length(years), length(files_ecm), length(maps_names));
    for mp = 1:length(maps_names)
        for yr = 1:length(years)
            A_map = readgeoraster(sprintf('%s/%s_%s_%s_%d.tif', path_maps, VIs{vi}, product, maps_names{mp}, years(yr)));
            
            % Create a composite index for the parfor loop
            num_combinations = length(files_park) * length(files_ecm);
            temp_maps = cell(num_combinations, 1);

            parfor idx = 1:num_combinations
                [pk, ec] = ind2sub([length(files_park), length(files_ecm)], idx);

                A_park = readgeoraster(sprintf('%s/%s', files_park(pk).folder, files_park(pk).name));
                A_park(A_park == 0) = nan;
                A_park(~isnan(A_park)) = 1;

                A_ecm = readgeoraster(sprintf('%s/%s', files_ecm(ec).folder, files_ecm(ec).name));
                A_ecm(A_ecm == 0) = nan;
                A_ecm(~isnan(A_ecm)) = 1;

                % Apply the park and ecm masks and calculate the mean
                combined_mask = A_park .* A_ecm;
                temp_maps{idx} = mean(A_map .* combined_mask, 'all', 'omitnan');
            end
            
            % Assign the results from temp_maps to the correct slice of maps
            for idx = 1:num_combinations
                [pk, ec] = ind2sub([length(files_park), length(files_ecm)], idx);
                maps(pk, yr, ec, mp) = temp_maps{idx};
            end
        
            fprintf('VI: %s, Year: %d DONE \n', VIs{vi}, years(yr))
        end
    end
    maps_vi{vi} = maps;
end

save('maps_parks_vi_pixel_level.mat','maps_vi')
