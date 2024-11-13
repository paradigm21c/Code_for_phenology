
% Pixel-level phenology
product = 'SDC_full';
VIs = {'EVI2', 'NDVI', 'NIRv'};
names_daymet = {'tmax','tmin','prcp','drad'};
name_drad = {'srad','dayl'};

path_datacal = 'D:/SetoLab/Phenology/data_cal/Parks_polygon/';
path_maps = 'D:/SetoLab/Phenology/map/SDC_full/';
path_park_mask = 'D:/SetoLab/Phenology/mask/parks_arc_buf_SDC_full/';
path_ecm_mask = 'D:/SetoLab/Phenology/mask/SDC_full/';
path_daymet = sprintf('C:/Temp/Temperature/DayMet/%s',product);
file_dates = importdata(sprintf('%s/file_dates_%s.mat',path_datacal,product));
file_dates.Format='yyyyMMdd';
file_name = string(file_dates);

files_park = dir([path_park_mask, '*.tif']);
files_ecm = dir([path_ecm_mask, '*.tif']);


for dy =1:length(names_daymet)
    pixels = nan(length(file_name),length(files_park), length(files_ecm));
    for fn = 1:length(file_name)
        if dy ~= 4
            [A_daymet, ~]=readgeoraster(sprintf('%s/%s/%d/%s_%s_daymet.tif',path_daymet,names_daymet{dy},file_dates(fn).Year,...
                names_daymet{dy},file_name(fn)));
            A_daymet = double(A_daymet);
        elseif dy ==4
            [A_srad, ~]=readgeoraster(sprintf('%s/%s/%d/%s_%s_daymet.tif',path_daymet,names_drad{1},file_dates(fn).Year,...
                names_drad{1},file_name(fn)));
            A_srad = double(A_srad);

            [A_dayl, ~]=readgeoraster(sprintf('%s/%s/%d/%s_%s_daymet.tif',path_daymet,names_drad{2},file_dates(fn).Year,...
                names_drad{2},file_name(fn)));
            A_dayl = double(A_dayl);
            A_daymet=((A_srad.* A_dayl)./1000000);
        end

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
            temp_maps{idx} = mean(A_daymet .* combined_mask, 'all', 'omitnan');
        end
        
        % Assign the results from temp_maps to the correct slice of maps
        for idx = 1:num_combinations
            [pk, ec] = ind2sub([length(files_park), length(files_ecm)], idx);
            pixels(fn, pk, ec) = temp_maps{idx};
        end
    
        fprintf('Daymet: %s, Year: %s DONE \n', daymet_names{dy}, file_name(fn)) 
        toc;
    end
    save(sprintf('daymet_%s_parks_pixel_ecm.mat',pixels_daymet{dy}),'pixels') %pixels_daymet{dy} = pixels;
end

