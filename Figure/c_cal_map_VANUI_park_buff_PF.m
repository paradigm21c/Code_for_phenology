if isempty(gcp('nocreate'))
    parpool('Processes',12); % Adjust the number of workers as needed
end

%%
product ='PF';

np=2;
if np ==1
    product_ntl = 'VNP46A4';
else
    product_ntl = 'VIIRS_like';
end
name_product = {'VNP46A4','LongNTL'};

path_VANUI = sprintf('C:/Temp/NTL/%s/VANUI',product_ntl);
path_mask = 'D:/SetoLab/Phenology/mask';
path_phn = 'D:/SetoLab/Phenology/map';
parks_files = [100,200,300];

name_phn= {'greenup_date', 'max_date','senescence_date','dormancy_date','max_value'};
name_vis = {'NDVI','NIRv','EVI2','EVI'}; %,'EVI2'};
years = 2018:2022;

files_park_OG = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));
files_park_OG_buff100 = dir(sprintf('%s/parks_arc_%s_DS_buffer%d/*.tif',path_mask,product,parks_files(1)));
files_park_OG_buff200 = dir(sprintf('%s/parks_arc_%s_DS_buffer%d/*.tif',path_mask,product,parks_files(2)));
files_park_OG_buff300 = dir(sprintf('%s/parks_arc_%s_DS_buffer%d/*.tif',path_mask,product,parks_files(3)));


clear temp
for pn = 1:length(files_park_OG)
    temp(pn,1) = string(files_park_OG(pn).name);
end
temp =split(temp,'_');
park_num_OG =sort(str2double(temp(:,2)));


% clear temp
% for pn = 1:length(files_park_OG_buffer)
%     temp(pn,1) = string(files_park_OG_buffer(pn).name);
% end
% temp =split(temp,'_');
% park_num_OG_buffer =sort(str2double(temp(:,2)));

for num_threshold=150 %100:100:500
    for vi =1:length(name_vis)
        for yr =1:length(years)
            [map_vanui, ~] = readgeoraster(sprintf('%s_thre%d/VANUI_%s_%d_%s.tif',path_VANUI,num_threshold,product,years(yr),lower(name_vis{vi})));
            map_vanui = double(map_vanui);
        
            parfor pn = 1:length(files_park_OG)
                [map_park_OG, R] = readgeoraster(sprintf('%s/park_%d_mask.tif',files_park_OG(pn).folder,park_num_OG(pn)));
                [map_park_buffer100, ~] = readgeoraster(sprintf('%s/park_%d_mask.tif',files_park_OG_buff100(pn).folder,park_num_OG(pn)));
                [map_park_buffer200, ~] = readgeoraster(sprintf('%s/park_%d_mask.tif',files_park_OG_buff200(pn).folder,park_num_OG(pn)));
                [map_park_buffer300, ~] = readgeoraster(sprintf('%s/park_%d_mask.tif',files_park_OG_buff300(pn).folder,park_num_OG(pn)));
         
                map_park_OG = double(map_park_OG);
                map_park_OG(map_park_OG~=0)=1;
                map_park_OG(map_park_OG==0)=nan;
        
                map_park_buffer100 = double(map_park_buffer100);
                map_park_buffer100(map_park_buffer100~=0)=1;
                map_park_buffer100(map_park_buffer100==0)=nan;
        
                map_park_buffer200 = double(map_park_buffer200);
                map_park_buffer200(map_park_buffer200~=0)=1;
                map_park_buffer200(map_park_buffer200==0)=nan;
        
                map_park_buffer300 = double(map_park_buffer300);
                map_park_buffer300(map_park_buffer300~=0)=1;
                map_park_buffer300(map_park_buffer300==0)=nan;
           
                % Calculate map_LST_park and mean value
                map_VANUI_park = map_vanui .* map_park_OG;
                VANUI_park(pn, yr) = mean(map_VANUI_park, "all", "omitmissing");
                
                % Process for 100m buffer
                temp_buff = map_vanui .* map_park_buffer100;
                VANUI_buffer100(pn, yr) = mean(temp_buff, "all", "omitmissing");
        
                temp_buff(~isnan(map_VANUI_park)) = nan;
                VANUI_diff100(pn, yr) = mean(temp_buff, "all", "omitmissing");
                
                % Process for 200m buffer
                temp_buff = map_vanui .* map_park_buffer200;
                VANUI_buffer200(pn, yr) = mean(temp_buff, "all", "omitmissing");
        
                temp_buff(~isnan(map_VANUI_park)) = nan;
                VANUI_diff200(pn, yr) = mean(temp_buff, "all", "omitmissing");
                
                % Process for 300m buffer
                temp_buff = map_vanui .* map_park_buffer300;
                VANUI_buffer300(pn, yr) = mean(temp_buff, "all", "omitmissing");
        
                temp_buff(~isnan(map_VANUI_park)) = nan;
                VANUI_diff300(pn, yr) = mean(temp_buff, "all", "omitmissing");
            end 
        
            fprintf('%d year %s, %d thresholddone\n', years(yr),name_vis{vi},num_threshold)
        end
    
    save(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_park_%s_DS_%s_%s.mat',num_threshold,product,name_vis{vi},name_product{np}),'VANUI_park')
    save(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}),'VANUI_buffer100')
    save(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}),'VANUI_buffer200')
    save(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}),'VANUI_buffer300')
    save(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}),'VANUI_diff100')
    save(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}),'VANUI_diff200')
    save(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}),'VANUI_diff300')
    
    end
end