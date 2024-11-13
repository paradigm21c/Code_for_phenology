if isempty(gcp('nocreate'))
    parpool('Processes',12); % Adjust the number of workers as needed
end

%%
product ='PF';

path_LST = sprintf('D:/SetoLab/Landsat_ARD/Annual_%s_Mean',product);
path_mask = 'D:/SetoLab/Phenology/mask';
path_phn = 'D:/SetoLab/Phenology/map';
parks_files = [100,200,300];

name_phn= {'greenup_date', 'max_date','senescence_date','dormancy_date','max_value'};
name_vis = {'EVI2','NIRv'};
years = 2000:2022;

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

vi=1;
% for vi=1:length(name_vis)
    for yr =1:length(years)
        [map_LST, ~] = readgeoraster(sprintf('%s/mean_%d.tif',path_LST,years(yr)));
        map_LST = double(map_LST);
        map_LST = map_LST - 273.15;
        map_LST = kron(map_LST,ones(10));
        
        for phn = 1%:length(name_phn)
            %      [map_phn, ~] = readgeoraster(sprintf('%s/%s/%s_%s_%s_%d.tif',path_phn,product,name_vis{vi},product,name_phn{phn},years(yr))); % EVI2_SDC_full_greenup_date_2018
            % map_phn = double(map_phn);
            % 
            % if phn==1
            %     map_phn(map_phn<30 | map_phn>180)=nan;
            % elseif phn==2
            %     map_phn(map_phn<150 | map_phn>240)=nan;
            % elseif phn==3
            %     map_phn(map_phn<180 | map_phn>300)=nan;
            % elseif phn==4
            %     map_phn(map_phn<240 | map_phn>350)=nan;
            % elseif phn==5
            %     map_phn(map_phn<0 | map_phn>1)=nan;
            % end

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

                % map_phn_park = map_phn.*map_park_OG;
                % phn_park(pn,yr,phn) = mean(map_phn_park,"all","omitmissing");

                
                % map_LST_park = map_LST.*map_park_OG;
                % 
                % map_LST_buffer100 = map_LST.*map_park_buffer100;
                % map_LST_diff100 = map_LST_buffer100;
                % map_LST_diff100(~isnan(map_LST_park))=nan;
                % 
                % map_LST_buffer200 = map_LST.*map_park_buffer200;
                % map_LST_diff200 = map_LST_buffer200;
                % map_LST_diff200(~isnan(map_LST_park))=nan;
                % 
                % map_LST_buffer300 = map_LST.*map_park_buffer300;
                % map_LST_diff300 = map_LST_buffer300;
                % map_LST_diff300(~isnan(map_LST_park))=nan;

                % phn_park(pn,yr,phn) = mean(map_phn_park,"all","omitmissing");
                % LST_park(pn,yr,phn) = mean(map_LST_park,"all","omitmissing");
                % 
                % LST_buffer100(pn,yr,phn) = mean(map_LST_buffer100,"all","omitmissing");
                % LST_diff100(pn,yr,phn) = mean(map_LST_diff100,"all","omitmissing");
                % 
                % LST_buffer200(pn,yr,phn) = mean(map_LST_buffer200,"all","omitmissing");
                % LST_diff200(pn,yr,phn) = mean(map_LST_diff200,"all","omitmissing");
                % 
                % LST_buffer300(pn,yr,phn) = mean(map_LST_buffer300,"all","omitmissing");
                % LST_diff300(pn,yr,phn) = mean(map_LST_diff300,"all","omitmissing");


                % Calculate map_LST_park and mean value
                map_LST_park = map_LST .* map_park_OG;
                LST_park(pn, yr) = mean(map_LST_park, "all", "omitmissing");
                
                % Process for 100m buffer
                temp_buff = map_LST .* map_park_buffer100;
                LST_buffer100(pn, yr) = mean(temp_buff, "all", "omitmissing");

                temp_buff(~isnan(map_LST_park)) = nan;
                LST_diff100(pn, yr) = mean(temp_buff, "all", "omitmissing");
                
                % Process for 200m buffer
                temp_buff = map_LST .* map_park_buffer200;
                LST_buffer200(pn, yr) = mean(temp_buff, "all", "omitmissing");

                temp_buff(~isnan(map_LST_park)) = nan;
                LST_diff200(pn, yr) = mean(temp_buff, "all", "omitmissing");
                
                % Process for 300m buffer
                temp_buff = map_LST .* map_park_buffer300;
                LST_buffer300(pn, yr) = mean(temp_buff, "all", "omitmissing");

                temp_buff(~isnan(map_LST_park)) = nan;
                LST_diff300(pn, yr) = mean(temp_buff, "all", "omitmissing");

            end 
        end
        fprintf('%d done\n', years(yr))
    end
    % vi_all{vi,1} = phn_park;
% end

% save(sprintf('D:/SetoLab/Phenology/data_cal/LST/vi_all_DS.mat'),'vi_all')
save(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_park_%s_DS.mat',product),'LST_park')
save(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(1)),'LST_buffer100')
save(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(2)),'LST_buffer200')
save(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(3)),'LST_buffer300')
save(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(1)),'LST_diff100')
save(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(2)),'LST_diff200')
save(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(3)),'LST_diff300')
