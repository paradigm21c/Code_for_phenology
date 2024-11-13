%
product = 'SDC_full';

newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
num_bin = 40;

years =2000:2022;
parks_files = [100,200,300];
% Path
path_all = 'D:/SetoLab/Phenology/';
path_map = sprintf('%s/map/%s/',path_all,product);
path_buff = sprintf('%s/mask/%s_buffer/',path_all,product);
path_park = sprintf('%s/mask/parks_arc_%s_DS/',path_all,product);
path_ecm = sprintf('%s/mask/%s_temp/',path_all,product);
path_figure = sprintf('%s/figure/Preseason_%s_DS/all_years',path_all,product);

% Files
file_park = dir(sprintf('%s/*.tif',path_park));
file_buff = dir(sprintf('%s/*.tif',path_buff));
file_ecm =  dir(sprintf('%s/tree_mask_%s.tif',path_ecm,product));

park_num= importdata(sprintf("park_num_%s_DS.mat",product));
ind_park_all=1:length(park_num);
ind_park_vs = find(park_num<=11); % 0.009 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>11 & park_num<=50); % 0.99 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>50 & park_num<=200); % 4.50 to 18 hectares
ind_park_l = find(park_num>200); % 18.09 to more hectares

park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};
phenology_names = {'SOS', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
name_vis = {'EVI2','NIRv'};

vi =2;
% Park names
park_name = nan(length(file_park),1);
for i=1:length(file_park)
    temp = split(file_park(i).name,'_');
    park_name(i) = str2double(temp{2});
end
park_name = sort(park_name);

% ECM
A_ecm = readgeoraster(sprintf('%s/%s',file_ecm.folder,file_ecm.name));
A_ecm = double(A_ecm);
A_ecm(A_ecm==0)=nan;
A_ecm = gpuArray(A_ecm);


% Buff
for bf =1:length(file_buff)
    A_temp = readgeoraster(sprintf('%s/%s',file_buff(bf).folder,file_buff(bf).name));
    A_buff(:,:,bf)=double(A_temp);
end
% A_buff = gpuArray(A_buff);

for yr =1:length(years)
    A_phn = readgeoraster(sprintf('%s/%s_%s_greenup_date_%d.tif',path_map,name_vis{vi},product,years(yr)));  
    A_phn = double(A_phn);
    A_phn(A_phn<30 | A_phn>180)=nan;
    % A_phn = gpuArray(A_phn);

    for pk = 1:length(file_park)
        A_park = readgeoraster(sprintf('%s/park_%d_mask.tif',file_park(pk).folder,park_name(pk)));
        A_park = double(A_park);
        A_park(A_park~=0)=1;
        A_park(A_park==0)=nan; 
        % A_park = gpuArray(A_park);

        A_temp = A_park .* A_ecm; % imagesc(A_temp_park)
        A_temp(A_temp==0)=nan;
        % A_temp=gather(A_temp);
        park_pixel_full(pk,1) = sum(~isnan(A_temp),'all','omitmissing');
        pheno_date_full(yr,pk,1) = mean(A_phn(~isnan(A_temp)),'all','omitmissing');

        % 10 m interval
        for bf =1:length(file_buff)-1
            A_temp=A_buff(:,:,bf) + A_buff(:,:,bf+1);
            A_temp(isnan(A_park))=nan;

            A_temp = A_temp .* A_ecm;
            % A_temp=gather(A_temp);
            pheno_date_bd_int10(yr,pk,bf) = mean(A_phn(A_temp==1),'all','omitmissing');
            pheno_date_in_int10(yr,pk,bf) = mean(A_phn(A_temp==2),'all','omitmissing');

            park_pixel_bd_int10(yr,pk,bf) = sum((A_temp==1),'all','omitmissing');
            park_pixel_in_int10(yr,pk,bf) = sum((A_temp==2),'all','omitmissing');
        end


        % 20 m interval
        A_temp=(A_park + A_buff(:,:,3)).* A_ecm;
        A_temp(isnan(A_park))=nan;  
        % A_temp=gather(A_temp);
        pheno_date_bd_int20(yr,pk,1) = mean(A_phn(A_temp==1),'all','omitmissing');
        park_pixel_bd_int20(yr,pk,1) = sum((A_temp==1),'all','omitmissing');

        A_temp=(A_buff(:,:,3) + A_buff(:,:,5)).* A_ecm;
        A_temp(isnan(A_park))=nan;  
        % A_temp=gather(A_temp);
        pheno_date_bd_int20(yr,pk,2) = mean(A_phn(A_temp==1),'all','omitmissing');
        park_pixel_bd_int20(yr,pk,2) = sum((A_temp==1),'all','omitmissing');

        % 30 m interval
        A_temp=(A_park + A_buff(:,:,4)).* A_ecm;
        A_temp(isnan(A_park))=nan;  
        % A_temp=gather(A_temp);
        pheno_date_bd_int30(yr,pk) = mean(A_phn(A_temp==1),'all','omitmissing');
        park_pixel_bd_int30(yr,pk) = sum((A_temp==1),'all','omitmissing');

        % 40 m interval
        A_temp=(A_park + A_buff(:,:,5)).* A_ecm;
        A_temp(isnan(A_park))=nan;
        % A_temp=gather(A_temp);
        pheno_date_bd_int40(yr,pk) = mean(A_phn(A_temp==1),'all','omitmissing');
        park_pixel_bd_int40(yr,pk) = sum((A_temp==1),'all','omitmissing');

        % 50 m interval
        A_temp=(A_park + A_buff(:,:,6)).* A_ecm;
        A_temp(isnan(A_park))=nan;  
        % A_temp=gather(A_temp);
        pheno_date_bd_int50(yr,pk) = mean(A_phn(A_temp==1),'all','omitmissing');
        park_pixel_bd_int50(yr,pk) = sum((A_temp==1),'all','omitmissing');
    end
    fprintf('%d done \n',years(yr))
end


park_pixel_full= gather(park_pixel_full);
park_pixel_in_int10 = gather(park_pixel_in_int10);
pheno_date_full = gather(pheno_date_full);
pheno_date_in_int10 = gather(pheno_date_in_int10);

park_pixel_bd_int10 = gather(park_pixel_bd_int10);
park_pixel_bd_int20 = gather(park_pixel_bd_int20);
park_pixel_bd_int30 = gather(park_pixel_bd_int30);
park_pixel_bd_int40 = gather(park_pixel_bd_int40);
park_pixel_bd_int50 = gather(park_pixel_bd_int50);

pheno_date_bd_int10 = gather(pheno_date_bd_int10);
pheno_date_bd_int20 = gather(pheno_date_bd_int20);
pheno_date_bd_int30 = gather(pheno_date_bd_int30);
pheno_date_bd_int40 = gather(pheno_date_bd_int40);
pheno_date_bd_int50 = gather(pheno_date_bd_int50);
