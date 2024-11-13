% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
    parpool('Processes',18); % Adjust the number of workers as needed
end


product = 'PF';
product = 'SDC_full';

path_mask = 'D:/SetoLab/Phenology/mask';
mask_park = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));
mask_ecm = dir(sprintf('%s/%s/*.tif',path_mask,product));

clear temp
for pr = 1:length(mask_park)
    temp(pr,1) = string(mask_park(pr).name);
end
temp =split(temp,'_');
park_name =sort(str2double(temp(:,2)));


park_eco_ratio = nan(length(mask_park),length(mask_ecm));
park_eco_num = nan(length(mask_park),length(mask_ecm));
park_num = nan(length(mask_park),1); 

for ec = 1:length(mask_ecm)
    [A_temp, ~] = readgeoraster(sprintf('%s/%s',mask_ecm(ec).folder,mask_ecm(ec).name));
    A_temp = double(A_temp); A_temp(A_temp==0)=nan;
    A_ecm(:,:,ec)=A_temp;
end

for pk = 1:length(mask_park)
    [A_park, ~] = readgeoraster(sprintf('%s/park_%d_mask.tif',mask_park(pk).folder,park_name(pk)));
    A_park = double(A_park); A_park(A_park==0)=nan;

    num_park = sum(~isnan(A_park),'all');
    park_num(pk,1)=num_park;

    parfor ec = 1:length(mask_ecm)        
        A_ecm_park = A_park .* A_ecm(:,:,ec);
        num_ecm_park = sum(~isnan(A_ecm_park),'all');
        park_eco_num(pk,ec) = num_ecm_park;
        park_eco_ratio(pk,ec) = num_ecm_park./num_park;
    end
    fprintf('%d park done \n',park_name(pk))
end

park_eco_ratio=gather(park_eco_ratio);
park_eco_num=gather(park_eco_num);
park_num=gather(park_num);

save(sprintf("park_eco_ratio_%s_DS.mat",product),'park_eco_ratio');
save(sprintf("park_eco_num_%s_DS.mat",product),'park_eco_num');
save(sprintf("park_num_%s_DS.mat",product),'park_num');