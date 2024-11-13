
% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
   parpool('Processes',24); % Adjust the number of workers as needed
end

product = 'ST';

raster_files_directory = 'D:/SetoLab/Landsat_ARD/Pre_Fusion_Clear_ST/';
path_mask = sprintf('D:/SetoLab/Phenology/mask');

% files_mask = dir(sprintf('%s/new_parks/*.tif',path_mask));
files_mask = dir(sprintf('%s/parks_arc_buf/*.tif',path_mask));
mkn = length(files_mask);

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask'};
% ecotype = {'combined_mask'};
ecn = length(ecotype);

st = 1;

sdc_files = dir(sprintf('%s/*.tif',raster_files_directory));

for fn =1:length(sdc_files)
    temp = split(sdc_files(fn).name,'_');
    temp = temp{2};
    file_dates(fn,1) = datetime(str2num(temp(1:4)),str2num(temp(5:6)),str2num(temp(7:8)),"Format",'yyyyMMdd');
end

ST_all = nan(length(sdc_files),length(ecotype),length(files_mask));

sdcn = length(sdc_files);

parfor fn=1:length(sdc_files)
    [A, ~] = readgeoraster(sprintf('%s/%s',sdc_files(fn).folder,sdc_files(fn).name));
    A(A>61441 | A<291)=nan;
    A=double(A); A=A*0.00341802 + 149.0;
    
    for mk = 1:mkn
        [A_mask2, ~] = readgeoraster(sprintf('%s/%s',files_mask(mk).folder,files_mask(mk).name));
        A_mask2= double(A_mask2);
        A_mask2(A_mask2==0)=nan; A_mask2(~isnan(A_mask2))=1; 
   
        ST_area = A .* A_mask2;
    
        for ec = 1:ecn
            [A_mask, ~] = readgeoraster(sprintf('%s/%s/%s_%s.tif',path_mask,product,ecotype{ec},product));
            A_mask= double(A_mask);
            A_mask(A_mask==0)=nan;
            
            ST_eco = ST_area .*A_mask;

            ST_all(fn,ec,mk) = mean(ST_eco(~isnan(ST_eco)),'omitnan');

        end
    end
    fprintf('%d / out of %d\n',fn,sdcn)
end

save(sprintf('ST_all_arc3_%s.mat',product),'ST_all')
save(sprintf('file_dates_%s.mat',product),'file_dates')
