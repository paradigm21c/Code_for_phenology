
% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
   parpool('Processes',24); % Adjust the number of workers as needed
end

% product = 'MODIS';
product = 'SDC_full';

% raster_files_directory = 'C:/Temp/MODIS/SDC500/Pre_Fusion_Crop/';
raster_files_directory = 'C:/Temp/MODIS/SDC500/Fusion_500m_full/';

path_mask = 'D:/SetoLab/Phenology/mask';
path_daymet = sprintf('C:/Temp/Temperature/DayMet/%s',product);

% files_mask = dir(sprintf('%s/new_parks/*.tif',path_mask));
% files_mask = dir(sprintf('%s/parks_arc_buf/*.tif',path_mask));
files_mask = dir(sprintf('%s/parks_arc_buf_%s/*.tif',path_mask,product));
mkn = length(files_mask);

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask'};
% ecotype = {'combined_mask'};
ecn = length(ecotype);


sdc_files = dir(sprintf('%s/*.tif',raster_files_directory));

for fn =1:length(sdc_files)
    temp = split(sdc_files(fn).name,'_');
    temp = temp{2};
    file_dates(fn,1) = datetime(str2num(temp(1:4)),str2num(temp(5:6)),str2num(temp(7:8)),"Format",'yyyyMMdd');
end

evi2_all = nan(length(sdc_files),length(ecotype),length(files_mask));
tmax_all = nan(length(sdc_files),length(ecotype),length(files_mask));

% sdcn = 2800; %length(sdc_files);

for fn=1:length(sdc_files)
    tic;
    [A, ~] = readgeoraster(sprintf('%s/%s',sdc_files(fn).folder,sdc_files(fn).name));
    A=double(A); A=A./10000;
    red = A(:,:,3); red(red>1 | red<0)=nan;
    nir = A(:,:,4); nir(nir>1 | nir<0)=nan;
    ndvi = (nir-red)./(nir+red);
    ndvi(ndvi>1 | ndvi<0)=nan;

    nirv = ndvi .* nir;
    nirv(nirv>1 | nirv<0)=nan;

    evi2 = 2.5 * (nir - red) ./ (nir + 2.4 * red + 1); %(nir-red)./(nir+red);
    evi2(evi2>1 | evi2<0)=nan;
    temp_date = file_dates(fn);

    A_tmax = nan(size(A,1),size(A,2));
    A_tmin = nan(size(A,1),size(A,2));
    A_srad = nan(size(A,1),size(A,2));
    A_dayl = nan(size(A,1),size(A,2));
    A_prcp = nan(size(A,1),size(A,2));
    % A_vp = nan(size(A,1),size(A,2));

    if ~(temp_date.Month==12 && temp_date.Day== 31)
        [A_tmax, ~] = readgeoraster(sprintf('%s/tmax/%s/tmax_%s%02d%02d_daymet.tif',path_daymet, num2str(temp_date.Year),...
            num2str(temp_date.Year),temp_date.Month,temp_date.Day));
        A_tmax=double(A_tmax);

        [A_tmin, ~] = readgeoraster(sprintf('%s/tmin/%s/tmin_%s%02d%02d_daymet.tif',path_daymet, num2str(temp_date.Year),...
            num2str(temp_date.Year),temp_date.Month,temp_date.Day));
        A_tmin=double(A_tmin);

        [A_srad, ~] = readgeoraster(sprintf('%s/srad/%s/srad_%s%02d%02d_daymet.tif',path_daymet, num2str(temp_date.Year),...
            num2str(temp_date.Year),temp_date.Month,temp_date.Day));
        A_srad=double(A_srad);

        [A_dayl, ~] = readgeoraster(sprintf('%s/dayl/%s/dayl_%s%02d%02d_daymet.tif',path_daymet, num2str(temp_date.Year),...
            num2str(temp_date.Year),temp_date.Month,temp_date.Day));
        A_dayl=double(A_dayl);

        [A_prcp, ~] = readgeoraster(sprintf('%s/prcp/%s/prcp_%s%02d%02d_daymet.tif',path_daymet, num2str(temp_date.Year),...
            num2str(temp_date.Year),temp_date.Month,temp_date.Day));
        A_prcp=double(A_prcp);

        % [A_vp, ~] = readgeoraster(sprintf('%s/vp/%s/vp_%s%02d%02d_daymet.tif',path_daymet, num2str(temp_date.Year),...
        %     num2str(temp_date.Year),temp_date.Month,temp_date.Day));
        % A_vp=double(A_vp); 

    end

    % Daily total radiation (MJ/m2/day) can be calculated as follows: ((srad (W/m2) * dayl (s/day)) / l,000,000)
    for mk = 1:mkn
        [A_mask2, ~] = readgeoraster(sprintf('%s/%s',files_mask(mk).folder,files_mask(mk).name));
        A_mask2= double(A_mask2);
        A_mask2(A_mask2==0)=nan; A_mask2(~isnan(A_mask2))=1; 
   
        evi2_area = evi2 .* A_mask2;
        ndvi_area = ndvi .* A_mask2;
        nirv_area = nirv .* A_mask2;

        tmax_area = A_tmax .* A_mask2;
        tmin_area = A_tmin .* A_mask2;
        prcp_area = A_prcp .* A_mask2;
        % vp_area = A_vp .* A_mask2;
        dayl_area = A_dayl .* A_mask2;
        drad_area = ((A_srad.* A_dayl)./1000000) .* A_mask2;
    
        for ec = 1:ecn
            [A_mask, ~] = readgeoraster(sprintf('%s/%s/%s_%s.tif',path_mask,product,ecotype{ec},product));
            A_mask= double(A_mask);
            A_mask(A_mask==0)=nan;
            
            evi2_eco = evi2_area .*A_mask;
            ndvi_eco = ndvi_area .*A_mask;
            nirv_eco = nirv_area .*A_mask;

            tmax_eco = tmax_area .*A_mask;
            tmin_eco = tmin_area .*A_mask;
            prcp_eco = prcp_area .*A_mask;
            % vp_eco = vp_area .*A_mask;
            drad_eco = drad_area .*A_mask;
            dayl_eco = dayl_area .*A_mask;

            evi2_all(fn,ec,mk) = mean(evi2_eco(~isnan(evi2_eco)),'omitnan');
            ndvi_all(fn,ec,mk) = mean(ndvi_eco(~isnan(ndvi_eco)),'omitnan');
            nirv_all(fn,ec,mk) = mean(nirv_eco(~isnan(nirv_eco)),'omitnan');

            tmax_all(fn,ec,mk) = mean(tmax_eco(~isnan(tmax_eco)),'omitnan');
            tmin_all(fn,ec,mk) = mean(tmin_eco(~isnan(tmin_eco)),'omitnan');
            prcp_all(fn,ec,mk) = mean(prcp_eco(~isnan(prcp_eco)),'omitnan');
            % vp_all(fn,ec,mk) = mean(vp_eco(~isnan(vp_eco)),'omitnan');
            drad_all(fn,ec,mk) = mean(drad_eco(~isnan(drad_eco)),'omitnan');
            dayl_all(fn,ec,mk) = mean(dayl_eco(~isnan(dayl_eco)),'omitnan');
        end
    end
    fprintf('%d / out of %d\n',fn,length(sdc_files))
    toc;
end

save(sprintf('evi2_all_arc_%s.mat',product),'evi2_all')
save(sprintf('ndvi_all_arc_%s.mat',product),'ndvi_all')
save(sprintf('nirv_all_arc_%s.mat',product),'nirv_all')

save(sprintf('tmax_all_arc_%s.mat',product),'tmax_all')
save(sprintf('tmin_all_arc_%s.mat',product),'tmin_all')
save(sprintf('prcp_all_arc_%s.mat',product),'prcp_all')
% save(sprintf('vp_all_arc_%s.mat',product),'vp_all')
save(sprintf('drad_all_arc_%s.mat',product),'drad_all')
save(sprintf('dayl_all_arc_%s.mat',product),'dayl_all')

save(sprintf('file_dates_%s.mat',product),'file_dates')
