% Comparing pixel-based SOS and polygon-based SOS in park size-level

product = 'SDC_full';

path_figure = 'D:/SetoLab/Phenology/figure/Preseason_SDC_full_DS/all_years/';
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

path_cal = 'D:/SetoLab/Phenology/data_cal/Parks_polygon/ParkDS';
path_mask = sprintf('D:/SetoLab/Phenology/mask/');
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));

park_eco_ratio = importdata(sprintf('park_eco_ratio_%s_DS.mat',product));
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);
park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};

date_all=importdata(sprintf('%s/date_daymet_DS_%s.mat',path_cal,product));  %(fn,ec,mk)
tmax_all=importdata(sprintf('%s/tmax_daymet_DS_%s.mat',path_cal,product));  %(fn,ec,mk)
tmin_all=importdata(sprintf('%s/tmin_daymet_DS_%s.mat',path_cal,product));  %(fn,ec,mk)
tmean_all = (tmax_all + tmin_all)./2;

park_num= importdata("park_num_DS.mat");
ind_park_all = 1:length(park_num);
ind_park_vs = find(park_num<=11); % 0.09 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>=12 & park_num<=50); % 0.09 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>=51 & park_num<=200); % 4.59 to 18 hectares
ind_park_l = find(park_num>=201); % 18.09 to more hectares

vi_all=importdata(sprintf('%s/vi_range_threshold_%s.mat',path_cal,product));  %(yr,ec,mk,phn)
vi_all=vi_all{1,1}; % ps == 1

phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};
vin = {'EVI2','NDVI','NIRv'};
years=2000:2022;
mkn = length(files_mask);
etn={'Tree'};
tree = [1 3 4 5];
etype = tree;

vi =1; phn = 1; % Greenup

filter_ranges = [30, 180; 150, 240; 180, 300; 240, 350];
filter_min = filter_ranges(phn , 1);
filter_max = filter_ranges(phn , 2);

pheno_date = vi_all{vi,1}; % vi =1 , 1 means date // %(yr,ec,mk,phn)
% pheno_date(13,:,:,:)=NaT; % year 2012
% pheno_date(14,:,:,:)=NaT; % year 2013

pheno_date.Format='DDD';
pheno_date = str2double(string(pheno_date));
pheno_date(pheno_date<filter_min | pheno_date> filter_max)=nan;

temp = [];
 % (pk, yr, ec)
for i = 1:length(etype)
    temp(:,:,i) = permute(squeeze(pheno_date(:, etype(i), :, phn)), [2, 1]) .* park_veg_ratio(:, i);
end

% Sum ignoring NaNs
temp =  sum(temp,3,'omitnan'); 
temp(temp==0)=nan;
data_polygon  = permute(temp,[2,1]); % (yr,mk)
data_polygon(data_polygon<filter_min | data_polygon> filter_max)=nan;
   

% load(sprintf('file_%s_pheno_date_inNbound_cpu_%s.mat',vin{vi},product)) % (yr,pk)
% data_pixel = pheno_date_full;

vi_all=importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/vi_all_%s_DS.mat',product));
if vi ==3
    pheno_date = vi_all{vi-1,1};
elseif vi ==1
    pheno_date = vi_all{vi,1};
end

pheno_date= pheno_date(:,:,phn);
data_pixel = permute(pheno_date,[2, 1]);

%% Phenodate: Polygon vs Pixel
% X = reshape(data_ploygon(13:end,:),[],1); % NTL data
X = reshape(data_polygon,[],1);
Y = reshape(data_pixel,[],1);

Y(isnan(X))=nan;
X(isnan(Y))=nan;

X(isnan(X))=[];
Y(isnan(Y))=[];

figure;
h=dscatter(X,Y,'plottype','scatter');
h.Colormap = newmap;

ylim([50 150])
xlim([50 150])


hr =refline(1,0);
hr.LineWidth = 2;
hr.Color = 'k';

mdl = fitlm(X,Y);

hr =refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineStyle ='--';
hr.LineWidth = 2;
hr.Color = 'r';

xlabel("Polygon-based SOS")
ylabel("Pixel-based SOS")

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;

