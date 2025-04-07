
path_figure = 'D:/SetoLab/Phenology/figure/Preseason_SDC_full_DS/all_years/HUMID/';
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

product='SDC_full';
path_cal = 'D:/SetoLab/Phenology/data_cal/Parks_polygon/ParkDS/HUMID';
path_mask = sprintf('D:/SetoLab/Phenology/mask/');
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));

park_eco_ratio = importdata(sprintf('park_eco_ratio_DS.mat'));
tree = [4, 8, 10, 15];
park_veg_ratio=park_eco_ratio(:,tree);
park_veg_ratio=park_veg_ratio./sum(park_veg_ratio,2);
park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};

date_all=importdata(sprintf('%s/date_HUMID_DS_%s_parkBD.mat',path_cal,product));  %(fn,mk)
tmax_all=importdata(sprintf('%s/tmax_HUMID_DS_%s_parkBD.mat',path_cal,product));  %(fn,mk)
tmin_all=importdata(sprintf('%s/tmin_HUMID_DS_%s_parkBD.mat',path_cal,product));  %(fn,mk)
tmean_all = (tmax_all + tmin_all)./2;

park_num= importdata("park_num_SDC_full_DS.mat");
ind_park_all = 1:length(park_num);
ind_park_vs = find(park_num<=11); % 0.09 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>=12 & park_num<=50); % 0.09 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>=51 & park_num<=200); % 4.59 to 18 hectares
ind_park_l = find(park_num>=201); % 18.09 to more hectares

% Polygon-based
cal_flag =1;
vi_all=importdata(sprintf('%s/vi_range_threshold_%s.mat',path_cal,product));  %(yr,ec,mk,phn)
vi_all=vi_all{1,1}; % ps == 1

% % Pixel-based
% cal_flag =2;
% vi_all=importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/vi_all_%s_DS.mat',product));  %phn_park(pn,yr,phn);


phenology_names = {'Greenup', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};
vin = {'EVI2','NDVI','NIRv'};
years=1999:2018;
mkn = length(files_mask);
etn={'Tree'};
tree = [1 3 4 5];
etype = tree;

dayn = 1; phn = 1; % Greenup

filter_ranges = [30, 180; 150, 240; 180, 300; 240, 350];
filter_min = filter_ranges(phn , 1);
filter_max = filter_ranges(phn , 2);

vi =3;%:length(vin) 

 % vi =1 , 1 means date // %(yr,ec,mk,phn)
% pheno_date(13,:,:,:)=NaT; % year 2012
% pheno_date(14,:,:,:)=NaT; % year 2013

if cal_flag ==1 % pheno_date(yr,ec,mk,phn)
    pheno_date = vi_all{vi,1};
    pheno_date.Format='DDD';
    pheno_date = str2double(string(pheno_date));
    pheno_date(pheno_date<filter_min | pheno_date> filter_max)=nan;

    temp = [];
    tmean = [];
     % (pk, yr, ec)
    for i = 1:length(etype)
        temp(:,:,i) = permute(squeeze(pheno_date(:, etype(i), :, phn)), [2, 1]) .* park_veg_ratio(:, i);
    end

    % Sum ignoring NaNs
    temp =  sum(temp,3,'omitnan'); 
    temp(temp==0)=nan;
    vi_pheno_date = permute(temp,[2,1]); % (yr,mk)

elseif cal_flag ==2 %pheno_date(pn,yr,phn)
    if vi ==3
        pheno_date = vi_all{vi-1,1};
    elseif vi ==1
        pheno_date = vi_all{vi,1};
    end

    pheno_date= pheno_date(:,:,phn);
    vi_pheno_date= permute(pheno_date,[2, 1]);
end

vi_pheno_date(vi_pheno_date<filter_min | vi_pheno_date> filter_max)=nan;

CA=nan(length(years)-1,mkn);
HR=nan(length(years)-1,mkn);
WT_mean=nan(length(years)-1,mkn);
ST_mean=nan(length(years)-1,mkn);
Annual=nan(length(years)-1,mkn);

for yr = 2:length(years)
    date_chilling = datetime(years(yr-1),11,01);
    date_limit = datetime(years(yr),06,01);
    date_heat = datetime(years(yr),01,01);

    date_d1 = datetime(years(yr),01,01);
    date_dend = datetime(years(yr),12,30);

    date_mar = datetime(years(yr),03,1);
    date_apr = datetime(years(yr),04,30);

    date_dec = datetime(years(yr-1),12,01);
    date_feb = datetime(years(yr),03,01);

    ind_chilling = find(date_all == date_chilling);
    ind_heat = find(date_all == date_heat);

    ind_dec = find(date_all == date_dec);
    ind_feb = find(date_all == date_feb);

    ind_mar = find(date_all == date_mar);
    ind_apr = find(date_all == date_apr);
    

    ind_d1 = find(date_all == date_d1);
    ind_dend = find(date_all == date_dend);

    for mk=1:mkn
        date_sos = vi_pheno_date(yr-1,mk);
        if ~isnan(date_sos)
            [month, day]= doy2day(years(yr),round(date_sos));
            temp_date = datetime(years(yr),month,day);
            ind_pheno = find(date_all == temp_date);
        
            temp = tmean_all(ind_chilling:ind_pheno,mk);
            temp = temp <= 5 & temp > -10;
            CA(yr-1,mk) = sum(temp,'all','omitmissing');

            temp = tmax_all(ind_heat:ind_pheno,mk);
            temp_HR = temp >= 5;
            HR(yr-1,mk) = sum(temp(temp_HR),'all','omitmissing'); 
        end

        temp = tmean_all(ind_dec:ind_feb,mk);
        WT_mean(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmean_all(ind_mar:ind_apr,mk);
        ST_mean(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmax_all(ind_dec:ind_feb,mk);
        WT_max(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmax_all(ind_mar:ind_apr,mk);
        ST_max(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmin_all(ind_dec:ind_feb,mk);
        WT_min(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmin_all(ind_mar:ind_apr,mk);
        ST_min(yr-1,mk) = mean(temp,"all","omitmissing");

        temp = tmean_all(ind_d1:ind_dend,mk);
        Annual(yr-1,mk) = mean(temp,"all","omitmissing");
    end
end


%%
% X-axis: Deciduous ratio
% 
% Y-axis: SOS
% 
% Color: Park size


%% Barplot SOS dates in park size.
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% park_year = {'2000-209','2013-2022','2000-2022'};
park_year = {'2000-2009','2009-2018','2000-2018'};

temp_park = ind_park_vs;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:10,temp_park); % MV =10
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(10:end,temp_park); % MV =10
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_pocket = [X;Y;Z];
Park_pocket = repmat(string(park_categor(1)),length(SOS_pocket),1);
Years_pocket = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_pocket = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_s;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:10,temp_park); % MV =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(10:end,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_neighbor = [X;Y;Z];
Park_neighbor = repmat(string(park_categor(2)),length(SOS_neighbor),1);
Years_neighbor = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_neighbor = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_m;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:10,temp_park); % MV =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(10:end,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_community = [X;Y;Z];
Park_community = repmat(string(park_categor(3)),length(SOS_community),1);
Years_community = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_community = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_l;

% X = vi_pheno_date(2:11,temp_park); % MV =10
X = vi_pheno_date(1:10,temp_park); % MV =12
X=mean(X,1,'omitmissing');
% Y = vi_pheno_date(14:23,temp_park); % MV =10
Y = vi_pheno_date(10:end,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = vi_pheno_date(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

SOS_regional = [X;Y;Z];
Park_regional = repmat(string(park_categor(4)),length(SOS_regional),1);
Years_regional = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];

Mean_regional = [X_temp;Y_temp;Z_temp];

tbl_sos = table([SOS_pocket;SOS_neighbor;SOS_community; SOS_regional],...
    [Park_pocket;Park_neighbor;Park_community;Park_regional],...
    [Years_pocket;Years_neighbor;Years_community;Years_regional],...
    'VariableNames',["SOS","Park","MV"]);

tbl_sos.Park = categorical(tbl_sos.Park,park_categor);


bar_mean =[];
bar_std=[];
n = [];

for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(tbl_sos.SOS(tbl_sos.Park == park_categor{ps} & tbl_sos.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(tbl_sos.SOS(tbl_sos.Park == park_categor{ps} & tbl_sos.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(tbl_sos.SOS(tbl_sos.Park == park_categor{ps} & tbl_sos.MV == park_year{mv}));
        % histogram(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'Normalization','probability','NumBins',50);
    end
end

bar_se = bar_std./ sqrt(n);

x_label = categorical(park_categor,park_categor);
bar_data =[];
bar_err=[];

bar_data(1,:) = bar_mean(1,:);
bar_data(2,:) = bar_mean(3,:);
bar_data(3,:) = bar_mean(2,:);

bar_err(1,:) = bar_std(1,:);
bar_err(2,:) = bar_std(3,:);
bar_err(3,:) = bar_std(2,:);

bar_err(1,:) = bar_se(1,:);
bar_err(2,:) = bar_se(3,:);
bar_err(3,:) = bar_se(2,:);

model_series = permute(bar_data,[2,1]);
model_error = permute(bar_err,[2,1]); 

figure('Position',[600,200,500,200]);
hold on
b = bar(x_label ,model_series, 'grouped');
b(1).FaceColor = newmap(1,:);
b(1).EdgeColor = newmap(1,:);
b(1).LineWidth = 1;
b(1).FaceAlpha = 0.3;

b(2).FaceColor = newmap(end/2,:);
b(2).EdgeColor = newmap(end/2,:);
b(2).LineWidth = 1;
b(2).FaceAlpha = 0.3;

b(3).FaceColor = newmap(end,:);
b(3).EdgeColor = newmap(end,:);
b(3).LineWidth = 1;
b(3).FaceAlpha = 0.3;

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.YAxisLocation=fig_y_axis_loc;

ylim([90,110])

if vi ==1
    ylim([80,100])
end

yline(ax.YLim, 'Color', 'k', 'LineWidth', 0.5);
xline(ax.XLim(1), 'Color', 'k', 'LineWidth', 0.5);
ax.Box = 'off';                % Turn off top and right axis lines
ax.TickDir = 'out';            % Ticks pointing outward   

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_SOS_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_SOS_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,vin{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_SOS_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,vin{vi}),'BackgroundColor','none','ContentType','vector')
