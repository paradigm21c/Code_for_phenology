
path_figure = 'D:/SetoLab/Phenology/figure/Preseason_SDC_full_DS/all_years/HUMID';
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


vi =3;%:length(vin) 
dayn = 1; phn = 1; % Greenup

filter_ranges = [30, 180; 150, 240; 180, 300; 240, 350];
filter_min = filter_ranges(phn , 1);
filter_max = filter_ranges(phn , 2);


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
            temp = temp <= 5 ;
            CA(yr-1,mk) = sum(temp,'all','omitmissing');
            CA(CA==0)=nan;

            temp = tmax_all(ind_heat:ind_pheno,mk);
            temp_HR = temp > 5;
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


%% Histogram - CA : Parks
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
park_year = {'2000-2009','2000-2018','2009-2018'};

CA_all =[];
mv_all =[];
park_all = [];

for ps =1:4
    if ps ==1
        temp_park = ind_park_vs;
    elseif ps ==2
        temp_park = ind_park_s;
    elseif ps ==3
        temp_park = ind_park_m;
    elseif ps ==4
        temp_park = ind_park_l;
    elseif ps ==5
        temp_park = ind_park_all;
    end

    A = [];
    An = [];

    for dp =1:3

        if dp == 1
            mv = 1:10;
            cp = 1;
        elseif dp == 3
            mv = 10:19;
            cp = length(newmap);
        elseif dp == 2
            mv =1:19;
            cp = length(newmap)/2;
        end

        X = reshape(CA(mv,temp_park),[],1);
        % X(X<30)=nan;
        Y = reshape(HR(mv,temp_park),[],1);
        % Y(isnan(X))=nan;
        % X(isnan(Y))=nan;
        % 
        % X(isnan(X))=[];
        % Y(isnan(Y))=[];
        % 
        temp_n = repmat(string(park_year{dp}),length(X),1);
        A  = [A;X];
        An = [An;temp_n];
    end

    temp_park = repmat(string(park_categor{ps}),length(A),1);

    CA_all =[CA_all;A];
    mv_all =[mv_all;An];
    park_all = [park_all;temp_park];

end

table_ca = table(CA_all,mv_all,park_all,'VariableNames',{'CA','MV','Park'});
table_ca.MV = categorical(table_ca.MV ,park_year);
table_ca.Park = categorical(table_ca.Park ,park_categor);

% BARPLOT
for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}));
        % histogram(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'Normalization','probability','NumBins',50);
    end
end

bar_se = bar_std./ sqrt(n);

x_label = categorical(park_categor,park_categor);
bar_data =[];
bar_err=[];

bar_data = bar_mean;
bar_err= bar_std;
bar_err = bar_se;


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
ylim([75,90])
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

yline(ax.YLim, 'Color', 'k', 'LineWidth', 0.5);
xline(ax.XLim(1), 'Color', 'k', 'LineWidth', 0.5);
ax.Box = 'off';                % Turn off top and right axis lines
ax.TickDir = 'out';            % Ticks pointing outward   

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,vin{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_CA_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,vin{vi}),'BackgroundColor','none','ContentType','vector')



%% Histogram - HR : Parks
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
park_year = {'2000-2009','2000-2018','2009-2018'};

HR_all =[];
mv_all =[];
park_all = [];

for ps =1:4
    if ps ==1
        temp_park = ind_park_vs;
    elseif ps ==2
        temp_park = ind_park_s;
    elseif ps ==3
        temp_park = ind_park_m;
    elseif ps ==4
        temp_park = ind_park_l;
    elseif ps ==5
        temp_park = ind_park_all;
    end
    
    A = [];
    An = [];

    for dp =1:3
    
        if dp == 1
            mv = 1:10;
            cp = 1;
        elseif dp == 3
            mv = 10:19;
            cp = length(newmap);
        elseif dp == 2
            mv =1:19;
            cp = length(newmap)/2;
        end
        
        X = reshape(CA(mv,temp_park),[],1);
        % X(X<30)=nan;
        Y = reshape(HR(mv,temp_park),[],1);
 
        % 
        temp_n = repmat(string(park_year{dp}),length(Y),1);
        A  = [A;Y];
        An = [An;temp_n];
    end
    
    temp_park = repmat(string(park_categor{ps}),length(A),1);

    HR_all =[HR_all;A];
    mv_all =[mv_all;An];
    park_all = [park_all;temp_park];

end

table_ca = table(HR_all,mv_all,park_all,'VariableNames',{'CA','MV','Park'});
table_ca.MV = categorical(table_ca.MV ,park_year);
table_ca.Park = categorical(table_ca.Park ,park_categor);

% BARPLOT
for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}));
        % histogram(table_ca.CA(table_ca.Park == park_categor{ps} & table_ca.MV == park_year{mv}),'Normalization','probability','NumBins',50);
    end
end

bar_se = bar_std./ sqrt(n);

x_label = categorical(park_categor,park_categor);
bar_data =[];
bar_err=[];

bar_data = bar_mean;
bar_err= bar_std;
bar_err = bar_se;


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

ylim([500,1100])
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

yline(ax.YLim, 'Color', 'k', 'LineWidth', 0.5);
xline(ax.XLim(1), 'Color', 'k', 'LineWidth', 0.5);
ax.Box = 'off';                % Turn off top and right axis lines
ax.TickDir = 'out';            % Ticks pointing outward   


hold off


savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_HR_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,vin{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_HR_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,vin{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_HR_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,vin{vi}),'BackgroundColor','none','ContentType','vector')
