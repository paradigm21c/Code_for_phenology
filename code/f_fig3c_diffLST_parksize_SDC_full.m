% 1. Calculating LST - VANUI relationship
% 2. Calculating buffer VANUI and Phenodate edge

product = 'SDC_full';

newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

np = 2;
if np ==1
    years = 2012:2022;
    fig_maker = 'square';
elseif np ==2
    years = 2000:2022;
    fig_maker = 'o';
end
parks_files = [100,200,300];
% Path
path_all = 'D:/SetoLab/Phenology/';
path_map = sprintf('%s/map/%s/',path_all,product);
path_buff = sprintf('%s/mask/%s_buffer/',path_all,product);
path_park = sprintf('%s/mask/parks_arc_%s_DS/',path_all,product);
path_figure = sprintf('%s/figure/Preseason_%s_DS/all_years/HUMID/',path_all,product);

% Files
file_park = dir(sprintf('%s/*.tif',path_park));

park_num= importdata(sprintf("park_num_%s_DS.mat",product));
ind_park_all=1:length(park_num);
ind_park_vs = find(park_num<=11); % 0.09 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>=11 & park_num<=50); % 0.9 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>=50 & park_num<=200); % 4.50 to 18 hectares
ind_park_l = find(park_num>=200); % 18.09 to more hectares

park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};
park_size_name = {'Pkt','Nbrhd','Cmnty','Rgnl','All'};
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl','All'};

phenology_names = {'SOS', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
name_vis = {'EVI2','NIRv','NDVI','EVI'};
name_product = {'VNP46A4','LongNTL'};
name_length = {'','_full_length'};

% Park names
park_name = nan(length(file_park),1);
for i=1:length(file_park)
    temp = split(file_park(i).name,'_');
    park_name(i) = str2double(temp{2});
end
park_name = sort(park_name);

vi=1;

if vi ==1
    load(sprintf('file_EVI2_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np})) % (yr,pk) 
elseif vi ==2
    load(sprintf('file_NIRv_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np}))  % (yr,pk)
end

% LST -> (pn,yr,phn)
LST_park = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_park_%s_DS.mat',product));
LST_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(1)));
LST_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(2)));
LST_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(3)));
LST_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(1)));
LST_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(2)));
LST_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(3)));

vi=1;

%% Bar LST dates in park size. .. test 2
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% park_year = {'2000-209','2013-2022','2000-2022'};
park_year = {'2000-2011','2011-2022','2000-2022'};


% LST
temp_lst = permute(LST_diff200-LST_park, [2 1]);

temp_park = ind_park_vs;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_pocket = [X;Y;Z];
Park_pocket = repmat(string(park_categor(1)),length(LST_pocket),1);
Years_pocket = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_pocket = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_s;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_neighbor = [X;Y;Z];
Park_neighbor = repmat(string(park_categor(2)),length(LST_neighbor),1);
Years_neighbor = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_neighbor = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_m;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_community = [X;Y;Z];
Park_community = repmat(string(park_categor(3)),length(LST_community),1);
Years_community = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];
Mean_community = [X_temp;Y_temp;Z_temp];

temp_park = ind_park_l;

X = temp_lst(1:12,temp_park); % MV =12
X=mean(X,1,'omitmissing');
Y = temp_lst(12:23,temp_park); % MV =12
Y=mean(Y,1,'omitmissing');
Z = temp_lst(:,temp_park);
Z=mean(Z,1,'omitmissing');

X = reshape(X,[],1);
Y = reshape(Y,[],1);
Z = reshape(Z,[],1);

X_temp = mean(X,'all','omitmissing');
Y_temp = mean(Y,'all','omitmissing');
Z_temp = mean(Z,'all','omitmissing');

LST_regional = [X;Y;Z];
Park_regional = repmat(string(park_categor(4)),length(LST_regional),1);
Years_regional = [repmat(string(park_year(1)),length(X),1);
    repmat(string(park_year(2)),length(Y),1);
    repmat(string(park_year(3)),length(Z),1)];

Mean_regional = [X_temp;Y_temp;Z_temp];

tbl_lst = table([LST_pocket;LST_neighbor;LST_community; LST_regional],...
    [Park_pocket;Park_neighbor;Park_community;Park_regional],...
    [Years_pocket;Years_neighbor;Years_community;Years_regional],...
    'VariableNames',["LST","Park","MV"]);

tbl_lst.Park = categorical(tbl_lst.Park,park_categor);


bar_mean =[];
bar_std=[];
n = [];

for ps = 1:length(park_categor)
    for mv = 1:length(park_year)
        bar_mean(mv,ps) = mean(tbl_lst.LST(tbl_lst.Park == park_categor{ps} & tbl_lst.MV == park_year{mv}),'all','omitmissing');
        bar_std(mv,ps) = std(tbl_lst.LST(tbl_lst.Park == park_categor{ps} & tbl_lst.MV == park_year{mv}),'omitmissing');
        n(mv,ps) = length(tbl_lst.LST(tbl_lst.Park == park_categor{ps} & tbl_lst.MV == park_year{mv}));
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
% ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
ylim([-0.3,1.3])

yline(ax.YLim, 'Color', 'k', 'LineWidth', 0.5);
xline(ax.XLim(1), 'Color', 'k', 'LineWidth', 0.5);
ax.Box = 'off';                % Turn off top and right axis lines
ax.TickDir = 'out';            % Ticks pointing outward   

hold off

savefig(gcf,sprintf('%s/Barplot_T_Parksize_X_LST_Y_frequency_G_interval_%s_mv_full_%s.fig',path_figure,product,name_vis{vi}))
saveas(gcf,sprintf('%s/Barplot_T_Parksize_X_LST_Y_frequency_G_interval_%s_mv_full_%s.png',path_figure,product,name_vis{vi}))
exportgraphics(gcf,sprintf('%s/Barplot_T_Parksize_X_LST_Y_frequency_G_interval_%s_mv_full_%s.eps',path_figure,product,name_vis{vi}),'BackgroundColor','none','ContentType','vector')
