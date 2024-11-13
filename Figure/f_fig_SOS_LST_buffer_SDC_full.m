%
product = 'SDC_full';
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

parks_files = [100,200,300];
% Path
path_all = 'D:/SetoLab/Phenology/';
path_map = sprintf('%s/map/%s/',path_all,product);
path_buff = sprintf('%s/mask/%s_buffer/',path_all,product);
path_park = sprintf('%s/mask/parks_arc_%s_DS/',path_all,product);
path_ecm = sprintf('%s/mask/PF_temp/',path_all);
path_figure = sprintf('%s/figure/Preseason_%s_DS/',path_all,product);


% Files
file_park = dir(sprintf('%s/*.tif',path_park));

park_num= importdata("park_num_PF_DS.mat");
ind_park_all=1:length(park_num);
ind_park_vs = find(park_num<=1111); % 0.009 to 1 hectares (1 pixel = 0.09 hectares)
ind_park_s = find(park_num>=1112 & park_num<=5000); % 0.99 to 4.5 hectares (1 pixel = 0.09 hectares)
ind_park_m = find(park_num>=5001 & park_num<=20000); % 4.50 to 18 hectares
ind_park_l = find(park_num>=20100); % 18.09 to more hectares

park_size_name = {"Pocket","Neighborhood","Community","Regional","All"};
park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl','All'};
phenology_names = {'SOS', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};
name_vis = {'EVI2','NIRv'};
name_product = {'VNP46A4','LongNTL'};
name_length = {'','_full_length'};

np = 2;

if np ==1
    years = 2012:2022;
elseif np ==2
    years = 2000:2022;
end


% Park names
park_name = nan(length(file_park),1);
for i=1:length(file_park)
    temp = split(file_park(i).name,'_');
    park_name(i) = str2double(temp{2});
end
park_name = sort(park_name);

sos_lim = [90 115];
sos_lim_tick = [90, 105,115];

phn=1; vi =2;


if vi ==1
    load(sprintf('file_EVI2_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np})) % (yr,pk) 
elseif vi ==2
    load(sprintf('file_NIRv_pheno_date_inNbound_cpu_%s%s.mat',product,name_length{np}))  % (yr,pk)
end

% LST -> (pn,yr,phn)
LST_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(1)));
LST_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(2)));
LST_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_buffer_%s_DS_buff%d.mat',product,parks_files(3)));
LST_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(1)));
LST_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(2)));
LST_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/LST/LST_diff_%s_DS_buff%d.mat',product,parks_files(3)));

num_threshold = 150;
% VANUI -> (pn,yr,phn)
VANUI_buffer100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
VANUI_buffer200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
VANUI_buffer300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_buffer_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));
VANUI_diff100 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(1),name_vis{vi},name_product{np}));
VANUI_diff200 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(2),name_vis{vi},name_product{np}));
VANUI_diff300 = importdata(sprintf('D:/SetoLab/Phenology/data_cal/VANUI/VANUI_thre%d_diff_%s_DS_buff%d_%s_%s.mat',num_threshold,product,parks_files(3),name_vis{vi},name_product{np}));

%% Gradient of phenodate in 10 m interval with different urban park size

for yr = 1:length(years)
    temp_park = [];
    temp_interval = [];
    temp_sos = [];
    temp_mean=[];
    
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
       
        Abd = reshape(pheno_date_bd_int10(yr,ind_temp,1),[],1);
        Bbd = reshape(pheno_date_bd_int10(yr,ind_temp,2),[],1);
        Cbd = reshape(pheno_date_bd_int10(yr,ind_temp,3),[],1);    
        Dbd = reshape(pheno_date_bd_int10(yr,ind_temp,4),[],1);
        Ebd = reshape(pheno_date_bd_int10(yr,ind_temp,5),[],1);
    
    
        Cin = reshape(pheno_date_in_int10(yr,ind_temp,3),[],1);
    
        temp_mean(ps,1) =  mean(Abd,"all","omitmissing");
        temp_mean(ps,2) =  mean(Bbd,"all","omitmissing");
        temp_mean(ps,3) =  mean(Cbd,"all","omitmissing");
        temp_mean(ps,4) =  mean(Dbd,"all","omitmissing");
        temp_mean(ps,5) =  mean(Ebd,"all","omitmissing");
        
        An =repmat("10m",length(Abd),1);
        Bn =repmat("20m",length(Bbd),1);
        Cn =repmat("30m",length(Cbd),1);
        Dn =repmat("40m",length(Dbd),1);
        En =repmat("50m",length(Ebd),1);
    
        temp_park = [temp_park;repmat(string(park_size_name{ps}),length(Abd).*5,1)];
        temp_interval = [temp_interval;[An;Bn;Cn;Dn;En]];
        temp_sos = [temp_sos;[Abd;Bbd;Cbd;Dbd;Ebd]]; 
    
    end
    
    table_all = table(temp_park,temp_interval,temp_sos,...
        'VariableNames',{'Park','Interval','SOS'});
    % ans = table_all.SOS(table_all.Park=='Pocket' & table_all.Interval == '40m');
    
    interval = {'10m','20m','30m','40m','50m'};
    park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
    
    table_all.Interval = categorical(table_all.Interval,interval);
    table_all.Park = categorical(table_all.Park,park_size_name);
    
    ind_temp = table_all.Interval =='40m' | table_all.Interval =='50m';
    
    %
    color_pick = [1, 256/2 , 256];
    
    figure('Position',[600,200,500,200]);
    hold on
    h = boxchart(table_all.Park(~ind_temp),table_all.SOS(~ind_temp),'groupby',table_all.Interval(~ind_temp));
    for ps = 1:3 %length(park_size_name)
        h(ps).JitterOutliers = 'on';
        h(ps).MarkerStyle = '.';
        h(ps).BoxFaceColor = newmap(color_pick(ps),:);
        h(ps).MarkerColor = newmap(color_pick(ps),:);
    end
    ylabel('SOS')
    ylim([60 140])
    title(sprintf('%d',years(yr)))
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = 'bold';
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    % ax.YTick ='';
    hold off
    
    savefig(gcf,sprintf('%s/all_years/Boxchart_T_Year_%d_X_Park_size_Y_%s_%s_%s_10m_int.fig',path_figure,years(yr), phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Boxchart_T_Year_%d_X_Park_size_Y_%s_%s_%s_10m_int.png',path_figure,years(yr), phenology_names{phn},product,name_vis{vi}))
    
    % Radar chart
    % figure('Position',[600,200,600,600]);
    figure('Position',[600,200,350,350]);
    
    X=permute(temp_mean(:,1:3),[2 1]);
    RC=radarChart(X,'Type','Patch');
    % RC.RLim=[95,100]; % EVI2
    RC.RLim=sos_lim; %round([min(X,[],'all')-10,max(X,[],'all')+3]); % NIRv
    RC.RTick=sos_lim_tick; %round([min(X,[],'all')-10,mean(X,"all"),max(X,[],'all')+3]);
    RC.PropName={'Pocket','Neighborhood','Community','Regional','All'};
    RC.ClassName={'10m','20m','30m'};
    RC=RC.draw();
    % RC=RC.legend();
    % RC.setType('Patch');
    
    color_pick = [1, 256/2 , 256];
    for n=1:RC.ClassNum
        RC.setPatchN(n,'FaceColor',newmap(color_pick(n),:),'EdgeColor',newmap(color_pick(n),:));
    end
    
    title(sprintf('%d',years(yr)))
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    % ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Radarchart_T_Year_%d_X_Park_size_Y_%s_%s_%s_10m_int.fig',path_figure,years(yr), phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Radarchart_T_Year_%d_X_Park_size_Y_%s_%s_%s_10m_int.png',path_figure,years(yr), phenology_names{phn},product,name_vis{vi}))
end

%% Boxchart Phenodate annual mean for each park size - 10m interval

temp_park = [];
temp_interval = [];
temp_sos = [];
temp_mean=[];
temp_year=[];
for yr = 1:length(years)
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
       
        Abd = reshape(pheno_date_bd_int10(yr,ind_temp,1),[],1);
        Bbd = reshape(pheno_date_bd_int10(yr,ind_temp,2),[],1);
        Cbd = reshape(pheno_date_bd_int10(yr,ind_temp,3),[],1);    
        Dbd = reshape(pheno_date_bd_int10(yr,ind_temp,4),[],1);
        Ebd = reshape(pheno_date_bd_int10(yr,ind_temp,5),[],1);
    
    
        Cin = reshape(pheno_date_in_int10(yr,ind_temp,3),[],1);
    
        temp_mean(yr,ps,1) =  mean(Abd,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(Bbd,"all","omitmissing");
        temp_mean(yr,ps,3) =  mean(Cbd,"all","omitmissing");
        temp_mean(yr,ps,4) =  mean(Dbd,"all","omitmissing");
        temp_mean(yr,ps,5) =  mean(Ebd,"all","omitmissing");
        
        An =repmat("10m",length(Abd),1);
        Bn =repmat("20m",length(Bbd),1);
        Cn =repmat("30m",length(Cbd),1);
        Dn =repmat("40m",length(Dbd),1);
        En =repmat("50m",length(Ebd),1);
    
        temp_park = [temp_park;repmat(string(park_size_name{ps}),length(Abd).*5,1)];
        temp_year = [temp_year;repmat(string(years(yr)),length(Abd).*5,1)];
        temp_interval = [temp_interval;[An;Bn;Cn;Dn;En]];
        temp_sos = [temp_sos;[Abd;Bbd;Cbd;Dbd;Ebd]]; 
    end
end 

table_years = table(temp_year,temp_park,temp_interval,temp_sos,...
    'VariableNames',{'Year','Park','Interval','SOS'});

interval = {'10m','20m','30m','40m','50m'};
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
year_str = string(years);

table_years.Interval = categorical(table_years.Interval,interval);
table_years.Park = categorical(table_years.Park,park_size_name);
table_years.Year = categorical(table_years.Year,year_str);

color_pick = [1, 256/2 , 256];

for ps = 1:length(park_size_name)
    ind_temp = table_years.Interval =='40m' | table_years.Interval =='50m';
    ind_temp = ~ind_temp & table_years.Park == park_size_name{ps};
    figure;
    hold on
    h = boxchart(table_years.Year(ind_temp),table_years.SOS(ind_temp),'groupby',table_years.Interval(ind_temp));
    for int = 1:3 %length(park_size_name)
        h(int).JitterOutliers = 'on';
        h(int).MarkerStyle = '.';
        h(int).BoxFaceColor = newmap(color_pick(int),:);
        h(int).MarkerColor = newmap(color_pick(int),:);
    end
    ylabel('SOS')
    title(sprintf('%s',park_size_name{ps}))
    ylim([40 140])

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Boxchart_T_Park_%s_X_Year_Y_%s_%s_%s_10m_int_10m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Boxchart_T_Park_%s_X_Year_Y_%s_%s_%s_10m_int_10m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))

    X = permute(squeeze(temp_mean(:,ps,1:3)),[2 1]);
     % Radar chart
    figure('Position',[600,200,350,350]);
    RC=radarChart(X,'Type','Patch');
    % RC.RLim=[95,100]; % EVI2
    RC.RLim=sos_lim; %round([min(X,[],'all')-10,max(X,[],'all')+3]); % NIRv
    RC.RTick=sos_lim_tick; %round([min(X,[],'all')-10,mean(X,"all"),max(X,[],'all')+3]);
    RC.PropName=string(years);
    RC.ClassName={'10m','20m','30m'};
    RC=RC.draw();
    % RC=RC.legend();
    % RC.setType('Patch');
    
    color_pick = [1, 256/2 , 256];
    for n=1:RC.ClassNum
        RC.setPatchN(n,'FaceColor',newmap(color_pick(n),:),'EdgeColor',newmap(color_pick(n),:));
    end
    
    title(sprintf('%s',park_size_name{ps}))
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    % ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Radarchart_T_Park_%s_X_Year_Y_%s_%s_%s_10m_int_10m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Radarchart_T_Park_%s_X_Year_Y_%s_%s_%s_10m_int_10m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    
end

% Mean plot
for ps = 1:length(park_size_name)
    figure;
    hold on
    
    for int = 1:3 %length(park_size_name)
        h = plot(years,squeeze(temp_mean(:,ps,int)));
        h.Color= newmap(color_pick(int),:);
        h.LineWidth = 2;
    end
    
    title(sprintf('%s',park_size_name{ps}))
    ylim([90 120])

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight =fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Plot_mean_T_Park_%s_X_Year_Y_%s_%s_%s_10m_int_10m_width.fig',path_figure, park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Plot_mean_T_Park_%s_X_Year_Y_%s_%s_%s_10m_int_10m_width.png',path_figure, park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    

end

close all


%% Boxchart Phenodate annual mean for each park size - 20m interval

temp_park = [];
temp_interval = [];
temp_sos = [];
temp_mean=[];
temp_year=[];
for yr = 1:length(years)
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
       
        Fbd = reshape(pheno_date_bd_int20(yr,ind_temp,1),[],1);
        Gbd = reshape(pheno_date_bd_int20(yr,ind_temp,2),[],1);

        temp_mean(yr,ps,1) =  mean(Fbd,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(Gbd,"all","omitmissing");

        
        Fn =repmat("20m",length(Fbd),1);
        Gn =repmat("40m",length(Gbd),1);

        temp_park = [temp_park;repmat(string(park_size_name{ps}),length(Fbd).*2,1)];
        temp_year = [temp_year;repmat(string(years(yr)),length(Fbd).*2,1)];
        temp_interval = [temp_interval;[Fn;Gn]];
        temp_sos = [temp_sos;[Fbd;Gbd]]; 
    end
end 

table_20m = table(temp_year,temp_park,temp_interval,temp_sos,...
    'VariableNames',{'Year','Park','Interval','SOS'});

interval = {'20m','40m'};
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
year_str = string(years);

table_20m.Interval = categorical(table_20m.Interval,interval);
table_20m.Park = categorical(table_20m.Park,park_size_name);
table_20m.Year = categorical(table_20m.Year,year_str);

color_pick = [1, 256];

for ps = 1:length(park_size_name)
    ind_temp = table_20m.Park == park_size_name{ps};
    figure;
    hold on
    h = boxchart(table_20m.Year(ind_temp),table_20m.SOS(ind_temp),'groupby',table_20m.Interval(ind_temp));
    for int = 1:2 %length(park_size_name)
        h(int).JitterOutliers = 'on';
        h(int).MarkerStyle = '.';
        h(int).BoxFaceColor = newmap(color_pick(int),:);
        h(int).MarkerColor = newmap(color_pick(int),:);
    end
    ylabel('SOS')
    title(sprintf('%s',park_size_name{ps}))
    ylim([40 140])

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight =fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Boxchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_20m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Boxchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_20m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    

    X = permute(squeeze(temp_mean(:,ps,1:2)),[2 1]);
    % Radar chart
    figure('Position',[600,200,350,350]);
    RC=radarChart(X,'Type','Patch');
    % RC.RLim=[95,100]; % EVI2
    RC.RLim=sos_lim; %round([min(X,[],'all')-10,max(X,[],'all')+3]); % NIRv
    RC.RTick=sos_lim_tick; %round([min(X,[],'all')-10,mean(X,"all"),max(X,[],'all')+3]);
    RC.PropName=string(years);
    RC.ClassName={'20m','40m'};
    RC=RC.draw();
    % RC=RC.legend();
    % RC.setType('Patch');
    
    color_pick = [1, 256];
    for n=1:RC.ClassNum
        RC.setPatchN(n,'FaceColor',newmap(color_pick(n),:),'EdgeColor',newmap(color_pick(n),:));
    end
    
    title(sprintf('%s',park_size_name{ps}))
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight =fig_font_weight;
    ax.FontSize = fig_font_size;
    % ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Radarchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_20m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Radarchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_20m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))


end

color_pick = [1, 256];
% Mean plot
for ps = 1:length(park_size_name)
    figure;
    hold on
    
    for int = 1:2 %length(park_size_name)
        h = plot(years,squeeze(temp_mean(:,ps,int)));
        h.Color= newmap(color_pick(int),:);
        h.LineWidth = 2;
    end
    
    title(sprintf('%s',park_size_name{ps}))
    ylim([90 110])

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Plot_mean_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_20m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Plot_mean_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_20m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))

end


close all
%% Boxchart Phenodate annual mean for each park size - 20m interval- 10m width

temp_park = [];
temp_interval = [];
temp_sos = [];
temp_mean=[];
temp_year=[];
for yr = 1:length(years)
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
       
        Abd = reshape(pheno_date_bd_int10(yr,ind_temp,1),[],1);
        Cbd = reshape(pheno_date_bd_int10(yr,ind_temp,3),[],1);
        Ebd = reshape(pheno_date_bd_int10(yr,ind_temp,5),[],1);

        temp_mean(yr,ps,1) =  mean(Abd,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(Cbd,"all","omitmissing");
        temp_mean(yr,ps,3) =  mean(Ebd,"all","omitmissing");

        
        An =repmat("10m",length(Abd),1);
        Cn =repmat("30m",length(Cbd),1);
        En =repmat("50m",length(Ebd),1);

        temp_park = [temp_park;repmat(string(park_size_name{ps}),length(Abd).*3,1)];
        temp_year = [temp_year;repmat(string(years(yr)),length(Abd).*3,1)];
        temp_interval = [temp_interval;[An;Cn;En]];
        temp_sos = [temp_sos;[Abd;Cbd;Ebd]]; 
    end
end 

table_20m = table(temp_year,temp_park,temp_interval,temp_sos,...
    'VariableNames',{'Year','Park','Interval','SOS'});

interval = {'10m','30m','50m'};
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
year_str = string(years);

table_20m.Interval = categorical(table_20m.Interval,interval);
table_20m.Park = categorical(table_20m.Park,park_size_name);
table_20m.Year = categorical(table_20m.Year,year_str);

color_pick = [1, 256/2, 256];

for ps = 1:length(park_size_name)
    ind_temp = table_20m.Park == park_size_name{ps};
    figure;
    hold on
    h = boxchart(table_20m.Year(ind_temp),table_20m.SOS(ind_temp),'groupby',table_20m.Interval(ind_temp));
    for int = 1:3 %length(park_size_name)
        h(int).JitterOutliers = 'on';
        h(int).MarkerStyle = '.';
        h(int).BoxFaceColor = newmap(color_pick(int),:);
        h(int).MarkerColor = newmap(color_pick(int),:);
    end
    ylabel('SOS')
    title(sprintf('%s',park_size_name{ps}))
    ylim([40 140])

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Boxchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_10m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Boxchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_10m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))

    X = permute(squeeze(temp_mean(:,ps,1:3)),[2 1]);
    % Radar chart
    figure('Position',[600,200,350,350]);
    RC=radarChart(X,'Type','Patch');

    RC.RLim=sos_lim; %round([min(X,[],'all')-10,max(X,[],'all')+3]); % NIRv
    RC.RTick=sos_lim_tick; %round([min(X,[],'all')-10,mean(X,"all"),max(X,[],'all')+3]);
    % RC.RLim=round([90,115]); % NIRv
    % RC.RTick=round([90,105,115]);
    RC.PropName=string(years);
    RC.ClassName={'10m','30m','50m'};
    RC=RC.draw();
    % RC=RC.legend();
    % RC.setType('Patch');
    
    color_pick = [1, 256/2 , 256];
    for n=1:RC.ClassNum
        RC.setPatchN(n,'FaceColor',newmap(color_pick(n),:),'EdgeColor',newmap(color_pick(n),:));
    end
    
    % title(sprintf('%s',park_size_name{ps}))
    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    % ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/Radarchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_10m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/Radarchart_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_10m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
end

% Mean plot
for ps = 1:length(park_size_name)
    figure;
    hold on
    
    for int = 1:3 %length(park_size_name)
        h = plot(years,squeeze(temp_mean(:,ps,int)));
        h.Color= newmap(color_pick(int),:);
        h.LineWidth = 2;
    end
    
    title(sprintf('%s',park_size_name{ps}))
    ylim([90 110])

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    % ax.YTick ='';
    hold off

    savefig(gcf,sprintf('%s/all_years/Plot_mean_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_10m_width.fig',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))
    saveas(gcf,sprintf('%s/all_years/Plot_mean_T_Park_%s_X_Year_Y_%s_%s_%s_20m_int_10m_width.png',path_figure,park_size_name{ps}, phenology_names{phn},product,name_vis{vi}))

end

close all

%% Dscatter: LST and Edge SOS

park_categor = {'Pocket','Neighborhood', 'Community' , 'Regional'};
park_categor = {'Pkt','Nbrhd','Cmnty','Rgnl'};
% LST_temp = LST_diff300;
% LST = permute(LST_temp(:,19:23,phn),[2 1]);

temp_park = [];
temp_interval = [];
temp_sos = [];
temp_mean=[];
temp_year=[];
temp_lst = [];

for yr = 1:length(years)
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
       
        Abd = reshape(pheno_date_bd_int10(yr,ind_temp,1),[],1);
        Fbd = reshape(pheno_date_bd_int20(yr,ind_temp,1),[],1);
        Hbd = reshape(pheno_date_bd_int30(yr,ind_temp),[],1);
        Ibd = reshape(pheno_date_bd_int40(yr,ind_temp),[],1);
        Jbd = reshape(pheno_date_bd_int50(yr,ind_temp),[],1);
        
        
        An =repmat("10m",length(Abd),1);
        Fn =repmat("20m",length(Fbd),1);
        Hn =repmat("30m",length(Hbd),1);
        In =repmat("40m",length(Ibd),1);
        Jn =repmat("50m",length(Jbd),1);
        
        if np ==1
            adding=18;
        else
            adding =0;
        end
        LST_100 = reshape(permute(LST_diff100(ind_temp,yr+adding,phn),[2 1]),[],1);
        LST_200 = reshape(permute(LST_diff200(ind_temp,yr+adding,phn),[2 1]),[],1);
        LST_300 = reshape(permute(LST_diff300(ind_temp,yr+adding,phn),[2 1]),[],1);

        temp_mean(yr,ps,1) =  mean(Abd,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(Fbd,"all","omitmissing");
        temp_mean(yr,ps,3) =  mean(Hbd,"all","omitmissing");
        temp_mean(yr,ps,4) =  mean(Ibd,"all","omitmissing");
        temp_mean(yr,ps,5) =  mean(Jbd,"all","omitmissing");

        temp_park = [temp_park;repmat(string(park_size_name{ps}),length(Abd).*5,1)];
        temp_year = [temp_year;repmat(string(years(yr)),length(Abd).*5,1)];
        temp_interval = [temp_interval;[An;Fn;Hn;In;Jn]];
        temp_sos = [temp_sos;[Abd;Fbd;Hbd;Ibd;Jbd]]; 
        temp_lst = [temp_lst;repmat([LST_100,LST_200,LST_300],5,1)];

    end
end 

table_edges = table(temp_year,temp_park,temp_interval,temp_sos,temp_lst(:,1),temp_lst(:,2),temp_lst(:,3),...
    'VariableNames',{'Year','Park','Interval','SOS','LST100','LST200','LST300'});
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
table_edges.Park = categorical(table_edges.Park ,park_size_name);

rho_all=nan(length(years),1);
pval_all=nan(length(years),1);

int=30;
for yr = 1:length(years)
    if yr ~= 13
    ind_temp = table_edges.Year == sprintf('%d',years(yr)) & table_edges.Interval == sprintf('%dm',int);
    
    X = table_edges.LST200(ind_temp);
    Y = table_edges.SOS(ind_temp);
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];

    [rho_all(yr), pval_all(yr)] =corr(X,Y);

    end
end


X_corr_sig = rho_all;
X_corr_nonsig = rho_all;

X_corr_sig(pval_all>=0.05)=nan;
X_corr_nonsig(pval_all<0.05)=nan;

figure('Position',[600,200,500,200]);
hold on


h = plot(years,rho_all);
h.Color = 'k';

h = scatter(years,X_corr_sig);
h.Marker = "o";
h.MarkerFaceColor = 'k';
h.MarkerEdgeColor = 'k';

h = scatter(years,X_corr_nonsig);
h.Marker = "o";
h.MarkerFaceColor = 'none';
h.MarkerEdgeColor = 'k';

mdl = fitlm(years,rho_all);

hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
hr.LineWidth =2;
hr.LineStyle = '--';
hr.Color = 'r';

hr = yline(0);
hr.Color = 'k';
hr.LineWidth = 1;

ylim([-0.3 0.3])
ax = gca;
ax.YColor = 'k';

xlim([years(1) years(end)])
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;


savefig(gcf,sprintf('%s/all_years/Plot_X_year_Y_R_SOS_LST_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
saveas(gcf,sprintf('%s/all_years/Plot_X_year_Y_R_SOS_LST_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))

% Dscatter= All
for int = 10:10:50
    ind_temp = table_edges.Interval == sprintf('%dm',int);
    
    X = table_edges.LST200(ind_temp);
    Y = table_edges.SOS(ind_temp);
    Y(isnan(X))=nan;
    X(isnan(Y))=nan;
    
    X(isnan(X))=[];
    Y(isnan(Y))=[];
    
    figure('Position',[600,200,500,250]);
    hold on 
    h=dscatter(X,Y,'plottype','scatter');
    h.Colormap = newmap;
    
    xlabel('Annual mean LST')
    ylabel('SOS')
    xlim([10 35])
    ylim([30 140])
    
    mdl = fitlm(X,Y);
    hr = refline(mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
    hr.LineWidth =2;
    hr.Color = 'k';

    corr(X,Y)

    ax = gca;
    ax.FontName = fig_font;
    ax.FontWeight = fig_font_weight;
    ax.FontSize = fig_font_size;
    ax.Box = 'on';
    ax.YAxisLocation=fig_y_axis_loc;
    
    savefig(gcf,sprintf('%s/all_years/Dscatter_X_annual_LST_Y_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
    saveas(gcf,sprintf('%s/all_years/Dscatter_X_annual_LST_Y_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))
end

close all


%% Boxchart: LST and Edge SOS 
temp_park = [];
temp_interval = [];
temp_mean=[];
temp_year=[];
temp_lst = [];

for yr = 1:length(years)
    for ps = 1:length(park_size_name)
    
        if ps ==1
            ind_temp = ind_park_vs;        
        elseif ps ==2
            ind_temp = ind_park_s;
        elseif ps ==3
            ind_temp = ind_park_m;
        elseif ps ==4
            ind_temp = ind_park_l;
        elseif ps ==5
            ind_temp = ind_park_all;
        end
        
       if np ==1
            adding=18;
        else
            adding =0;
        end
        LST_100 = reshape(permute(LST_diff100(ind_temp,yr+adding,phn),[2 1]),[],1);
        LST_200 = reshape(permute(LST_diff200(ind_temp,yr+adding,phn),[2 1]),[],1);
        LST_300 = reshape(permute(LST_diff300(ind_temp,yr+adding,phn),[2 1]),[],1);

        LST_100n =repmat("100m",length(LST_100),1);
        LST_200n =repmat("200m",length(LST_200),1);
        LST_300n =repmat("300m",length(LST_300),1);


        temp_mean(yr,ps,1) =  mean(LST_100,"all","omitmissing");
        temp_mean(yr,ps,2) =  mean(LST_200,"all","omitmissing");
        temp_mean(yr,ps,3) =  mean(LST_300,"all","omitmissing");
        
        temp_lst = [temp_lst;[LST_100;LST_200;LST_300]];
        temp_interval = [temp_interval;[LST_100n;LST_200n;LST_300n]];
        temp_park = [temp_park;repmat(string(park_size_name{ps}),length([LST_100;LST_200;LST_300]),1)];
        temp_year = [temp_year;repmat(string(years(yr)),length([LST_100;LST_200;LST_300]),1)];
        
    end
end 

table_lst = table(temp_year,temp_park,temp_interval,temp_lst,...
    'VariableNames',{'Year','Park','Interval','LST',});
park_size_name = {'Pocket','Neighborhood','Community','Regional','All'};
year_str = string(years);

table_lst.Park = categorical(table_lst.Park ,park_size_name);
table_lst.Year = categorical(table_lst.Year ,year_str);


ind_temp = table_lst.Interval == "200m" & table_lst.Park ~= "All";

figure('Position',[600,200,500,250]);
hold on 
h=boxchart(table_lst.Park(ind_temp),table_lst.LST(ind_temp)); %,'GroupByColor',table_lst.Year(ind_temp),);
h(1).MarkerStyle = '.';
h(1).BoxFaceColor = newmap(1,:);
h(1).MarkerColor = newmap(1,:);


% xlabel('Annual mean LST')
ylabel('LST')
% xlim([10 35])
ylim([0 35])

ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
ax.Box = 'on';
ax.YAxisLocation=fig_y_axis_loc;
% ax.YTick ='';
hold off

savefig(gcf,sprintf('%s/all_years/Boxchart_X_annual_LST_Y_%s_%s_%s_%dm_int_edge.fig',path_figure, phenology_names{phn},product,name_vis{vi},int))
saveas(gcf,sprintf('%s/all_years/Boxchart_X_annual_LST_Y_%s_%s_%s_%dm_int_edge.png',path_figure, phenology_names{phn},product,name_vis{vi},int))



%
X = permute(squeeze(temp_mean(:,ps,1:3)),[2 1]);
% Radar chart
figure('Position',[600,200,350,350]);
RC=radarChart(X,'Type','Patch');
RC.RLim=round([min(X,[],'all')-10,max(X,[],'all')+3]); % NIRv
RC.RTick=round([min(X,[],'all')-10,mean(X,"all"),max(X,[],'all')+3]);% RC.RLim=round([90,115]); % NIRv

% RC.RLim=sos_lim; %round([min(X,[],'all')-10,max(X,[],'all')+3]); % NIRv
% RC.RTick=sos_lim_tick; %round([min(X,[],'all')-10,mean(X,"all"),max(X,[],'all')+3]);% RC.RLim=round([90,115]); % NIRv
% RC.RTick=round([90,105,115]);
RC.PropName=string(years);
RC.ClassName={'100m','200m','300m'};
RC=RC.draw();
% RC=RC.legend();
% RC.setType('Patch');

color_pick = [1, 256/2 , 256];
for n=1:RC.ClassNum
    RC.setPatchN(n,'FaceColor',newmap(color_pick(n),:),'EdgeColor',newmap(color_pick(n),:));
end

% title(sprintf('%s',park_size_name{ps}))
ax = gca;
ax.FontName = fig_font;
ax.FontWeight = fig_font_weight;
ax.FontSize = fig_font_size;
% ax.Box = 'on';
% ax.YTick ='';
hold off
%%



%
%% Boxchart Phenodate Change from Edge to Interior 
% figure;
% hold on
% boxchart([Abd,Fbd,Hbd,Ibd,Jbd])
% ylim([90 120])
% 
% temp(1) = mean(Abd,"all","omitmissing");
% temp(2) = mean(Bbd,"all","omitmissing");
% temp(3) = mean(Cbd,"all","omitmissing");
% temp(4) = mean(Dbd,"all","omitmissing");
% temp(5) = mean(Ebd,"all","omitmissing");
% 
% figure;
% plot(temp)


%% Histogram Phenodate Change from Edge to Interior 
% Abd(isnan(Abd))=[];
% Bbd(isnan(Bbd))=[];
% Cbd(isnan(Cbd))=[];
% Dbd(isnan(Dbd))=[];
% Ebd(isnan(Ebd))=[];
% 
% figure;
% hold on
% histogram(Abd,40,Normalization="probability")
% histogram(Bbd,40,Normalization="probability")
% histogram(Cbd,40,Normalization="probability")
% histogram(Dbd,40,Normalization="probability")
% histogram(Ebd,40,Normalization="probability")
% 
% 
% Fbd(isnan(Fbd))=[];
% Gbd(isnan(Gbd))=[];
% 
% temp = Abd-Dbd;
% 
% temp(isnan(temp))=[];
% 
% figure;
% hold on
% % histogram(Fbd,40,Normalization="probability")
% % histogram(Gbd,40,Normalization="probability")
% 
% hf = histogram(temp,'Normalization','probability');
% hf.NumBins=num_bin;
% hf.FaceColor =newmap(end/2,:);
% hf.EdgeColor = 'none';
% hf.FaceAlpha = 0.2;
% 
% h = xline(0);
% h.Color = 'k';
% h.LineWidth=2;
% 
% h = xline(mean(temp));
% h.Color = 'r';
% h.LineWidth=2;
% 
% pd = fitdist(temp,'Normal');
% % ci95 = paramci(pd);
% x_values = -40:1:40;
% y = pdf(pd,x_values);%./num_bin;
% h=plot(x_values,y);
% h.Color=newmap(end/2,:);
% h.LineWidth=2;

%%
% figure;
% hold on
% histfit(Adiff)
% % boxchart([Abd,Bbd,Cbd,Dbd,Ebd])
% % boxchart([Abd,Fbd,Hbd,Ibd,Jbd])
% % boxchart([Abd,Bbd,Cbd,Dbd,Ebd])
% boxchart([Ain,Bin,Cin,Din,Ein])
% 
% boxchart([Fbd,Gbd])
% 
% 
% figure;
% hold on
% h=histfit(Adiff);
% mean(Adiff,'all','omitmissing')
% 
% h=xline(0);
% h.LineWidth=2;
% h.Color='k';


