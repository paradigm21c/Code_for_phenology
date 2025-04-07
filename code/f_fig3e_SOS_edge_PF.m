%
product = 'PF';
newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;
val_alpha=0.4;
fig_font_weight = "normal"; % "bold"
fig_y_axis_loc = 'right';

years =2018:2022;
parks_files = [100,200,300];
% Path
path_all = 'D:/SetoLab/Phenology/';
path_map = sprintf('%s/map/PF/',path_all);
path_buff = sprintf('%s/mask/PF_buffer/',path_all);
path_park = sprintf('%s/mask/parks_arc_PF_DS/',path_all);
path_ecm = sprintf('%s/mask/PF_temp/',path_all);
path_figure = sprintf('%s/figure/Preseason_%s_DS/all_years/HUMID',path_all,product);


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

% Park names
park_name = nan(length(file_park),1);
for i=1:length(file_park)
    temp = split(file_park(i).name,'_');
    park_name(i) = str2double(temp{2});
end
park_name = sort(park_name);



phn=1; vi =2;

if vi ==1
    load('file_EVI2_pheno_date_inNbound_cpu_PF.mat') % (yr,pk) 
elseif vi ==2
    load('file_NIRv_pheno_date_inNbound_cpu_PF.mat') % (yr,pk)
end

if vi==2
    sos_lim = [95 115];
    sos_lim_tick = [95, 105,115];
elseif vi ==1
    sos_lim = [90 110];
    sos_lim_tick = [90, 100,110];
end

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
    X = permute(squeeze(temp_mean(:,ps,1:3)),[2 1]);
    % Radar chart
    figure('Position',[600,200,350,350]);
    set(gcf,'renderer','Painters')
    RC=radarChart(X,'Type','Patch');

    RC.RLim=sos_lim; %round([min(X,[],'all')-10,max(X,[],'all')+3]); % NIRv
    RC.RTick=sos_lim_tick; %round([min(X,[],'all')-10,mean(X,"all"),max(X,[],'all')+3]);
    % RC.RLim=round([90,115]); % NIRv
    % RC.RTick=round([90,105,115]);
    RC.PropName=string(years);
    RC.ClassName={'0-10m','20-30m','40-50m'};
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

