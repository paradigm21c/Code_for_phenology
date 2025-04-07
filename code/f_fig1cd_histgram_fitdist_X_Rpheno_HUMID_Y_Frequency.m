%%
product='SDC_full';

ecotype = {'Forested_Wetland','Maintained_Lawn_Shrub','Maritime_Forest',...
    'Other_Tree_Canopy','Upland_Forest','Upland_Grass_Shrub', 'combined_mask',...
    'Bare_Soil','Building','Freshwater_Aquatic_Vegetation','Freshwater_Wetland',...
    'Other_Paved_Surface','Other_Water','Road_Railroad','Saltwater_Aquatic_Vegetation',...
    'Tidal_Wetland','All_layer'};

phenology_names = {'SOS', 'MidGreenUp', 'Maturity', 'Peak', 'Senescence', 'MidGreenDown', 'Dormancy'};

daymet_names = {'Tmax','Tmin','Precipitation','Daily total radiation','Vapor pressure'};

path_mask = sprintf('D:/SetoLab/Phenology/mask/');
files_mask = dir(sprintf('%s/parks_arc_%s_DS/*.tif',path_mask,product));
mkn = length(files_mask);

vin = {'EVI2','NDVI','NIRv'};
etn={'Tree'};

path_cal = 'D:/SetoLab/Phenology/data_cal/Parks_polygon/ParkDS/HUMID';
path_figure = sprintf('D:/SetoLab/Phenology/figure/Preseason_%s_DS/all_years',product);

newmap = importdata('D:/SetoLab/Colors/batlow/batlow.mat');
fig_font = 'Arial';
fig_font_size = 14;

%
shrub =[2 6];
tree = [1 3 4 5];

winsize =12;
windows = [10, 12 ,15];

vi=3; etype = tree; phn=1; dayn=1; 

rho_stack_full=importdata(sprintf('%s/Rho_SDC_full_DS_mv_full.mat',path_cal));
pval_stack_full=importdata(sprintf('%s/Pval_SDC_full_DS_mv_full.mat',path_cal));
num_bin =40;

for mv =2%:length(windows)
    winsize=windows(mv);

    color_interval = floor(length(newmap)./(18-winsize));
    percent_save_neg= nan(size(1:18-winsize,2),1);
    percent_save_pos= nan(size(1:18-winsize,2),1);

    rho_stack=importdata(sprintf('%s/Rho_SDC_full_DS_mv%d.mat',path_cal,windows(mv)));
    pval_stack=importdata(sprintf('%s/Pval_SDC_full_DS_mv%d.mat',path_cal,windows(mv)));

    for dayn =1 %[1,2,4]
        for vi =3 %1:length(vin)
            for phn=1 %[1 4 5 7]
                rho_max_full=  rho_stack_full{dayn,vi,phn};
                pval_max_full = pval_stack_full{dayn,vi,phn};
    
                figure;
                hold on
                hf = histogram(rho_max_full(:,phn),'Normalization','probability');
                hf.NumBins=num_bin;
                hf.FaceColor =newmap(end/2,:);
                hf.EdgeColor = 'none';
                hf.FaceAlpha = 0.2;
    
                pd = fitdist(rho_max_full(:,phn),'Normal');
                % ci95 = paramci(pd);
                x_values = -1:0.01:1;
                y = pdf(pd,x_values); %./num_bin;
                y_scaled = y* hf.BinWidth; % Scale by the bin width

                h=plot(x_values,y_scaled);
                h.Color=newmap(end/2,:);
                h.LineWidth=2;
      
                for i = [1,18-winsize]
                    rho_max=  rho_stack{dayn,vi,phn,i};
                    pval_max = pval_stack{dayn,vi,phn,i};
                    
                    hf = histogram(rho_max(:,phn),'Normalization','probability');
                    hf.NumBins=num_bin;
                    if i==1
                        hf.FaceColor=newmap(1,:);
                    else
                        hf.FaceColor=newmap(1+color_interval*(i-1),:);
                    end

                    % hf.FaceColor = newmap(1+color_interval*(i-1),:);
                    hf.EdgeColor = 'none';
                    hf.FaceAlpha = 0.2;
    
                    pd = fitdist(rho_max(:,phn),'Normal');
                    
                    x_values = -1:0.01:1;
                    y = pdf(pd,x_values); %./num_bin;
                    y_scaled = y* hf.BinWidth; % Scale by the bin width

                    h=plot(x_values,y_scaled);
                    if i==1
                        h.Color=newmap(1,:);
                    else
                        h.Color=newmap(1+color_interval*(i-1),:);
                    end
                    h.LineWidth=2;
                
                    percent_pos=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) > 0)./sum(~isnan(rho_max(:,phn)));
                    percent_neg=sum(pval_max(:,phn) < 0.05 & rho_max(:,phn) < 0)./sum(~isnan(rho_max(:,phn)));
                
                    percent_save_pos(i) = percent_pos;
                    percent_save_neg(i) = percent_neg;
                
                    % annotation(gcf,'textbox',...
                    % [0.655357142857142 0.534238094284421 0.166785717836448 0.0638095247631981],...
                    % 'String',{sprintf('%.3g (%.3g)',round(percent_pos,3),round(percent_neg,3))});
                
                    ylabel('Frequency (%)');
                    xlabel(sprintf('R_{(%s - %s)}',phenology_names{phn}, daymet_names{dayn}))
                
                    hy=xline(0);
                    hy.LineWidth =2;
                    hy.Color = 'k';
                
                    ax = gca;
                    ax.FontName = fig_font;
                    ax.FontWeight = 'bold';
                    ax.FontSize = fig_font_size;
                    ax.Box = 'on';
                
                    xlim([-1 1])
                    % ylim([0 35])
                
                    % rho_stack(:,:,i) = rho_max;

                    % Change y-axis to percentage
                    % Get current y-ticks
                    yticks = get(gca, 'YTick');
                    
                    % Multiply y-ticks by 100 to convert to percentage
                    new_yticks = yticks * 100;

                    % Update y-tick labels with percentage values
                    set(gca, 'YTickLabel', strcat(num2str(new_yticks')));
                end
    
                % savefig(gcf,sprintf('%s/Histfit_X_%s_%s_tree_Fusion_%s_mv_%d.fig',path_figure, phenology_names{phn},daymet_names{dayn},vin{vi},winsize))
                % saveas(gcf,sprintf('%s/Histfit_X_%s_%s_tree_Fusion_%s_mv_%d.png',path_figure, phenology_names{phn},daymet_names{dayn},vin{vi},winsize))
    
                % close all
            end
        end
    end
end