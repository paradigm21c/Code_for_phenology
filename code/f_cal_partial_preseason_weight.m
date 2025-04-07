function [rho, pval] = f_cal_partial_preseason_weight(vi_daymet, daymet_names, vi, phn, dayn, mkn, etype, w_s, w_e,park_veg_ratio)

% w_s=i;
% w_e=i+winsize(i)-1;


    temp_1 = vi_daymet{vi,1}; % (yr,ec,mk,phn)
    % temp_1(13,:,:,:)=NaT; % year 2012
    % temp_1(14,:,:,:)=NaT; % year 2013
    for td = 1:length(daymet_names)
        temp_2 = vi_daymet{vi,td+1}; % Tmax
        daymet_temp(:,:,:,td) = temp_2(:,:,:,phn); % (yr,ec,mk,phn)
    end
    % td_range = 1:length(daymet_names);
    % td_range = td_range(td_range ~= dayn);
    % td_range =[3,4,5];
    if dayn ==1 || dayn ==2
        td_range =[3,4];
    elseif dayn ==3
        td_range =[1,4];
    elseif dayn ==4
        td_range =[1,3];
    end

    temp_1.Format='DDD';
    temp_1 = str2double(string(temp_1));
    temp = [];
     % (pk, yr, ec)
    for i = 1:length(etype)
        temp(:,:,i) = permute(squeeze(temp_1(:, etype(i), :, phn)), [2, 1]) .* park_veg_ratio(:, i);  
    end

    % Sum ignoring NaNs
    temp =  sum(temp,3,'omitnan'); 
    temp(temp==0)=nan;
    vi_pheno_date = permute(temp,[2,1]); % (yr,ec,mk,phn)

    temp = [];
     % (pk, yr, ec)
    for i =1:length(etype)
        temp(:,:,i) = permute(squeeze(daymet_temp(:,etype(i),:,dayn)),[2, 1]).*park_veg_ratio(:,i);
    end
    % Sum ignoring NaNs
    temp =  sum(temp,3,'omitnan');
    temp(temp==0)=nan;
    Daymet_pheno_date =permute(temp,[2,1]);% (yr,ec,mk,td)

    rho = nan(mkn,1);
    pval = nan(mkn,1);
                  
    for mk = 1:mkn
        pheno_date = reshape(squeeze(vi_pheno_date(w_s:w_e,mk)),[],1);
        % pheno_date.Format='DDD';
        % pheno_date = str2double(string(pheno_date));
        daymet_main =  reshape(squeeze(Daymet_pheno_date(w_s:w_e,mk)),[],1);
        daymet_remain=[];
        
        for on = 1:length(td_range)
            temp = []; % (yr,ec,mk,td)
            for i =1:length(etype) 
                temp(:,:,i) = permute(squeeze(daymet_temp(w_s:w_e,etype(i),mk,td_range(on))),[2, 1]).*park_veg_ratio(mk,i);
            end
            temp =  sum(temp,3,'omitnan');
            % temp(temp==0)=nan;
            daymet_temp2 = reshape(squeeze(temp),[],1);
            daymet_remain = [daymet_remain,daymet_temp2];
        end

        idx = isnan(daymet_main);
        pheno_date(idx,:)=[];
        daymet_main(idx,:)=[];
        daymet_remain(idx,:)=[];

        if ~isempty(pheno_date) && sum(~isnan(pheno_date)) >= round((w_e-w_s)./2)  && ~isempty(daymet_main)&& sum(~isnan(daymet_main)) >= round((w_e-w_s)./2)
            [rho(mk), pval(mk)] = partialcorr(pheno_date, daymet_main, daymet_remain,'Rows','complete');
        end
    end
end