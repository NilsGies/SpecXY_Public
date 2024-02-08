function [corrected_intensity,baseline] = BaselineFun(~,x_signal,y_signal,signal_split_n,value,value2,method,interp_method)


if signal_split_n>1
    vq_in = reshape([y_signal;nan(mod(-numel(y_signal),signal_split_n),1)],[],signal_split_n);
    wn_in = reshape([x_signal;nan(mod(-numel(x_signal),signal_split_n),1)],[],signal_split_n);
    corrected_intensity=[];
else
    vq_in = y_signal;
    wn_in = x_signal;
    corrected_intensity=[];
end

% figure
% hold on
for m=1:signal_split_n

    intensity_reshaped=vq_in(:,m);
    intensity_reshaped=intensity_reshaped(~isnan(intensity_reshaped),1);
    wavenumber_reshaped=wn_in(~isnan(intensity_reshaped),m);


    if strcmp(method,'msbackadj')
        % interval=round(length(intensity_reshaped)/50);
        if not(exist("value2","var"))
            value2=100;
        end
        interval=value2;
        wf = @(mz) interval + value .* mz;
        if wavenumber_reshaped(1)<wavenumber_reshaped(end)
            resultpartout=msbackadj(wavenumber_reshaped,intensity_reshaped, 'StepSize',wf,'ShowPlot', false);
        else
            resultpartout=flipud(msbackadj(flipud(wavenumber_reshaped),flipud(intensity_reshaped), 'StepSize',wf,'ShowPlot', false));
        end


    elseif strcmp(method,'interp1') %% ATH add possibility for UI points and find minimum
        %         dataPoints=[intensity_reshaped(1) intensity_reshaped(end)];
        %         wavePoints=[wavenumber_reshaped(1) wavenumber_reshaped(end)];
        %         baseline=interp1(wavePoints,dataPoints,wavenumber_reshaped,interp_method);
        %        mask=wavenumber_reshaped>3800 & wavenumber_reshaped<5800;
        mask=true(size(wavenumber_reshaped));
        baseline=interp1(wavenumber_reshaped(mask),intensity_reshaped(mask),wavenumber_reshaped,interp_method,'extrap');

        resultpartout=intensity_reshaped-baseline;

    elseif strcmp(method,'Rubberband') || strcmp(method,'Rubberband boundary (slow but better)')


        intensity_reshaped(isnan(intensity_reshaped))=0;
        if unique(intensity_reshaped)==0
            resultpartout=zeros(size(intensity_reshaped));
            corrected_intensity=[corrected_intensity;resultpartout];
            continue
        end


        if strcmp(method,'Rubberband boundary (slow but better)')
            k = boundary(wavenumber_reshaped,intensity_reshaped,value); %slow
        else
            k = convhull(wavenumber_reshaped,intensity_reshaped);
        end
        % but better

        [~, locs]=findpeaks(k);
        longest_data=max([length(x_signal(k(1:locs))) length(x_signal(k(locs:end)))]);

        k1=interp1(x_signal(k(1:locs)),y_signal(k(1:locs)),linspace(x_signal(1),x_signal(end),longest_data),"linear");
        k2=interp1(x_signal(k(locs:end)),y_signal(k(locs:end)),linspace(x_signal(1),x_signal(end),longest_data),"linear");

        if sum(k1)>sum(k2)%trapz(intensity(k(1:locs)))<trapz(flipud(intensity(k(locs:end)))) %trapz(intensity(k(1:locs)))>trapz(intensity(k(locs:end)))
            k=(k(locs:end));else;k=(k(1:locs));
        end
        resultpartout=intensity_reshaped-(interp1(wavenumber_reshaped(k),intensity_reshaped(k),wavenumber_reshaped,interp_method));

    elseif strcmp(method,'Rubberband NG')

        resultpartout=rubberband_NG(wavenumber_reshaped,intensity_reshaped,value2,value,interp_method);

    elseif strcmp(method,'remcon')
        resultpartout=remCont2(wavenumber_reshaped,intensity_reshaped);
    elseif strcmp(method,'Polynominal')

        [~, sorted_indices] = sort(intensity_reshaped);
baseline_indices = sorted_indices(1:numel(intensity_reshaped));

% Fit a polynomial to the baseline points
baseline_x = (1:numel(intensity_reshaped))';
baseline_y = intensity_reshaped(baseline_indices);
p = polyfit(baseline_x, baseline_y, value2);

% Subtract the polynomial fit from the original signal
baseline = polyval(p, 1:length(intensity_reshaped))';
resultpartout = intensity_reshaped - baseline;

    elseif strcmp(method,'sgoley')
        resultpartout = intensity_reshaped-smoothdata(intensity_reshaped,1,'sgolay',value,'Degree',value2);
    elseif strcmp(method,'linCorrect2')
        resultpartout=linCorrect2(intensity_reshaped);
    elseif strcmp(method,'linCorrect')
        resultpartout=linCorrect2(intensity_reshaped);
    elseif strcmp(method,'linear')
        signal=lin_BL_corr([wavenumber_reshaped,intensity_reshaped]);
        resultpartout=signal(:,2:end);
    elseif strcmp(method,'linear2')
        resultpartout=linCorrect(intensity_reshaped);
   
    elseif strcmp(method,'Auto Raman Baseline')
        [Base, resultpartout]=baseline_raman(wavenumber_reshaped);
    end
    corrected_intensity=[corrected_intensity;resultpartout];
end
baseline=y_signal-corrected_intensity;
end

