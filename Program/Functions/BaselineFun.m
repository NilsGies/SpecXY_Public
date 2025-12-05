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
    % plot(wavenumber_reshaped,intensity_reshaped)

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
        %figure;plot(resultpartout)

    elseif strcmp(method,'interp1') %% ATH add possibility for UI points and find minimum



%         dataPoints=[intensity_reshaped(1) intensity_reshaped(end)];
%         wavePoints=[wavenumber_reshaped(1) wavenumber_reshaped(end)];
%         baseline=interp1(wavePoints,dataPoints,wavenumber_reshaped,interp_method);
%         
        
        mask=wavenumber_reshaped>3800 & wavenumber_reshaped<5800;
       





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

        if sum(k1)>sum(k2)%trapz(y_signal(k(1:locs)))<trapz(flipud(y_signal(k(locs:end)))) %trapz(y_signal(k(1:locs)))>trapz(y_signal(k(locs:end)))
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
    %    resultpartout = intensity_reshaped-polyval(polyfit(wavenumber_reshaped,intensity_reshaped,value2),wavenumber_reshaped);
% 
%         opol = value2;
%         t=wavenumber_reshaped(wavenumber_reshaped>4000)
%         ecgnl=intensity_reshaped(wavenumber_reshaped>4000);
% 
%         [p,s,mu] = polyfit(t,ecgnl,opol);
%         f_y = polyval(p,wavenumber_reshaped,[],mu);
% 
%         resultpartout = intensity_reshaped - f_y;


% mask=wavenumber_reshaped>4000 ;%| (wavenumber_reshaped<2700 & wavenumber_reshaped>2300);
% cfit(wavenumber_reshaped(mask,:),intensity_reshaped(mask,:),'poly2')
% 
%         yhat = cf(wavenumber_reshaped(wavenumber_reshaped>3000&wavenumber_reshaped<3720 ))
% 
% 
%   cf = fit(movmean(wavenumber_reshaped(mask,:),4),intensity_reshaped(mask,:),'smoothingspline');
% 
%        prediction_range=wavenumber_reshaped(wavenumber_reshaped>3000 )
%         yhat = cf(prediction_range)
% 
%         plot(prediction_range,yhat,'-x')
%         hold on
%         plot(wavenumber_reshaped,intensity_reshaped,'-')
% 
%    %     Tbl = table(wavenumber_reshaped,intensity_reshaped);
% % 
%         mask=wavenumber_reshaped>4000 ;%| (wavenumber_reshaped<2700 & wavenumber_reshaped>2300);
%         Tbl_=Tbl(mask,:);%={0};
%         Mdl1 = fitrensemble(Tbl_,'intensity_reshaped','Method','Bag'); %,'OptimizeHyperparameters',{'NumLearningCycles','MaxNumSplits'});
% 
%         disp('done');
% 
% 
%         ExtrapolationX = (3000:3800)';
% ExtrapolationY = feval(Mdl1, ExtrapolationX);
% 
% 
%         resultpartout = predict(Mdl1,wavenumber_reshaped);
% 
% 

        
        %resultpartout = intensity_reshaped-smoothdata(intensity_reshaped,1,'movmedian',value);

        %     [refractive_index, corrected_spectrum] = thin_film_compensation_auto(wavenumber_reshaped, intensity_reshaped, 165);
        %    resultpartout= mean([real(corrected_spectrum) intensity_reshaped],2);
        %%
%         Y_full = fft(intensity_reshaped);
%         ir=intensity_reshaped(wavenumber_reshaped>2415 &wavenumber_reshaped<2825);
%         Y = fft(ir);
% 
%         f = 1./1/4*linspace(0,1,length(wavenumber_reshaped));
% 
%         power = abs(Y).^2/ length(ir);
%         [~,idx]=sort(power,'descend');
% 
%         %frequencies = f(idx)'; % Extract the frequencies with high power
%       % selected_freqs = f(idx([1:value ceil(numel(f)/value2):numel(idx)]))'; % Extract the frequencies with high power
%      %   selected_freqs = f(idx([1:value numel(f)-value2:numel(idx)]))'; % Extract the frequencies with high power
%      %  selected_freqs = f(idx([1:value ]))'; % Extract the frequencies with high power
%       numel(f)
%      selected_freqs = f(idx(1:10))'; % Extract the frequencies with high power
% 
%         % Remove selected frequencies from spectrum
%         mask = ones(size(f));
%         for i = 1:length(selected_freqs)
%             [~, idx] = min(abs(f - selected_freqs(i)));
%             mask(idx) = 0;
%         end
%         Y_filtered = Y_full .* mask';
%         resultpartout=real(ifft(Y_filtered));
        %%
    elseif strcmp(method,'sgoley')
        resultpartout = intensity_reshaped-smoothdata(intensity_reshaped,1,'sgolay',value,'Degree',value2);
    elseif strcmp(method,'linCorrect2')
        resultpartout=(intensity_reshaped)-min(intensity_reshaped);
    elseif strcmp(method,'linCorrect')
        resultpartout=linCorrect2(intensity_reshaped);
    elseif strcmp(method,'linear')
        signal=lin_BL_corr([wavenumber_reshaped,intensity_reshaped]);
        resultpartout=signal(:,2:end);
    elseif strcmp(method,'linear2')
        resultpartout=linCorrect(intensity_reshaped);
        %         a=(intensity_reshaped(1,:));
        %         b=(intensity_reshaped(end,:));
        %         steps = size(intensity_reshaped,1);                      %// number of steps
        %         correction_line = bsxfun(@plus,((b(:)-a(:))./(steps-1))*[0:steps-1],a(:));
        %         resultpartout =         intensity_reshaped-correction_line';
    elseif strcmp(method,'Auto Raman Baseline')
        [Base, resultpartout]=baseline_raman(wavenumber_reshaped);
    end
    corrected_intensity=[corrected_intensity;resultpartout];
end
baseline=y_signal-corrected_intensity;
end


% Rubberband
% Rubberband boundary (slow but better)
% msbackadj
% interp1
% Polynominal
% linear
% linCorrect2
% linCorrect
% linear2
% Auto Raman Baseline


%% archive
% 
%  function [corrected_signal]=norm_signal_wavenumber(~,signal)
%             if findpeaks(signal(:,1))>2
%                 [ ~, peak]=findpeaks(signal(:,1));
%                 signal=signal(1:peak(1),:);
%             end
% 
%             x=signal(:,1);
%             region_min=floor(min(x));
%             region_max=ceil(max(x));
% 
%             xq=linspace(region_min,region_max,region_max-region_min+1);
% 
% 
%             corrected_signal=zeros(region_max-region_min+1,size(signal,2));
%             corrected_signal(:,1)= xq;
%             for k=2:size(signal,2)
%                 v=signal(:,k);
%                 vq = interp1(x,v,xq);
%                 vq(isnan(vq))=0;
%                 corrected_signal(:,k)=vq;
%             end
%         end
% 
%         function [signal_corrected]=lin_BG_corr(app,signal)
%             signal=norm_signal_wavenumber(app,signal);
% 
%             if findpeaks(signal(:,1))>2
%                 [ ~, peak]=findpeaks(signal(:,1));
%                 signal=signal(1:peak(1),:);
%             end
% 
% 
%             min_signal_y=round(signal(1,2));
%             max_signal_y=round(signal(end,2));
%             min_signal_x=round(signal(1,1));
%             max_signal_x=round(signal(end,1));
% 
% 
% 
%             signal_corrected=zeros(max_signal_x-min_signal_x+1,size(signal,2));
%             signal_corrected(:,1)=signal(:,1);
%             for k=2:size(signal,2)
% 
%                 correction_line=linspace(min_signal_y,max_signal_y,(max_signal_x-min_signal_x+1));
%                 signal_corrected(:,k)=signal(:,k)-correction_line';
%                 signal_corrected((signal_corrected(:,k)<0),k)=0;
%             end
% 
%         end
% 
%         function[signal_corrected]=spline_BG_cor(app,signal)
%             signal_n=norm_signal_wavenumber(app,signal);
% 
%             min_signal_x=round(signal_n(1,1));
%             max_signal_x=round(signal_n(end,1));
%             signal_corrected=zeros(max_signal_x-min_signal_x,size(signal_n,2));
%             % size_cs1=size(signal_corrected)
%             signal_corrected(:,1)=signal_n(:,1);
%             %  size_cs2=size(signal_n(:,1))
% 
%             correction=signal_corrected;
%             xq=signal_n(:,1);
% 
%             for k=2:size(signal_n,2)
% 
%                 vq=signal_n(:,k);
%                 wf = @(mz) app.interval + 0.01 .* mz;
%                 signal_BL_corr=msbackadj(xq,vq, 'StepSize',wf,'ShowPlot', false);
% 
%                 vq=signal_BL_corr;
%                 %vq=vq-min(vq);
% 
%                 signal_corrected(:,k)=vq;
%                 %signal_corrected((signal_corrected(:,k)<0),k)=0
%                 correction(:,k)=signal_n(:,k)-vq;
%             end
%         end
