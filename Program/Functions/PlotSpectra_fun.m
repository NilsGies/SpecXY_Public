function PlotSpectra_fun(axis2plot,data,settings)

if settings(1).reset==true
    try
        cla(axis2plot)
    catch

    end
end

if  isfield(settings,'XLabel')
    axis2plot.XLabel.String=settings(1).XLabel;
else
    axis2plot.XLabel.String='Wavenumbers [cm^{-1}]';
end
if isfield(settings,'YLabel')
    axis2plot.YLabel.String=settings(1).YLabel;
else
    axis2plot.YLabel.String='Absorbance [cm^{-1}]';
end

if isfield(settings,'ZLabel') && isfield(axis2plot,'ZLabel')
    axis2plot.ZLabel.String=settings(1).ZLabel;
elseif isfield(axis2plot,'ZLabel')
    axis2plot.ZLabel.String='';

end

try
    if  isfield(settings,'x_limits') && not(isempty(settings(1).x_limits)) &&not(isequal(settings(1).x_limits(1),settings(1).x_limits(2)))
        xlim(axis2plot,sort(settings(1).x_limits))
    else
        xlim(axis2plot,'auto')
    end
    xlim(axis2plot,settings(1).x_limits)
    if isfield(settings,'y_limits') && not(isempty(settings(1).y_limits)) &&not(isequal(settings(1).y_limits(1),settings(1).y_limits(2)))
        ylim(axis2plot,sort(settings(1).y_limits))
    else
        ylim(axis2plot,'auto')
    end
catch
    axis2plot=axis2plot.Children;
    hold on
    % xlim(axis2plot.Children,settings(1).x_limits)
    if isfield(settings,'x_limits') && not(isempty(settings(1).x_limits)) &&not(isequal(settings(1).x_limits(1),settings(1).x_limits(2)))
        xlim(axis2plot,settings(1).x_limits)
    else
        xlim(axis2plot,'auto')
    end

    if isfield(settings,'y_limits') && not(isempty(settings(1).y_limits)) &&not(isequal(settings(1).y_limits(1),settings(1).y_limits(2)))
        ylim(axis2plot,settings(1).y_limits)
    else
        ylim(axis2plot,'auto')
    end

end
if not(isfield(settings(1),'Annotation_interpreter'))
    settings(1).Annotation_interpreter='none';
end
if strcmp(settings(1).type,'Default') || strcmp(settings(1).type,'3D')
    offset=zeros(1,numel(data));
elseif strcmp(settings(1).type,'Stacked')
    for n=1:numel(data)
        y1=data(n).y_signal;
        x1=data(n).x_signal;
        y1=y1(x1>min(settings(1).x_limits) & x1 < max(settings(1).x_limits),:);
        offset(n)=max(mean(y1,2,"omitnan"),[],"omitnan"); %#ok
    end

    offset=[0 cumsum(offset)];
    if isfield(settings(1),'y_offset')
        y_offset=settings(1).y_offset;

    else
        y_offset=1.05;
    end

    %  offset(2:end)+([1:numel(offset(2:end))]* repmat(offset(2)*(y_offset-1),size( offset(2:end))))
    %   offset(2:end)=offset(2:end)+(1:numel(offset(2:end))*repmat(offset(2)*(y_offset-1),size( offset(2:end))));
    offset(2:end)=offset(2:end)+((1:numel(offset(2:end))).*(offset(2)*(y_offset-1))); %offset=offset*y_offset;

end
if strcmp(settings(1).type,'Stacked') || strcmp(settings(1).type,'Default')

    yticks(axis2plot,'auto')
    yticklabels(axis2plot,'auto');
    axis2plot.TickLabelInterpreter='tex';
    
    view(axis2plot,0,90);



    baseline=0;

    for n=1:numel(data)
        if isfield(data,'baseline') && not(isempty(data(n).baseline)) && size(data(n).baseline,1)==size(data(n).x_signal,1)
            if startsWith(settings(n).LineStyle.Baseline{:},'A')
                LineStyle=erase(settings(n).LineStyle.Baseline{:},'A');
                y_bl=data(n).baseline;
                y_bl(end)=min([data(n).baseline(:); 0]);
                y_bl(1)=min([data(n).baseline(:); 0]);
                fill(axis2plot,data(n).x_signal,y_bl+offset(n),settings(n).Color.Baseline,'LineWidth',settings(n).LineWidth.Baseline,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Baseline,'FaceColor',settings(n).Color.Baseline)
            else
                plot(axis2plot,data(n).x_signal,data(n).baseline+offset(n),LineWidth=settings(n).LineWidth.Baseline,Color=settings(n).Color.Baseline,LineStyle=settings(n).LineStyle.Baseline);%# ok

            end
        end

        if settings(n).PlotSpectra==true
            if startsWith(settings(n).LineStyle.Spectra{:},'A')
                LineStyle=erase(settings(n).LineStyle.Spectra{:},'A');
                
                fill(axis2plot,data(n).x_signal,data(n).y_signal+offset(n),settings(n).Color.Spectra,'LineWidth',settings(n).LineWidth.Spectra,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Spectra,'FaceColor',settings(n).Color.Spectra)
            else
                %plt(n)=
                plot(axis2plot,data(n).x_signal,data(n).y_signal+offset(n),LineWidth=settings(n).LineWidth.Spectra,Color=settings(n).Color.Spectra,LineStyle=settings(n).LineStyle.Spectra);
            end
        end

        if settings(n).PlotRaw==true && isfield(data,'y_signal_raw') && not(isempty(data(n).y_signal_raw)) && isfield(data,'x_signal_raw') && not(isempty(data(n).x_signal_raw))
            if startsWith(settings(n).LineStyle.Raw{:},'A')
                LineStyle=erase(settings(n).LineStyle.Raw{:},'A');
                fill(axis2plot,data(n).x_signal_raw,data(n).y_signal_raw+offset(n),settings(n).Color.Raw,'LineWidth',settings(n).LineWidth.Raw,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Raw,'FaceColor',settings(n).Color.Raw)
            else
                plot(axis2plot,data(n).x_signal_raw,data(n).y_signal_raw+offset(n),LineWidth=settings(n).LineWidth.Raw,Color=settings(n).Color.Raw,LineStyle=settings(n).LineStyle.Raw);%# ok
            end
        end

        if  settings(n).PlotMean==true
            if startsWith(settings(n).LineStyle.Mean{:},'A')
                LineStyle=erase(settings(n).LineStyle.Mean{:},'A');
                fill(axis2plot,data(n).x_signal,mean(data(n).y_signal,2,'omitnan')+offset(n),settings(n).Color.Mean,'LineWidth',settings(n).LineWidth.Mean,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Mean,'FaceColor',settings(n).Color.Mean)
            else
                plt(n)=plot(axis2plot,data(n).x_signal,mean(data(n).y_signal,2,'omitnan')+offset(n),LineWidth=settings(n).LineWidth.Mean,Color=settings(n).Color.Mean,LineStyle=settings(n).LineStyle.Mean);
            end
        end

        if  settings(n).PlotSum==true
        
            if startsWith(settings(n).LineStyle.Sum{:},'A')
                LineStyle=erase(settings(n).LineStyle.Sum{:},'A');
                fill(axis2plot,data(n).x_signal,sum(data(n).y_signal,2,'omitnan')+offset(n),settings(n).Color.Sum,'LineWidth',settings(n).LineWidth.Sum,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Sum,'FaceColor',settings(n).Color.Sum)
            else
            plt(n)=plot(axis2plot,data(n).x_signal,sum(data(n).y_signal,2,'omitnan')+offset(n),LineWidth=settings(n).LineWidth.Sum,Color=settings(n).Color.Sum,LineStyle=settings(n).LineStyle.Sum);
            end

        end


        if isfield(data,'DeconData') && not(isempty(data(n).DeconData)) && isfield(data(n).DeconData,'DeconResults') && not(isempty(data(n).DeconData(1).DeconResults)) && not(isempty(data(n).DeconData(1).DeconResults)) && (settings(n).PlotDeconDataSum ==true || any(settings(n).PlotDeconDataSingle==true))
            DeconID=data(n).DeconData.ActiveDeconGroup;
            y_calc=DeconPlotFun(data(n).DeconData.PeakFit(DeconID).peakshape,data(n).DeconData.PeakFit(DeconID).signal(:,1)',data(n).DeconData.DeconResults(DeconID).pos',data(n).DeconData.DeconResults(DeconID).width',data(n).DeconData.DeconResults(DeconID).hight',data(n).DeconData.PeakFit(DeconID).extra',data(n).DeconData.PeakFit(DeconID).DELTA');

            if ~isequal(baseline,0) || isequal(size(baseline),data(n).DeconData.PeakFit(DeconID).signal(:,1))
                baseline_dc=baseline(data(n).x_signal>=min(data(n).DeconData.PeakFit(DeconID).signal(:,1)) & data(n).x_signal<=(max(data(n).DeconData.PeakFit(DeconID).signal(:,1))),:);
            else
                baseline_dc=baseline;
            end

            for m=1:size(data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color,1)
                if data(n).DeconData.PlotSettings(DeconID).DeconSingle.Active(m)==true
                    if not(isempty(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m})) && startsWith(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m},'A')
                        LineStyle=erase(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m},'A');
                        X=[data(n).DeconData.PeakFit(DeconID).signal(:,1) data(n).DeconData.PeakFit(DeconID).signal(:,1) ];
                        Y=[offset(n)+baseline_dc.*ones(size(y_calc(m,:)))' y_calc(m,:)'+offset(n)+baseline_dc];
                        patch(axis2plot,X,Y,data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color(m,:),'LineWidth',data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineWidth(m),'LineStyle',LineStyle,'FaceAlpha',data(n).DeconData.PlotSettings(DeconID).DeconSingle.FaceAlpha(m))
                    end
                end

            end
            %do i need the next plot decon input?
            if settings(n).PlotSpectra==true
                plt(n)=plot(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),data(n).DeconData.PeakFit(DeconID).signal(:,2:end)+offset(n)+baseline_dc,'LineWidth',data(n).DeconData.PlotSettings(DeconID).Signal.LineWidth,'Color',data(n).DeconData.PlotSettings(DeconID).Signal.Color,'LineStyle',erase(data(n).DeconData.PlotSettings(DeconID).Signal.LineStyle,'A'));
            end
            if    settings(n).PlotDeconDataSum ==true
                try
                    plot(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),sum(y_calc(logical(data(n).DeconData.PlotSettings(DeconID).DeconSingle.Active),:))'+offset(n)+baseline_dc,'LineWidth',settings.LineWidth.DeconSum,'Color',settings.Color.DeconSum,'LineStyle',settings.LineStyle.DeconSum)

                catch
                    plot(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),sum(y_calc(logical(data(n).DeconData.PlotSettings(DeconID).DeconSingle.Active),:))'+offset(n)+baseline_dc,'LineWidth',data(n).DeconData.PlotSettings(DeconID).DeconTotal.LineWidth,'Color',data(n).DeconData.PlotSettings(DeconID).DeconTotal.Color,'LineStyle',data(n).DeconData.PlotSettings(DeconID).DeconTotal.LineStyle)

                end

            end

            for m=1:size(data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color,1)
                if data(n).DeconData.PlotSettings(DeconID).DeconSingle.Active(m)==true
                    LineStyle=erase(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m},'A');
                    plot(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),y_calc(m,:)'+offset(n)+baseline_dc,'LineWidth',data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineWidth(m),'Color',data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color(m,:),'LineStyle',LineStyle)
                end
            end
        end
        %%
        if settings(n).PlotPeaks==true
            if isfield(data,'DeconData') && not(isempty(data(n).DeconData)) && isfield(data(n).DeconData,'DeconResults') && not(isempty(data(n).DeconData(1).DeconResults)) && not(isempty(data(n).DeconData(1).DeconResults)) && (settings(n).PlotDeconDataSum ==true || any(settings(n).PlotDeconDataSingle==true))
                DeconID=data(n).DeconData.ActiveDeconGroup;


                label=join([ string(round((data(n).DeconData.DeconResults(DeconID).pos)')) data(n).DeconData.DeconInput(DeconID).Peaks.PeakDB.PeakNames]);
                linecolor=data(n).DeconData.PlotSettings(DeconID).PeakPos.Color;
                %old      line_settings=data(n).DeconData.PlotSettings(DeconID).PeakPos.LineType;
                line_settings=settings(n).LineStyle.Peaks;


                y_pks=data(n).y_signal+offset(n);
                y_pks=y_pks(data(n).x_signal>min(settings(1).x_limits) & data(n).x_signal < max(settings(1).x_limits),:);
                x_pks=data(n).x_signal(data(n).x_signal>min(settings(1).x_limits) & data(n).x_signal < max(settings(1).x_limits),:);

                PlotPeaksFun(axis2plot,x_pks,y_pks,data(n).DeconData.DeconResults(DeconID).pos,line_settings,label,linecolor)
            end
        end

        if isfield(data,'baseline') && not(isempty(data(n).baseline)) && size(data(n).baseline,1)==size(data(n).x_signal,1)
            plot(axis2plot,data(n).x_signal,data(n).baseline+offset(n),'--k',LineWidth=1);
        end

        %         if settings(1).PlotSpecNames==true && isfield(data(n),'Annotation')
        %             if n==1
        %                 t=text(axis2plot,min(settings(1).x_limits),offset(n),data(n).SampleName{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','none');
        %             elseif n==numel(data)
        %                 text(axis2plot,min(settings(1).x_limits),offset(n),data(n).SampleName{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','none')
        %
        %                 ylim(axis2plot,[t.Extent(2)*1.5 axis2plot.YLim(2)])
        %             else
        %                 text(axis2plot,min(settings(1).x_limits),offset(n),data(n).SampleName{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','none')
        %             end
        %         end
        if settings(1).PlotSpecNames==true && isfield(data(n),'Annotation')

            if n==1
                t(n)=text(axis2plot,min(settings(1).x_limits),offset(n),data(n).Annotation{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter',settings(1).Annotation_interpreter);
            elseif n==numel(data)

                t(n)=text(axis2plot,min(settings(1).x_limits),offset(n),data(n).Annotation{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter',settings(1).Annotation_interpreter);
                try
                    ylim(axis2plot,sort([t(1).Extent(2)*1.5 axis2plot.YLim(2)]))
                catch
                    t(n)=text(axis2plot,min(settings(1).x_limits),offset(1),data(1).Annotation{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter',settings(1).Annotation_interpreter);
                    ylim(axis2plot,sort([t(1).Extent(2)*1.5 axis2plot.YLim(2)]))

                end
            else
                t(n)=text(axis2plot,min(settings(1).x_limits),offset(n),data(n).Annotation{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter',settings(1).Annotation_interpreter);
            end
            if isfield(settings(1),'AnnotationFontSize')
                t(n).FontSize=settings(1).AnnotationFontSize;
            end


        end


    end
    if settings(1).PlotLegend==true
        legend(axis2plot,flip(plt),flip(data.SampleName),'Location','eastoutside',Interpreter='none');
    end


elseif strcmp(settings(1).type,'3D')

    [az,el] = view(axis2plot);
    if az==0 && el==90
        %         az=25;
        %         el=25;
        az=-10;
        el=10;
    end
    view(axis2plot,az,el)

    for n=1:numel(data)
        baseline=0;

        if isfield(data,'baseline') && not(isempty(data(n).baseline)) && size(data(n).baseline,1)==size(data(n).x_signal,1)

            if startsWith(settings(n).LineStyle.Baseline{:},'A')
                LineStyle=erase(settings(n).LineStyle.Baseline{:},'A');

                y_bl=data(n).baseline;
                y_bl(end)=min([data(n).baseline(:); 0]);
                y_bl(1)=min([data(n).baseline(:); 0]);

                fill3(axis2plot,data(n).x_signal   ,ones(size(data(n).x_signal)).*n,y_bl,settings(n).Color.Baseline,'LineWidth',settings(n).LineWidth.Baseline,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Baseline,'FaceColor',settings(n).Color.Baseline)
            else
                plot3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,data(n).baseline+offset(n),LineWidth=settings(n).LineWidth.Baseline,Color=settings(n).Color.Baseline,LineStyle=settings(n).LineStyle.Baseline);
            end

        else

        end

        if settings(n).PlotSpectra==true
            if startsWith(settings(n).LineStyle.Spectra{:},'A')
                LineStyle=erase(settings(n).LineStyle.Spectra{:},'A');
                fill3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,data(n).y_signal+offset(n),settings(n).Color.Spectra,'LineWidth',settings(n).LineWidth.Spectra,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Spectra,'FaceColor',settings(n).Color.Spectra)
            else
                plot3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,data(n).y_signal+offset(n),LineWidth=settings(n).LineWidth.Spectra,Color=settings(n).Color.Spectra,LineStyle=settings(n).LineStyle.Spectra);
            end
        end
        if settings(n).PlotRaw==true && isfield(data,'y_signal_raw') && not(isempty(data(n).y_signal_raw)) && isfield(data,'x_signal_raw') && not(isempty(data(n).x_signal_raw))
            if startsWith(settings(n).LineStyle.Raw{:},'A')
                LineStyle=erase(settings(n).LineStyle.Raw{:},'A');
                fill3(axis2plot,data(n).x_signal_raw,ones(size(data(n).x_signal_raw)).*n,data(n).y_signal_raw+offset(n),settings(n).Color.Raw,'LineWidth',settings(n).LineWidth.Raw,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Raw,'FaceColor',settings(n).Color.Raw)
            else
                plot3(axis2plot,data(n).x_signal_raw,ones(size(data(n).x_signal_raw)).*n,data(n).y_signal_raw+offset(n),LineWidth=settings(n).LineWidth.Raw,Color=settings(n).Color.Raw,LineStyle=settings(n).LineStyle.Raw);
            end
        end

        if  settings(n).PlotMean==true
            %%       plt(n)=plot3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,mean(data(n).y_signal,2,'omitnan')+offset(n),LineWidth=2);%# ok;%# ok

            if startsWith(settings(n).LineStyle.Mean{:},'A')
                LineStyle=erase(settings(n).LineStyle.Mean{:},'A');
                fill3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,mean(data(n).y_signal,2,'omitnan')+offset(n),settings(n).Color.Mean,'LineWidth',settings(n).LineWidth.Mean,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Mean,'FaceColor',settings(n).Color.Mean)
            else
                plot3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,mean(data(n).y_signal,2,'omitnan')+offset(n),LineWidth=settings(n).LineWidth.Mean,Color=settings(n).Color.Mean,LineStyle=settings(n).LineStyle.Mean);
            end
        end

        if  settings(n).PlotSum==true
            if startsWith(settings(n).LineStyle.Sum{:},'A')
                LineStyle=erase(settings(n).LineStyle.Sum{:},'A');
                fill3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,sum(data(n).y_signal,2,'omitnan')+offset(n),settings(n).Color.Sum,'LineWidth',settings(n).LineWidth.Sum,'LineStyle',LineStyle,'FaceAlpha',settings(n).FaceAlpha.Sum,'FaceColor',settings(n).Color.Sum)
            else
                plot3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,sum(data(n).y_signal,2,'omitnan')+offset(n),LineWidth=settings(n).LineWidth.Sum,Color=settings(n).Color.Sum,LineStyle=settings(n).LineStyle.Sum);
            end
        end


        if isfield(data,'DeconData') && not(isempty(data(n).DeconData)) && isfield(data(n).DeconData,'DeconResults') && not(isempty(data(n).DeconData(1).DeconResults)) && not(isempty(data(n).DeconData(1).DeconResults)) && (settings(n).PlotDeconDataSum ==true || any(settings(n).PlotDeconDataSingle==true))
            DeconID=data(n).DeconData.ActiveDeconGroup;
            y_calc=DeconPlotFun(data(n).DeconData.PeakFit(DeconID).peakshape,data(n).DeconData.PeakFit(DeconID).signal(:,1)',data(n).DeconData.DeconResults(DeconID).pos',data(n).DeconData.DeconResults(DeconID).width',data(n).DeconData.DeconResults(DeconID).hight',data(n).DeconData.PeakFit(DeconID).extra',data(n).DeconData.PeakFit(DeconID).DELTA');

            if ~isequal(baseline,0) || isequal(size(baseline),data(n).DeconData.PeakFit(DeconID).signal(:,1))
                baseline_dc=baseline(data(n).x_signal>=min(data(n).DeconData.PeakFit(DeconID).signal(:,1)) & data(n).x_signal<=(max(data(n).DeconData.PeakFit(DeconID).signal(:,1))),:);
            else
                baseline_dc=baseline;
            end

            for m=1:size(data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color,1)
                if data(n).DeconData.PlotSettings(DeconID).DeconSingle.Active(m)==true
                    if not(isempty(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m})) && startsWith(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m},'A')
                        LineStyle=erase(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m},'A');
                        fill3(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),ones(size(data(n).x_signal)).*n,offset(n)+baseline_dc.*ones(size(y_calc(m,:)))+offset(n)...
                            ,data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color(m,:),'LineWidth',data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineWidth(m),'LineStyle',LineStyle,'FaceAlpha',data(n).DeconData.PlotSettings(DeconID).DeconSingle.FaceAlpha(m),'FaceColor',data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color(m,:))

                    end
                end

            end
            %do i need the next plot3 decon input?
            if settings(n).PlotSpectra==true
                plot3(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),ones(size(data(n).DeconData.PeakFit(DeconID).signal(:,1))).*n,data(n).DeconData.PeakFit(DeconID).signal(:,2:end)+offset(n)+baseline_dc,'LineWidth',data(n).DeconData.PlotSettings(DeconID).Signal.LineWidth,'Color',data(n).DeconData.PlotSettings(DeconID).Signal.Color,'LineStyle',erase(data(n).DeconData.PlotSettings(DeconID).Signal.LineStyle,'A'));
            end
            if    settings(n).PlotDeconDataSum ==true
                plot3(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),ones(size(data(n).DeconData.PeakFit(DeconID).signal(:,1))).*n,sum(y_calc(logical(data(n).DeconData.PlotSettings(DeconID).DeconSingle.Active),:))'+offset(n)+baseline_dc,'LineWidth',data(n).DeconData.PlotSettings(DeconID).DeconTotal.LineWidth,'Color',data(n).DeconData.PlotSettings(DeconID).DeconTotal.Color,'LineStyle',data(n).DeconData.PlotSettings(DeconID).DeconTotal.LineStyle)
            end

            for m=1:size(data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color,1)
                if data(n).DeconData.PlotSettings(DeconID).DeconSingle.Active(m)==true
                    LineStyle=erase(data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineStyle{m},'A');
                    plot3(axis2plot,data(n).DeconData.PeakFit(DeconID).signal(:,1),ones(size(data(n).DeconData.PeakFit(DeconID).signal(:,1))).*n,y_calc(m,:)'+offset(n)+baseline_dc,'LineWidth',data(n).DeconData.PlotSettings(DeconID).DeconSingle.LineWidth(m),'Color',data(n).DeconData.PlotSettings(DeconID).DeconSingle.Color(m,:),'LineStyle',LineStyle)
                end
            end
        end
        %%
        %         if settings(n).PlotPeaks==true
        %             label=join([ string(round((data(n).DeconData.DeconResults.pos)')) data(n).DeconData.DeconInput.Peaks.PeakDB.PeakNames]);
        %             linecolor=data(n).DeconData.PlotSettings.PeakPos.Color;
        %             line_settings=data(n).DeconData.PlotSettings.PeakPos.LineType;
        %             PlotPeaksFun(axis2plot,data(n).x_signal,data(n).y_signal+offset(n),data(n).DeconData.DeconResults.pos,line_settings,label,linecolor)
        %         end

        if isfield(data,'baseline') && not(isempty(data(n).baseline)) && size(data(n).baseline,1)==size(data(n).x_signal,1)
            plot3(axis2plot,data(n).x_signal,ones(size(data(n).x_signal)).*n,data(n).baseline+offset(n),'--k',LineWidth=1);
        end

        %         if settings(1).PlotSpecNames==true && isfield(data(n),'Annotation')
        %             t(n)=text(axis2plot,min(settings(1).x_limits),n,offset(n),n,data(n).Annotation{:},'HorizontalAlignment','right','VerticalAlignment','top','Interpreter',settings(1).Annotation_interpreter);
        %         if isfield(settings(1),'AnnotationFontSize')
        %             t(n).FontSize=settings(1).AnnotationFontSize;
        %         end
        %         end


    end

    if settings(1).PlotSpecNames==true && isfield(data(n),'Annotation')
        yticks(axis2plot,1:numel([data.Annotation]))
        yticklabels(axis2plot,[data.Annotation]);
        axis2plot.TickLabelInterpreter=settings(1).Annotation_interpreter;
    else
        yticklabels(axis2plot,'auto');
    end


    if settings(1).PlotLegend==true
        legend(axis2plot,plt,data.SampleName,Interpreter='none')
    end
    %
    % end
    axis2plot.XLabel.String='Wavenumbers [cm^{-1}]';
    axis2plot.YLabel.String='';
    axis2plot.ZLabel.String='Absorbance [cm^{-1}]';

if  isfield(settings,'XLabel')
    axis2plot.XLabel.String=settings(1).XLabel;
else
    axis2plot.XLabel.String='Wavenumbers [cm^{-1}]';
end

if isfield(settings,'YLabel')
    axis2plot.YLabel.String=settings(1).ZLabel;
else
    axis2plot.ZLabel.String='';
end
if isfield(settings,'ZLabel')
    axis2plot.ZLabel.String=settings(1).YLabel;
else
    axis2plot.YLabel.String='Absorbance [cm^{-1}]';
end

end

