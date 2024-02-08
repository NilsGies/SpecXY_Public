function PlotPeaksFun(axis2plot,wn,in,PeakPosition,line_settings,label,linecolor,peaklablesize,peaklinewidth)
if not(exist('linecolor','var')) || not(size(PeakPosition,1)==numel(PeakPosition))
    linecolor=zeros(numel(PeakPosition),3);
end
if not(exist('peaklablesize','var')) || isnan(peaklablesize)
    peaklablesize=15;
end

if not(exist('peaklinewidth','var'))
    peaklinewidth=2;
end


if not(isempty(PeakPosition)) %&& app.Peaks_CheckBox.Value==true || strcmp(app.PeakPreviewSwitch.Value,'On')

    if strcmpi(line_settings,'line')
        Peak_marker_offset_value=0.25; %ATH add settings
        offset=max(in,[],'all')*Peak_marker_offset_value*1.1;

        PeakIntensity=nan(size(PeakPosition));

        for m=1:length(PeakPosition)
            [ ~, ix ] = min( abs(wn-PeakPosition(m) ) );
            PeakIntensity(m)=max(in(ix,:));
        end
        x=[PeakPosition',PeakPosition'];
        y=[zeros(size(PeakIntensity'))+(min(in,[],'all')),PeakIntensity'];

        plot(axis2plot,x',y' ,'--k','linewidth',2);
        if not(strcmp(label,''))
            plot(axis2plot,PeakPosition,PeakIntensity+offset*0.25,"Marker",'v','MarkerFaceColor','k','MarkerSize',10,'LineStyle','none')
            try
                text(axis2plot,PeakPosition,PeakIntensity+offset*0.5,cellstr([num2str(PeakPosition) repmat([' ' label],size(PeakPosition))]),"HorizontalAlignment","center","Rotation",0)
            catch
                text(axis2plot,PeakPosition,PeakIntensity+offset*0.5,label,"HorizontalAlignment","center","Rotation",0)
            end% text(axis2plot,PeakPosition,PeakIntensity+offset,cellstr([num2str(PeakPosition) repmat([' ' label],size(PeakPosition))]),"VerticalAlignment","middle")
        end

    elseif strcmpi(line_settings,'XLine')
        try
       %     xln=xline(axis2plot,PeakPosition,'--',cellstr([num2str(PeakPosition) repmat([' ' label],size(PeakPosition))]));
        %    xln=xline(axis2plot,PeakPosition,'--',cellstr([num2str(PeakPosition) string(repmat(' ',size(PeakPosition))),label]));
        xln=xline(axis2plot,PeakPosition,'--',cellstr(join([num2str(PeakPosition) string(repmat('-',size(PeakPosition))) string(label) ])));
        catch
            xln=xline(axis2plot,PeakPosition,'--',label);
        end
        [xln.LineWidth]=deal(peaklinewidth);
        [xln.FontSize]=deal(peaklablesize);

       
        for m=1:length(PeakPosition)
        if size(linecolor,1)==1
                        [xln(m).Color]=deal(linecolor(1,:));
        else
            [xln(m).Color]=deal(linecolor(m,:));
        end
        end
    end



end