function [signal_corrected]=lin_BL_corr(signal,plot_on)

        min_signal_y=(signal(1,2));
        max_signal_y=(signal(end,2));
        min_signal_x=min(signal(:,1));
        max_signal_x=max(signal(:,1));
        
        signal_corrected=zeros(size(signal));

        for k=2:size(signal,2)
            correction_line=linspace(min_signal_y,max_signal_y,size(signal,1));

            signal_corrected(:,k)=signal(:,k)-correction_line';
            signal_corrected((signal_corrected(:,k)<0),k)=0;
            
            if exist("plot_on","var") && plot_on==1
            figure (30000+k)
            hold on 
             
            plot(signal_corrected(:,1),signal_corrected(:,k))
            plot(signal(:,1),signal(:,k))

            legend('corrected signal','signal')
            set(gca,'XDir','reverse') 
            xlabel('Wavenumbers (cm-1)');
            ylabel('Absorbance');

            hold off

            end
         
        end
