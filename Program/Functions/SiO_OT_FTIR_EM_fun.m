function SiO_OT_FTIR_EM_fun(y_signal,x_signal,Mineral,num_EM,cmap)

if strcmp(Mineral,'CPX') && num_EM==6

    % here sort
    % beta max(2200)
    % alpha max (1973) & min(1582)
    % gamma min(1973)

    [ ~, ix ] = min(abs(x_signal-2200));
    [ ~, pos_b ]=max(y_signal(ix,:));

    [ ~, ix ] = min(abs(x_signal-1973));
    [ ~, pos_a ]=max(y_signal(ix,:));

    [ ~, ix ] = min(abs(x_signal-1973));
    [ ~, pos_c ]=min(y_signal(ix,:));

    % X= max 1641rea
    % Z= max 1865

    %X  1965/1765
    %Z  1565/1657

    [ ~, ix ] = min(abs(x_signal-1965));
    [ ~, ix2 ] = min(abs(x_signal-1765));
    [ val, pos_X ]=max(y_signal(ix,:)./y_signal(ix2,:));

    pos_X=pos_X(not(ismember(pos_X,[(pos_a),(pos_b),(pos_c)])));
    pos_X=pos_X(end);

    [ ~, ix ] = min(abs(x_signal-1565));
    [ ~, ix2 ] = min(abs(x_signal-1657));

    [ val, pos_Z ]=max(y_signal(ix,:)./y_signal(ix2,:));

    figure
    hold on
    plot(x_signal,y_signal(:,1:end-1),'-','markersize',1,'color', [0.4 0.4 0.4 0.15],'handlevisibility','off');
    plot(x_signal,y_signal(:,end),'-','markersize',1,'color', [0.4 0.4 0.4 0.15]);

    plot(x_signal,y_signal(:,pos_a),'linewidth',2,'Color',cmap(1,:))
    plot(x_signal,y_signal(:,pos_b),'linewidth',2,'Color',cmap(2,:))
    plot(x_signal,y_signal(:,pos_c),'linewidth',2,'Color',cmap(3,:))
    plot(x_signal,y_signal(:,pos_X),'linewidth',2,'Color',cmap(4,:))
    plot(x_signal,y_signal(:,pos_Z),'linewidth',2,'Color',cmap(5,:))

    legend({'All  Spectra','\alpha','\beta','\gamma','X','Z'})

    figure

    %% subplots for each selected band and then iterate through all subplots
%     nexttile
%     [ ~, ix ] = min(abs(x_signal-1965));
% 
%     [ ~, ix2 ] = min(abs(x_signal-x_signal(1)));
%     s0bJ= scatter(y_signal(ix,:),y_signal(ix2,:),'.k')
%     hold on
%     s0bJ_a= scatter(y_signal(ix,pos_a),y_signal(ix2,pos_a),100,'filled','MarkerFaceColor',cmap(1,:))
%     s0bJ_b= scatter(y_signal(ix,pos_b),y_signal(ix2,pos_b),100,'filled','MarkerFaceColor',cmap(2,:))
%     s0bJ_c= scatter(y_signal(ix,pos_c),y_signal(ix2,pos_c),100,'filled','MarkerFaceColor',cmap(3,:))
% 
% 
%     for n=1:size(x_signal,1)
%         [ ~, ix2 ] = min(abs(x_signal-x_signal(n)));
%         s0bJ.YData=y_signal(ix2,:);
%         s0bJ_a.YData=y_signal(ix2,pos_a);
%         s0bJ_b.YData=y_signal(ix2,pos_b);
%         s0bJ_c.YData=y_signal(ix2,pos_c);
%         ylabel(num2str(x_signal(ix)))
%         xlabel(num2str(1965))
%         pause(0.1)
%     end

elseif strcmp(Mineral,'OPX') && num_EM==3
    % alpha max(1944)
    % beta max(2146)
    % gamma min(1944)

elseif strcmp(Mineral,'Olivine') && num_EM==6
    % alpha max(2035)
    % beta max(1696)
    % gamma max(1787)

end

end