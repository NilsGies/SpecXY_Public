function [fig] = CorPlotFun(DataIn)
try
    ScreenSize = get(0,'ScreenSize');

    
DataMatrix=DataIn.DataMat;
notfill={'+','*','.','x', '_','|'};


if iscell(DataMatrix)
    ngroup2plot=numel(DataMatrix);
    [~,n] = size(DataMatrix{1});   % # of rows and columns of the input matrix

else
    %% This could maybe go?
    [~,n] = size(DataMatrix);   % # of rows and columns of the input matrix
    if n ~= length(VarNames)
        disp('Warning: # of variables & # of columns not matched.'), pause
    end
    Symbol=DataIn.Symbol;
    ngroup2plot=1;
    DataMatrix={DataMatrix};
    if not(iscell(Symbol))
        Symbol{1:size(cmap,1)}=Symbol;
    end
end
% Define a figure window position and size (full screen)
fig=figure('Position',[.05*ScreenSize(3) .05*ScreenSize(4) .9*ScreenSize(3) .87*ScreenSize(4)]);
% Adjust marker and font sizes according to the # of columns
    two_only=false;
if n ==2
     gap='tight';
    symbolsize = 8.909-0.4545*n;      % Marker size
    labelsize = 12.909-.5*n;         % Font size for text and ticklabel
       t = tiledlayout(1,1,'TileSpacing','none');
    disp_set='half';
    start=2;
    two_only=true;
elseif n < 8 && not(isfield(DataIn,'disp_set') && not(isempty(DataIn.disp_set)) && strcmp(DataIn.disp_set,'half'))
    gap='tight';
    symbolsize = 8.909-0.4545*n;      % Marker size
    labelsize = 12.909-.5*n;         % Font size for text and ticklabel
  %  t = tiledlayout(n,n+1,'TileSpacing','none');
    t = tiledlayout(n,n+1,'TileSpacing','tight');
    disp_set='full';
    start=1;
else     % when # of variables >= 14
    gap='none';
    gap='tight';

    symbolsize = 3;
    labelsize = 10;
    %t = tiledlayout(n-1,n+1,'TileSpacing','none');
    t = tiledlayout(n-1,n+1,'TileSpacing','tight');
    disp_set='half';
    start=2;
end
if isfield(DataIn,'labelsize') && not(isempty(DataIn.labelsize))
    labelsize=DataIn.labelsize;
end
caption='';
if isfield(DataIn,'link') && not(isempty(DataIn.link))
    link_set=DataIn.link;
else
    link_set=   false;
end

if isfield(DataIn,'ExcludeEmpty') && not(isempty(DataIn.ExcludeEmpty))
    ExcludeEmpty=DataIn.ExcludeEmpty;
else
    ExcludeEmpty=true;
end






%% Create subplot axes %%

for rows2plot = start:n      % Loop for columns
    for cols2plot = 1:n+1      % Loop for rows

       
        if rows2plot <= n
            if cols2plot <=n
                    if not(n==2)
                tc{rows2plot-(start-1),cols2plot}= nexttile;
               
                    end
                                            ax = gca;

                    if rows2plot ==start &&  cols2plot == 1
              ax.Toolbar.Visible=false;
        else
                   ax.Toolbar.Visible=true;
                    end

                if strcmp(disp_set,'half')
                    if    cols2plot > rows2plot-(start-1) && not(cols2plot==n+1)&& not(n==2)

                        ax.Visible=false;
                        continue
                    elseif n==2  && cols2plot>=2
                    continue
                    end
                end

                disp(['row: ' num2str(rows2plot) ' - col: ' num2str(cols2plot)])
                % title(['row: ' num2str(rows2plot) ' - col: ' num2str(col2plot)])
                hold on

                for k=1:ngroup2plot
                    % for h=1:nmat2plot

                    A=DataMatrix{k};
                    try
                    rows2keep=not((isnan(A(:,cols2plot)) | isnan(A(:,rows2plot))));
                    catch ME
                        0
                        continue
                    end
                    % X = rand(1,1000);
                    % X(2,:) = NaN;
                    % Y = rand(1,1000);
                    % Y(2,:) = NaN;
                    % H = plot(X,Y,'*');
                    X=[];
                    Y=[];
%                     X(1,:) = A(rows2keep,cols2plot);
%                    % X(2,:) = NaN;
%                     Y(1,:) = A(rows2keep,rows2plot);       % Data to be plotted
%                   %  X(2,:) = NaN;
                   % X(2,:) = NaN;
                  %  X(2,:) = NaN
                    VarNames=DataIn.VarNames;
                    classes=DataIn.classNumber{k};
                    ClassNames=DataIn.ClassNames{k};
                    cmap=DataIn.cmap{k};
                    Symbol=DataIn.Symbol{k};
                    SymbolSize=DataIn.SymbolSize{k};

                    X(:,1) = A(rows2keep,cols2plot);
                    Y(:,1) = A(rows2keep,rows2plot);       % Data to be plotted
                    nclasses_=unique(classes(not(isnan(classes))));
                    for j=1:numel(nclasses_)
                        nclasses=nclasses_(j);
                        if DataIn.Active{k}(j)==true
                            if ismember(Symbol{j},notfill)
                                plot(X(classes(rows2keep)==nclasses,:),Y(classes(rows2keep)==nclasses,:),Symbol{j},'MarkerFaceColor',cmap(j,:),'MarkerSize',SymbolSize(j),'MarkerEdgeColor',cmap(j,:),'LineWidth',SymbolSize(j)/10);
                            else
                                plot(X(classes(rows2keep)==nclasses,:),Y(classes(rows2keep,:)==nclasses),Symbol{j},'MarkerFaceColor',cmap(j,:),'MarkerSize',SymbolSize(j),'MarkerEdgeColor','k','LineWidth',SymbolSize(j)/10);
                            end
                        end

                    end
                end

                % Add ylabels to the subplots in the first column %
                if cols2plot == 1
                    %ylabel(regexprep(VarNames(rows2plot), '_', '-'),'FontSize',labelsize+2);
                    string_skalar = split_label_lines(VarNames(rows2plot),20);
                    % ylabel(regexprep(string_skalar, '_', '-'),'FontSize',labelsize+2,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
                    ylabel(string_skalar,'Interpreter','none','FontSize',labelsize+2,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
                    ax = gca;
                    ax.FontSize = labelsize;     % Ticklabel size
                else
                    set(gca, 'YTickLabel','');
                end

                % Add xlabels to the subplots in the last row %
                if rows2plot == n
                  %  xlabel(regexprep(VarNames(cols2plot), '_', '-'),'FontSize',labelsize+2);

                 
string_skalar = split_label_lines(VarNames(cols2plot),20);
                   xlabel(string_skalar,'Interpreter','none','FontSize',labelsize+2,'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle');

                    ax = gca;
                    ax.FontSize = labelsize;
                else
                    set(gca, 'XTickLabel','');
                end
                %
%                 xl =xlim;
%                 xlim([xl(1)*0.95 xl(2)*1.05])
                %
                clear X Y* Po

            else
                if rows2plot==start

                    if not(n==2)
                    nexttile([n-start+1,1])

                    end
                    xl=xlim;
                    yl=ylim;
                    hold on
                    ClassNames=[];

                    for k=1:ngroup2plot

                        classes=DataIn.classNumber{k};
                        cmap=DataIn.cmap{k};
                        Symbol=DataIn.Symbol{k};
                        SymbolSize=DataIn.SymbolSize{k};
                        Active=DataIn.Active{k};
                        %           ClassNames=[ClassNames;string(DataIn.ClassNames{k})];
                        nclasses_=unique(classes(not(isnan(classes))));
                        for j=1:numel(nclasses_)
                            nclasses=nclasses_(j);

                            testmat=DataMatrix{k}(classes==nclasses,:);

                           % if ExcludeEmpty==true && (Active(j)==false || any(size(testmat,1) == sum(isnan(testmat),1) ))
                            if ExcludeEmpty==true && (Active(j)==false || all(size(testmat,1) == sum(isnan(testmat),1)))
                                continue
                            else

                                ClassNames=[ClassNames;string(DataIn.ClassNames{k}(j))];

                                if ismember(Symbol{j},notfill)
                                    plot(-10^30,-10^30,Symbol{j},'MarkerFaceColor',cmap(j,:),'MarkerSize',SymbolSize(j)*10,'MarkerEdgeColor',cmap(j,:),'LineWidth',SymbolSize(j)/10);

                                else
                                    plot(-10^30,-10^30,Symbol{j},'MarkerFaceColor',cmap(j,:),'MarkerSize',SymbolSize(j),'MarkerEdgeColor','k','LineWidth',SymbolSize(j)/10);
                                end
                            end
                        end
                    end

                     xlim(xl)
                     ylim(yl)
                    if not(isempty(ClassNames))
                        if not(n==2)

                            legend(replace(ClassNames,""," "),'interpreter','none')
                        else
                            legend(replace(ClassNames,""," "),'interpreter','none',Location='eastoutside')
                        end
                    end

                    if exist('caption','var')

                        ts=sgtitle(caption);
                        set(ts, 'Interpreter', 'none')
                    end
                    if not(n==2)
                        ax = gca;
                        ax.Visible=false;
                    end
                end
            end   % if j > counter+1

        end    % j-loop for columns
    end     % if i > counter

end      % i-loop for rows
if link_set==true
    for k=1:n
        if k<=size(tc,1)
            linkaxes([tc{k,:}],'y')
        end
        if k<=size(tc,2)
            linkaxes([tc{:,k}],'x')
        end

    end
end

catch ME
0
end