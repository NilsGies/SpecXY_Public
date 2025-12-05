function ax2plot=SpecGifFun(x,y,options)
not(exist('options','var'))

if not(isfield(options,'outfile'))
    f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
    [file,path] = uiputfile('sweep.gif');
    options.outfile=fullfile(path,file);
    delete(f)
end


hfig=figure;
hfig.Position(3)=hfig.Position(3)*2;
hfig.Position(4)=hfig.Position(4)*1.5;
centerfig(hfig)
ax2plot=nexttile;
ax2plot.FontSize=22;
grid on

if not(isfield(options,'xlim'))
    options.xlim=[min(x) max(x)];
end

if not(isfield(options,'xlabel'))
    options.xlabel='x';
end

if not(isfield(options,'xlabel'))
    options.xlabel='y';
end


xlabel(options.xlabel)
ylabel(options.ylabel)
xlim(options.xlim)
ylim([min(y(x(:,1)>min(options.xlim) & x(:,1)< max(options.xlim),:),[],'all')*0.95 max(y(x(:,1)>min(options.xlim) & x(:,1)< max(options.xlim),:),[],'all')*1.05])
hold on
n=1;
plt_0=plot(x(:,1:n),y(:,1:n),'Color',[0 0 0,.3],'LineWidth',.5);

plt_1=plot(x(:,n),y(:,n),'Color',[0 0.4470 0.7410],'LineWidth',2);

for n=1:size(y,2)
    delete(plt_0)
    delete(plt_1)
    if (isfield(options,'title'))
        title([char(options.title(n)) ' (' num2str(n) '/' num2str(size(y,2)) ')'],'Interpreter','none')
    end


    if n>1
        plt_0=plot(x(:,1:n),y(:,1:n),'Color',[0.8 0.8 0.8 0.8],'LineWidth',.2);
    end
    plt_1=plot(x(:,n),y(:,n),'Color',[0 0.4470 0.7410],'LineWidth',3);


    drawnow;
    frame = getframe(1);
    im = frame2im(frame);

    [imind,cm] = rgb2ind(im,256);

    % On the first loop, create the file. In subsequent loops, append.
    if n==1
        t=0;
        imwrite(imind,cm ,options.outfile,'gif','DelayTime',0,'loopcount',inf);
    else
        imwrite(imind,cm ,options.outfile,'gif','DelayTime',.1,'writemode','append');
    end



end
delete(plt_1)
plt_0=plot(x(:,1),y(:,1:n),'Color',[0.8 0.8 0.8 0.8],'LineWidth',.5);
imwrite(imind,cm ,options.outfile,'gif','DelayTime',.1,'writemode','append');
imwrite(imind,cm ,options.outfile,'gif','DelayTime',.1,'writemode','append');
imwrite(imind,cm ,options.outfile,'gif','DelayTime',.1,'writemode','append');
plt_1=plot(x,y,'Color',[0 0.4470 0.7410],'LineWidth',1);
imwrite(imind,cm ,options.outfile,'gif','DelayTime',1,'writemode','append');
imwrite(imind,cm ,options.outfile,'gif','DelayTime',1,'writemode','append');

end