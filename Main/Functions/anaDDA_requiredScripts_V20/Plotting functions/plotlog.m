function plotlog(framenr,parameters,D, input,rangeD,bootstrapparamstd,frametimerange,singlefigure)

load(fullfile(fileparts(mfilename('fullpath')),'layoutparameters.mat'))
layoutparameters.range = -2.5:0.05:1.5;
maxindex = numel(rangeD);

D = D(1,D(2,:)==framenr);
logrange = 10.^(min(layoutparameters.range):0.01:min([max(rangeD) max(layoutparameters.range)]));% Converted to x values
plotrange = min(layoutparameters.range):layoutparameters.plotstep:max(layoutparameters.range)-0.01;
histrange = min(layoutparameters.range):layoutparameters.histstep:max(layoutparameters.range);
histplotratio = layoutparameters.histstep/layoutparameters.plotstep;
histogram(log10(D), histrange,'Normalization','probability','Facecolor', layoutparameters.facecolor, 'Edgecolor',layoutparameters.edgecolor,'HandleVisibility','off');
%D is histogram data
totallogpdf = zeros(numel(logrange)-1,1);

[fx,fy] = Generateconfinedfunction(0:0.05:5,rangeD,input);
fx = fx';
fy = fy';
if input.trackingwindow < 100
      maxD = (input.trackingwindow*input.pixelsize)^2/(4*input.frametime);
      maxDindtracking = round(maxD./(rangeD(2)-rangeD(1)));
else
    maxDindtracking = 0;
end

color = [0.3 0.8 0.3];
for ii = 1:size(parameters)
koff = parameters(ii,2);
kon = parameters(ii,3);
Dfree = parameters(ii,4);
c = parameters(ii,1);
framescombined = DDistributiongenerator(koff,kon,Dfree,rangeD,input.dist(frametimerange).locerrorpdfcorrected,maxindex,fx,fy,maxDindtracking,input,1);
framescombined = framescombined./sum(framescombined);
framescombined = c*framescombined(:,framenr);

func = @(x) interp1(framescombined,x,'spline');
lograngetrue = (logrange-rangeD(1))*framenr./(rangeD(3)-rangeD(2))+1;

for i = 1:numel(lograngetrue)-1
    %so here JV actually defines the function for whatever fucking reason
logpdf(i) = integral(func,lograngetrue(i),lograngetrue(i+1));
end
color = circshift(color,1);
% framescombined = interp1(framescombined,1:1/precisionfactor:numel(framescombined),'spline');
% maximumplotDvalue = max(logrange);% Max x value for your log plot
% logrange2 = round(logrange*(precisionfactor*framenr)/(rangeDStracy(3)-rangeDStracy(2))); % edges for each bin 
% % logpdf converted from linear precisepdf
% logpdf = zeros(numel(logrange)-1,1);
% 
% for i = 1:numel(logrange)-1
%     logpdf(i) = sum(framescombined(logrange2(i):logrange2(i+1))); 
% end
totallogpdf = totallogpdf + logpdf';
if size(parameters)>1
plot(plotrange,logpdf*histplotratio,'LineWidth',layoutparameters.linewidth,'Color',color)
end
legendinfo(ii) = strcat({'species'},{num2str(ii)},{' '},{num2str(round(100*parameters(ii,1)))},{'\pm'},{num2str(round(100*bootstrapparamstd(ii,1)));}, {' %'});
end
% keyboard
plot(plotrange,histplotratio*totallogpdf,'LineWidth',layoutparameters.linewidth,'Color','Black')
if singlefigure == true
defaultlayout(layoutparameters,legendinfo)
end
function defaultlayout(layoutparameters,legendinfo)
%% Plots figures with the defaultlayout
set(gca,'fontsize',layoutparameters.fontsize)
axis([-2 1 0 0.035])

legend(legendinfo,'FontSize',layoutparameters.fontsizelegend,'Location','northwest')
xticks([-2 -1 0 1])
xticklabels({'0.01','0.1','1','10'})
%xlabel(['Diffusion Coefficient  ' 956 'm^2/s'])
yticks([0:0.01:0.05])
yticklabels({'0','1','2','3','4','5'})
if layoutparameters.halflength == 1
legend(legendinfo,'FontSize',layoutparameters.fontsizelegend*2,'Location','northwest')
pos = get(gcf, 'Position');
 pos(3) = pos(3)*2;
 set(gcf, 'Position', pos)
 layoutparameters.fontsize = layoutparameters.fontsize*2;
 set(gca,'fontsize',layoutparameters.fontsize);
%  yticks([0.0:0.02:0.04])
% yticklabels({'0','2','4'})
 yticks([0.0:0.015:0.03])
yticklabels({'0','2','4'})
end
%ylabel('Tracks (%)')
set(gca,'TickDir','out')
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gca,'ZColor','k')
set(gca,'linewidth',3)
set(gca,'box','off')