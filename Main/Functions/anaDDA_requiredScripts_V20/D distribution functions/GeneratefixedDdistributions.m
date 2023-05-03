function [Dfixed, fitspecies, fixedspecies] = GeneratefixedDdistributions(numberofspecies, fixedparameters, rangeD,input)
%% Generate fixed distributions of already predetermined parameters that you want to include in the fitting as a separate species
fixedspecies = [];
numberoffixedspecies = 0;
maxindex = length(rangeD);
for j = 1:numel(input.frametimerange)
input.frametime = input.frametimerange(j);
[fx(:,j),fy(:,j)] = Generateconfinedfunction(0:0.05:5,rangeD,input);

if input.trackingwindow < 100
      maxD = (input.trackingwindow*input.pixelsize)^2/(4*input.frametime);
      maxDindtracking = round(maxD./(rangeD(2)-rangeD(1)));
else
    maxDindtracking = 0;
end
%% Generate fixed distributions
for i = 1:numberofspecies
if sum(fixedparameters(i,2:4)>-1) == 3
    Dfixed(i,j).dist = DDistributiongenerator(fixedparameters(i,2),fixedparameters(i,3),fixedparameters(i,4),rangeD,input.dist(j).locerrorpdfcorrected,maxindex,fx,fy,maxDindtracking,input,j)';
    %Dfixed(i).dist = Dfixed(i,:)./(sum(Dfixed(i,:)));
    numberoffixedspecies = numberoffixedspecies + 1;
    fixedspecies(numberoffixedspecies) = i;
else
Dfixed(i,j).dist = zeros(8,maxindex);
end
end
end
fixedspecies = unique(fixedspecies);
fitspecies = setdiff(1:numberofspecies,fixedspecies);