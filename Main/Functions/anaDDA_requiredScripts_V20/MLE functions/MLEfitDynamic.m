function [parametersbest,bootstrapparammean,bootstrapparamstd] = MLEfitDynamic(Dlistdata, numberofspecies,fixedparameters,rangeD, Dfixed,fitspecies,fixedspecies,input)
%% Used to extract parameters based on MLE estimation.
% Can estimate parameters for up to three separate species, each with a koff, kon, Dfree and abundance c.
% Each species adds 4 parameters to be estimated if not restricted. However
% two degrees of freedom are removed due to the following equations:

% meanD of dataset
meanD = mean(Dlistdata(1,:));
fixedparameters(1,1) = 1;
% Number of degrees of freedom (all the ones that are 0)
% depend on species and which parameters are fixed
fixedparameterstemp = fixedparameters(1:numberofspecies,:);

startparameters = [0.3 50 30 1;0.3 1 1 1; 0.3 1 1 1];
indexfittingparameters = fixedparameterstemp==-1;
indexfittingparameters = [indexfittingparameters; zeros(3-numberofspecies,4)];
indexfittingparameters = logical(indexfittingparameters);

startparam = startparameters(indexfittingparameters);
lowerbound = zeros(length(startparam),1)+0.000001;
upperbound = [ones(sum(fixedparameterstemp(:,1)==-1),1) ;100*input.upperstartkon*input.upperstartkoff*ones(sum(sum(indexfittingparameters(:,2:3)==1)),1);50*ones(sum(indexfittingparameters(:,4)==1),1)];

Frametimelist = Dlistdata(3,:);
Numberofframes = Dlistdata(2,:);
Dlistdata = Dlistdata(1,:);
[Numberofframes, sortind]=sort(Numberofframes);
Dlistdata = Dlistdata(sortind);
Frametimelist = Frametimelist(sortind);
[Frametimelist, sortind]=sort(Frametimelist);
Dlistdata = Dlistdata(sortind);
Numberofframes = Numberofframes(sortind);

for j = 1:numel(input.frametimerange)
    table = tabulate(Numberofframes(Frametimelist==input.frametimerange(j)));
    frequency(j,:) = table(:,2);
    input.frametime = input.frametimerange(j);
    [fx(:,j),fy(:,j)] = Generateconfinedfunction(0:0.05:input.upperDfree,rangeD,input);
end

if input.compensatetracking == true
    maxD = (input.trackingwindow*input.pixelsize)^2/(4*input.frametime);
    maxDindtracking = round(maxD./(rangeD(2)-rangeD(1)));
else
    maxDindtracking = 0;
end

% Fitting function
custompdf = @(Dlistdata,varargin) pdfDvaluesMLE(Dlistdata, Numberofframes,Frametimelist,input,rangeD,fitspecies, fixedspecies, Dfixed, fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,frequency,varargin);


%% Start values for MLE
mincyclenumber = input.cyclenumber;
lowerstartkoff = input.lowerstartkoff;
upperstartkoff = input.upperstartkoff;

lowerstartkon = input.lowerstartkon;
upperstartkon = input.upperstartkon;
%startc = ;
%  startDfree = [1];
nllbest = 10^99;
% parametersbest = [];
maxDfree = input.upperDfree;
%parameters = zeros(sum(indexfittingparameters(:)),cyclenumber);
i = 1;
parametersbestsofar = 1000;
pass = 0;
while pass == 0
    
    disp(['Running fitting cycle ' num2str(i) ' of at least ' num2str(mincyclenumber)])
    for j = 1:3
        koffstart(j) = 10^((log10(upperstartkoff)-log10(lowerstartkoff))*rand()+log10(lowerstartkoff));
        konstart(j) = (upperstartkon-lowerstartkon)*rand+lowerstartkon;
    end
    startparameters = [rand koffstart(1) konstart(1)*koffstart(1) maxDfree*rand;rand koffstart(2) konstart(2)*koffstart(2) maxDfree*rand; rand koffstart(3) konstart(3)*koffstart(3) maxDfree];
    startparam = startparameters(indexfittingparameters);
    tableout = maketable(startparam,fixedparameters, indexfittingparameters,numberofspecies);
    disp('Starting parameters are:')
    disp(tableout)
    try
        parameters(:,i) = mle(Dlistdata, 'pdf',custompdf,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound);
        nll(:,i) = -sum(log(custompdf(Dlistdata',parameters(:,i))));
        disp('Found parameters are:')
        tableout = maketable(parameters(:,i),fixedparameters, indexfittingparameters,numberofspecies);
        disp(tableout)
    catch
        keyboard
        warning('error in fit: no result for this cycle')
        parameters(:,i) = 0;
        nll(:,i) = 1e10;
    end
    if i >= mincyclenumber
        nllrank = sort(nll);
        bestrun = find(nll==nllrank(1));
        secondbestrun = find(nll==nllrank(2));
        parametersbest = parameters(:,bestrun(1));
        parameterssecondbest = parameters(:,secondbestrun(1));
        if any(abs(parametersbest./parameterssecondbest-1)>0.05)
            pass = 0;
        else
            pass = 1;
        end
    end
    i =  i + 1;
end

parameters(:,nll==0)=[];
nll(:,nll==0)=[];
nllbest = min(nll);
bestrun = find(nll==nllbest);
parametersbest1 = parameters(:,bestrun(1));
disp('End of run, best parameters were:')
tableout = maketable(parametersbest1,fixedparameters, indexfittingparameters,numberofspecies);
disp(tableout)
%   if nllbest >  -sum(log(custompdf(Dlistdata',input.koff1_A,input.kon1_A,input.Dfree_A)))
%         sprintf('not enough cyclenumbers')
%   end

%% Add fitted parameters to already fixed parameters
parametersbest = fixedparameters;
try
    parametersbest(indexfittingparameters) = parametersbest1;
catch
    parametersbest(indexfittingparameters) = startparam;
    sprintf('no parameters found')
end
startparam = parametersbest1;

c = parametersbest(:,1);
koff = parametersbest(:,2);
kon = parametersbest(:,3);
Dfree = parametersbest(:,4);
if input.numberofspecies<3
    c(input.numberofspecies+1:end) = 0;
end
c(1) = 1 - c(2) - c(3);
parametersbest(1,1) = c(1);

Dmean2 = c(2)* Dfree(2)*(1-kon(2)/(koff(2)+kon(2)));
Dmean3 = c(3)* Dfree(3)*(1-kon(3)/(koff(3)+kon(3)));


% if fixedparameters(1,4) == 0 % If Dfree is not fixed
% parametersbest(1,4) = (meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*(1-(kon(1)/(koff(1)+kon(1)))));
% elseif fixedparameters(1,3) == 0 % If kon is not fixed
%    a = 1-(meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*Dfree(1));
%     parametersbest(1,3) = a*koff(1)/(1-a);
% else
% end

parametersbest = parametersbest(1:numberofspecies,:);
paramsize = size(parametersbest);


if input.bootstrapping == true
    disp('Initialising Bootstrapping')
    AmountofSubsamples = input.numberofbootstraps;
    bootstrapparameters=zeros(paramsize(1),paramsize(2),AmountofSubsamples);
    for i = 1:AmountofSubsamples
        disp(['Running Bootstrap ' num2str(i) ' of ' num2str(AmountofSubsamples)])
        Dbootstrap = [];
        for j = input.frametimerange
            for k = input.framerange
                index = Frametimelist == j & Numberofframes == k;
                Dlistdatatemp = Dlistdata(index);
                randomindex = randi(numel(Dlistdatatemp),numel(Dlistdatatemp),1);
                Dbootstrap(index) = Dlistdatatemp(randomindex);
            end
        end
        [bootstraptemp] = mle(Dbootstrap, 'pdf',custompdf,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound);
        
        bootstrapparametersbest = fixedparameters;
        bootstrapparametersbest(indexfittingparameters) = bootstraptemp;
        
        c = bootstrapparametersbest(:,1);
        % koff = bootstrapparametersbest(:,2);
        % kon = bootstrapparametersbest(:,3);
        % Dfree = bootstrapparametersbest(:,4);
        if input.numberofspecies<3
            c(input.numberofspecies+1:end) = 0;
        end
        c(1) = 1 - c(2) - c(3);
        bootstrapparametersbest(1,1) = c(1);
        
        % Dmean2 = c(2)* Dfree(2)*(1-kon(2)/(koff(2)+kon(2)));
        % Dmean3 = c(3)* Dfree(3)*(1-kon(3)/(koff(3)+kon(3)));
        %
        %
        % if fixedparameters(1,4) == 0 % If Dfree is not fixed
        % bootstrapparametersbest(1,4) = (meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*(1-(kon(1)/(koff(1)+kon(1)))));
        % elseif fixedparameters(1,3) == 0 % If kon is not fixed
        %    a = 1-(meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*Dfree(1));
        %     bootstrapparametersbest(1,3) = a*koff(1)/(1-a);
        % else
        % end
        
        bootstrapparametersbest = bootstrapparametersbest(1:numberofspecies,:);
        %paramsize = size(bootstrapparametersbest);
        bootstrapparameters(:,:,i)=bootstrapparametersbest;
        
    end
    
    bootstrapparamstd = std(bootstrapparameters,0,3);
    bootstrapparammean = mean(bootstrapparameters,3);
else
    paramsize = size(parametersbest);
    bootstrapparamstd = zeros(paramsize(1),paramsize(2));
    bootstrapparammean= parametersbest;
end

function [output] = pdfDvaluesMLE(x, Numberofframes,Frametimelist,input,rangeD, fitspecies, fixedspecies, Dfixed,fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,frequency,varargin)
parameters = fixedparameters;
try
    parameters(indexfittingparameters) = cell2mat(varargin{1});
catch
    keyboard
end
%% One degree of freedom less because meanD is linked to fonDNA and Dfree of all species
c = parameters(:,1);
koff = parameters(:,2);
kon = parameters(:,3);
Dfree = parameters(:,4);
if input.numberofspecies<3
    c(input.numberofspecies+1:end) = 0;
end
c(1) = 1 - c(2) - c(3);
%% If one of the kinetic parameters is not fixed it can be directly calculated from the meanD
Dmean2 = c(2)* Dfree(2)*(1-kon(2)/(koff(2)+kon(2)));
Dmean3 = c(3)* Dfree(3)*(1-kon(3)/(koff(3)+kon(3)));
% if fixedparameters(1,4) == 0 % If Dfree is not fixed
%      Dfree(1) = (meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*(1-(kon(1)/(koff(1)+kon(1)))));
% elseif fixedparameters(1,3) == 0 % If kon is not fixed
%     a = 1-(meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*Dfree(1));
%     kon(1) = a*koff(1)/(1-a);
% else
% end

%% Used to be able to use the PDA/Stracy distribution together with MLE in MATLAB. Returns a probability for each point in the distribution.
maxindex = numel(rangeD);
% maxDfree = max(Dfree(fitspecies))+input.sigmaerror^2/min(input.frametimerange);
% maxDfree = 10;
% [~,maxindex] = min(abs(-log(maxDfree*1e-10)*maxDfree-rangeD));
output = zeros(numel(Numberofframes),1);
dataind = 1;
for j = 1:numel(input.frametimerange)
    input.frametime = input.frametimerange(j);
    %locerrorpdf = input.dist(j).locerrorpdf(1:maxindex,:);
    locerrorpdfcorrected = input.dist(j).locerrorpdfcorrected;
    %locerrorpdfcorrected = locerrorpdfcorrected-locerrorpdf;
    combinedpdf = zeros(maxindex,8);
    for i = fitspecies
        [pdf]= DDistributiongenerator(koff(i),kon(i),Dfree(i),rangeD,locerrorpdfcorrected,maxindex,fx,fy,maxDindtracking,input,j);
        pdf = pdf./sum(pdf);
        combinedpdf = combinedpdf + c(i)*pdf;
    end
    if ~isempty(fixedspecies)
        for i = fixedspecies
            combinedpdf = combinedpdf + c(i)*Dfixed(i,j).dist';
        end
    end
    for i = input.framerange
        numberofdata = frequency(j,i);
        ind = i*(x(dataind:dataind+numberofdata-1)-rangeD(1))/(rangeD(2)-rangeD(1))+1;
        
        %KM adding: Removing the outliers from ind
        prctileCutOff = 95;
        cutoff = prctile(ind,prctileCutOff);
        ind(ind>=cutoff) = cutoff;
        
        %rangeDx = rangeD(1):(rangeD(2)-rangeD(1))/i:max(rangeD);
        %rangeDx = rangeDx(1:numel(combinedpdf(1:round(max(ind))+100,i)));
        outputtemp = interp1(combinedpdf(1:round(max(ind))+100,i),ind,'spline');
        %outputtemp = interp1(rangeDx,combinedpdf(1:round(max(ind))+100,i),x(dataind:dataind+numberofdata-1),'spline');
        output(dataind:dataind+numberofdata-1) = outputtemp;
        dataind = dataind+numberofdata;
    end
end

output = max(1e-99,output);
