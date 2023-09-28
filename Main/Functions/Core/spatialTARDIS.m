%% Exemplary script on how to perform spatialTARDIS routine
%The core concept: use a small part of the FoV as 'start point', but use the whole FoV
%as 'lookup'.
%Only works if the full FoV data is known

%Note: takes hours to run!

%First, we load the data
clear all
load('.\Example_data\spatialTARDIS_2Populations_confined.mat');
poslist = pos;

%We obtain settings
settingsTARDIS = GenerateTARDISsettings();
%User-define the most important settings
settingsTARDIS.maxdist = 3e-6; %Maximum dist (m)
settingsTARDIS.startpointBG = 2e-6;
settingsTARDIS.populations = 1; %number of populations to fit
settingsTARDIS.loc_unc = 30*1e-9; %Localization uncertainty in m
settingsTARDIS.fitWithBleach = 0; %Bleach fit or not
settingsTARDIS.performestimationfit = 0; %Estimation fit or not
settingsTARDIS.linorlogvis = 'lin'; %Linear or logartihmic fit+visualisation
settingsTARDIS.visualisationMLEIntFit = 0;
settingsTARDIS.bgbinningnr = 100;
settingsTARDIS.frame_dist_BG = [20:50];
settingsTARDIS.dt_arr = [1:3];

%Define a sub-area size and sub-area IDs. SpatialTARDIS will run in every
%subarea specified by this size.
subareasize = .25e-6;
subareaid_x = [1:round(max(poslist(:,2))/subareasize)];
subareaid_y = [1:round(max(poslist(:,3))/subareasize)];

%Looping over all subareas
for sx = subareaid_x
    for sy = subareaid_y
        disp(sx)
        disp(sy)
        %Define the start, end position (in m) of the sub-area
        subArea = [(sx-1)*subareasize,sx*subareasize; (sy-1)*subareasize,sy*subareasize];
        
        %Perform 3x TARDIS: 1x with restricted to slow diffusion, 1x with
        %restricted to fast diffusion, 1x with restricted to
        %double-diffusion.
        %Later, we check via the Olkin-Pratt estimator which of the three
        %(or just noise-fit) is the best.

        %P = 1 --> Single population, slow diffusion
        %P = 2 --> Single population, fast diffusion
        %P = 3 --> Double populations
        for p = 1:3
            if p < 3
                settingsTARDIS.populations = 1;
                if p == 1
                    %Set slow diffusion fit settings, lb, ub
                    settingsTARDIS.start_1pop= '[0.45]';
                    settingsTARDIS.lb_1pop= '[0.225]';
                    settingsTARDIS.ub_1pop= '[0.68]';
                elseif p == 2
                    %Set fast diffusion fit settings, lb, ub
                    settingsTARDIS.start_1pop= '[3.1]';
                    settingsTARDIS.lb_1pop= '[1.55]';
                    settingsTARDIS.ub_1pop= '[4.65]';
                end
            else
                %Set 2-population fit settings - quite constrained
                settingsTARDIS.populations = 2;
                settingsTARDIS.start_2pop= '[0.45 3.1 1]';
                settingsTARDIS.lb_2pop= '[0.225 1.55 0.5]';
                settingsTARDIS.ub_2pop= '[0.68 4.65 2]';
            end

            %Sub-area as TARDIS start area, defined as poslistsub
            poslistsub = pos;
            poslistsub(poslistsub(:,2)<subArea(1,1),:) = [];
            poslistsub(poslistsub(:,2)>subArea(1,2),:) = [];
            poslistsub(poslistsub(:,3)<subArea(2,1),:) = [];
            poslistsub(poslistsub(:,3)>subArea(2,2),:) = [];
            
            %Whole coordinates as lookup
            settingsTARDIS.AlternativeLookupPosList = poslist;
            
            %Run TARDIS as normal
            [parameters{p}, parametersCI{p}, paramEsts, HOfitCI, tottime, time, bgarr, truthoffsetpartial, ...
                Visualisation_HO_outputCell{p},Visualisation_FF_outputCell{p},anaDDAvisInfoHO,anaDDAvisInfoFF,...
                JDonlydata,SwiftParameters] =...
                TARDIS(poslistsub,settingsTARDIS);
            
            %Getting the R2 fit values for Olkin-Pratt estimator
            ydata = [];
            if p < 3
                for dt = 1:max(settingsTARDIS.dt_arr)
                    %Get visualisation info of the intra-/inter-emitter
                    %species
                    v1 = Visualisation_FF_outputCell{p}.extraoutput.main{dt}{1};
                    v2 = Visualisation_FF_outputCell{p}.extraoutput.main{dt}{2};
                    v2(isnan(v2)) = 0;
                    %Calculate the ydata as being the intra-/inter-emitter
                    %species with certain fractions
                    ydata = [ydata;...
                        v1/sum(v1)*parameters{p}(end-(dt-1))+...
                        v2/sum(v2)*(1-parameters{p}(end-(dt-1)))];
                end
            elseif p == 3
                for dt = 1:max(settingsTARDIS.dt_arr)
                    %Get visualisation info of the intra-/inter-emitter
                    %species
                    v1 = Visualisation_FF_outputCell{p}.extraoutput.main{dt}{1};
                    v2 = Visualisation_FF_outputCell{p}.extraoutput.main{dt}{2};
                    v2(isnan(v2)) = 0;
                    v3 = Visualisation_FF_outputCell{p}.extraoutput.main{dt}{3};
                    v3(isnan(v3)) = 0;
                    %Determine fractions between the populations
                    fraction_pops = (1-parameters{p}(end-(dt-1)));
                    fraction_betweenpops = [1/(1+parameters{p}(3)), 1-1/(1+parameters{p}(3))];
                    %Calculate the ydata as being the intra-/inter-emitter
                    %species with certain fractions
                    ydata = [ydata;...
                        v1/sum(v1)*parameters{p}(end-(dt-1))+...
                        v2/sum(v2)*fraction_pops*fraction_betweenpops(1)+...
                        v3/sum(v3)*fraction_pops*fraction_betweenpops(2)];
                end
            end
            
            %Calculate the R-squared value
            %Basically, it's R-squared of the raw JD data compared to the
            %TARDIS fit
            nrbins_forR2 = settingsTARDIS.bgbinningnr+1;
            %binval is histogram of inputdata
            hist_binEdges = [linspace(0,settingsTARDIS.maxdist,nrbins_forR2+1)]';
            binval = zeros(nrbins_forR2*max(settingsTARDIS.dt_arr),1);
            hist_R2_binmids = zeros(nrbins_forR2*max(settingsTARDIS.dt_arr),1);
            for dt = 1:max(settingsTARDIS.dt_arr)
                hist_R2 = histcounts(Visualisation_FF_outputCell{p}.JDarrSignalCell{1,dt},hist_binEdges,'Normalization','Probability');
                hist_R2_binmids(1+(dt-1)*nrbins_forR2:(dt)*nrbins_forR2) = hist_binEdges(1:end-1)+hist_binEdges(2)/2;
                binval(1+(dt-1)*nrbins_forR2:(dt)*nrbins_forR2) = hist_R2;
            end
            modelData = ydata;
            modelData = modelData./sum(modelData)*max(settingsTARDIS.dt_arr);
            SStot = sum((binval-mean(binval)).^2);
            SSres = sum((binval-modelData).^2);
            R2 = 1-SSres/SStot;
            
            %Obtain the adjusted R2 values
            R2_adj_arr(p+1) = Olkin_Pratt_Estimator(R2,max(size(hist_R2_binmids)),size(parameters{p},2),10);
        end
        
        %Determine the adjusted R2 value for noise-only
        modelData_noPop = [];
        for dt =  1:max(settingsTARDIS.dt_arr)
            modelData_noPop_dt{dt} = Visualisation_FF_outputCell{1}.extraoutput.main{dt}{1}./sum(Visualisation_FF_outputCell{1}.extraoutput.main{dt}{1});
            modelData_noPop = [modelData_noPop;modelData_noPop_dt{dt}];
        end
        
        SStot_nopop = sum((binval-mean(binval)).^2);
        SSres_nopop = sum((binval-modelData_noPop).^2);
        
        R2_nopop = 1-SSres_nopop/SStot_nopop;
        R2_adj_arr(1) = Olkin_Pratt_Estimator(R2_nopop,max(size(hist_R2_binmids)),0,10);

        %Store this as information - later this can be recalled and used
        %for visualisation
        R2adjarr_cell{sx,sy} = R2_adj_arr;
        params_1pop{sx,sy} = parameters{1};
        params_2pop{sx,sy} = parameters{2};
        params_1popCI{sx,sy} = parametersCI{1};
        params_2popCI{sx,sy} = parametersCI{2};
    end
end