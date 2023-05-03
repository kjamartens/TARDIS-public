%Create and store trueMinDistJD

% Darr = [0.1 1 10]; %vertical
% densarr = [0.1 0.5 1 2.5 5 10 15];

% for diff = 1:size(Darr,2)
%     for densa = 1:size(densarr,2)
%         dens = densarr(densa);
%         Dfree = Darr(diff);
%         
%         storagefolder = 'G:\My Drive\CMU Employment\Articles\2020_relativeDisplacementTracking\ManuscriptFigureData\Figure2\RawData\';
%         mainstoragestring = ['singlePop_densMult10_' num2str(dens*10) '_locum2_DfreeMult10_' num2str(Dfree*10) '_um2s.mat'];
%         %Load data
%         load([storagefolder mainstoragestring])
        
inputdir = 'G:\My Drive\CMU Employment\Articles\2020_relativeDisplacementTracking\ManuscriptFigureData\Figure2_noFoVcutoff\SWIFT\'
S = dir(fullfile(inputdir,'*.mat'));
for k = 1:numel(S)
    fnm = fullfile(inputdir,S(k).name)
    load(fnm);
    
    filename = [inputdir '/' S(k).name(1:end-4) '_SwiftInput.csv'];
    fileID = fopen(filename,'w');
    fprintf(fileID,'id,"frame","x [nm]","y [nm]","z [nm]"\r\n');
    for i = 1:size(tracks,1)
        id = i;
        frame = tracks(i,3);
        xpos = tracks(i,1)*1e9;
        ypos = tracks(i,2)*1e9;
        zpos = 0;
        str = [num2str(id) ',' num2str(frame) ',' num2str(xpos) ',' num2str(ypos) ',' num2str(zpos) '\r\n'];
        fprintf(fileID,str);
    end
    fclose(fileID);
        %????? don't use below line And save the swift data
%         save([inputdir 'SWIFT_input/' mainstoragestring(1:end-4) '_SwiftInput.csv'],'TrueMinDisttracklist');
end
        %     end
% end