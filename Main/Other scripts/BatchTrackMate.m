

inputdir = 'G:\My Drive\CMU Employment\Articles\2020_relativeDisplacementTracking\ManuscriptFigureData\Figure2_noFoVcutoff\'
S = dir(fullfile(inputdir,'*.mat'));
for k = 1:numel(S)
    fnm = fullfile(inputdir,S(k).name);
    load(fnm);
    
    loclist = tracks(:,[3 1 2]);
    loclist(:,[2 3]) = loclist(:,[2 3])*1e6;
    savename = [inputdir 'TrackMate\' S(k).name(1:end-4) '_TrackMate.xml']
    locArrToTrackMate(loclist,savename,0.1);
end