function trackmatrix = TrueMinDistTracking(poslist,maxdist)
%% True minimal jump distance analysis for tracking
%goal: fix localizations without blinking from 1 frame to the other, based
%solely on local minimal distance 'optimization'.
%Input: loclist, maxdist
%Output: tracklist

% maxdist = 2e-6;

%Get frame info
frinfo = accumarray(poslist(:,1),1);

%Loop over frames and start connecting
firstframe = min(poslist(:,1));
% curframe = 2;
clear trackmatrix
for curframe = firstframe:max(poslist(:,1))-1
    try
    if mod(curframe,100) == 0
        fprintf('Frame %d of %d\n',curframe,max(poslist(:,1)));
    end
    if curframe == firstframe
        trackmatrix = [poslist((sum(frinfo(1:curframe-1))+1):(sum(frinfo(1:curframe))),[1 2 3]) ...
            [1:size(poslist((sum(frinfo(1:curframe-1))+1):(sum(frinfo(1:curframe))),[1 2 3]),1)]'];
    end
    
    trinfo = accumarray(trackmatrix(:,1),1);
    stuck = 0;
    loclistcurframe = poslist((sum(frinfo(1:curframe-1))+1):(sum(frinfo(1:curframe))),[1 2 3]);
    loclistnextframe = poslist((sum(frinfo(1:curframe))+1):(sum(frinfo(1:curframe+1))),[1 2 3]);
    
    while (stuck == 0) && (size(loclistcurframe,1) > 0 && size(loclistnextframe,1) > 0)
        %create 2d dist matrix for all locs in curframe to curframe+1
        distMatrix = pdist2(loclistcurframe(:,[2 3]),loclistnextframe(:,[2 3]));
        %if min dist is larger than maxdist, then stop analysing
        if min(distMatrix(:)) > maxdist
            stuck = 1;
        else
            %find min distMatrix
            [mincol,minrow] = find(distMatrix == min(distMatrix(:)));
            %give the mincol entry in the loclistnextframe the same trackID as the
            %minrow entry in the loclistcurframe, and remove both from current list
    %         positionInTrackMatrixMinRow = sum(trinfo(1:curframe-1))+1+minrow-1;
            positionInTrackMatrixMinRow = find(sum(trackmatrix(:,[2 3]),2) == sum(loclistcurframe(mincol,[2 3]),2));
            trackmatrix = [trackmatrix; [loclistnextframe(minrow,:) trackmatrix(positionInTrackMatrixMinRow,4)]];
            loclistcurframe(mincol,:) = [];
            loclistnextframe(minrow,:) = [];
        end
    end
    %if there is still info in loclistnextframe, give that a new trackid
    trackmatrix = [trackmatrix; loclistnextframe max(trackmatrix(:,4))+...
        [1:size(loclistnextframe,1)]'];
    catch
        keyboard
    end
end
end