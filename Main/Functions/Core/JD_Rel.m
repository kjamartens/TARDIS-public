%% Function JD_rel
%Calculate the relative Jumping distance between sets of localizations
%Input parameters:
%inputarr: full frame-x-y array of whole data
%framedist: distance in frames to be calculated (should be integer)
%maxdist: maximum distance than linking should be reported - same units as
%inputarr
%verbosemin: 1 if verbose is wanted, 0 otherwise.
%Output parameters:
%totfinalarr: 1D array of all relative jumping distances
%Required scripts: calculaterelativepos
function totfinalarr = JD_Rel(inputarr,framedist,maxdist,verbosemin)
    %Get frame list for further processing
    framelist = accumarray(inputarr(:,1),1);
    % Now we have a list of localisations belonging to tracks, ordered by frame
    %lets do some relative distance calculation
    startfr = 1;
    endfr = max(inputarr(:,1))-framedist;
    totsize = 0;
    verbose = verbosemin;
    
    %Remove 0-values only if we're comparing same frame
    if framedist == 0
        removediagonal = 1;
    else
        removediagonal = 0;
    end
    for fr = startfr:endfr
        startid = sum(framelist(1:fr-1))+1;
        frame1list = inputarr(startid:startid+framelist(fr)-1,:);
        startid = sum(framelist(1:fr-1+framedist))+1;
        frame2list = inputarr(startid:startid+framelist(fr+framedist)-1,:);
        finalarrcell{fr} = calculaterelativepos(frame1list,frame2list,maxdist,verbose,removediagonal);
        totsize = totsize+size(finalarrcell{fr},1);
        if (mod(fr,100) == 0) && (verbosemin)
            disp(['Frame ' num2str(fr) ' of ' num2str(endfr)])
        end
    end
    totfinalarr = zeros(totsize,1);
    curpos = 1;
    for fr = startfr:endfr
        if ~isempty(finalarrcell{fr})
            totfinalarr(curpos:curpos+size(finalarrcell{fr},1)-1) = finalarrcell{fr};
            curpos = curpos+size(finalarrcell{fr},1);
        end
    end
end

