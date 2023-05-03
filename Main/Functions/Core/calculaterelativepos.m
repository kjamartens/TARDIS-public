%% Function calculaterelativepos
%Calculates the relative position between two lists of localizations
%input parameters:
%frame1list, frame2list: list of localizations where relative position
%should be calculated for
%maxdist: Maximum distance to be linked
%verbose: whether verbose output is wanted
%remove diagonal: if same dataset is compared with itself, half the data
%should be removed (along with the diagonal). Putting this to 1 does that,
%leaving it at 0 keeps all data
%Output parameter:
%Finalarr: 1D arr with relative distances
function finalarr = calculaterelativepos(frame1list,frame2list,maxdist,verbose,removediagonal)
%Calculates relative positions between all localizations in frame 1 and
%frame 2 - currently assuming 2d arrays for simplicity
finalarr = pdist2(frame1list(:,[2 3]),frame2list(:,[2 3]));
%If we're removing the diagonal, we assume same size in x and y (i.e. same
%number of localizations) - this is needed for 0-frame 0 removal
if removediagonal
    for xx = 1:size(frame1list,1)
        for yy = 1:size(frame1list,1)
            if xx > yy
                finalarr(xx,yy) = 0;
            end
        end
    end
end
%Make matrix 1D
finalarr = finalarr(:);
%finish remove diagonal option
if removediagonal
    finalarr(finalarr == 0) = [];
end
a=size(finalarr,1);
finalarr(finalarr>maxdist) = [];
b=size(finalarr,1);
if verbose
disp(['Pre filt: ' num2str(a) '  Post filt: ' num2str(b)])
end
end