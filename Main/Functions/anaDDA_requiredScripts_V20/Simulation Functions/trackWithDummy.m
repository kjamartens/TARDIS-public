function tracks = trackWithDummy(pos,trackParams)
% add a dummy track outside the field of view to ensure continuous spacing in time
dummyTrack(:,1:size(pos,2)-1) = 10000*ones(int32(max(pos(:,end))),size(pos,2)-1);
dummyTrack(:,size(pos,2)) = 1:int32(max(pos(:,end))); %continuous frame counter
posAndDummy = [pos; dummyTrack];

posAndDummy = sortrows(posAndDummy,size(pos,2)); % sort tracks by ascending frame number

tracks = track(posAndDummy,trackParams.maxDisp,trackParams);
% delete dummy track
tracks(tracks(:,1) == 10000,:) = [];

end