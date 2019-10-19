function rawData = prep(~,rawData)

dimdata = 1;
dimtime = 2;


nAcqs = length(rawData);

% deal with data badness (set bad voxels to be always all 0)
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),rawData,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for ii=1:nAcqs
    rawData{ii}(repmat(bad,[ones(1,dimdata) size(rawData{ii},dimtime)])) = 0;
end


end
