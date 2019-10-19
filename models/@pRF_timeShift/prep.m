function data = prep(~,data)

dimdata = 1;
dimtime = 2;
nAcqs = length(data);

% Set "bad" voxels to have uniform zero values
bad = cellfun(@(x) any(~isfinite(x),dimtime) | all(x==0,dimtime),data,'UniformOutput',0);  % if non-finite or all 0
bad = any(cat(dimtime,bad{:}),dimtime);  % badness in ANY run
for ii=1:nAcqs
    data{ii}(repmat(bad,[ones(1,dimdata) size(data{ii},dimtime)])) = 0;
end


end
