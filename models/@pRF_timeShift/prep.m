function signal = prep(obj,signal)

signal = signal - nanmean(signal,1);

end
