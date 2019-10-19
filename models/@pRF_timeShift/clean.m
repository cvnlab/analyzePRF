function signal = clean(obj,signal)

signal = signal - nanmean(signal,1);

end
