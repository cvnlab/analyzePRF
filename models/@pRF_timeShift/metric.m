function metric = metric(obj, signal, x)

    metric = calccorrelation(signal, obj.forward(x));

end

