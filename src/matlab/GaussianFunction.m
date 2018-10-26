function f1 = GaussianFunction(H,x0,sigma,x)
    f1 = H*exp(-(x-x0).^2/sigma^2);
