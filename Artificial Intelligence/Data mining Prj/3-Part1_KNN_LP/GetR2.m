function R2=GetR2(y,yhat)

    R2=1-mean((y-yhat).^2)/var(y,1);

end