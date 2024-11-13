% RMSE: Computes the root-mean squared difference between two vectors
%                stored in vector x and y 
%
%              x,y = 1-D data vector
%              usage: z= RMSE(x,y)
function z=RMSE(x,y)
a=(x-y).^2;
b=nanmean(a);
z=b.^0.5;
