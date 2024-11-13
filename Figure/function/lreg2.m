% This code is written for linear regression analysis.
% Reference: Zar JH. 1999. Biostatistical analysis 4th Edition. Prentice-Hall, Inc. New Jersey, USA
% 
% /////////////////////////////////////////////////////////////
% Contact information:
% Youngryel Ryu, Ph.D. Student
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% Email: yryu@nature.berkeley.edu
% Homepage: http://nature.berkeley.edu/~yryug
% /////////////////////////////////////////////////////////////

% opt=0: force to pass the origin (zero intercept)
% opt=1: simple linear regression

% 2016 07 15: Included RMSE and bias

function [coe]=lreg2(x, y, opt)
usex=~isnan(x);
usey=~isnan(y);
use=find(usex.*usey==1);
x=x(use);
y=y(use);

if opt==0
    int=0;
    slo=sum(x.*y)./sum(power(x, 2)); % x\y
    SStot=sum(y.^2);
    SSreg=power(sum(x.*y), 2)/sum(power(x, 2));
    SSerr=sum(power(y-x*slo, 2));
    SSy=sum(power(y-mean(y), 2));
    MSE=(SStot-SSreg)/(length(use)-1);
    SEslo=sqrt(MSE/sum(x.^2));
    r2=1-SSerr/SSy;        
    t=slo/SEslo;
    p=2*(1-tcdf(abs(t), length(x)-1));
end

if opt==1
    coef=polyfit(x, y,1);
    slo=coef(1);
    int=coef(2);
    r=corrcoef(x, y);
    r2=r(1, 2)^2;
    SEslo=abs(slo/sqrt(length(use)-2)*sqrt(1/r2-1));
    t=slo/SEslo;
    p=2*(1-tcdf(abs(t), length(x)-1));
end
    
rmse=RMSE(x(:), y(:));
rmseR=rmse./nanmean(x(:));
bias=nanmean(y(:)-x(:));
biasR=bias./nanmean(x(:));

coe(1)=slo;
coe(2)=int;
coe(3)=r2;
coe(4)=SEslo;
coe(5)=t;
coe(6)=p;
coe(7)=rmse;
coe(8)=rmseR;
coe(9)=bias;
coe(10)=biasR;
coe(11)=length(x);


