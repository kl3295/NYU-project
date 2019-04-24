
function AsianOption(S0,sigma,T,r,K,m,n)

path=zeros(n,1+m);
X=zeros(n,1);
Y=zeros(n,1);
path(:,1)=S0;
dt=T/m;
nudt=(r-sigma*sigma/2)*(dt);
sidt=sigma*sqrt(dt);
rng('default');
for i=1:n
    for j=1:m
        path(i,j+1)=path(i,j)*exp(nudt+sidt*randn);
    end
    X(i,1)=exp(-r*T)*max(power(prod(path(i,2:m+1)),1/m)-K,0);
    Y(i,1)=exp(-r*T)*max((sum(path(i,2:m+1))/m-K),0);
end
C=cov(X,Y)
b=C(1,2)/var(X);
Targ=(m+1)/2*dt;
sigma2=sigma*sigma/m/Targ*((m+1)*(2*m+1)/6)*dt;
delta=sigma*sigma/2-sigma2/2;
d=(log(S0/K)+(r-delta+sigma2/2)*Targ)/sqrt(sigma2*Targ);
EX=normcdf(d,0,1)*S0*exp(-delta*Targ-r*(T-Targ))-normcdf(d-sqrt(sigma2*Targ),0,1)*K*exp(-r*T);

optionprice=mean(Y-b*(X-EX))

error=var(Y)-2*b*C(1,2)+b*b*var(X)
coef=C(1,2)/sqrt(var(X)*var(Y))
end

