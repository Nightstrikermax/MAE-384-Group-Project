%%Problem 3 setup & eq's
% Y = Gamma
%Part 1
clc
h = 1;
t = linspace(1,30,1);
B = .3;
Y = .1;
I =@(t,I) I(0)*exp(k*t)



%Part 2
t = linspace(1,30,1);
N = 1000;
S = 990;
Y = .1;

%k = ((B*S)/N) - Y;
%log(I) = log(I(0)) + kt;

%least squares setup
for i=1:length(t)
k = ((length(t)*(sum(t(i)*I(t)))-(sum(t(i))*sum(I(t)))))/((length(t(i))*sum((t(i)).^2))-(sum(t(i)).^2))
end

%Part 3
t = linspace(1,10,1);


