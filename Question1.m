clear
clc
load('dataset1.1.mat')
%%
% Parameters Y and duration buildup
Y = zeros(15,85);
for i=1:15  % from day 1 to day 15
    for j=1:85  % from proj 1 to day 85
        Y(date_ready(j):deadline(j),j)=1;
    end
end
% We duplicate both Y and duration 7 times to account for each crew.
duration_temp = [];
Y_temp = [];
for i=1:8
    duration_temp = [duration_temp; duration];
    Y_temp = [Y_temp Y];
end
%%
% Objective function: Max No of Proj. 
f = [zeros(1,85*15*8) ones(1,85*8)];
%%
% Objective function: Solving more proj. with higher WPIS
f = [zeros(1,85*8)];
f_weight = [];

for i=1:85
    len = deadline(i)-date_ready(i)+1;
    proj_weight = Y(:,i)';
    for j=date_ready(i):deadline(i)
        proj_weight(j) = WPIS(i);
    end
    f_weight = [f_weight proj_weight];
end
for c=1:8
    f = [f_weight f];
end
%%
% Model Constraints Setup
intcon = 1:length(f);  % integer constraint applied to all variables

% inequality constraints are set up below:
A1 = zeros(15*8,length(f));
for c=1:8
    for i=1:15
        for j=0:85-1
            A1(i+(c-1)*15,(c-1)*85*15+i+j*15)=1;
        end
    end
end
b1 = ones(15*8,1);

A2 = zeros(85,length(f));
for i=1:85
    for j=0:8-1
        A2(i,15*85*8+i+j*85)=1;
    end
end
b2 = ones(85,1);

A3 = zeros(8,length(f));
for c=1:8
    A3(c,15*85*(c-1)+1:15*85*c) = 1;
end
b3 = availability';

A = [A1;A2;A3];
b = [b1;b2;b3];

% equation constraints are set up below:
Aeq1 = zeros(85*8,length(f));
for i=1:85*8
    Aeq1(i,1+(i-1)*15:i*15)=1;  % column sum
    Aeq1(i,15*85*8+i)=-duration_temp(i); 
end
beq1 = zeros(85*8,1);

Aeq2 = zeros(85*8,length(f));
for i=1:85*8
    Aeq2(i,1+(i-1)*15:i*15)=Y_temp(:,i)';  % x_(ij) * y_(ij)
    Aeq2(i,15*85*8+i)=-duration_temp(i);
end
beq2 = zeros(85*8,1);

Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];

% lower and upper bounds are set up below:
lb = zeros(85*15*8+85*8,1);
ub = ones(85*15*8+85*8,1);

%%
% Utilizing problem solver
% if the objective function is the first:
x = intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub);
% if the objective function is the second:
options = optimoptions('intlinprog','MaxNodes', 2e4);
x = intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub,options);

%%
% Transformation: into a explicit way
X = zeros(15,85*8);
G = x(15*85*8+1:end);
for i=1:85*8
    X(1:15,i) = x(1+(i-1)*15:i*15);
end
X = round(X);  % X consists 8 X_i matraces
G = round(G);  % G consists 8 G_i matraces

result = zeros(15,8);
for k=1:8
    X_indiv = X(1:15,1+85*(k-1):85*k);
    for i=1:15
        for j=1:85
            if X_indiv(i,j)==1
                result(i,k) = j;
            end
        end
    end
end
result  % result is the work schedule of each crew during 15 working days
sum(G)  % the number of performed projects