clear
clc
load('dataset1.1.mat')
%%
Y = zeros(15,85);
for i=1:15  % from workday 1 to workday 15
    for j=1:85  % for each project j
        Y(date_ready(j):deadline(j),j)=1;
    end
end

f = [zeros(1,86*8)];
f_weight = [];

for i=1:85
    len = deadline(i)-date_ready(i)+1;
    proj_weight = Y(:,i)';
    for j=date_ready(i):deadline(i)
        proj_weight(j) = (WPIS(i));
    end
    f_weight = [f_weight proj_weight];
end
index=[];
for i=1:85
    index = [index 15*(i-1)+1:15*(i-1)+5];
end
f_weight(index)=[];
f_weight = [f_weight 10 10 0 0 0 0 0 0 0 0];
f = zeros(1,86*8);
for c=1:8
    f = [f_weight f];
end

Y = Y(6:15,:);
Y = [Y [1 1 0 0 0 0 0 0 0 0]'];

Y_scheduled = [1 2 3 6 7 8 9 10 11 14 15 16 17 21 23 24 26 32 33 34 38 40 44 46 47 50 51 54 55 57 61 63 64 65 67 71 75 78 80 85];
Y(:,Y_scheduled)=0;

availability=[0 10 0 0 10 0 0 0];
duration = [duration; 2];

duration_temp = [];
Y_temp = [];
for i=1:8
    duration_temp = [duration_temp; duration];
    Y_temp = [Y_temp Y];
end

%%

intcon = 1:length(f);

A1 = zeros(10*8,length(f));
for c=1:8
    for i=1:10
        for j=0:86-1
            A1(i+(c-1)*10,(c-1)*86*10+i+j*10)=1;
        end
    end
end
b1 = ones(10*8,1);

A2 = zeros(86,length(f));
for i=1:86
    for j=0:8-1
        A2(i,10*86*8+i+j*86)=1;
    end
end
b2 = ones(86,1);

A3 = zeros(8,length(f));
for c=1:8
    A3(c,10*86*(c-1)+1:10*86*c) = 1;
end
b3 = availability';

A = [A1;A2;A3];
b = [b1;b2;b3];

Aeq1 = zeros(86*8,length(f));
for i=1:86*8
    Aeq1(i,1+(i-1)*10:i*10)=1;  % column sum
    Aeq1(i,10*86*8+i)=-duration_temp(i); 
end
beq1 = zeros(86*8,1);

Aeq2 = zeros(86*8,length(f));
for i=1:86*8
    Aeq2(i,1+(i-1)*10:i*10)=Y_temp(:,i)';  % x_(ij) * y_(ij)
    Aeq2(i,10*86*8+i)=-duration_temp(i);
end
beq2 = zeros(86*8,1);

Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];

lb = zeros(86*10*8+86*8,1);
ub = ones(86*10*8+86*8,1);

x = intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub);
X = zeros(10,86*8);
G = x(10*86*8+1:end);
for i=1:86*8
    X(1:10,i) = x(1+(i-1)*10:i*10);
end

X = round(X)
G = round(G)

result = zeros(10,8);
for k=1:8
    X_indiv = X(1:10,1+86*(k-1):86*k);
    for i=1:10
        for j=1:86
            if X_indiv(i,j)==1
                result(i,k) = j;
            end
            
        end
    end
end
result
