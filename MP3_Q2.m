clc; 
clear;
close all;

x=zeros (1, n);
dataset_1=zeros (n, 7);
x(1,1:31)=1.3+0.2*rand;

for k=31:n-1
x (1, k+1)=0.2* ((x(1, k-30))/ (1+x (1, k-30)^10))+0.9*x(1, k);
dataset_1 (k, 2:6)= [x(1, k-3) x(1, k-2) x(1, k-1) x(1, k) x(1, k+1)];
end
dataset (1:600, 2:6)=dataset_1 (201: 800, 2:6);
t=1:600;

figurel = figure ('Color', [1 1 1]); plot (t,x (201:800), 'Linewidth', 2)
[Number_training, ~]=size (dataset);
Rul=zeros (Number_training/2,6);
Rules_total=zeros (Number_training/2, 6);
﻿

for s=1:2
    switch s
        case 1
        num_membership_functions=7; c=linspace (0.5, 1.3,5);
        h=0.2;
        membership_functions=cell(num_membership_functions, 2);
        for k=1:num_membership_functions
            if k==1
                membership_functions {k, 1}= [0, 0, 0.3, 0.5];
                membership_functions {k, 2}='trapmf'; 
            elseif k==num_membership_functions 
                membership_functions{k, 1}=[1.3, 1.5, 1.8, 1.8];
                membership_functions {k, 2}='trapmf';
            else
                membership_functions {k, 1}=[c(k-1)-h, c(k-1), c(k-1)+h];
                membership_functions {k, 2}='trimf';
            end
        end
        case 2
        num_membership_functions=15; 
        c=linspace(0.3,1.5, 13); 
        h=0.1;
        membership_functions=cell(num_membership_functions, 2);
        for k=1:num_membership_functions
            if k==1
                membership_functions{k, 1}=[0, 0, 0.2, 0.3];
                membership_functions{k, 2}='trapmf';
            elseif k==num_membership_functions 
                membership_functions{k, 1}=[1.5, 1.6, 1.8, 1.8];
                membership_functions{k,2}='trapmf';
            else
                membership_functions{k, 1}=[c(k-1)-h, c(k-1), c(k-1)+h];
                membership_functions{k,2}='trimf';
            end
        end
    end


    vec_x=zeros (1, num_membership_functions); 
    vec=zeros (1,5);
    for t=1: Number_training
        dataset(t, 1)=t;
            for i=2:6
                x=dataset(t, i);
                    for j=1:num_membership_functions
                        if j==1
                        vec_x (1, j) = trapmf(x, membership_functions {1,1});
                        elseif j==num_membership_functions
                        vec_x (1, j)=trapmf (x, membership_functions{num_membership_functions, 1});
                        else
                        vec_x (1, j) = trimf (x, membership_functions {j,1});
                        end
                    end
                [valu_x, column_x]=max(vec_x);
                vec (1, i-1)=max (vec_x);
                Rules(t, i-1)=column_x;
                Rules(t, 6) =prod(vec);
                dataset (t,7) =prod(vec);
            end
     end
Rules_total(1, 1:6)=Rules(1,1:6);
i=1;
for t=2:Number_training
    m=zeros (1,1);
    for j=1:i
        m(1, j)=isequal(Rules(t, 1:4), Rules_total(j, 1:4));
        if m(1,j)==1 && Rules(t, 6)>=Rules_total (j,6) 
            Rules_total(j, 1:6)=Rules (t, 1:6);
        end
    end
    if sum (m)==0
        Rules_total(i+1, 1:6)=Rules(t, 1:6);
        i=i+1;

    end
end

disp('******************************')
disp(['Final rules for ', num2str(num_membership_functions),' membership functions for each input variables']) 
final_Rules=Rules_total(1:1, :);
Fisname='Prediction controller';
Fistype='mamdani';
Andmethod='prod';
Ormethod='max';
Impmethod='prod';
Aggmethod='max';
Defuzzmethod='centroid';
fis=newfis(Fisname, Fistype, Andmethod, Ormethod, Impmethod, Aggmethod, Defuzzmethod);
for num_input = 1:4
    fis = addInput(fis, [0.1 1.7], "Name", ['x', num2str(num_input)]);
end
fis = addOutput(fis,[0.1, 1.7],  'Name', 'x5');
for num_input = 1:4
    for input_Rul = 1:num_membership_functions 
        fis = addMF(fis, ['x', num2str(num_input)], membership_functions{input_Rul,2},membership_functions{input_Rul,1}, 'Name', ['A', num2str(input_Rul)]);
    end
end
for input_Rul = 1:num_membership_functions 
    fis = addMF(fis, 'x5',membership_functions{input_Rul, 2}, membership_functions{input_Rul, 1}, 'Name', ['MF_', num2str(input_Rul)]);
end
non_zero_rows = any(Rules_total(:, 1:5), 2); 
fis_Rules = ones(sum(non_zero_rows), 7);
fis_Rules(:, 1:6) = Rules_total(non_zero_rows, 1:6);
fis = addrule(fis, fis_Rules);
jadval_prediction=zeros(300,2);
f=1;
for i=301:600
    input=dataset(i, 2:6);
    output1=dataset(i, 6);
    x5=evalfis([input(1, 1); input(1, 2); input(1,3); input(1,4)], fis);
    jadval_prediction(f, :)= [f, x5];
    f=f+1;
end
figure;
plot(jadval_prediction(:,1),jadval_prediction(:,2), 'r-.', 'Linewidth', 2);
hold on;
plot(jadval_prediction(:,1),dataset(301: 600, 6), 'b', 'Linewidth', 2); 
legend('estimate value', 'real value')
end
inputVariableIndex = 1;

figure;
plotmf(fis, 'output', inputVariableIndex);
title(['Membership Functions for Input Variable ', num2str(inputVariableIndex)]);
