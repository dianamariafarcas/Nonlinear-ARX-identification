%%
clc
clear
close all
load('iddata-13.mat')

idU = id.InputData;  % Identification input data
idY = id.OutputData;  % Identification output data
valU = val.InputData;  % Validation input data
valY = val.OutputData;  % Validation output data

%%
% figure
% subplot(211)
% plot(idU);
% subtitle('Identification Input');
% grid on
% 
% subplot(212)
% plot(idY);
% subtitle('Identification Output');
% grid on
% 
% figure
% subplot(211)
% plot(valU);
% subtitle('Validation Input');
% grid on
% 
% subplot(212)
% plot(valY);
% subtitle('Validation Output');
% grid on




%% 

na_max = 10; % Maximum model order
nb_max = 10; % Maximum model order
m_max = 3; % Maximum allowed power
MSEs_val = [];
MSE_min_pred_val = 10;  % Minimum MSE for prediction on validation data
MSE_min_sim_val = 10;  % Minimum MSE for simulation on validation data
MSE_min_sim_id = 10;  % Minimum MSE for simulation on identification data
MSE_min_pred_id = 10;  % Minimum MSE for prediction on identification data

for m = 1 : m_max
    k = 1;
    for na = 1 : na_max
        for nb = 1 : nb_max

            global array;  % Global variable to store powers
            current = zeros(1, na + nb); % Current power vector initialization
            array = [];  % Reset the array for each configuration
            backtrack(1, 0, na + nb, m, current);  % Generation of combinations of powers
            powers = array; % Store all power combinations
                
            clear global array;

            % Prediction
            delayed_outputs_id = delayed(na, nb, idU, idY);  % Create delayed input-output data for identification
                
            R_id = calculateR(na, nb, powers, delayed_outputs_id);  % Compute the regression matrix for identification data
                
            theta = R_id \ idY;

            delayed_outputs_val = delayed(na, nb, valU, valY);  % Create delayed input-output data for validation
    
            R_val = calculateR(na, nb, powers, delayed_outputs_val);  % Compute the regression matrix for validation data

            % Prediction
    
            Yh_pred_id = R_id * theta;  % Predicted output for identification data
            Yh_pred_val = R_val * theta;  % Predicted output for validation data
            
            % Simulation

            Yh_sim_id = Yh_sim(na, nb, idU, theta, powers);  % Simulated output for identification data            
            Yh_sim_val = Yh_sim(na, nb, valU, theta, powers);  % Simulated output for validation data
            
            % MSE for prediction and simulation on identification data 
           
            MSEs_id_pred(k, m) = sum((idY - Yh_pred_id) .^ 2) / length(idY);  % MSE for prediction on identification data
            MSEs_id_sim(k, m) = sum((idY - Yh_sim_id) .^ 2) / length(idY);  % MSE for simulation on identification data
            
            if MSEs_id_pred(k,m) < MSE_min_pred_id
                MSE_min_pred_id = MSEs_id_pred(k,m);
                m_min_pred_id = m;
                na_min_pred_id = na;
                nb_min_pred_id = nb;
                Yh_pred_optim_id = Yh_pred_id;              
            end

            if MSEs_id_sim(k,m) < MSE_min_sim_id
                MSE_min_sim_id = MSEs_id_sim(k,m);
                m_min_sim_id = m;
                na_min_sim_id = na;
                nb_min_sim_id = nb;
                Yh_sim_optim_id = Yh_sim_id;
            end

            % MSE for prediction and simulation on validation data 

            MSEs_val_pred(k, m) = sum((valY - Yh_pred_val) .^ 2) / length(valY);
            MSEs_val_sim(k, m) = sum((valY - Yh_sim_val) .^ 2) / length(valY);


            if MSEs_val_pred(k,m) < MSE_min_pred_val
                MSE_min_pred_val = MSEs_val_pred(k,m);
                m_min_pred_val = m;
                na_min_pred_val = na;
                nb_min_pred_val = nb;
                Yh_pred_optim_val = Yh_pred_val;
                
            end

            if MSEs_val_sim(k,m) < MSE_min_sim_val
                MSE_min_sim_val = MSEs_val_sim(k,m);
                m_min_sim_val = m;
                na_min_sim_val = na;
                nb_min_sim_val = nb;
                Yh_sim_optim_val = Yh_sim_val;
            end

            k = k+1; 
        end
    end
end  



%%

% Graphical representation on identification

yPredId = iddata(Yh_pred_optim_id, idU, id.Ts);
ySimId = iddata(Yh_sim_optim_id, idU, id.Ts);
IdY1= iddata(idY, idU, id.Ts);
figure;
compare(IdY1, yPredId, ySimId);
title('Comparison between the original, predicted and simulated outputs on identification');


% Graphical representation on validation

yPredVal = iddata(Yh_pred_optim_val, valU, val.Ts);
ySimVal = iddata(Yh_sim_optim_val, valU, val.Ts);
valY1 = iddata(valY, valU, val.Ts);
figure;
compare(valY1, yPredVal, ySimVal);
title('Comparison between the original, predicted and simulated outputs on validation');

%% MSE plots for prediction on identification

clc
MSE1 = convert(MSEs_id_pred(:, 1), na_max, nb_max);
figure
subplot(131)
mesh(1 : 10, 1 : 10, MSE1);
title('MSEs_{id}_{ pred} for m = 1');

MSE2 = convert(MSEs_id_pred(:,2), na_max, nb_max);
subplot(132)
mesh(1 : 10, 1 : 10, MSE2);
title('MSEs_{id}_{ pred} for m = 2');

MSE3 = convert(MSEs_id_pred(:,3), na_max, nb_max);
subplot(133)
mesh(1 : 10, 1 : 10, MSE3);
title('MSEs_{id}_{ pred} for m = 3');

%% MSE plots for simulation on identification

clc
MSE4 = convert(MSEs_id_sim(:, 1), na_max, nb_max);
figure
subplot(131)
mesh(1 : 10, 1 : 10, MSE4);
title('MSEs_{id}_{ sim} for m = 1');

MSE5 = convert(MSEs_id_sim(:,2), na_max, nb_max);
subplot(132)
mesh(1 : 10, 1 : 10, MSE5);
title('MSEs_{id}_{ sim} for m = 2');

MSE6 = convert(MSEs_id_sim(:,3), na_max, nb_max);
subplot(133)
mesh(1 : 10, 1 : 10, MSE6);
title('MSEs_{id}_{ sim} for m = 3');

%% MSE plots for prediction on validation

clc
MSE7 = convert(MSEs_val_pred(:, 1), na_max, nb_max);
figure
subplot(131)
mesh(1 : 10, 1 : 10, MSE7);
title('MSEs_{val}_{ pred} for m = 1');

MSE8 = convert(MSEs_val_pred(:,2), na_max, nb_max);
subplot(132)
mesh(1 : 10, 1 : 10, MSE8);
title('MSEs_{val}_{ pred} for m = 2');

MSE9 = convert(MSEs_val_pred(:,3), na_max, nb_max);
subplot(133)
mesh(1 : 10, 1 : 10, MSE9);
title('MSEs_{val}_{ pred} for m = 3');

%% MSE plots for simulation on validation

clc
MSE10 = convert(MSEs_val_sim(:, 1), na_max, nb_max);
figure
subplot(131)
mesh(1 : 10, 1 : 10, MSE10);
title('MSEs_{val}_{ sim} for m = 1');

MSE11 = convert(MSEs_val_sim(:, 2), na_max, nb_max);
subplot(132)
mesh(1 : 10, 1 : 10, MSE11);
title('MSEs_{val}_{ sim} for m = 2');

MSE12 = convert(MSEs_val_sim(:, 3), na_max, nb_max);
subplot(133)
mesh(1 : 10, 1 : 10, MSE12);
title('MSEs_{val}_{ sim} for m = 3');

%% Plots used in the presentation

figure
subplot(121)
mse1 = convert(MSEs_id_pred(:, m_min_pred_id), na_max, nb_max);
mesh(1 : 10, 1 : 10, mse1)
grid on;
title(['MSEs_{id}_{ pred} for m = ' num2str(m_min_pred_id)]);

subplot(122)
mse2 = convert(MSEs_val_pred(:, m_min_pred_val), na_max, nb_max);
mesh(1 : 10, 1 : 10, mse2)
grid on;
title(['MSEs_{val}_{ pred} for m = ' num2str(m_min_pred_val)]);

figure
subplot(121)
mse3 = convert(MSEs_id_sim(:, m_min_sim_id), na_max, nb_max);
mesh(1 : 10, 1 : 10, mse3) 
grid on;
title(['MSEs_{id}_{ sim} for m = ' num2str(m_min_sim_id)]);

subplot(122)
mse4 = convert(MSEs_val_sim(:, m_min_pred_val), na_max, nb_max);
mesh(1 : 10, 1 : 10, mse4)
grid on;
title(['MSEs_{val}_{ sim} for m = ' num2str(m_min_sim_val)]);

%% Generation of all combinations of powers that sum up to m

function backtrack(index, currentSum, length, m, current)
    global array; % Global variable used to store all the valid power combinations
    
    if currentSum > m
        return;
    end
        
    if index > length % Verify if all the elements in the power vector are set
        array(end + 1,:) = current; % Append the current vector to the array
        return;
    end 

    for digit = 0 : m
        if currentSum + digit <= m
            current(index) = digit; % Assign the power for the current index
            backtrack(index + 1, currentSum + digit, length, m, current);
        end
    end
end

%% Computation of the regressors matrix

function [R] = calculateR(na, nb, powers, delayed)
    Nid = size(delayed, 1);
    R = zeros(Nid, size(powers,1));

    for i = 1 : Nid % Loop over the identification data 
        for j = 1 : size(powers,1) % Loop over each row of powers
            sum = 1;
            for k = 1 : na + nb % Loop over all delayed input-output terms
                if powers(j, k) ~= 0
                    sum = sum * (delayed(i, k) ^ powers(j, k));
                end
            end

            R(i, j) = sum;

        end    
    end
end

%% Generation of the delayed outputs

function delayed = delayed(na, nb, input, output)
    Nid = size(input, 1);
    delayed = zeros(Nid, na + nb);

    for i = 1 : Nid
        for j = 1 : na
            if i - j <= 0
                delayed(i, j) = 0;
            else
                delayed(i, j) = input(i - j);
            end
        end 
    
        for j = 1 : nb
            if i - j <= 0
                delayed(i, j + na) = 0;
            else
                delayed(i, j + na) = output(i - j);
            end
        end
    end
end

%% Computation of the simulated outputs

function Yh_sim = Yh_sim(na, nb, input, theta, powers)
    Yh_sim = zeros(size(input,1),1);
    delayed = zeros(1, na+nb);
    R = zeros(1,size(powers,1));
    Nid = size(input,1);

    for i = 1 : Nid
        for j = 1 : na
            if i - j <= 0
                delayed(j) = 0;
            else
                delayed(j) = input(i - j);
            end
        end
    
        for j = 1 : nb
            if i - j <= 0
                delayed(j + na) = 0;
            else
                delayed(j + na) = Yh_sim(i - j);
            end
        end

        for j = 1 : size(powers,1)
            sum = 1;
            for k = 1 : na + nb
                if powers(j, k) ~= 0
                    sum = sum * (delayed(k) ^ powers(j, k));
                end
            end

            R(j) = sum;
        end
        Yh_sim(i) = R * theta;
    end
end

%% Function for plotting the MSEs

function MSE = convert(mse, na, nb)
    MSE = zeros(na, nb);
    k = 1;
    mse;
    for i = 1 : na
        for j = 1 : nb
            MSE(i, j) = mse(k);
            k = k + 1;
        end
    end
end