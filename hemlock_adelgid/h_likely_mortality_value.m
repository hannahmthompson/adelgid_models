%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T. canadensis time series data analysis: 
% Determine threshold of likely hemlock mortality by 
% determining proportion of trees recorded at or below a 
% particular percent tips alive that are recorded to have
% died
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in time series data
hem = readtable('ha_data.xlsx','Sheet','H_time_series','Range','B2:OI12');
% determine and print total number of trees in data
n = width(hem);
fprintf('Total number of trees: %3.0f\n', n)

% determine and print total number of trees that died
l = 0;
for i = 1:n
    if any(hem.(i) == -5)
        l = l + 1;
    end
end
fprintf('Number of trees that died: %3.0f\n', l);

% determine and print number and proportion of trees
% that reach 0 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <=  0)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 5 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 5)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.05 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.05: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 10 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 10)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.1 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.1: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 15 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 15)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.15 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.15: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 20 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 20)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.2 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.2: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 25 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 25)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.25 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.25: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 30 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 30)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.30 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.30: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 35 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 35)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.35 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.35: %3.0f\n', k);

% determine and print number and proportion of trees
% that reach 40 percent tips alive that died
j = 0;
k = 0;
for i = 1:n
    if any(hem.(i) <= 40)
        k = k + 1;
        if any(hem.(i) == -5)
            j = j + 1;
        end
    end
end
fprintf('Proportion of trees that reach 0.40 that died: %1.4f\n', j/k);
fprintf('Number of trees that reached 0.40: %3.0f\n', k);
