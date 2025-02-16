% Sample nested loop
x = 1:10;
figure; % Create a new figure
hold on; % Allow multiple plots on the same axes

for i = 1:3 % Outer loop
    for j = 1:2 % Inner loop
        % Simulate some data
        y = (i + j) * x; 
        
        % Generate a label for this simulation
        simLabel = sprintf('Simulation i=%d, j=%d', i, j);
        
        % Plot the data and set the DisplayName for the legend
        plot(x, y, 'DisplayName', simLabel);
    end
end

% Automatically generate the legend
legend('show');
hold off;
