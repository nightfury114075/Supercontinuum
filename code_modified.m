clc
close all
clear all
format long
tic

%%%%%%%% Data for DL
w_clad = 7;
h_clad = 5;
w_core = 4;
h_core = 0.5;
power = 200;         		      	% peak power of input [W]2
lambdapump = 1550;           	    % pump wavelength [nm]
FWHM = 75e-3;                          % pulse width in FWHM [ps] 



Aeff = 2.39405580745888E-12; 	  	% effective mode area [m^2]

%%Dispersion Engineering                      
%%==============X===================================================================================================
N_betas = 10;                     	% number of higer order dispersion 
N_pump = 12;                      	% index number of pumpwavelength1
c = 3e-7; 			 	         % Unit of velocity (c) of light is km/ps
del_lambda = 50*1e-12;
lambda = (1000*1e-12:del_lambda:5000*1e-12);	
% Unit of wavelength (lambda) in km
N = length(lambda);

n_eff =[
3.069313760736980
3.064372188544280
3.059612561078730
3.054988636635420
3.050465030481140
3.046014414403410
3.041615524206640
3.037251717878630
3.032909916294810
3.028579814255680
3.024253285549110
3.019923929245370
3.015586720128080
3.011237736819420
3.006873948503980
3.002493046293850
2.998093308918330
2.993673495033750
2.989232756343810
2.984770567110300
2.980286666662740
2.975781012284430
2.971253740431620
2.966705134683170
2.962135599154830
2.957545636372570
2.952935828801300
2.948306823382770
2.943659318560950
2.938994053371010
2.934311798246010
2.929613347258130
2.924899511561030
2.920171113840610
2.915428983614390
2.910673953246180
2.905906854565110
2.901128515995730
2.896339760121030
2.891541401612470
2.886734245471610
2.881919085536140
2.877096703210870
2.872267866389870
2.867433328541260
2.862593827930740
2.857750086963130
2.852902811624940
2.848052691013220
2.843200396938540
2.838346583591830
2.833491887266350
2.828636926127920
2.823782300027170
2.818928590349320
2.814076359897250
2.809226152804910
2.804378494478290
2.799533891562260
2.794692831931410
2.789855784703830
2.785023200277010
2.780195510384970
2.775373128176360
2.770556448312910
2.765745847088120
2.760941682565890
2.756144294738780
2.751354005705860
2.746571119869650
2.741795924152130
2.737028688229150
2.732269664783050
2.727519089772950
2.722777182721960
2.718044147020890
2.713320170247520
2.708605424500670
2.703900066748220
2.699204239188010
2.694518069620630

];

[betas,D]=dispersion(lambda,n_eff,N,del_lambda,N_pump);
lambda=lambda*1e3;
betas=betas(1:N_betas);

%%%%%Figure-1
figure(1)
plot (lambda(2:N-1)/1e-6,D(2:N-1),'linewidth',2,'Color','r');
xlabel ('Wavelength [\mum]');
ylabel ('D [ps/nm/km]');
%xlim([0.9,1.6]);
grid on


n2 = 2.8e-17;                           % nonlinear refrective index [m^2/W] for Si2N
loss = 600;                               % loss [dB/m]
gamma = 2*pi*n2/(1e-9*lambdapump*Aeff);	% nonlinear coefficient [1/W/m]
nt = 2^13;                   			    % number of grid points (FFT)
%dt = 2.76e-3;                              % time step [ps]
T = 20;                                     % Time window [ps]
dt = T/nt;                                 	% time step
%T = nt*dt;                                 % width of time window [ps]
c = 3e8*1e9/1e12;                           % speed of light [nm/ps]
 
f0 = c/lambdapump;                          % pump frequency [THz] [pump frequency = 193.55 THz at 1550 nm]
w0 = 2*pi*f0;                               % angular pump frequency
t = linspace(-T/2, T/2, nt); 				% time grid
V = 2*pi*(-nt/2:nt/2-1)'/(nt*dt);   		% frequency grid

%[W_length]=Waveguide_length();
%dlength = W_length(1);                 	% device length [m]
dlength = 30e-3;

% === input pulse and parameters	

chirp = 0;
mshape = 0;                									
T0 = FWHM/(2*acosh(sqrt(2)));

if mshape == 0
    A =sqrt(power)*sech(t/T0).*exp(-1i*chirp*t.^2/(2*T0.^2));	% input field [W^(1/2)]
else
    A =exp(-0.5*(1+1i*chirp).*(t/T0).^(2*mshape));
end

LD = T0^2/abs(betas(2));				% dispersion length [m] 
LNL = 1/(gamma*power);					% non-linear length [m]
N = sqrt(LD/LNL);						% soliton order
Lfiss = (LD/N);							% soliton fission length [m]  
                				

fr = 0.013;						 % fractional Raman contribution forGeAsSe
tau1 = 21.3e-3;						% duration [ps]
tau2 = 195e-3;						    % duration [ps]
RT = ((tau1^2+tau2^2)/(tau1*tau2^2)).*sin(t/tau1).*exp(-t/tau2);
RT(t<0) = 0;							% heaviside step function
RT = RT/trapz(t,RT);  					% normalise RT to unit integral
Tr = 0;          						

% === simulation parameters
nplot = 240;							% number of length steps to save field at 
   						
%----propagation field
[Z, AT, AW, W] = gnlse(t, V, A, w0, gamma, betas, loss, fr, Tr, dlength, nplot);

IA = abs(A).^2;
IAf = abs(fftshift(ifft(A)).*(nt*dt)/sqrt(2*pi)).^2;
mIAf = max(max(IAf));
lIAf = 10*log10(IAf/mIAf);



IW = abs(AW).^2; 						% spectral intensity
mIW = max(max(IW));        				% max value, for scaling plot
lIW = 10*log10(IW/mIW);

IT = abs(AT).^2; 						% temporal intensity
mIT = max(max(IT));       				% max value, for scaling plot
lIT = 10*log10(IT/mIT);

WL = 2*pi*c./W;							% wavelength grid 
iis = (WL>300 & WL<20000);

% === plot input and final pulse spectral shape
figure(2)
plot(WL(iis)/1000,lIW(240,iis),'-b','linewidth',3);
xlabel('Wavelength [\mum]','FontSize',16);
ylabel('Spectral Power [dB]','FontSize',16);
set(gca,'FontSize',16);
xlim([0.6,7]);
ylim([-40,0]);
grid

% % === plot log scale spectral and temporal density
figure(3)
pcolor(WL(iis)/1000,Z*1000,lIW(:,iis));			 		% plot as pseudocolor map
%surf(WL(iis),Z,lIW(:,iis)),'MeshStyle', 'col', 'EdgeColor', 'none');
set(gca,'CLim',[-40 0]);
xlim([0.75,7]); 
shading interp;
colormap(jet.^3)
%title('Spectral density','FontSize',16);
xlabel({'Wavelength [\mum]',''},'FontSize',26); 
ylabel('Distance [mm]','FontSize',26);
set(gca,'FontSize',22);
view(2);

figure(4)
pcolor(t,Z*1000,lIT);						% plot as pseudocolor map
%surf(t,Z,lIT),'MeshStyle', 'col', 'EdgeColor', 'none');         
set(gca,'CLim',[-40 0]);
xlim([-.5,2.2]);   
shading interp;
colormap(jet.^3)
%title('Temporal Density');
xlabel({'Delay [ps]',''},'FontSize',26); 
ylabel('Distance [mm]','FontSize',26);
set(gca,'FontSize',22);
view(2);

toc;
cputime=toc;
disp('CPU time:'), disp(cputime);
X1=WL(iis)/1000;
Y1= lIW(240,iis);
Y1=Y1';

AT_output_1=AT;
AW_1 = AW;


% Define the matrix variable (example: a 240x8192 matrix)
data = lIW; % Replace this with your actual data

% Extract the last row (240th row)
lastRow = data(240, :);

% Combine single values and last row into a single row vector
dataRow = [w_clad, h_clad, w_core, h_core, power, lambdapump, FWHM, lastRow];

% Define the CSV file name
csvFileName = 'DL_data.csv';

% Check if the CSV file already exists
if isfile(csvFileName)
    % Append the new row to the existing file
    fileID = fopen(csvFileName, 'a'); % Open in append mode
else
    % Create a new file and write the header
    fileID = fopen(csvFileName, 'w');
    
    % Write the header row
    header = ['width_clad, height_clad ,width, height, Power,LambdaPump,FWHM,' sprintf('Col%d,', 1:8191) 'Col8192\n'];
    fprintf(fileID, header);
end

% Write the data row to the CSV file
fprintf(fileID, ['%f,%f,%f,%f,%f,%f,%f,' repmat('%f,', 1, size(lastRow, 2) - 1) '%f\n'], dataRow);

% Close the file
fclose(fileID);

disp('Simulation data saved successfully.');


% Define the Excel file name for D values
dispersionFileName = 'Dispersion_data.csv';

% Create a dynamic header based on the simulation variables
header = sprintf('wcore_%g_hcore_%g_wclad_%g_hclad_%g', w_core, h_core, w_clad, h_clad);

% Check if the file already exists
if isfile(dispersionFileName)
    % Read the existing data and preserve the original headers
    opts = detectImportOptions(dispersionFileName);
    opts.VariableNamesLine = 1; % Ensure the headers are read correctly
    existingData = readtable(dispersionFileName, opts);
else
    % Create an empty table if the file doesn't exist
    existingData = table();
end

% Ensure D is a column vector
D = D(:);

% Create a new table with the current D values and dynamic header
newTable = array2table(D, 'VariableNames', {header});

% Append the new column to the existing table
updatedTable = [existingData, newTable];

% Write the updated table to the CSV file
writetable(updatedTable, dispersionFileName);

disp('Dispersion data saved successfully to CSV with dynamic header.');
