% Helheim version of Hewitt et al. 2013 with some lines used in Stevens et al. 2018.
% Morlighem Helheim BedMachine3 80km up from terminus 
% with Helheim RACMO input 2007 DOY 1 to 365
% loop over the 2007 melt season. Starting condition is 2007 DOY 1. 
% LAS 25 March 2019 - this run is the initial test for Helheim geometry
% LAS 2 DEC 2019 - Add in GPS locations and save N, p_w at GPS + lake locations
% LAS 17 DEC 2019 - Add in Helheim RACMO daily runoff forcing
% LAS 03 MAR 2020 - Spatially-varying Ub based on MeASUreS July 2007 vels
% LAS 19 MAR 2020 - Ramp up lake forcing in 9 hrs, four input spots (correct), ub=spatially varying
% LAS 24 MAR 2020 - Keep DOY 227 RACMO on for entire time
% LAS 25 MAR 2020 - Start a flood on DOY67 (winter model)

format compact;
clear
oo.root = '';          % filename root
oo.fn = 'nevis_h22222_ubspatial_R67_lakerampM_4tiles_Ks100_s1e6H_repo';  % filename
oo.code = './nevis';   % code directory
addpath(oo.code);

%% parameters
% default parameters
[pd,oo] = nevis_defaults([],oo);    
% alter default parmaeters 
pd.c_e_reg2 = 0.01/1e3/9.81;        % elastic sheet thickness [m/Pa]
pd.N_reg2 = 1e4; % 1e3              % regularisation pressure for elastic sheet thickness 
pd.u_b = 500/pd.ty;                 % sliding speed [m/day] <-- this gets overwriten below with MeAsUrEs obs (western Margin = 100 m/day)
pd.sigma = 1e-6;                    % englacial void fraction default = 1e-3. loose = 1e-2. tight = 1e-4.
pd.h_r = 0.1; % .1                  % roughness height [m]
pd.l_r = 10; % 10                   % roughness length [m]
pd.l_c = 1000; % 10                 % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 1;                         % sheet permeability constant
pd.tau_b = 60e3; %60 kPa            % driving stress [Pa]
%pd.melt = (pd.G+(0))/pd.rho_w/pd.L;%   % geothermal heat only [m/s] 0.0059 m/yr
pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  % geothermal heat + frictional heating derived basal melt [m/s] 0.0262 m/yr
pd.meltinterior = ((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; % flux of basal melt up to the ~icedivide (200 km) [m2/s]
%pd.meltinterior = pd.melt.*2;      % flux of basal melt for one node outside fo the domain coming into the domain
% pd.gamma_cc = 7.8e-8;		    % beta, pressure dependent freezing pt.
% non-dimensionalise
[ps,pp] = nevis_nondimension(pd);   

%% grid and geometry input in m with [0,0] at the terminus centerline 2007 DOY 236
% Origin: lat = -38.1757246659987, lon = 66.3661134769539
load('morlighem_helheim_for_nevis_30km'); % load bedmap
dd = morlighem_helheim_for_nevis_30km; dd.S_30km = double(dd.S_30km);
dd.skip = 1; % [this uses lower resolution, though would be better to do 
% some spatial averaging to get lower resolution maps]
gg = nevis_grid(dd.X_30km(1:dd.skip:end,1)/ps.x,dd.Y_30km(1,1:dd.skip:end)/ps.x,oo); 
b = reshape(dd.B_30km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
s = reshape(dd.S_30km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);

%% mask with minimum ice thickness and seaward of calving front
H = max(s-b,0);
Hmin = 0/ps.z; 
X_km = reshape(gg.nx, gg.nIJ, 1);
Y_km = reshape(gg.ny, gg.nIJ, 1);
seaward_x = 0.00;
nout = find((H<=Hmin) | (X_km>=seaward_x) | (Y_km<=-0.5) | ((X_km>=-0.4) & (Y_km>=1.4))); % minimum ice thickness and calving front and a geometry control
% nout = find(H<=Hmin); % minimum ice thickness
% nout = unique([nout; gg.bdy.nlbdy; gg.bdy.nrbdy; gg.bdy.ntopbdy; gg.bdy.nbotbdy]); 
gg = nevis_mask(gg,nout); 
gg.n1m = gg.n1;             % label all edge nodes as boundary nodes for pressure

%% label boundary nodes
gg = nevis_label(gg,gg.n1m);
oo.adjust_boundaries = 1;   % enable option of changing conditions

%% plot grid
%nevis_plot_grid(gg); return;  % check to see what grid looks like

%% initialize variables
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);      % default initialisation
pd.k_f=0.9;                                    % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);  % initial pressure  k_f*phi_0
N=aa.phi_0-vv.phi;                             % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)
%vv.hs = pp.c16*aa.Ub./(pp.c16*aa.Ub*pp.c17+pp.c18*abs(aa.phi_0-vv.phi).^(pp.n_Glen-1).*(aa.phi_0-vv.phi)); % initial cavity size (equilibrium for given pressure)

% set Ub to a spatially varying field based on Measures data
Ub_load = load('measures_july_2007_nevis_30km.mat'); %previously colated w/ a min speed = 10 m/yr
aa.Ub = (Ub_load.measures_july_2007_nevis_30km./pd.ty)./ps.u_b; % need pd.ty and ps.ub because have already initialized

%% boundary conditions
% aa.phi = aa.phi_a(gg.nbdy)+k_factor*(aa.phi_0(gg.nbdy)-aa.phi_a(gg.nbdy));    % prescribed boundary pressure
% aa.phi_b = aa.phi_0;             % prescribed boundary pressure at overburden
aa.phi_b = max(aa.phi_0,aa.phi_a); % prescribed boundary pressure at overburden or atmospheric LAS 18 Nov. 2015

% define boundary nodes with flux pp.meltinterior for a large inland area
%aa.m(gg.nbdy) =  pp.meltinterior; % this is the ice margin boundary
%aa.m(gg.nout) =  pp.meltinterior; % this is the tundra
% find nodes where gg.nx <= -2.95 (the interior-most line of model domain inland boundary) 
X_km = reshape(gg.nx, gg.nIJ, 1); [row,col] = find(X_km <= -2.95);
aa.m(row) = pp.melt.*50;          % 50 times the basal melt flux (integrate over 50 more km inland of model domain)

%% alter initial channels
pd.S_0 = 0*0.1;                      % to keep and save in pd
S_0 = pd.S_0/ps.S;                  % initial channel cross-section (m^2) set to 0.01 originally
vv.Sx = S_0*vv.Sx.^0;
vv.Sy = S_0*vv.Sy.^0;
vv.Sr = S_0*vv.Ss.^0;
vv.Ss = S_0*vv.Sr.^0;

%% surface meltwater input
% % Used in Stevens et al. 2018
% % to include distributed runoff from eg RACMO, define function
% % runoff(t,gg) to return runoff (m/s) at time t (s), at each point on the
% % grid (ie runoff(t,gg) should return a vector of size gg.nIJ-by-1), then
% % include as:
% % pp.runoff_function = @(t) runoff(ps.t*t,gg)/ps.m; oo.runoff_function = 1;
load('runoff_2007_nevis_h30km_lake.mat');
% make DOY 67 everyday in the run 
runoff_2007_nevis_h30km(68,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(69,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(70,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(71,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(72,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(73,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(74,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(75,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(76,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(77,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(78,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(79,:) = runoff_2007_nevis_h30km(67,:);
runoff_2007_nevis_h30km(80,:) = runoff_2007_nevis_h30km(67,:);
% change this racmo forcing to be at a 0.1 day timestep for 228 through 240
runoff_here = runoff_2007_nevis_h30km(68:80,:);
runoff_here_20 = vertcat(repmat(runoff_here(1,:),20,1), repmat(runoff_here(2,:),20,1),...
                 repmat(runoff_here(3,:),20,1),repmat(runoff_here(4,:),20,1),...
                 repmat(runoff_here(5,:),20,1),repmat(runoff_here(6,:),20,1),...
                 repmat(runoff_here(7,:),20,1),repmat(runoff_here(8,:),20,1),...
                 repmat(runoff_here(9,:),20,1),repmat(runoff_here(10,:),20,1),...
                 repmat(runoff_here(11,:),20,1),repmat(runoff_here(12,:),20,1),...
                 repmat(runoff_here(13,:),20,1));
runoff_time = [68:0.05:80.95];

% RAMP up the lake drainage 
% full drainage over 9.4 hours 69.5 to 69.9
volume = 8.8888e+05; % spread out 10^10 mm w.e. over 150m x 150m tile over 12 hrs.
drainage_rate9_flat = volume/(0.4*4*2); % over 0.4 days and 4 tiles to get to Q = 290m3/s for the lake
% over 4 tiles evenly
runoff_here_20(31:32, 25658) =  drainage_rate9_flat.*0.25;
runoff_here_20(31:32, 25452) =  drainage_rate9_flat.*0.25;
runoff_here_20(31:32, 25657) =  drainage_rate9_flat.*0.25;
runoff_here_20(31:32, 25451) =  drainage_rate9_flat.*0.25;

runoff_here_20(33:34, 25658) =  drainage_rate9_flat.*0.75;
runoff_here_20(33:34, 25452) =  drainage_rate9_flat.*0.75;
runoff_here_20(33:34, 25657) =  drainage_rate9_flat.*0.75;
runoff_here_20(33:34, 25451) =  drainage_rate9_flat.*0.75;

runoff_here_20(35:36, 25658) =  drainage_rate9_flat.*1.125;
runoff_here_20(35:36, 25452) =  drainage_rate9_flat.*1.125;
runoff_here_20(35:36, 25657) =  drainage_rate9_flat.*1.125;
runoff_here_20(35:36, 25451) =  drainage_rate9_flat.*1.125;

runoff_here_20(37:38, 25658) =  drainage_rate9_flat.*1.5;
runoff_here_20(37:38, 25452) =  drainage_rate9_flat.*1.5;
runoff_here_20(37:38, 25657) =  drainage_rate9_flat.*1.5;
runoff_here_20(37:38, 25451) =  drainage_rate9_flat.*1.5;

% RACMO distributed input
oo.distributed_input = 1;              % turn on distributed input
oo.runoff_function = 1;                % turn on racmo distributed input
pp.runoff_function = @(t) runoff_lake(((t*ps.t)/pd.td),runoff_here_20,runoff_time)./ps.m;  % distributed input (m/sec)
% 
% % RACMO moulin input
oo.input_function = 0;          % turn on moulin input (m3/sec)
% pp.input_function = @(t) runoff_moulins(((t*ps.t)/pd.td),runoff_2009_nevis200,pp.sum_m,gg.Dx(1))./ps.m; % RACMO moulin input (m3/sec)

% % Input the lake drainage event at location: (-38.46, 66.46)
[pp.ni_lake,pp.ni_sum_lake] = nevis_moulins(-12900./ps.x,8694./ps.x,gg,oo); % lake location 25658
[pp.ni_lake2,pp.ni_sum_lake2] = nevis_moulins(-13050./ps.x,8694./ps.x,gg,oo); % lake location 25657
[pp.ni_lake3,pp.ni_sum_lake3] = nevis_moulins(-13050./ps.x,8544./ps.x,gg,oo); % lake location 25451
[pp.ni_lake4,pp.ni_sum_lake4] = nevis_moulins(-12750./ps.x,8544./ps.x,gg,oo); % lake location 25452

% GPS points as points of interest to save (treat them as moulins, but don't put water down them!)
load gps_2007_flood_nevis.mat
[pp.ni_gps,pp.ni_sum_gps] = nevis_moulins(gps_2007_flood_nevis(:,1)./ps.x,gps_2007_flood_nevis(:,2)./ps.x,gg,oo);

oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = [pp.ni_lake;pp.ni_gps]; % save first lake and GPS

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
% % for a lake drainage scenario
load('nevis_h22222_ubspatial_R67_lakerampM_4tiles_Ks100_s1e5H/0281.mat','vv','tt') % load a former timestep as an initial
% run at a 0.1 day timestep
[tt,vv,info] = nevis_timesteps([68:0.025:75]*pd.td/ps.t,vv,aa,pp,gg,oo,pd,ps);

% % for an annual run
% load('nevis_h22222_gps_racmo_restart/0365.mat','vv','tt') % load a former timestep as an initial
% run at a 1 day timestep
% [tt,vv,info] = nevis_timesteps([1:1:365]*pd.td/ps.t,vv,aa,pp,gg,oo,pd,ps);
