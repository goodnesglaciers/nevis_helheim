% Helheim version of Hewitt et al. 2013 with some lines used in Stevens et al. 2018.
% Morlighem Helheim BedMachine3 80km up from terminus 
% with Helheim RACMO input 2007 DOY 1 to 365
% loop over the 2007 melt season. Starting condition is 2007 DOY 1. 
% LAS 25 March 2019 - this run is the initial test for Helheim geometry
% LAS 2 DEC 2019 - Add in GPS locations and save N, p_w at GPS + lake locations
% LAS 03 MAR 2020 - Spatially-varying Ub based on MeASUreS July 2007 vels

format compact;
clear
oo.root = '';           % filename root
oo.fn = 'nevis_h22222_ubspatial_R365';  % filename
oo.code = './nevis';   % code directory
addpath(oo.code);

%% parameters
% default parameters
[pd,oo] = nevis_defaults([],oo);    
% alter default parmaeters 
pd.c_e_reg2 = 0.01/1e3/9.81;        % elastic sheet thickness [m/Pa]
pd.N_reg2 = 1e4; % 1e3              % regularisation pressure for elastic sheet thickness 
pd.u_b = 500/pd.ty;                 % sliding speed [m/day] (western Margin = 100 m/day)
pd.sigma = 1e-4;                    % englacial void fraction
pd.h_r = 0.1; % .1                  % roughness height [m]
pd.l_r = 10; % 10                   % roughness length [m]
pd.l_c = 1000; % 10                 % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 1e-2;                      % sheet permeability constant
pd.tau_b = 60e3; %60 kPa            % driving stress [Pa]
%pd.melt = (pd.G+(0))/pd.rho_w/pd.L;%   % geothermal heat only [m/s] 0.0059 m/yr
pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  % geothermal heat + frictional heating derived basal melt [m/s] ~0.023 m/yr
pd.meltinterior = ((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; % flux of basal melt up to the ~icedivide (200 km) [m2/s]
%pd.meltinterior = pd.melt.*2;      % flux of basal melt for one node outside fo the domain coming into the domain
% pd.gamma_cc = 7.8e-8;		        % beta, pressure dependent freezing pt.
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
load('runoff_2007_nevis_h30km_lake.mat'); % load data (previously collated) in mm w.e. per day
runoff_2007_nevis_h30km(229,25658) = 30; % revert day 229 back to normal racmo level at position of lake
% RACMO distributed input
oo.distributed_input = 1;              % turn on distributed input
oo.runoff_function = 1;                % turn on racmo distributed input
pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2007_nevis_h30km)./ps.m;  % distributed input (m/sec)

% RACMO moulin input
oo.input_function = 0;          % turn on moulin input (m3/sec)
% pp.input_function = @(t) runoff_moulins(((t*ps.t)/pd.td),runoff_2009_nevis200,pp.sum_m,gg.Dx(1))./ps.m; % RACMO moulin input (m3/sec)

% GPS points as points of interest to save (treat them as moulins, but don't put water down them)
[pp.ni_lake,pp.ni_sum_lake] = nevis_moulins(-12900./ps.x,8694./ps.x,gg,oo); % lake location
load gps_2007_flood_nevis.mat
[pp.ni_gps,pp.ni_sum_gps] = nevis_moulins(gps_2007_flood_nevis(:,1)./ps.x,gps_2007_flood_nevis(:,2)./ps.x,gg,oo);

oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = [pp.ni_lake;pp.ni_gps]; % save both lake and GPS

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
% for an annual run
%load('nevis_h22222_ubspatial_sheetperm2_Ks1e2/0233.mat','vv','tt') % load a former timestep as an initial
%condition
[tt,vv,info] = nevis_timesteps([1:1:365]*pd.td/ps.t,vv,aa,pp,gg,oo,pd,ps);
