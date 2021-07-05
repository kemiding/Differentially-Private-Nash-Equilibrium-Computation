%%% construct linear quadratic game
clear 
close all
clc

% Call gsp_BOX
gsp_start
gsp_reset_seed(1)

% General parameters
N =30;    % number of nodes
dim = 2; % embedding dimentionality
paramplot.show_edges = 1; % parameter for graph plot
m0 =2;  
m =2;

%%%% construct the Community graph
param_Com.Nc = 5;
G_Com = gsp_community(N,param_Com);
subplot(1,3,3)
gsp_plot_graph(G_Com,paramplot);
title('Social influence Community graph','interpreter','latex','FontSize', 15)

%%% construct the scale-free graph 
G_SF = gsp_barabasi_albert( N, m0, m );
%G_SF = gsp_barabasi_albert( N);
%G_SF = gsp_barabasi_albert(N, m0, m);
%SF_c = gsp_compute_coordinates(G_SF , dim);
%%% randomly generate the coordination for graph plot
paramgraph.distribute = 1;
[ XCoords, YCoords] = create_coords(N,paramgraph.distribute);
%G.coords =[ XCoords, YCoords];

G_SF2 = gsp_update_coordinates(G_SF, [ XCoords, YCoords]);
subplot(1,3,2)
gsp_plot_graph(G_SF2,paramplot);
title('Social influence Scale-free graph','interpreter','latex','FontSize', 15)

%%%% ER graph
p = 0.5;
param.connected = 1;
G_ER = gsp_erdos_renyi(N,p,param);
%ER_c = gsp_compute_coordinates(G_ER , dim, 'isomap');
paramgraph.distribute = 1;
[ XCoords, YCoords] = create_coords(N,paramgraph.distribute);
G_ER2 = gsp_update_coordinates(G_ER, [ XCoords, YCoords]);
subplot(1,3,1)
gsp_plot_graph(G_ER2,paramplot);
title('Social influence ER graph','interpreter','latex','FontSize', 15)

%%%%% Normalized the graph weight matrix to satify Assumption 1
W = G_Com.W/max(max(G_Com.L));
G_Com2 = gsp_update_weights( G_Com, W );
clear W
W= G_SF.W/max(max(G_SF.L));
G_SF2 = gsp_update_weights( G_SF2, W );
clear W
W = G_ER.W/max(max(G_ER.L));
G_ER2 = gsp_update_weights( G_ER2, W );
clear W

%%% generate the marginal payoff and compute NE result
b = rand(N,1);
H_ER = (eye(N)-G_ER2.W)';
a_star_ER = (eye(N)- G_ER2.W)^(-1)*b;

H_Com = (eye(N)-G_Com2.W)';
a_star_Com = (eye(N)- G_Com2.W)^(-1)*b;

H_SF = (eye(N)-G_SF2.W)';
a_star_SF = (eye(N)- G_SF2.W)^(-1)*b;

%%% save the workspace for the simulation of Algorithm 1-3 
save('matlab.m')

%[rho_m,h_M,alpha_1,alpha_2,s_m,alpha] = parameter_check1(H_SF,G_Com2.L,N, s);

%%%% create coords function

function [ XCoords, YCoords] = create_coords(N, distribute)

% TODO: VECTORIZE!!!!
XCoords = zeros(N,1);
YCoords = zeros(N,1);
if distribute
    mdim = ceil(sqrt(N));
    ind = 1;
    for ii = 0:mdim-1
        for jj=0:mdim-1
            if ind<=N
                XCoords(ind) = 1/mdim*rand(1)+ii*1/mdim;
                YCoords(ind) = 1/mdim*rand(1)+jj*1/mdim;
            end
            ind = ind+1;
        end
    end
else
    % take random coordinates in a 1 by 1 square
    XCoords = rand(N,1);
    YCoords = rand(N,1);
end


end


