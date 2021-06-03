%%%%%privacy_marginal_three_graph
%%%% ER_graph BA_graph Community_graph and ER_graph is also the
%%%% Communication graph

clear 
clc
load('matlab.mat')

clear T s rho_m h_M alpha_1 alpha_2 alpha s_m

T = 16001;
%s  = 0.3;
s  = 0.3;
test_no = 10;
sigma =0.01; % noise variance
%%%% check 
%[rho_m_ER,h_M_ER,alpha_1_ER,alpha_2_ER,s_m_ER,alpha_ER] = parameter_check1(H_ER,G_Com2.L,N, s);

%[rho_m,h_M,alpha_1,alpha_2,s_m,alpha] = parameter_check1(H_Com,G_Com2.L,N, s);

%[rho_m,h_M,alpha_1,alpha_2,s_m,alpha] = parameter_check1(H_SF,G_Com2.L,N, s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot ER graph related results
[mean_z_ER,mean_z2_ER,upp_bound_ER] = Algorithm3(T,s,G_ER2,a_star_ER,H_ER,test_no,sigma,N,b);

[mean_z_SF,mean_z2_SF,upp_bound_SF] = Algorithm3(T,s,G_ER2,a_star_SF,H_SF,test_no,sigma,N,b);

[mean_z_Com,mean_z2_Com,upp_bound_Com] = Algorithm3(T,s,G_ER2,a_star_Com,H_Com,test_no,sigma,N,b);

%%%% plot the convergence result of node 1
figure(3)
i=1;
xx = 1:T;
yy = reshape(mean_z_ER(:,i,:),[N,T]);
%p1 = plot(xx,vecnorm(yy)/norm(a_star));
%plot(xx,vecnorm(yy),'r');
semilogy(xx,vecnorm(yy)/norm(a_star_ER),'r','LineWidth',2);
 hold on

yyy = reshape(mean_z_SF(:,i,:),[N,T]);
%p1 = plot(xx,vecnorm(yy)/norm(a_star));
 %plot(xx,vecnorm(yyy),'g');
 semilogy(xx,vecnorm(yyy)/norm(a_star_SF),'g','LineWidth',2);
 
hold on
yyyy = reshape(mean_z_Com(:,i,:),[N,T]);
%p1 = plot(xx,vecnorm(yy)/norm(a_star));
semilogy(xx,vecnorm(yyyy)/norm(a_star_Com),'b','LineWidth',2);
% plot(xx,vecnorm(yyyy));

xlabel('Iterations: $t$','interpreter','latex','FontSize',15)
ylabel('Normalized error: $||E[\textbf{x}_i(t)]-\textbf{a}^{\star}||/|\textbf{a}^{\star}|$','interpreter','latex','FontSize',15)
%ylabel('Normalized error: $||E[\textbf{x}_i(t)]-\textbf{a}^{\star}||$','interpreter','latex','font',15)
legend('Erdos-Renyi graph','Scale-free graph','Community graph','Location','northeast','interpreter','latex','FontSize',15)
grid on 
clear xx yy  yyy  yyyy 

%%%% plot the convergence result with average over all nodes
figure(4)
xx = 1:T;
yy = reshape(mean(mean_z2_ER,2),[1,T]);
yyy= reshape(mean(mean_z2_SF,2),[1,T]);
yyyy = reshape(mean(mean_z2_Com,2),[1,T]);

%yyy = vecnorm(yy).^2;
semilogy(xx,yy,'r','LineWidth',2);
hold on
semilogy(xx,yyy,'g','LineWidth',2);
hold on 
semilogy(xx,yyyy,'b','LineWidth',2);
p1 = semilogy(xx, upp_bound_ER*ones(1,T),'Color','r','LineStyle',':','LineWidth',2);
p2 = semilogy(xx, upp_bound_SF*ones(1,T),'Color','g','LineStyle','--','LineWidth',2);
p3 = semilogy(xx, upp_bound_Com*ones(1,T),'Color','b','LineStyle','-.','LineWidth',2);

xlabel('iterations: $t$','interpreter','latex','FontSize',15)
ylabel('Variance: $E[||\textbf{x}_i(t)-\textbf{a}^{\star}||^2]$','interpreter','latex','FontSize',15)
legend('Erdos-Renyi graph','Scale-free graph','Community graph','Location','northeast','interpreter','latex','FontSize',15)
grid on 
%clear xx yy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%