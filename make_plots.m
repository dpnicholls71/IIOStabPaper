% make_plots.m
%
% MATLAB script to plot relative error versus perturbation order.

clear all; close all;

SavePlots = 0;

% Interface: 1->U, 2->V_u, 3->V_ell, 4->W
Interface = 1;
% SumType: 1->Taylor, 2->Pade
SumType = 1;
M = SumType;

num_files = input('Number of datafiles? ');

Eps_full = zeros(num_files,1);
nplot = [];
relerr = [];

for j=1:num_files
  % three.mat
  filename = input('Datafile name? ','s');
  load(filename);
  Eps_full(j) = Eps;
  
  if(j==1)
    line_str = 'b-o';
  elseif(j==2)
    line_str = 'g-*';
  elseif(j==3)
    line_str = 'r-<';
  elseif(j==4)
    line_str = 'c->';
  elseif(j==5)
    line_str = 'y-d';
  else
    line_str = 'k-^';
  end
  disp_name = sprintf('$\\varepsilon = %g$',Eps);

  if(Interface==1)
    nplot = [nplot nplot_U];
    relerr = [relerr relerr_U(:,M)];
  elseif(Interface==2)
    nplot = [nplot nplot_V_u];
    relerr = [relerr relerr_V_u(:,M)];
  elseif(Interface==3)
    nplot = [nplot nplot_V_ell];
    relerr = [relerr relerr_V_ell(:,M)];
  else
    nplot = [nplot nplot_W];
    relerr = [relerr relerr_W(:,M)];
  end
  figure(101); hh101 = gca;
  semilogy(nplot(:,j),relerr(:,j),line_str,...
      'LineWidth',2,'MarkerSize',8,...
      'DisplayName',disp_name);
  axis([min(nplot(:,j)) max(nplot(:,j)) 1e-16 max(relerr(:,j))]);
  title('Relative Error versus $N$','Interpreter','latex');
  xlabel('$N$','Interpreter','latex');
  ylabel('Relative Error','Interpreter','latex');
  set(hh101,'fontsize',16);
  hold on;
end

for n=0:2:N
  if(n==0)
    line_str = 'b-o';
  elseif(n==2)
    line_str = 'g-*';
  elseif(n==4)
    line_str = 'r-<';
  elseif(n==6)
    line_str = 'c->';
  elseif(n==8)
    line_str = 'y-d';
  else
    line_str = 'k-^';
  end
  disp_name = sprintf('$N = %d$',n);
  
  figure(102); hh102 = gca;
  loglog(Eps_full,relerr(n+1,:),line_str,...
      'LineWidth',2,'MarkerSize',8,...
      'DisplayName',disp_name);
  axis([min(Eps_full) max(Eps_full) 1e-16 max(relerr(:,j))])
  title('Relative Error versus $\\\varepsilon$','interpreter','latex');
  xlabel('$\\\varepsilon$','interpreter','latex');
  ylabel('Relative Error','interpreter','latex');
  set(hh102,'fontsize',16);
  hold on;
end

figure(101);
ll = legend('show','Location','northeast');
set(ll,'interpreter','latex');
figure(102);
ll = legend('show','Location','northeast');
set(ll,'interpreter','latex');

if(SavePlots==1)
  saveas(hh101,'conv_N','epsc');
  saveas(hh102,'conv_Eps','epsc');
end