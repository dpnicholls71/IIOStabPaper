function [U_n,Vu_n,Vell_n,W_n] ...
  = threelayer_iio_tfe_helmholtz(zeta_n,psi_n,theta_n,mu_n,...
  hbar,eta,fu,fell,tau2,sigma2,p,alphap,gamma_up,gamma_vp,gamma_wp,...
  eep,eem,a,b,Nx,Nz,N)

[Dz,z] = cheb(Nz);

% hbar is the "half-height" of the layer

U_n = zeros(Nx,N+1);
Vu_n = zeros(Nx,N+1);
Vell_n = zeros(Nx,N+1);
W_n = zeros(Nx,N+1);
Q_U_rs = zeros(Nx,N+1,N+1);
Ru_V_rs = zeros(Nx,N+1,N+1);
Rell_V_rs = zeros(Nx,N+1,N+1);
S_W_rs = zeros(Nx,N+1,N+1);

Q0 = (-1i*gamma_up+1i*eta)./(-1i*gamma_up-1i*eta);
S0 = (-1i*gamma_wp+1i*eta)./(-1i*gamma_wp-1i*eta);
sh = sinh(1i*gamma_vp*hbar)./(1i*gamma_vp);
ch = cosh(1i*gamma_vp*hbar);
aa = -(gamma_vp.^2).*sh - 1i*eta*ch;
bb = ch - 1i*eta*sh;
Ruu0 = (1.0/2.0)*( conj(aa)./aa + conj(bb)./bb );
Ruell0 = (1.0/2.0)*( conj(aa)./aa - conj(bb)./bb );
Rellu0 = Ruell0;
Rellell0 = Ruu0;

for n=0:N
  zetahat = fft(eem.*zeta_n(:,n+1));
  psihat = fft(eem.*psi_n(:,n+1));
  thetahat = fft(eem.*theta_n(:,n+1));
  muhat = fft(eem.*mu_n(:,n+1));
  Y_Uhat = -2*1i*eta*zetahat;
  Y_Vuhat = -2*psihat;
  Y_Vellhat = -2*1i*eta*thetahat;
  Y_What = -2*muhat;

  for m=0:n-1
    YY = -Q_U_rs(:,n-m+1,m+1) + Ru_V_rs(:,n-m+1,m+1);
    Y_Uhat = Y_Uhat - fft(eem.*YY);
    YY = Q_U_rs(:,n-m+1,m+1) + tau2*Ru_V_rs(:,n-m+1,m+1);
    Y_Vuhat = Y_Vuhat - fft(eem.*YY);
    YY = -Rell_V_rs(:,n-m+1,m+1) + S_W_rs(:,n-m+1,m+1);
    Y_Vellhat = Y_Vellhat - fft(eem.*YY);
    YY = Rell_V_rs(:,n-m+1,m+1) + sigma2*S_W_rs(:,n-m+1,m+1);
    Y_What = Y_What - fft(eem.*YY);
  end

  Uhat = zeros(Nx,1);
  Vuhat = zeros(Nx,1);
  Vellhat = zeros(Nx,1);
  What = zeros(Nx,1);

  for j=1:Nx
    MM = [1.0-Q0(j) -1.0+Ruu0(j) Ruell0(j) 0.0;...
        1.0+Q0(j) tau2*(1.0+Ruu0(j)) tau2*Ruell0(j) 0.0;...
        0.0 -Rellu0(j) 1.0-Rellell0(j) -1.0+S0(j);...
        0.0 Rellu0(j) 1.0+Rellell0(j) sigma2*(1.0+S0(j))];
    bb = [Y_Uhat(j);Y_Vuhat(j);Y_Vellhat(j);Y_What(j)];
    xx = MM\bb;
    Uhat(j) = xx(1);
    Vuhat(j) = xx(2);
    Vellhat(j) = xx(3);
    What(j) = xx(4);
  end

  U_n(:,n+1) = eep.*ifft(Uhat);
  Vu_n(:,n+1) = eep.*ifft(Vuhat);
  Vell_n(:,n+1) = eep.*ifft(Vellhat);
  W_n(:,n+1) = eep.*ifft(What);

  % Compute and store Q_r[U_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  xi(:,0+1) = U_n(:,s+1);
  un = field_tfe_helmholtz_upper(xi,eta,fu,...
      p,alphap,gamma_up,eep,eem,Dz,a,Nx,Nz,N-s);
  Qn = iio_tfe_helmholtz_upper(un,eta,fu,...
      p,alphap,gamma_up,eep,eem,Dz,a,Nx,Nz,N-s);
  for r=0:N-s
    Q_U_rs(:,r+1,s+1) = Qn(:,r+1);
  end

  % Compute and store Ru_r[V_s], Rell_r[V_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  zeta = zeros(Nx,N-s+1);
  xi(:,0+1) = Vu_n(:,s+1);
  zeta(:,0+1) = Vell_n(:,s+1);
  vn = field_tfe_helmholtz_middle(xi,zeta,hbar,eta,fu,fell,...
      p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
  [Run,Relln] = iio_tfe_helmholtz_middle(vn,hbar,eta,fu,fell,...
      p,alphap,gamma_vp,eep,eem,Dz,a,Nx,Nz,N-s);
  for r=0:N-s
    Ru_V_rs(:,r+1,s+1) = Run(:,r+1);
    Rell_V_rs(:,r+1,s+1) = Relln(:,r+1);
  end

  % Compute and store S_r[W_s]
  s = n;
  xi = zeros(Nx,N-s+1);
  xi(:,0+1) = W_n(:,s+1);
  wn = field_tfe_helmholtz_lower(xi,eta,fell,...
      p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
  Sn = iio_tfe_helmholtz_lower(wn,eta,fell,...
      p,alphap,gamma_wp,eep,eem,Dz,b,Nx,Nz,N-s);
  for r=0:N-s
    S_W_rs(:,r+1,s+1) = Sn(:,r+1);
  end
  
end

return;

%
% field_tfe_helmholtz_lower
%

function [wn] = field_tfe_helmholtz_lower(Wn,eta,f,p,alphap,gammap,...
    eep,eem,Dz,b,Nx,Nz,N)

wn = zeros(Nx,Nz+1,N+1);
Wnhat = zeros(Nx,N+1);

k2 = alphap(0+1)^2 + gammap(0+1)^2;

ell_top = 0 + 1;
ell_bottom = Nz + 1;
for n=0:N
  Wnhat(:,n+1) = fft(eem.*Wn(:,n+1));
end
f_x = real(ifft( (1i*p).*fft(f) ));

ll = [0:Nz]';
z_min = -b; z_max = 0.0;
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;

f_full = zeros(Nx,Nz+1);
f_x_full = zeros(Nx,Nz+1);
for ell=0:Nz
  f_full(:,ell+1) = f;
  f_x_full(:,ell+1) = f_x;
end
temphat = zeros(Nx,Nz+1);

z_plus_b_full = zeros(Nx,Nz+1);
for j=1:Nx
  z_plus_b_full(j,:) = z.' + b;
end

% Order zero

A0 = -1i*gammap - 1i*eta;
for ell=0:Nz
  wn(:,ell+1,0+1) = eep.*ifft( exp(-1i*gammap*z(ell+1)).*Wnhat(:,0+1)./A0 );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx,Nz+1);
  Hn = Wn(:,n+1);
  Jn = zeros(Nx,1);
  
  A1_xx = (2.0/b)*f_full;
  A1_xz = -(1.0/b)*(z_plus_b_full).*f_x_full;
  A1_zx = A1_xz;
  %A1_zz = 0;
  
  A2_xx = (1.0/b^2)*f_full.^2;
  A2_xz = -(1.0/b^2)*(z_plus_b_full).*(f_full.*f_x_full);
  A2_zx = A2_xz;
  A2_zz = (1.0/b^2)*((z_plus_b_full).^2).*(f_x_full.^2);
  
  B1_x = (1.0/b)*f_x_full;
  %B1_z = 0;
  
  B2_x = (1.0/b^2)*f_full.*f_x_full;
  B2_z = -(1.0/b^2).*(z_plus_b_full).*(f_x_full.^2);
  
  C1 = k2*(2.0/b)*f_full;
  C2 = k2*(1.0/b^2)*f_full.^2;
  
  if(n>=1)
    w_x = dxp(wn(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    temp = A1_xx.*w_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A1_zx.*w_x;
    Fn = Fn - dz(temp,Dz,b,Nx,Nz);
    temp = B1_x.*w_x;
    Fn = Fn + temp;
    
    w_z = dz(wn(:,:,n-1+1),Dz,b,Nx,Nz);
    temp = A1_xz.*w_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    %A1_zz = 0
    %B1_z = 0
    
    temp = C1.*wn(:,:,n-1+1);
    Fn = Fn - temp;
    
    Su = eep.*ifft( (-1i*gammap).*fft(eem.*wn(:,ell_bottom,n-1+1)) );
    Jn = Jn + (1.0/b)*f.*Su;
    Hn = Hn + (1.0/b)*f.*Wn(:,n-1+1)...
        + (1i*eta/b)*f.*wn(:,ell_top,n-1+1)...
        + f_x.*w_x(:,ell_top);
  end
  
  if(n>=2)
    w_x = dxp(wn(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    temp = A2_xx.*w_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zx.*w_x;
    Fn = Fn - dz(temp,Dz,b,Nx,Nz);
    temp = B2_x.*w_x;
    Fn = Fn + temp;
    
    w_z = dz(wn(:,:,n-2+1),Dz,b,Nx,Nz);
    temp = A2_xz.*w_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zz.*w_z;
    Fn = Fn - dz(temp,Dz,b,Nx,Nz);
    temp = B2_z.*w_z;
    Fn = Fn + temp;
    
    temp = C2.*wn(:,:,n-2+1);
    Fn = Fn - temp;
    
    Hn = Hn + (1.0/b)*f.*f_x.*w_x(:,ell_top)...
        - f_x.*f_x.*w_z(:,ell_top);
  end
  
  % Solve elliptic equation
  
  Fnhat = zeros(Nx,Nz+1);
  for ell=0:Nz
    Fnhat(:,ell+1) = fft(eem.*Fn(:,ell+1));
  end
  Hnhat = fft(eem.*Hn);
  Jnhat = fft(eem.*Jn);
  
  for j=1:Nx
    Fnhat_p = Fnhat(j,:).';
    alphaalpha = 1.0;
    betabeta = 0.0;
    gammagamma = k2 - (alphap(j))^2;
    d_a = -(-1i*gammap(j));
    n_a = 1.0;
    r_a = Jnhat(j);
    d_b = -1i*eta;
    n_b = 1.0;
    r_b = Hnhat(j);
    what_p = solvebvp_colloc(Fnhat_p,alphaalpha,betabeta,gammagamma,...
        (2.0/(z_max-z_min))*Dz,d_a,n_a,r_a,d_b,n_b,r_b);
    
    temphat(j,:) = what_p.';
  end
  
  for ell=0:Nz
    wn(:,ell+1,n+1) = eep.*ifft(temphat(:,ell+1));
  end

end

return;

%
% field_tfe_helmholtz_middle
%

function [un] = field_tfe_helmholtz_middle(Vun,Velln,hbar,eta,fu,fell,...
    p,alphap,gammap,eep,eem,Dz,a,Nx,Nz,N)

un = zeros(Nx,Nz+1,N+1);
Vunhat = zeros(Nx,N+1);
Vellnhat = zeros(Nx,N+1);

k2 = alphap(0+1)^2 + gammap(0+1)^2;

ell_top = 0 + 1;
ell_bottom = Nz + 1;
for n=0:N
  Vunhat(:,n+1) = fft(eem.*Vun(:,n+1));
  Vellnhat(:,n+1) = fft(eem.*Velln(:,n+1));
end
fu_x = real(ifft( (1i*p).*fft(fu) ));
fell_x = real(ifft( (1i*p).*fft(fell) ));

ll = [0:Nz]';
z_min = -hbar; z_max = hbar;
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;

fu_full = zeros(Nx,Nz+1);
fell_full = zeros(Nx,Nz+1);
fu_x_full = zeros(Nx,Nz+1);
fell_x_full = zeros(Nx,Nz+1);
for ell=0:Nz
  fu_full(:,ell+1) = fu;
  fell_full(:,ell+1) = fell;
  fu_x_full(:,ell+1) = fu_x;
  fell_x_full(:,ell+1) = fell_x;
end
temphat = zeros(Nx,Nz+1);

Z_L = zeros(Nx,Nz+1);
Z_U = zeros(Nx,Nz+1);
for ell=0:Nz
  for j=1:Nx
    Z_L(j,ell+1) = (z(ell+1) - z_min)/(2.0*hbar);
    Z_U(j,ell+1) = (z_max - z(ell+1))/(2.0*hbar);
  end
end

% Order zero

sh = sinh(1i*gammap*hbar)./(1i*gammap);
ch = cosh(1i*gammap*hbar);
aa = -(gammap.^2).*sh - 1i*eta*ch;
bb = ch - 1i*eta*sh;
Bp = (Vunhat(:,0+1) + Vellnhat(:,0+1))./(2.0*aa);
Cp = (Vunhat(:,0+1) - Vellnhat(:,0+1))./(2.0*bb);
for ell=0:Nz
  un(:,ell+1,0+1) = eep.*ifft( Bp.*cosh(1i*gammap*z(ell+1))...
      + Cp.*sinh(1i*gammap*z(ell+1))./(1i*gammap) );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx,Nz+1);
  Hun = Vun(:,n+1);
  Helln = Velln(:,n+1);
  
  A10_xx = (-2.0/(2.0*hbar))*fell_full;
  A10_xz = -Z_U.*fell_x_full;
  A10_zx = A10_xz;
  %A10_zz = 0;
  
  A01_xx = (2.0/(2.0*hbar))*fu_full;
  A01_xz = -Z_L.*fu_x_full;
  A01_zx = A01_xz;
  %A01_zz = 0;
  
  A1_xx = A10_xx + A01_xx;
  A1_xz = A10_xz + A01_xz;
  A1_zx = A10_zx + A01_zx;
  %A1_zz = 0;
  
  A20_xx = (1.0/(2*hbar)^2)*fell_full.^2;
  A20_xz = (1.0/(2*hbar))*Z_U.*fell_full.*fell_x_full;
  A20_zx = A20_xz;
  A20_zz = (Z_U.^2).*(fell_x_full).^2;
  
  A11_xx = -(2.0/(2*hbar)^2)*fell_full.*fu_full;
  A11_xz = (1.0/(2*hbar))*(Z_L.*fell_full.*fu_x_full...
      -Z_U.*fu_full.*fell_x_full);
  A11_zx = A11_xz;
  A11_zz = 2*Z_U.*Z_L.*fell_x_full.*fu_x_full;
  
  A02_xx = (1.0/(2*hbar)^2)*fu_full.^2;
  A02_xz = -(1.0/(2*hbar))*Z_L.*fu_full.*fu_x_full;
  A02_zx = A02_xz;
  A02_zz = (Z_L.^2).*(fu_x_full).^2;
    
  A2_xx = A20_xx + A11_xx + A02_xx;
  A2_xz = A20_xz + A11_xz + A02_xz;
  A2_zx = A2_xz;
  A2_zz = A20_zz + A11_zz + A02_zz;
  
  B10_x = (1.0/(2.0*hbar))*fell_x_full;
  %B10_z = 0
  B01_x = -(1.0/(2.0*hbar))*fu_x_full;
  %B01_z = 0
  B20_x = -(1.0/(2.0*hbar)^2)*fell_full.*fell_x_full;
  B20_z = -(1.0/(2.0*hbar))*Z_U.*fell_x_full.^2;
  B11_x = (1.0/(2.0*hbar)^2)*(fu_full.*fell_x_full + fell_full.*fu_x_full);
  B11_z = (1.0/(2.0*hbar))*(Z_U-Z_L).*fu_x_full.*fell_x_full;
  B02_x = -(1.0/(2.0*hbar)^2)*fu_full.*fu_x_full;
  B02_z = (1.0/(2.0*hbar))*Z_L.*fu_x_full.^2;
  
  B1_x = B10_x + B01_x;
  %B1_z = 0
  B2_x = B20_x + B11_x + B02_x;
  B2_z = B20_z + B11_z + B02_z;
  
  C10 = -k2*(2.0/(2.0*hbar))*fell_full;
  C01 = k2*(2.0/(2.0*hbar))*fu_full;
  C20 = k2*(1.0/(2.0*hbar)^2)*fell_full.^2;
  C11 = -k2*(2.0/(2.0*hbar)^2)*fu_full.*fell_full;
  C02 = k2*(1.0/(2.0*hbar)^2)*fu_full.^2;
  
  C1 = C10 + C01;
  C2 = C20 + C11 + C02;
  
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    temp = A1_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A1_zx.*u_x;
    Fn = Fn - dz(temp,Dz,2*hbar,Nx,Nz);
    temp = B1_x.*u_x;
    Fn = Fn - temp;
    
    u_z = dz(un(:,:,n-1+1),Dz,2*hbar,Nx,Nz);
    temp = A1_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    %A1_zz = 0
    %B1_z = 0
    
    temp = C1.*un(:,:,n-1+1);
    Fn = Fn - temp;
    
    Hun = Hun ...
        + (1.0/(2*hbar))*fu.*Vun(:,n-1+1)...
        - (1.0/(2*hbar))*fell.*Vun(:,n-1+1)...
        + (1i*eta/(2*hbar))*fu.*un(:,ell_top,n-1+1)...
        - (1i*eta/(2*hbar))*fell.*un(:,ell_top,n-1+1)...
        + fu_x.*u_x(:,ell_top);
    
    Helln = Helln ...
        + (1.0/(2*hbar))*fu.*Velln(:,n-1+1)...
        - (1.0/(2*hbar))*fell.*Velln(:,n-1+1)...
        + (1i*eta/(2*hbar))*fu.*un(:,ell_bottom,n-1+1)...
        - (1i*eta/(2*hbar))*fell.*un(:,ell_bottom,n-1+1)...
        - fell_x.*u_x(:,ell_bottom);
  end
  
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    temp = A2_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zx.*u_x;
    Fn = Fn - dz(temp,Dz,2*hbar,Nx,Nz);
    temp = B2_x.*u_x;
    Fn = Fn - temp;
    
    u_z = dz(un(:,:,n-2+1),Dz,2*hbar,Nx,Nz);
    temp = A2_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zz.*u_z;
    Fn = Fn - dz(temp,Dz,2*hbar,Nx,Nz);
    temp = B2_z.*u_z;
    Fn = Fn - temp;
    
    temp = C2.*un(:,:,n-2+1);
    Fn = Fn - temp;
    
    Hun = Hun ...
        + (1.0/(2*hbar))*fu.*fu_x.*u_x(:,ell_top)...
        - (1.0/(2*hbar))*fell.*fu_x.*u_x(:,ell_top)...
        - fu_x.*fu_x.*u_z(:,ell_top);

    Helln = Helln ...
        - (1.0/(2*hbar))*fu.*fell_x.*u_x(:,ell_bottom)...
        + (1.0/(2*hbar))*fell.*fell_x.*u_x(:,ell_bottom)...
        + fell_x.*fell_x.*u_z(:,ell_bottom);
  end
  
  % Solve elliptic equation
  
  Fnhat = zeros(Nx,Nz+1);
  for ell=0:Nz
    Fnhat(:,ell+1) = fft(eem.*Fn(:,ell+1));
  end
  Hunhat = fft(eem.*Hun);
  Hellnhat = fft(eem.*Helln);
  
  for j=1:Nx
    Fnhat_p = Fnhat(j,:).';
    alphaalpha = 1.0;
    betabeta = 0.0;
    gammagamma = k2 - (alphap(j))^2;
    d_a = -1i*eta;
    n_a = -1.0;
    r_a = Hellnhat(j);
    d_b = -1i*eta;
    n_b = 1.0;
    r_b = Hunhat(j);
    uhat_p = solvebvp_colloc(Fnhat_p,alphaalpha,betabeta,gammagamma,...
        (2.0/(z_max-z_min))*Dz,d_a,n_a,r_a,d_b,n_b,r_b);
    temphat(j,:) = uhat_p.';
  end
  
  for ell=0:Nz
    un(:,ell+1,n+1) = eep.*ifft(temphat(:,ell+1));
  end

end

return;

%
% field_tfe_helmholtz_upper
%

function [un] = field_tfe_helmholtz_upper(Un,eta,f,p,alphap,gammap,...
    eep,eem,Dz,a,Nx,Nz,N)

un = zeros(Nx,Nz+1,N+1);
Unhat = zeros(Nx,N+1);

k2 = alphap(0+1)^2 + gammap(0+1)^2;

ell_top = 0 + 1;
ell_bottom = Nz + 1;
for n=0:N
  Unhat(:,n+1) = fft(eem.*Un(:,n+1));
end
f_x = real(ifft( (1i*p).*fft(f) ));

ll = [0:Nz]';
z_min = 0.0; z_max = a;
z = ((z_max-z_min)/2.0)*(cos(pi*ll/Nz) - 1.0) + z_max;

f_full = zeros(Nx,Nz+1);
f_x_full = zeros(Nx,Nz+1);
for ell=0:Nz
  f_full(:,ell+1) = f;
  f_x_full(:,ell+1) = f_x;
end
temphat = zeros(Nx,Nz+1);

a_minus_z_full = zeros(Nx,Nz+1);
for j=1:Nx
  a_minus_z_full(j,:) = a - z.';
end

% Order zero

A0 = -1i*gammap - 1i*eta;
for ell=0:Nz
  un(:,ell+1,0+1) = eep.*ifft( exp(1i*gammap*z(ell+1)).*Unhat(:,0+1)./A0 );
end

% Order n>0

for n=1:N
    
  % Form Fn, Jn
  
  Fn = zeros(Nx,Nz+1);
  Hn = Un(:,n+1);
  Jn = zeros(Nx,1);
  
  A1_xx = -(2.0/a)*f_full;
  A1_xz = -(1.0/a)*(a_minus_z_full).*f_x_full;
  A1_zx = A1_xz;
  %A1_zz = 0;
  
  A2_xx = (1.0/a^2)*f_full.^2;
  A2_xz = (1.0/a^2)*(a_minus_z_full).*(f_full.*f_x_full);
  A2_zx = A2_xz;
  A2_zz = (1.0/a^2)*((a_minus_z_full).^2).*(f_x_full.^2);
  
  B1_x = -(1.0/a)*f_x_full;
  %B1_z = 0;
  
  B2_x = (1.0/a^2)*f_full.*f_x_full;
  B2_z = (1.0/a^2).*(a_minus_z_full).*(f_x_full.^2);
  
  C1 = -k2*(2.0/a)*f_full;
  C2 = k2*(1.0/a^2)*f_full.^2;
  
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    temp = A1_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A1_zx.*u_x;
    Fn = Fn - dz(temp,Dz,a,Nx,Nz);
    temp = B1_x.*u_x;
    Fn = Fn + temp;
    
    u_z = dz(un(:,:,n-1+1),Dz,a,Nx,Nz);
    temp = A1_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    %A1_zz = 0
    %B1_z = 0
    
    temp = C1.*un(:,:,n-1+1);
    Fn = Fn - temp;
    
    Su = eep.*ifft( (1i*gammap).*fft(eem.*un(:,ell_top,n-1+1)) );
    Jn = Jn - (1.0/a)*f.*Su;
    Hn = Hn - (1.0/a)*f.*Un(:,n-1+1)...
        - (1i*eta/a)*f.*un(:,ell_bottom,n-1+1)...
        - f_x.*u_x(:,ell_bottom);
  end
  
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    temp = A2_xx.*u_x;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zx.*u_x;
    Fn = Fn - dz(temp,Dz,a,Nx,Nz);
    temp = B2_x.*u_x;
    Fn = Fn + temp;
    
    u_z = dz(un(:,:,n-2+1),Dz,a,Nx,Nz);
    temp = A2_xz.*u_z;
    Fn = Fn - dxp(temp,alphap,eep,eem,Nx,Nz);
    temp = A2_zz.*u_z;
    Fn = Fn - dz(temp,Dz,a,Nx,Nz);
    temp = B2_z.*u_z;
    Fn = Fn + temp;
    
    temp = C2.*un(:,:,n-2+1);
    Fn = Fn - temp;
    
    Hn = Hn + (1.0/a)*f.*f_x.*u_x(:,ell_bottom)...
        + f_x.*f_x.*u_z(:,ell_bottom);
  end
  
  % Solve elliptic equation
  
  Fnhat = zeros(Nx,Nz+1);
  for ell=0:Nz
    Fnhat(:,ell+1) = fft(eem.*Fn(:,ell+1));
  end
  Hnhat = fft(eem.*Hn);
  Jnhat = fft(eem.*Jn);
  
  for j=1:Nx
    Fnhat_p = Fnhat(j,:).';
    alphaalpha = 1.0;
    betabeta = 0.0;
    gammagamma = k2 - (alphap(j))^2;
    d_a = -1i*eta;
    n_a = -1.0;
    r_a = Hnhat(j);
    d_b = -1i*gammap(j);
    n_b = 1.0;
    r_b = Jnhat(j);
    uhat_p = solvebvp_colloc(Fnhat_p,alphaalpha,betabeta,gammagamma,...
        (2.0/(z_max-z_min))*Dz,d_a,n_a,r_a,d_b,n_b,r_b);
    temphat(j,:) = uhat_p.';
  end
  
  for ell=0:Nz
    un(:,ell+1,n+1) = eep.*ifft(temphat(:,ell+1));
  end

end

return;

%
% iio_tfe_helmholtz_lower
%

function [Sn] = iio_tfe_helmholtz_lower(wn,eta,f,p,alphap,gammap,...
    eep,eem,Dz,b,Nx,Nz,N)

Sn = zeros(Nx,N+1);

ell_top = 0 + 1;
f_x = ifft( (1i*p).*fft(f) );

for n=0:N
  w_z = dz(wn(:,:,n+1),Dz,b,Nx,Nz);
  Sn(:,n+1) = w_z(:,ell_top) + 1i*eta*wn(:,ell_top,n+1);
  if(n>=1)
    w_x = dxp(wn(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    Sn(:,n+1) = Sn(:,n+1) - f_x.*w_x(:,ell_top);
    
    Sn(:,n+1) = Sn(:,n+1) + (1i*eta/b)*f.*wn(:,ell_top,n-1+1);
    
    Sn(:,n+1) = Sn(:,n+1) - (1.0/b)*(f.*Sn(:,n-1+1));
  end
  if(n>=2)
    w_x = dxp(wn(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    Sn(:,n+1) = Sn(:,n+1) - (1.0/b)*(f.*(f_x.*w_x(:,ell_top)));

    w_z = dz(wn(:,:,n-2+1),Dz,b,Nx,Nz);
    Sn(:,n+1) = Sn(:,n+1) + f_x.*(f_x.*w_z(:,ell_top));
  end
end

return;

%
% iio_tfe_helmholtz_middle
%

function [Run,Relln] = iio_tfe_helmholtz_middle(un,hbar,eta,fu,fell,...
    p,alphap,gammap,eep,eem,Dz,a,Nx,Nz,N)

Run = zeros(Nx,N+1);
Relln = zeros(Nx,N+1);

ell_top = 0 + 1;
ell_bottom = Nz + 1;
fu_x = ifft( (1i*p).*fft(fu) );
fell_x = ifft( (1i*p).*fft(fell) );

for n=0:N
  u_z = dz(un(:,:,n+1),Dz,2*hbar,Nx,Nz);
  Run(:,n+1) = u_z(:,ell_top) + 1i*eta*un(:,ell_top,n+1);
  Relln(:,n+1) = -u_z(:,ell_bottom) + 1i*eta*un(:,ell_bottom,n+1);
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    
    Run(:,n+1) = Run(:,n+1) - fu_x.*u_x(:,ell_top);
    Run(:,n+1) = Run(:,n+1) + (1i*eta/(2*hbar))*fu.*un(:,ell_top,n-1+1);
    Run(:,n+1) = Run(:,n+1) - (1i*eta/(2*hbar))*fell.*un(:,ell_top,n-1+1);
    Run(:,n+1) = Run(:,n+1) - (1.0/(2*hbar))*fu.*Run(:,n-1+1);
    Run(:,n+1) = Run(:,n+1) + (1.0/(2*hbar))*fell.*Run(:,n-1+1);
    
    Relln(:,n+1) = Relln(:,n+1) + fell_x.*u_x(:,ell_bottom);
    Relln(:,n+1) = Relln(:,n+1) + (1i*eta/(2*hbar))*fu.*un(:,ell_bottom,n-1+1);
    Relln(:,n+1) = Relln(:,n+1) - (1i*eta/(2*hbar))*fell.*un(:,ell_bottom,n-1+1);
    Relln(:,n+1) = Relln(:,n+1) - (1.0/(2*hbar))*fu.*Relln(:,n-1+1);
    Relln(:,n+1) = Relln(:,n+1) + (1.0/(2*hbar))*fell.*Relln(:,n-1+1);
  end
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    u_z = dz(un(:,:,n-2+1),Dz,2*hbar,Nx,Nz);
    
    Run(:,n+1) = Run(:,n+1) - (1.0/(2*hbar))*fu.*(fu_x.*u_x(:,ell_top));
    Run(:,n+1) = Run(:,n+1) + (1.0/(2*hbar))*fell.*(fu_x.*u_x(:,ell_top));
    Run(:,n+1) = Run(:,n+1) + fu_x.*(fu_x.*u_z(:,ell_top));
    
    Relln(:,n+1) = Relln(:,n+1) + (1.0/(2*hbar))*fu.*(fell_x.*u_x(:,ell_bottom));
    Relln(:,n+1) = Relln(:,n+1) - (1.0/(2*hbar))*fell.*(fell_x.*u_x(:,ell_bottom));
    Relln(:,n+1) = Relln(:,n+1) - fell_x.*(fell_x.*u_z(:,ell_bottom));
  end
end

return;

%
% iio_tfe_helmholtz_upper
%

function [Qn] = iio_tfe_helmholtz_upper(un,eta,f,p,alphap,gammap,...
    eep,eem,Dz,a,Nx,Nz,N)

Qn = zeros(Nx,N+1);

ell_bottom = Nz + 1;
f_x = ifft( (1i*p).*fft(f) );

for n=0:N
  u_z = dz(un(:,:,n+1),Dz,a,Nx,Nz);
  Qn(:,n+1) = -u_z(:,ell_bottom) + 1i*eta*un(:,ell_bottom,n+1);
  if(n>=1)
    u_x = dxp(un(:,:,n-1+1),alphap,eep,eem,Nx,Nz);
    Qn(:,n+1) = Qn(:,n+1) + f_x.*u_x(:,ell_bottom);
    
    Qn(:,n+1) = Qn(:,n+1) - (1i*eta/a)*f.*un(:,ell_bottom,n-1+1);
    
    Qn(:,n+1) = Qn(:,n+1) + (1.0/a)*(f.*Qn(:,n-1+1));
  end
  if(n>=2)
    u_x = dxp(un(:,:,n-2+1),alphap,eep,eem,Nx,Nz);
    Qn(:,n+1) = Qn(:,n+1) - (1.0/a)*(f.*(f_x.*u_x(:,ell_bottom)));

    u_z = dz(un(:,:,n-2+1),Dz,a,Nx,Nz);
    Qn(:,n+1) = Qn(:,n+1) - f_x.*(f_x.*u_z(:,ell_bottom));
  end
end

return;

%
% dxp
%

function [u_x] = dxp(u,alphap,eep,eem,Nx,Ny)

u_x = zeros(Nx,Ny+1);

for ell=0:Ny
  f = u(:,ell+1);
  u_x(:,ell+1) = eep.*ifft( (1i*alphap).*fft(eem.*f) );
end

return;

%
% dz
%

function [u_z] = dz(u,Dz,b,Nx,Nz)

u_z = zeros(Nx,Nz+1);
g = zeros(Nz+1,1);

for j=1:Nx
  for ell=0:Nz
    g(ell+1) = u(j,ell+1);
  end
  u_z(j,:) = (2.0/b)*Dz*g;
end

return;

%
% solvebvp_colloc
%

function [utilde] = solvebvp_colloc(ftilde,alpha,beta,gamma,D,...
    d_a,n_a,r_a,d_b,n_b,r_b)

Ny = length(ftilde)-1;

D2 = D*D;
A = alpha*D2 + beta*D + gamma*eye(Ny+1);
b = ftilde;

A(Ny+1,:) = n_a*D(Ny+1,:);
A(Ny+1,Ny+1) = A(Ny+1,Ny+1) + d_a;
b(Ny+1) = r_a;

A(1,:) = n_b*D(1,:);
A(1,1) = A(1,1) + d_b;
b(1) = r_b;

utilde = A\b;

return;

%
% CHEB  compute D = differentiation matrix, x = Chebyshev grid
%
% From "Spectral Methods in MATLAB" by Lloyd N. Trefethen
%
% https://doi.org/10.1137/1.9780898719598
%

function [D,x] = cheb(N)
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)'; 
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';                  
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries

return;