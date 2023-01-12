clear all; clc; close all;
global chi0 chi1 gamma
tic;
%==========================================================================
%==========================================================================
                            %Set Parameters
%==========================================================================                         

%Production
A0=1;
alpha=0.33;
epsilon=10;
delta=0.07;

% Asset markets:
kappa=0.07;
rb_pos=0.03;
rb_neg=rb_pos + kappa;

%Preferences
rho = 0.06;
gamma=2;

%Percentage of income deposited in the iliquid asset account
xi=0.33;

%Government:
%Income tax
tau=0.25;
T=0.06;
Bgr=0.289;%1.11812;
Gr = tau*(1-alpha)*(epsilon-1)/epsilon + tau*(1/epsilon) - rb_pos*Bgr - T;

%Deposit Costs:
chi0 =0.03;
chi1=2;

%Productivity (for now, a simple poisson "jump" process)
%w = 4;
Nz = 2;
z      = [0.8,1.3];
la_mat = [-1/3, 1/3; 1/3, -1/3];
meanlabprod = (z(1)*la_mat(1,2) + z(2)*la_mat(2,1))/(la_mat(1,2) + la_mat(2,1));

%==========================================================================
                  %Parameters for the Steady State Convergence
%==========================================================================
Intr=40;
crit_S=10^(-5);

%==========================================================================
                                 %Grids
%==========================================================================
%Liquid Asset
L = 109;
min_b = -4;
%bmin = 0;
max_b =20;
b = linspace(min_b,max_b,L)';
db = (max_b-min_b)/(L-1);

%Iliquid Asset
J= 36;
min_a = 0;
max_a = 50;
a = linspace(min_a,max_a,J);
da = (max_a-min_a)/(J-1);

%==========================================================================
                    %Constructing the States' matrices
%==========================================================================

bb = b*ones(1,J);
aa = ones(L,1)*a;
zz = ones(J,1)*z;

bbb = zeros(L,J,Nz); 
aaa = zeros(L,J,Nz); 
zzz = zeros(L,J,Nz);

for nz=1:Nz
    bbb(:,:,nz) = bb;
    aaa(:,:,nz) = aa;
    zzz(:,:,nz) = z(nz);
end

Bswitch = [
    speye(L*J)*la_mat(1,1), speye(L*J)*la_mat(1,2);
    speye(L*J)*la_mat(2,1), speye(L*J)*la_mat(2,2)];

%==========================================================================
%==========================================================================
   % Step 1: Guess initial values for K and Z and rb
%==========================================================================
rk_max = 0.12;
rk_min = delta+0.02;
rk = 0.09;

Z=meanlabprod;

rb_min = 0.01;
rb_max = 0.99*rho;
rb = rb_pos;

%==========================================================================
%==========================================================================
    %Step 5: Solve the Households Problems using the Achdou et al(2020)
%==========================================================================


% Values necessary to the FD algorithm

maxit= 50;     % Maximum number of iterations allowed
crit = 10^(-6); % Criterium for convergence
Delta = 100;    %Step size


%==========================================================================
%==========================================================================
                              %HJB eq
%==========================================================================
	%Some matrices that will be useful at the end of the code
%==========================================================================
Vaf = zeros(L,J,Nz);
Vab = zeros(L,J,Nz);
Vbf = zeros(L,J,Nz);
Vbb = zeros(L,J,Nz);
c = zeros(L,J,Nz);
BBlowdiag=zeros(L*J,Nz);
BBUpdiag=zeros(L*J,Nz);
BBcentdiag=zeros(L*J,Nz);
AAi = cell(Nz,1);
BBi = cell(Nz,1);

for intr=1:Intr
 
 %ZD_v(intr) = Z;
 %ZD_minv(intr) = Z_min;
 %ZD_maxv(intr) = Z_max;
 
%==========================================================================
%==========================================================================
   % Step 2: Calculate Implied prices
%==========================================================================
%Marginal Cost
m = (epsilon-1)/epsilon;
%Capital interest rate
K = Z*((((alpha*A0*m))/rk)^(1/(1-alpha)));
%Aggregate Output
Y = A0*(K^(alpha))*(Z^(1-alpha));
%Government Supply of Bonds:
BgS(intr) = Bgr*Y;
%Iliquid Assets interest rate
ra=rk-delta;
%Wages
w=((1-alpha)*A0*m*((K/Z)^alpha));
%Aggregate profits
pi=(1-m)*Y;


rb_v(intr)=rb;
%rb_minv(intr)=rb_min;
%rb_maxv(intr)= rb_max;
KD_v(intr) = K;
%KD_minv(intr) = K_min;
%KD_maxv(intr) = K_max;

%==========================================================================
%==========================================================================
   % Step 3: Calculate Workers' bonuses
%==========================================================================
%cor=(1-tau)/(1-xi-tau);
profshare = (zzz/meanlabprod)*(pi);

%==========================================================================
%==========================================================================
   % Step 4: Fiscal Policy
%==========================================================================
Tr = (tau*(1-alpha)*(epsilon-1)/epsilon +tau*(1/epsilon) - rb*Bgr - Gr)*Y;

%==========================================================================

% First guess for V0

if intr>1
    v0 = V_r(:,:,:,intr-1);
else
    v0 = (((1-xi-tau)*(w*zzz + profshare) + (rb+kappa).*bbb + ra.*aaa + Tr).^(1-gamma))/(1-gamma)/rho;
end

v=v0;

%matrix of liquid returns
Rb = rb.*(bbb >0) + (rb+kappa).*(bbb <0);

%matrix of illiquid assets returns
tau2 = 9; raa = ra.*(1 - (1.35.*max_a./a).^(1-tau2));
%matrix of illiquid returns
Ra(:,:,1) = ones(L,1)*raa;
Ra(:,:,2) = ones(L,1)*raa;

for j = 1:maxit

    V = v;

    % Forward difference state b
    Vbf(1:L-1,:,:) = (V(2:L,:,:) - V(1:L-1,:,:))/db;
    Vbf(L,:,:) = ((1-xi-tau)*(w*zzz(L,:,:) + profshare(L,:,:)) + Rb(L,:,:).*max_b + Tr).^(-gamma);
    
	% Backward difference state b
    Vbb(2:L,:,:) = (V(2:L,:,:) - V(1:L-1,:,:))/db;
    Vbb(1,:,:) = ((1-xi-tau)*(w*zzz(1,:,:) + profshare(1,:,:)) + Rb(1,:,:).*min_b + Tr).^(-gamma);
   
     % Foward difference state a
    Vaf(:,1:J-1,:) = (V(:,2:J,:) - V(:,1:J-1,:))/da;
    % Backward difference state a
    Vab(:,2:J,:) = (V(:,2:J,:) - V(:,1:J-1,:))/da;
    
    cb = max(Vbb,10^(-6)).^(-1/gamma); %Consumption
	cf = max(Vbf,10^(-6)).^(-1/gamma); %Consumption
    
    dff = min(Vaf./Vbf -1 + chi0,0) .*(aaa/chi1) + max(Vaf./Vbf -1 - chi0,0).*(aaa/chi1);
	dfb = min(Vab./Vbf -1 + chi0,0) .*(aaa/chi1) + max(Vab./Vbf -1 - chi0,0).*(aaa/chi1);
	dbf = min(Vaf./Vbb -1 + chi0,0) .*(aaa/chi1) + max(Vaf./Vbb -1 - chi0,0).*(aaa/chi1);
	dbb = min(Vab./Vbb -1 + chi0,0) .*(aaa/chi1) + max(Vab./Vbb -1 - chi0,0).*(aaa/chi1);
    
    d_b = (dbf >0).*dbf + (dbb <0).*dbb;
    d_b(:,1,:) = (dbf(:,1,:) >10^(-12)).*dbf(:,1,:);
    d_b(:,J,:) = (dbb(:,J,:)<-10^(-12)).*dbb(:,J,:);
	d_b(1,1,:) = max(d_b(1,1,:),0);
    scb = (1-xi-tau)*(w*zzz + profshare) + Rb.*bbb - cb + Tr;
    % Laws of motion - Backward difference: upwind of split b
    sdb = - d_b - (chi0*abs(d_b) + (chi1.*d_b.^2/2)./max(aaa,10^(-5)));
    
    df= (dff>0).*dff + (dfb<0).*dfb;
    df(:,1,:) = (dff(:,1,:) >10^(-12)).*dff(:,1,:);
    df(:,J,:) = (dfb(:,J,:)<-10^(-12)).*dfb(:,J,:);
    scf = (1-xi-tau)*(w*zzz + profshare) + Rb.*bbb - cf + Tr;
    sdf = - df - (chi0*abs(df) + (chi1.*df.^2/2)./max(aaa,10^(-5)));
	sdf(L,:,:) = min(sdf(L,:,:),0);
    
    Icb = (scb <(-10^(-12)));
    Icf = (scf >(10^(-12))) .*(1 -Icb);
    Ic0 = 1 - Icf -Icb;

	Idf = (sdf >(10^(-12)));
    Idb = (sdb <(-10^(-12))) .*(1 -Idf);
    Idb(1,:,:) =0; Idb(L,:,:) = 1; Idf(L,:,:) = 0;
    Id0 = 1 -Idf -Idb;
    
    c0 = (1-xi-tau)*(w*zzz + profshare) + Rb.*bbb + Tr;


    c = cf.*Icf + cb.*Icb + c0.*Ic0;
    u = c.^(1-gamma)/(1-gamma);

    % Build the HJB transition Matrix regarding state b only

    XB = -Icb.*scb/db - Idb .*sdb/db;
	%XB[1,:,:].=0.0;

    ZB = Icf.*scf/db + Idf.*sdf/db;
	%ZB[L,:,:] .=0.0


	YB = (Icb.*scb - Icf.*scf)/db + (Idb.*sdb - Idf.*sdf)/db;
	%YB = -1*(XB + ZB);

	disp("\n Iteration number: ");
    disp(j);
    
    
    for i = 1:Nz
         BBcentdiag(:,i) = reshape(YB(:,:,i),L*J,1);
    end
    
    BBlowdiag(1:L-1,:) = XB(2:L,1,:);
    BBUpdiag(2:L,:) = ZB(1:L-1,1,:);
    for k = 2:J
        BBlowdiag(1:k*L,:) = [BBlowdiag(1:(k-1)*L,:);squeeze(XB(2:L,k,:));zeros(1,Nz)];
        BBUpdiag(1:k*L,:) = [BBUpdiag(1:(k-1)*L,:);zeros(1,Nz);squeeze(ZB(1:L-1,k,:))];
    end
    
    for nz=1:Nz
        BBi{nz}=spdiags(BBcentdiag(:,nz),0,L*J,L*J)+spdiags([BBUpdiag(:,nz);0],1,L*J,L*J)+spdiags([BBlowdiag(:,nz);0],-1,L*J,L*J);
    end

	BB = [BBi{1}, sparse(L*J,L*J);sparse(L*J,L*J), BBi{2}];

    if max(abs(sum(BB,2)))>10^(-6)
		disp("\n Sum of XB+YB+ZB:")
        disp(max(max(abs((XB+YB+ZB)))));
        disp("Improper matrix of transition");
        disp("Sum BB = ")
        disp(max(abs(sum(BB,2))));
        %break;
    end

    % Build the HJB transition Matrix regarding state a only

% Laws of motion: upwind a

    dB = Idb .*dbb + Idf .*dfb;
    dF = Idb .*dbf + Idf .*dff;
    adb = min(dB,0);  %Forward difference
	%adb(:,1,:) = 0.0;
    adf = max(dF,0) +Ra.*aaa + xi*(w*zzz); %Backward difference
	adb(:,J,:) = dB(:,J,:) + Ra(:,J,:).*max_a + xi*(w*zzz(:,J,:));
	adf(:,J,:) =0.0;
    
	XA = -adb/da;
	ZA = adf/da;
	YA = (adb-adf)/da;


	for nz = 1:Nz

		AAUpdiag = zeros(L,1);
		for k = 1:J
			AAUpdiag = [AAUpdiag;ZA(:,k,nz)];
		end

		AAcentdiag=YA(:,1,nz);
		for k = 2:J-1
			AAcentdiag = [AAcentdiag;YA(:,k,nz)];
		end
		AAcentdiag=[AAcentdiag;YA(:,J,nz)];

		AAlowdiag = XA(:,2,nz);
		for k = 3:J
			AAlowdiag = [AAlowdiag;XA(:,k,nz)];
		end

		%check_point=max(abs(sum(AAcentdiag+[zeros(L,2);AAlowdiag(1:end-L,:)]+AAUpdiag,2)));
		%disp("Checkpoint AA: ");
        %disp(check_point);

		AAi{nz} = spdiags(AAcentdiag,0,L*J,L*J) + spdiags(AAlowdiag,-L,L*J,L*J) + spdiags(AAUpdiag,L,L*J,L*J);
	end

	AA = [AAi{1}, sparse(L*J,L*J);sparse(L*J,L*J), AAi{2}];


    if max(abs(sum(AA,2))) > 10^(-9)
		disp("Sum of XA+YA+ZA:");
        disp(max(max(abs((XA+YA+ZA)))));
        disp("Improper matrix of transition");
        disp("Sum AA=");
        disp(max(abs(sum(AA,2))));
        %break;
    end

    A = AA + BB + Bswitch;

	if max(abs(sum(A,2))) > 10^(-9)
        disp("Improper matrix of transition")
        disp("Sum A=")
        disp(max(abs(sum(A,2))))
        %break;
    end

    B = (1/Delta + rho)*speye(L*J*Nz) -A;

    u_stacked = (reshape(u, L*J*Nz, 1));
    V_stacked = (reshape(V, L*J*Nz, 1));

    vec1 = u_stacked + V_stacked/Delta;

    V_stacked = B\vec1;

    V = reshape(V_stacked, L,J,Nz);

    Vchange = V-v;
    v = V;

	disp("Vchange:");
    disp(max(max(max(abs(Vchange)))));

	if max(max(max(abs(Vchange)))) < crit
        disp("\n Value function converged at: ")
        disp(j)
        break;
    end

end

if j==maxit
    disp('HJB did not converge')
   break;
end

d = Idb.*d_b + Idf.*df;
ad = d + xi*(w*zzz) + Ra.*aaa;
s = (1-xi-tau)*(w*zzz+profshare) + Rb.*bbb - d - (chi0*abs(d) + (chi1.*d.^2/2)./max(aaa,10^(-5))) - c + Tr;

sc = (1-xi-tau)*(w*zzz+profshare) + Rb.*bbb - c + Tr;
sd = - d - (chi0*abs(d) + (chi1.*d.^2/2)./max(aaa,10^(-5)));

%==========================================================================
%==========================================================================
                     %KF eq: (Possible) Stationary Distribution
%==========================================================================                        
AT = A';
% Fix one value so matrix isn't singular:
vec1=zeros(L*J*Nz,1);
iFix=1;
vec1(iFix)= 0.1;
AT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,L*J*Nz-iFix)];

% Solve system:
g_stacked = AT\vec1;
g_sum = g_stacked'*ones(L*J*Nz,1)*da*db;
g_stacked = g_stacked./g_sum;

g(:,:,1) = reshape(g_stacked(1:L*J),L,J);
g(:,:,2) = reshape(g_stacked(L*J+1:L*J*2),L,J);

g_v(:,:,:,intr)=g;
V_r(:,:,:,intr) = V;

%==========================================================================
 %Step 6: Aggegate individual policy functions and effective labour supply
%==========================================================================
AS_v(intr) = sum(sum(g(:,:,1)*a')*db)*da + sum(sum(g(:,:,2)*a'))*da*db;
%ZS_v(intr) 
meanlabprod = sum(sum(g(:,:,1)*z(1)))*da*db + sum(sum(g(:,:,2)*z(2)))*da*db;
BhD_v(intr) = sum(sum(g(:,:,1)'*b))*da*db + sum(sum(g(:,:,2)'*b))*da*db;

%==========================================================================
 %Step 7: Update Capital and Equity
%==========================================================================
q = xi*(1-m)*Y/ra;
if q<0
    disp('Negative value for q!')
    break;
end
if q > AS_v(intr)
    disp('q greater than A!')
    break;
end


KS_v(intr) = AS_v(intr) - q;

%==========================================================================
    %Step 8: Calculate excess demand
%==========================================================================

SK(intr) = (KS_v(intr) - KD_v(intr))/K;
%SZ(intr) = (ZS_v(intr) - ZD_v(intr))/Z;
SBhD(intr) = (BhD_v(intr) - BgS(intr))/BgS(intr);

S(intr) = abs(SK(intr)) + abs(SBhD(intr));% + abs(SZ(intr))

%==========================================================================

%UPDATE GUESSES
if SBhD(intr)<-crit_S*(10^-1)
    disp('Liquid Asset Excess Supply')
    rb_min = rb;
    rb = 0.5*(rb+rb_max);
elseif SBhD(intr)>crit_S*(10^-1)
    disp('Liquid Asset Excess Demand')
    rb_max = rb;
    rb = 0.5*(rb +rb_min);
end
    
if SK(intr)>crit_S*(10^-1)
    disp('Capital Excess Supply')
    rk_max = rk;
    rk = 0.5*rk+ 0.5*rk_min;
elseif SK(intr)<-crit_S*(10^-1)
    disp('Capital Excess Demand')
    rk_min = rk;
    rk = 0.5*rk+ 0.5*rk_max;
end
    
if abs(SK(intr))<crit_S*(10^-1) && abs(SBhD(intr))<crit_S*(10^-1)
    if abs(S(intr))<crit_S
        disp('Equilibrium Found, Interest rate =')
        disp(rb)
        disp('Equilibrium Found, Capital =')
        disp(K)
        disp('Accuracy Capital Stock =')
        disp(SK(intr))
        disp('Accuracy Liquid Asset =')
        disp(SBhD(intr))
        disp('Sum Accuracy =')
        disp(S(intr))
        break;
    end
end

end
toc;

%==========================================================================
                                    %Figures
%==========================================================================

figure(1)
%set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,c(:,:,1)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Consumption, Low Type')

subplot(1,2,2)
surf(b,a,c(:,:,2)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Consumption, High Type')
print -depsc consumption.eps

figure(2)
set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,d(:,:,1)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Deposits, Low Type')

subplot(1,2,2)
surf(b,a,d(:,:,2)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Deposits, High Type')
print -depsc deposits.eps


figure(3)
set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,s(:,:,1)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Liquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,s(:,:,2)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Liquid Savings, High Type')
print -depsc sav.eps

figure(4)
set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,ad(:,:,1)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Illiquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,ad(:,:,2)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b max_b])
ylim([min_a max_a])
title('Illiquid Savings, High Type')
print -depsc ill_sav.eps

figure(5)
icut = 25; jcut=36;
bcut=b(1:icut); acut=a(1:jcut);
gcut = g(1:icut,1:jcut,:);

set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(bcut,acut,gcut(:,:,1)')
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b b(icut)])
ylim([min_a a(jcut)])
title('Stationary, Low Type')

subplot(1,2,2)
surf(bcut,acut,gcut(:,:,2)')
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([min_b b(icut)])
ylim([min_a a(jcut)])
title('Stationary Distribution, High Type')
print -depsc stat_dist.eps