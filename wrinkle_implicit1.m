%Wrinkling of a thin film on viscous flow

% clear all

%load data from previous run

% Material parameters
Ey = 5e9;    %unit Pa
eta = 1170;  %unit Pa*s
v = 0.3;

% Dimensions
L = 20e-3;    %unit m
h = 10e-6;     
H0 = 400*h;

% Load
ed = -0.2;    %unit 1/s

% Normalization 
H0 = H0/h;
L = L/h;
beta = (1-v^2)*eta/Ey*ed;

% Numerical parameters
n = 1000;
dx = 2*L/n;

dt = dx*1000;

toleranceb = 5e-8;
toleranceR = 1e3;
ntotal = 50;           % total number of steps

redu = 1;



% % Initial condition
A0 = 0.0001;
xb = (-L+dx:dx:L-dx)';
w = A0*cos(pi/100*xb);
u = zeros(n-1,1);

w0 = A0*cos(pi/100*(-L));
wn = A0*cos(pi/100*(L));
u0 = 0;
un = 0;

N = zeros(n-1,1);
T0 = 0;
T = zeros(n-1,1);
Tn = 0;
p0 = 0;
p = zeros(n-1,1);
pn = 0;

ttotal = 0;


f = zeros(n-1,1);
g = zeros(n-1,1);
Ne = zeros(n-1,1);
Te = zeros(n-1,1);
pe = zeros(n-1,1);

e = zeros(5*n+3,5*n+3);

m = zeros(5*n+3,1);




for tt=1:ntotal
    
    tt
    
    ttotal = ttotal + dt;
    
wp = w;
wp0 = w0;
wpn = wn;

up = u;
up0 = u0;
upn = un;

bb = ones(5*n+3,1);
R = ones(5*n+3,1);
while(norm(bb)>toleranceb||norm(R)>toleranceR)
    

%Boundary condition    
    
    w_1 = -w(1) + 2*w0 ;
    w_2 = w(2) - 4*w(1) + 4*w0;
    u_1 = u(1) + (w(1) - w_1)^2/4/dx;
    
    wn1 = -w(n-1) + 2*wn;
    wn2 = w(n-2) - 4*w(n-1) + 4*wn;
    un1 = u(n-1) + (wn1 - w(n-1))^2/4/dx;
    
    
    N_1 = 0;
    N0 = 0;
    Nn = 0;
    Nn1 = 0;
    
    p_1 = 0;
    pn1 = 0;
    
% %Calculate forces    
%     
%     N_1 = 0;
%     N0 = 0;
%     N(1) = (u(2)-u0)/2/dx + 1/2*((w(2)-w0)/2/dx)^2;
%     for i = 2:n-2
%         N(i) = (u(i+1)-u(i-1))/2/dx + 1/2*((w(i+1)-w(i-1))/2/dx)^2;
%     end
%     N(n-1) = (un-u(n-2))/2/dx + 1/2*((wn-w(n-2))/2/dx)^2;
%     Nn = 0;
%     Nn1 = 0;
%     
%     
%     T0 = (N(1) - N_1)/2/dx;
%     T(1) = (N(2) - N0)/2/dx;
%     for i = 2:n-2
%         T(i) = (N(i+1) - N(i-1))/2/dx;
%     end
%     T(n-1) = (Nn - N(n-2))/2/dx;
%     Tn = (Nn1 - N(n-1))/2/dx;
% 
%     
%     p_1 = 0;
%     p0 = 1/12*(w(2)-4*w(1)+6*w0-4*w_1+w_2)/dx^4 - N0*(w(1)-2*w0+w_1)/dx^2 - T0*(w(1)-w_1)/2/dx;
%     p(1) =  1/12*(w(3)-4*w(2)+6*w(1)-4*w0+w_1)/dx^4 - N(1)*(w(2)-2*w(1)+w0)/dx^2 - T(1)*(w(2)-w0)/2/dx;
%     p(2) =  1/12*(w(4)-4*w(3)+6*w(2)-4*w(1)+w0)/dx^4 - N(2)*(w(3)-2*w(2)+w(1))/dx^2 - T(2)*(w(3)-w(1))/2/dx;
%     for i = 3:n-3 
%        p(i) =  1/12*(w(i+2)-4*w(i+1)+6*w(i)-4*w(i-1)+w(i-2))/dx^4 - N(i)*(w(i+1)-2*w(i)+w(i-1))/dx^2 - T(i)*(w(i+1)-w(i-1))/2/dx;
%     end
%     p(n-2) =  1/12*(wn-4*w(n-1)+6*w(n-2)-4*w(n-3)+w(n-4))/dx^4 - N(n-2)*(w(n-1)-2*w(n-2)+w(n-3))/dx^2 - T(n-2)*(w(n-1)-w(n-3))/2/dx;
%     p(n-1) =  1/12*(wn1-4*wn+6*w(n-1)-4*w(n-2)+w(n-3))/dx^4 - N(n-1)*(wn-2*w(n-1)+w(n-2))/dx^2 - T(n-1)*(wn-w(n-2))/2/dx;
%     pn =  1/12*(wn2-4*wn1+6*wn-4*w(n-1)+w(n-2))/dx^4 - Nn*(wn1-2*wn+w(n-1))/dx^2 - Tn*(wn1-w(n-1))/2/dx;
%     pn1 = 0;
    
    
    
%Eovlution equation

    f0 = w0 - wp0 - dt/dx^2*((H0+w0)^2/4*(w(1)-w_1)*(p(1)-p_1) + (H0+w0)^3/3*(p(1)-2*p0+p_1) - (H0+w0)^2/2*(N(1)-2*N0+N_1) - (H0+w0)/4*(w(1)-w_1)*(N(1)-N_1)) + ...
    dt*(beta*(H0+w0) + beta*(-L)*(w(1)-w_1)/2/dx);
    f(1) = w(1) - wp(1) - dt/dx^2*((H0+w(1))^2/4*(w(2)-w0)*(p(2)-p0) + (H0+w(1))^3/3*(p(2)-2*p(1)+p0) - (H0+w(1))^2/2*(N(2)-2*N(1)+N0) - (H0+w(1))/4*(w(2)-w0)*(N(2)-N0)) + ...
    dt*(beta*(H0+w(1)) + beta*xb(1)*(w(2)-w0)/2/dx);    
    
    for i = 2:n-2
        f(i) = w(i) - wp(i) - dt/dx^2*((H0+w(i))^2/4*(w(i+1)-w(i-1))*(p(i+1)-p(i-1)) + (H0+w(i))^3/3*(p(i+1)-2*p(i)+p(i-1)) - (H0+w(i))^2/2*(N(i+1)-2*N(i)+N(i-1)) - ...
        (H0+w(i))/4*(w(i+1)-w(i-1))*(N(i+1)-N(i-1))) + dt*(beta*(H0+w(i)) + beta*xb(i)*(w(i+1)-w(i-1))/2/dx);
    end
    
    f(n-1) = w(n-1) - wp(n-1) - dt/dx^2*((H0+w(n-1))^2/4*(wn-w(n-2))*(pn-p(n-2)) + (H0+w(n-1))^3/3*(pn-2*p(n-1)+p(n-2)) - (H0+w(n-1))^2/2*(Nn-2*N(n-1)+N(n-2)) - ...
    (H0+w(n-1))/4*(wn-w(n-2))*(Nn-N(n-2))) + dt*(beta*(H0+w(n-1)) + beta*xb(n-1)*(wn-w(n-2))/2/dx);
    fn = wn - wpn - dt/dx^2*((H0+wn)^2/4*(wn1-w(n-1))*(pn1-p(n-1)) + (H0+wn)^3/3*(pn1-2*pn+p(n-1)) - (H0+wn)^2/2*(Nn1-2*Nn+N(n-1)) - ...
    (H0+wn)/4*(wn1-w(n-1))*(Nn1-N(n-1))) + dt*(beta*(H0+wn) + beta*L*(wn1-w(n-1))/2/dx);
    
    
    g0 = u0 - up0 - dt/dx*(-(H0+w0)^2/4*(p(1)-p_1) + (H0+w0)/2*(N(1)-N_1)) - dt*beta*(-L);
    g(1) = u(1) - up(1) - dt/dx*(-(H0+w(1))^2/4*(p(2)-p0) + (H0+w(1))/2*(N(2)-N0)) - dt*beta*xb(1);
    for i = 2:n-2
       g(i) = u(i) - up(i) - dt/dx*(-(H0+w(i))^2/4*(p(i+1)-p(i-1)) + (H0+w(i))/2*(N(i+1)-N(i-1))) - dt*beta*xb(i); 
    end
    g(n-1) = u(n-1) - up(n-1) - dt/dx*(-(H0+w(n-1))^2/4*(pn-p(n-2)) + (H0+w(n-1))/2*(Nn-N(n-2))) - dt*beta*xb(n-1);
    gn = un - upn - dt/dx*(-(H0+wn)^2/4*(pn1-p(n-1)) + (H0+wn)/2*(Nn1-N(n-1))) - dt*beta*L;

    Ne(1) = N(1) - (u(2)-u0)/2/dx - 1/2*((w(2)-w0)/2/dx)^2;
    for i = 2:n-2
        Ne(i) = N(i) - (u(i+1)-u(i-1))/2/dx - 1/2*((w(i+1)-w(i-1))/2/dx)^2;
    end
    Ne(n-1) = N(n-1) - (un-u(n-2))/2/dx - 1/2*((wn-w(n-2))/2/dx)^2;
    
    Te0 = T0 - (N(1) - N_1)/2/dx;
    Te(1) = T(1) - (N(2) - N0)/2/dx;
    for i = 2:n-2
        Te(i) = T(i) - (N(i+1) - N(i-1))/2/dx;
    end
    Te(n-1) = T(n-1) - (Nn - N(n-2))/2/dx;
    Ten = Tn - (Nn1 - N(n-1))/2/dx;
    
    pe0 = p0 - 1/12*(w(2)-4*w(1)+6*w0-4*w_1+w_2)/dx^4 + N0*(w(1)-2*w0+w_1)/dx^2 + T0*(w(1)-w_1)/2/dx;
    pe(1) = p(1) - 1/12*(w(3)-4*w(2)+6*w(1)-4*w0+w_1)/dx^4 + N(1)*(w(2)-2*w(1)+w0)/dx^2 + T(1)*(w(2)-w0)/2/dx;
    pe(2) = p(2) - 1/12*(w(4)-4*w(3)+6*w(2)-4*w(1)+w0)/dx^4 + N(2)*(w(3)-2*w(2)+w(1))/dx^2 + T(2)*(w(3)-w(1))/2/dx;
    for i = 3:n-3 
       pe(i) = p(i) - 1/12*(w(i+2)-4*w(i+1)+6*w(i)-4*w(i-1)+w(i-2))/dx^4 + N(i)*(w(i+1)-2*w(i)+w(i-1))/dx^2 + T(i)*(w(i+1)-w(i-1))/2/dx;
    end
    pe(n-2) = p(n-2) - 1/12*(wn-4*w(n-1)+6*w(n-2)-4*w(n-3)+w(n-4))/dx^4 + N(n-2)*(w(n-1)-2*w(n-2)+w(n-3))/dx^2 + T(n-2)*(w(n-1)-w(n-3))/2/dx;
    pe(n-1) = p(n-1) - 1/12*(wn1-4*wn+6*w(n-1)-4*w(n-2)+w(n-3))/dx^4 + N(n-1)*(wn-2*w(n-1)+w(n-2))/dx^2 + T(n-1)*(wn-w(n-2))/2/dx;
    pen = pn - 1/12*(wn2-4*wn1+6*wn-4*w(n-1)+w(n-2))/dx^4 + Nn*(wn1-2*wn+w(n-1))/dx^2 + Tn*(wn1-w(n-1))/2/dx;
    
     
%Compose Jaccobi matrix

   e(1,1) = 1 - dt/dx^2*((H0+w0)/2*(w(1)-w_1)*(p(1)-p_1) + (H0+w0)^2/4*(-2)*(p(1)-p_1) + (H0+w0)^2*(p(1)-2*p0+p_1) - (H0+w0)*(N(1)-2*N0+N_1) - 1/4*(w(1)-w_1)*(N(1)-N_1)) - (H0+w0)/4*(-2)*(N(1)-N_1) + ...  
   dt*(beta*1 + beta*(-L)*(w(1)-2)/2/dx);
   e(1,2) = - dt/dx^2*((H0+w0)^2/4*(2)*(p(1)-p_1) - (H0+w0)/4*(2)*(N(1)-N_1)) + dt*(beta*(-L)*(2)/2/dx);
   e(1,2*n+3) = - dt/dx^2*(-(H0+w0)^2/2*(1) - (H0+w0)/4*(w(1)-w_1)*(1));
   e(1,4*n+3) = - dt/dx^2*(H0+w0)^3/3*(-2);
   e(1,4*n+4) = - dt/dx^2*((H0+w0)^2/4*(w(1)-w_1)*(1) + (H0+w0)^3/3*(1));
   
   e(2,1) = - dt/dx^2*((H0+w(1))^2/4*(-1)*(p(2)-p0) - (H0+w(1))/4*(-1)*(N(2)-N0)) + dt*(beta*xb(1)*(-1)/2/dx);
   e(2,2) = 1 - dt/dx^2*((H0+w(1))/2*(w(2)-w0)*(p(2)-p0) + (H0+w(1))^2*(p(2)-2*p(1)+p0) - (H0+w(1))*(N(2)-2*N(1)+N0) - 1/4*(w(2)-w0)*(N(2)-N0)) + dt*(beta); 
   e(2,3) = - dt/dx^2*((H0+w(1))^2/4*(1)*(p(2)-p0) - (H0+w(1))/4*(1)*(N(2)-N0)) + dt*(beta*xb(1)*(1)/2/dx);
   e(2,2*n+3) = - dt/dx^2*(-(H0+w(1))^2/2*(-2));
   e(2,2*n+4) = - dt/dx^2*(- (H0+w(1))^2/2*(1) - (H0+w(1))/4*(w(2)-w0)*(1));
   e(2,4*n+3) = - dt/dx^2*((H0+w(1))^2/4*(w(2)-w0)*(-1) + (H0+w(1))^3/3*(1));
   e(2,4*n+4) = - dt/dx^2*((H0+w(1))^3/3*(-2)); 
   e(2,4*n+5) = - dt/dx^2*((H0+w(1))^2/4*(w(2)-w0)*(1) + (H0+w(1))^3/3*(1)); 
   
    for i = 2:n-2
        e(i+1,i) = - dt/dx^2*((H0+w(i))^2/4*(-1)*(p(i+1)-p(i-1)) - (H0+w(i))/4*(-1)*(N(i+1)-N(i-1))) + dt*(beta*xb(i)*(-1)/2/dx);
        e(i+1,i+1) = 1 - dt/dx^2*((H0+w(i))/2*(w(i+1)-w(i-1))*(p(i+1)-p(i-1)) + (H0+w(i))^2*(p(i+1)-2*p(i)+p(i-1)) - (H0+w(i))*(N(i+1)-2*N(i)+N(i-1)) - 1/4*(w(i+1)-w(i-1))*(N(i+1)-N(i-1))) + dt*(beta);
        e(i+1,i+2) = - dt/dx^2*((H0+w(i))^2/4*(1)*(p(i+1)-p(i-1)) - (H0+w(i))/4*(1)*(N(i+1)-N(i-1))) + dt*(beta*xb(i)*(1)/2/dx);
        
        e(i+1,i+2*n+1) = - dt/dx^2*(- (H0+w(i))^2/2*(1) - (H0+w(i))/4*(w(i+1)-w(i-1))*(-1));
        e(i+1,i+2*n+2) = - dt/dx^2*(-(H0+w(i))^2/2*(-2));
        e(i+1,i+2*n+3) = - dt/dx^2*(- (H0+w(i))^2/2*(1) - (H0+w(i))/4*(w(i+1)-w(i-1))*(1));
        
        e(i+1,i+4*n+2) = - dt/dx^2*((H0+w(i))^2/4*(w(i+1)-w(i-1))*(-1) + (H0+w(i))^3/3*(1));
        e(i+1,i+4*n+3) = - dt/dx^2*((H0+w(i))^3/3*(-2));
        e(i+1,i+4*n+4) = - dt/dx^2*((H0+w(i))^2/4*(w(i+1)-w(i-1))*(1) + (H0+w(i))^3/3*(1));
    end
    
    e(n,n-1) = - dt/dx^2*((H0+w(n-1))^2/4*(-1)*(pn-p(n-2)) - (H0+w(n-1))/4*(-1)*(Nn-N(n-2))) + dt*(beta*xb(n-1)*(-1)/2/dx);
    e(n,n) = 1 - dt/dx^2*((H0+w(n-1))/2*(wn-w(n-2))*(pn-p(n-2)) + (H0+w(n-1))^2*(pn-2*p(n-1)+p(n-2)) - (H0+w(n-1))*(Nn-2*N(n-1)+N(n-2)) - 1/4*(wn-w(n-2))*(Nn-N(n-2))) + dt*(beta);
    e(n,n+1) = - dt/dx^2*((H0+w(n-1))^2/4*(1)*(pn-p(n-2)) - (H0+w(n-1))/4*(1)*(Nn-N(n-2))) + dt*(beta*xb(n-1)*(1)/2/dx);
    e(n,3*n) = - dt/dx^2*(- (H0+w(n-1))^2/2*(1) - (H0+w(n-1))/4*(wn-w(n-2))*(-1));
    e(n,3*n+1) = - dt/dx^2*(- (H0+w(n-1))^2/2*(-2));
    e(n,5*n+1) = - dt/dx^2*((H0+w(n-1))^2/4*(wn-w(n-2))*(-1) + (H0+w(n-1))^3/3*(1));
    e(n,5*n+2) = - dt/dx^2*((H0+w(n-1))^3/3*(-2));
    e(n,5*n+3) = - dt/dx^2*((H0+w(n-1))^2/4*(wn-w(n-2))*(1) + (H0+w(n-1))^3/3*(1));
    
    e(n+1,n) = - dt/dx^2*((H0+wn)^2/4*(-2)*(pn1-p(n-1))  - (H0+wn)/4*(-2)*(Nn1-N(n-1))) + dt*(beta*L*(-2)/2/dx);
    e(n+1,n+1) = 1 - dt/dx^2*((H0+wn)/2*(wn1-w(n-1))*(pn1-p(n-1)) + (H0+wn)^2/4*(2)*(pn1-p(n-1)) + (H0+wn)^2*(pn1-2*pn+p(n-1)) - (H0+wn)*(Nn1-2*Nn+N(n-1)) - ...
    1/4*(wn1-w(n-1))*(Nn1-N(n-1)) + (H0+wn)/4*(2)*(Nn1-N(n-1))) + dt*(beta + beta*L*(2)/2/dx);
    e(n+1,3*n+1) =  - dt/dx^2*(- (H0+wn)^2/2*(1) - (H0+wn)/4*(wn1-w(n-1))*(-1));
    e(n+1,5*n+2) = - dt/dx^2*((H0+wn)^2/4*(wn1-w(n-1))*(-1) + (H0+wn)^3/3*(1));
    e(n+1,5*n+3) = - dt/dx^2*((H0+wn)^3/3*(-2));
    
    
    e(n+2,1) = - dt/dx*(-(H0+w0)/2*(p(1)-p_1) + 1/2*(N(1)-N_1));
    e(n+2,n+2) = 1;
    e(n+2,2*n+3) = - dt/dx*((H0+w0)/2*(1));
    e(n+2,4*n+4) = - dt/dx*(-(H0+w0)^2/4*(1));
    
    e(n+3,2) = - dt/dx*(-(H0+w(1))/2*(p(2)-p0) + 1/2*(N(2)-N0));
    e(n+3,n+3) = 1;
    e(n+3,2*n+4) = - dt/dx*((H0+w(1))/2*(1));
    e(n+3,4*n+3) =  - dt/dx*(-(H0+w(1))^2/4*(-1));
    e(n+3,4*n+5) =  - dt/dx*(-(H0+w(1))^2/4*(1));

    for i = 2:n-2
       e(i+n+2,i+1) = - dt/dx*(-(H0+w(i))/2*(p(i+1)-p(i-1)) + 1/2*(N(i+1)-N(i-1)));
       e(i+n+2,i+n+2) = 1;
       e(i+n+2,i+2*n+1) = - dt/dx*((H0+w(i))/2*(-1));
       e(i+n+2,i+2*n+3) = - dt/dx*((H0+w(i))/2*(1));
       e(i+n+2,i+4*n+2) = - dt/dx*(-(H0+w(i))^2/4*(-1)); 
       e(i+n+2,i+4*n+4) = - dt/dx*(-(H0+w(i))^2/4*(1));
    end
    
    e(2*n+1,n) = dt/dx*(-(H0+w(n-1))/2*(pn-p(n-2)) + 1/2*(Nn-N(n-2)));
    e(2*n+1,2*n+1) = 1;
    e(2*n+1,3*n) = - dt/dx*((H0+w(i))/2*(-1));
    e(2*n+1,5*n+1) = - dt/dx*(-(H0+w(i))^2/4*(-1));
    e(2*n+1,5*n+3) = - dt/dx*(-(H0+w(i))^2/4*(1));
    
    e(2*n+2,n+1) = - dt/dx*(-(H0+wn)/2*(pn1-p(n-1)) + 1/2*(Nn1-N(n-1)));
    e(2*n+2,2*n+2) = 1;
    e(2*n+2,3*n+1) = - dt/dx*((H0+wn)/2*(-1));
    e(2*n+2,5*n+2) = - dt/dx*(-(H0+wn)^2/4*(-1));
    
    
    e(2*n+3,1) = - ((w(2)-w0)/2/dx)*(-1)/2/dx;
    e(2*n+3,1) = - ((w(2)-w0)/2/dx)*(1)/2/dx;
    e(2*n+3,n+2) = - (-1)/2/dx;
    e(2*n+3,n+4) = - (1)/2/dx;
    e(2*n+3,2*n+3) = 1;

    for i = 2:n-2
        e(i+2*n+2,i) = - ((w(i+1)-w(i-1))/2/dx)*(-1)/2/dx;
        e(i+2*n+2,i+2) = - ((w(i+1)-w(i-1))/2/dx)*(1)/2/dx;
        e(i+2*n+2,i+n+1) = - (-1)/2/dx;
        e(i+2*n+2,i+n+3) = - (1)/2/dx;
        e(i+2*n+2,i+2*n+2) = 1;
    end
    
    e(3*n+1,n-1) =  -((wn-w(n-2))/2/dx)*(-1)/2/dx;
    e(3*n+1,n+1) =  -((wn-w(n-2))/2/dx)*(1)/2/dx;
    e(3*n+1,2*n) = - (-1)/2/dx;
    e(3*n+1,2*n+2) = - (1)/2/dx;
    e(3*n+1,3*n+1) = 1;
    

    e(3*n+2,2*n+3) = - (1)/2/dx;
    e(3*n+2,3*n+2) = 1;

    e(3*n+3,2*n+4) = - (1)/2/dx;
    e(3*n+3,3*n+3) = 1;
    for i = 2:n-2
        e(i+3*n+2,i+2*n+1) = - (-1)/2/dx;
        e(i+3*n+2,i+2*n+3) = - (1)/2/dx;
        e(i+3*n+2,i+3*n+2) = 1;
    end
    e(4*n+1,3*n) = - (-1)/2/dx;
    e(4*n+1,4*n+1) = 1;
    
    e(4*n+2,3*n+1) = - (-1)/2/dx;
    e(4*n+2,4*n+2) = 1;

    
    e(4*n+3,1) = - 1/12*(6-4*2+4)/dx^4 + N0*(-2+2)/dx^2 + T0*(-2)/2/dx;
    e(4*n+3,2) = - 1/12*(-4-4*(-1)+(-4))/dx^4 + N0*(1+(-1))/dx^2 + T0*(1-(-1))/2/dx;
    e(4*n+3,3) = - 1/12*(1+1)/dx^4;
    e(4*n+3,3*n+2) = (w(1)-w_1)/2/dx;
    e(4*n+3,4*n+3) = 1;
    
    e(4*n+4,1) = - 1/12*(-4+2)/dx^4 + N(1)*(1)/dx^2 + T(1)*(-1)/2/dx;
    e(4*n+4,2) = - 1/12*(6+(-1))/dx^4 + N(1)*(-2)/dx^2;
    e(4*n+4,3) =  - 1/12*(-4)/dx^4 + N(1)*(1)/dx^2 + T(1)*(1)/2/dx;
    e(4*n+4,4) = - 1/12*(1)/dx^4;
    e(4*n+4,2*n+3) = (w(2)-2*w(1)+w0)/dx^2;
    e(4*n+4,3*n+3) = (w(2)-w0)/2/dx;
    e(4*n+4,4*n+4) = 1;
    
    e(4*n+5,1) = - 1/12*(1)/dx^4;
    e(4*n+5,2) = - 1/12*(-4)/dx^4 + N(2)*(1)/dx^2 + T(2)*(-1)/2/dx;
    e(4*n+5,3) = - 1/12*(6)/dx^4 + N(2)*(-2)/dx^2;
    e(4*n+5,4) = - 1/12*(-4)/dx^4 + N(2)*(1)/dx^2 + T(2)*(1)/2/dx;
    e(4*n+5,5) = - 1/12*(1)/dx^4;
    e(4*n+5,2*n+4) = (w(3)-2*w(2)+w(1))/dx^2;
    e(4*n+5,3*n+4) = (w(3)-w(1))/2/dx;
    e(4*n+5,4*n+5) = 1;
    
    
    
    for i = 3:n-3 
        e(i+4*n+3,i-1) = - 1/12*(1)/dx^4;
        e(i+4*n+3,i) = - 1/12*(-4)/dx^4 + N(i)*(1)/dx^2 + T(i)*(-1)/2/dx;
        e(i+4*n+3,i+1) = - 1/12*(6)/dx^4 + N(i)*(-2)/dx^2;
        e(i+4*n+3,i+2) = - 1/12*(-4)/dx^4 + N(i)*(1)/dx^2 + T(i)*(1)/2/dx;
        e(i+4*n+3,i+3) = - 1/12*(1)/dx^4;
        e(i+4*n+3,i+2*n+2) = (w(i+1)-2*w(i)+w(i-1))/dx^2;
        e(i+4*n+3,i+3*n+2) = (w(i+1)-w(i-1))/2/dx;
        e(i+4*n+3,i+4*n+3) = 1;
    end
    
    e(5*n+1,n-3) = - 1/12*(1)/dx^4;
    e(5*n+1,n-2) = - 1/12*(-4)/dx^4 + N(n-2)*(1)/dx^2 + T(n-2)*(-1)/2/dx;
    e(5*n+1,n-1) = - 1/12*(6)/dx^4 + N(n-2)*(-2)/dx^2;
    e(5*n+1,n) = - 1/12*(-4)/dx^4 + N(n-2)*(1)/dx^2 + T(n-2)*(1)/2/dx;
    e(5*n+1,n+1) = - 1/12*(1)/dx^4;
    e(5*n+1,3*n) = (w(n-1)-2*w(n-2)+w(n-3))/dx^2;
    e(5*n+1,4*n) = (w(n-1)-w(n-3))/2/dx;
    e(5*n+1,5*n+1) = 1;

    e(5*n+2,n-2) = - 1/12*(1)/dx^4;
    e(5*n+2,n-1) = - 1/12*(-4)/dx^4 + N(n-1)*(1)/dx^2 + T(n-1)*(-1)/2/dx;
    e(5*n+2,n) = - 1/12*(-1+6)/dx^4 + N(n-1)*(-2)/dx^2;
    e(5*n+2,n+1) = - 1/12*(2-4)/dx^4 + N(n-1)*(1)/dx^2 + T(n-1)*(1)/2/dx;
    e(5*n+2,3*n+1) = (wn-2*w(n-1)+w(n-2))/dx^2;
    e(5*n+2,4*n+1) = (wn-w(n-2))/2/dx;
    e(5*n+2,5*n+2) = 1;
    
    e(5*n+3,n-1) = - 1/12*(1+1)/dx^4;
    e(5*n+3,n) = - 1/12*(-4+4-4)/dx^4 + Nn*(-1+1)/dx^2 + Tn*(-1-1)/2/dx;
    e(5*n+3,n+1) = - 1/12*(4-4*2+6)/dx^4 + Nn*(2-2)/dx^2 + Tn*(2)/2/dx;
    e(5*n+3,4*n+2) = (wn1-w(n-1))/2/dx;
    e(5*n+3,5*n+3) = 1;
    
    
%Residual vector
    m(1) = w0;
    m(2:n) = w(1:n-1);
    m(n+1) = wn;
    m(n+2) = u0;
    m(n+3:2*n+1) = u(1:n-1);
    m(2*n+2) = un;
    m(2*n+3:3*n+1) = N(1:n-1);
    m(3*n+2) = T0;
    m(3*n+3:4*n+1) = T(1:n-1);
    m(4*n+2) = Tn;
    m(4*n+3) = p0;
    m(4*n+4:5*n+2) = p(1:n-1);
    m(5*n+3) = pn;


    bb(1) = f0;
    bb(2:n) = f(1:n-1);
    bb(n+1) = fn;
    bb(n+2) = g0;
    bb(n+3:2*n+1) = g(1:n-1);
    bb(2*n+2) = gn;
    bb(2*n+3:3*n+1) = Ne(1:n-1);
    bb(3*n+2) = Te0;
    bb(3*n+3:4*n+1) = Te(1:n-1);
    bb(4*n+2) = Ten;
    bb(4*n+3) = pe0;
    bb(4*n+4:5*n+2) = pe(1:n-1);
    bb(5*n+3) = pen;
    
    delta = -e\bb;
    
    R = delta./m;
    fprintf('norm R = %1.5e\n',norm(R))
    fprintf('norm bb = %1.5e\n',norm(bb))
    fprintf('\n') 
   
    if(norm(R)>1e0)
       delta = delta/redu;
    end
    
    m = m + delta;
    
    w0 = m(1);
    w(1:n-1) = m(2:n);
    wn = m(n+1);
    u0 = m(n+2);
    u(1:n-1) = m(n+3:2*n+1);
    un = m(2*n+2);
    N(1:n-1) = m(2*n+3:3*n+1);
    T0 = m(3*n+2);
    T(1:n-1) = m(3*n+3:4*n+1);
    Tn = m(4*n+2);
    p0 = m(4*n+3);
    p(1:n-1) = m(4*n+4:5*n+2);
    pn = m(5*n+3);
    
end

    
    
end

clear delta e m R
 
namemat = strcat('data(t=',num2str(ttotal),').mat');
save(namemat);
% namefig = strcat('profile(p=',num2str(p(mf)),').fig');
% plot(x,ddef)
% saveas(gcf,namefig)

namefig = strcat('profile_u(t=',num2str(ttotal),').fig');
plot(xb/L,u)
saveas(gcf,namefig)
namefig = strcat('profile_T(t=',num2str(ttotal),').fig');
plot(xb/L,T)
saveas(gcf,namefig)
namefig = strcat('profile_N(t=',num2str(ttotal),').fig');
plot(xb/L,N)
saveas(gcf,namefig)
namefig = strcat('profile_p(t=',num2str(ttotal),').fig');
plot(xb/L,p)
saveas(gcf,namefig)
namefig = strcat('profile_w(t=',num2str(ttotal),').fig');
plot(xb/L,w)
saveas(gcf,namefig)




