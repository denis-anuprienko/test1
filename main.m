clear,figure(1),clf,colormap(jet)
%physics
Lx          = 1;     % m
Ly          = 1;     % m
k_etaf      = 1;     % m^2/Pa/s
betaf       = 1;     % 1/Pa
lam         = 1e-2;  % ??
rhoCp       = 1;     % Pa/K
deltaT      = 1;     % K
w           = 1e-1;  % m
alphrhofg   = 1e2;   % Pa/m
%numerics
nx          = 100;
ny          = fix(nx*Ly/Lx);
nout        = 10;
nt          = 100*nout;
%preprocessing
dx          = Lx/(nx-1);
dy          = Ly/(ny-1);
dt_diff     = min(dx,dy)^2/(k_etaf/betaf)/4.1;
[x,y]       = ndgrid(-Lx/2:dx:Lx/2,-Ly/2:dy:Ly/2);
Raylegh     = alphrhofg*Ly*k_etaf/(lam/rhoCp)
%initial conditions
time        = 0;
T           = exp(-(x/w).^2-(y/w).^2);
Pf          = 0*T;
% boundary conditions
T(:,1)      =  deltaT/2;
T(:,end)    = -deltaT/2;
%action
for it = 1:nt   
    % Hydro
    qx                  = - k_eta * diff(Pf,1,1)/dx;
    qy                  = - k_eta * diff(Pf,1,1)/dx;
    
    qy(:,[1,end])       = 0;% boundary conditions

    Pf([1,end],:)       = Pf([2,(end-1)],:);% boundary conditions  

    % Thermo

    T([1,end],:)        = T([2,(end-1)],:);% boundary conditions  

    time                = time + dt;
 
    if mod(it,nout)==0%postprocessing
        pcolor(x,y, T),shading interp,colorbar,axis image,title(['T time = ',num2str(time)]),axis off,drawnow
    end
end