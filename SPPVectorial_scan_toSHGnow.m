
%%
%Program follows PRB67 165405 (until near-field and scan part)
clear all;
 
%Positions of scatterers.
%randxy04;	
%laserwrittingXY;
%Chains;
%laserwritting2_a;
laserwritting_C;
load xcc4;
load ycc4;
load Rscat;
load ww;
load x_start;		%nm
load x_end;
load y_start;
load y_end;
load Nscansteps;
load Rprobe;
load ALPHA;
load Limit;
load Rscat;

%file generating positions																		
%coord_y=[linspace(-1500,1500,5)]';%easy line of x scatters
%coord_x=coord_y*0;
%N=200;
if 1==1									%toggle load / don't load file with positions
coord_x=transpose(xcc4);
coord_y=transpose(ycc4);
N=length(coord_x);
									%Number of scatterers
end

%%%%%%%%%%%%%FH constants%%%%%%%%%%%%%%
%dielectric constant (lambda)
file=1;					2			%1 to use matrices in files, 2 to REcalculate these files;
lambda=800;								%Wavelength in nm
eps1=1;
filmmateriale=1;						%Gold 1, Silver 2
%EpsValue;								%find eps2
eps2=-23.11+1.4i;
%eps;								%value from interpolation in program EpsValue;

Rscat=Rscat;								%Radius of the scatterer in nm
Rprobe=Rprobe;
zscatindi=2*Rscat;					%The indirect distance in z from one scatter to another

zdi=Rscat+Rprobe;						%When thought of as spherical scatterers placed on a bulk substrate at z=0 (mirror plane)
														%which gives the indirect scattering
zindi=3*Rscat+Rprobe;				%The indirect distance in z (from probe to zms)

k0=2*pi/lambda;
k0_2=(2*pi/lambda)^2;
kappa=k0*sqrt(eps1*eps2/(eps1+eps2));					%SPP wavevector, kappa_p;
a_zz=kappa/2*i*(1/(sqrt(-eps1*eps2)*(1-(eps1^2/eps2^2))*((eps1+eps2)/(eps1*eps2)))); %seen from eq(25) PRB69-045422, eq(6) PRB67_165405 and a_zz= -3.2e-5+9.32e-4i at 750nm
k_1z=sqrt(k0_2-kappa^2);									%kappa_z
k_1z_kappa=k_1z/kappa;										%used (7*N^2+8*Nsize^2) times - gains 4 seconds with N=200
k_1z_kappa2=k_1z_kappa^2;									%used 3*N^2 + 4*Nsize^2 times
exp_ik_1z2Rscat=exp(i*k_1z*2*Rscat);					%used in Gspp - see below eq 6 in PRB67 165405 (p.2) (2003)
exp_ik_1z3Rscat=exp(i*k_1z*3*Rscat);								%height of source point plus height of observation point
alfa0=4*pi*(Rscat^3)*(eps2-1)/(eps2+2);				%Rayleigh polarizability;
renorm_xx_yy=1/(1-1/8*(eps2-1)/(eps2+2)*(eps2-1)/(eps2+1));
renorm_zz=1/(1-1/4*(eps2-1)/(eps2+2)*(eps2-1)/(eps2+1));
alfa_zz=alfa0*renorm_zz*k0_2;								%the k0^2 originates in from of the sum in eq(1) and (5) in PRB67 165405 
alfa_xx_yy=alfa0*renorm_xx_yy*k0_2;

nflimit=lambda/4;												%below this distance (nm) the nearfield propergator is used
a_nf= -1/(4*pi*k0_2);										%the front (-c^2/(4*pi*omega^2))-constant in D_nf (eq.3.8)from "wiley-paper"
e2pm1=(eps2-1)/(eps2+1);									%setting this to zero turns off indirect term in nearfield (when having no substrate)
n=sqrt(eps2);   												
reflec=(n-1)/(n+1);
reflec_1=reflec+1;											%for FH-incident field reflected from substrate

%%%%%%%%%%%%%FH-matrices%%%%%%%%%%%%%%%%

tic										%measure time from here (end with time=toc)

if 2==2,		%select true if theta should be recalculated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Theta(1:3*N,1:3*N)=0;					%For calculation of elements of green between scatterers	 
for I=1:N
   for J=1:N
      if I==J
      	Theta(I,J)=1;					%xx											
      	Theta(I+N,J+N)=1;				%yy
      	Theta(I+2*N, J+2*N)=1;		%zz
      	Theta(I,J+N)=0;				%xy
      	Theta(I,J+2*N)=0;				%xz
      	Theta(I+N,J)=0;				%yx
      	Theta(I+N,J+2*N)=0;			%yz
      	Theta(I+2*N,J)=0;				%zx	
      	Theta(I+2*N,J+N)=0;			%zy
   	else
      
      	x=(coord_x(I)-coord_x(J));
         y=(coord_y(I)-coord_y(J));
         
         rdi2=x^2+y^2;						%direct distances between scatterers
         rdi=sqrt(rdi2);
         
         if rdi<nflimit
            
         	rdi3=rdi2*rdi;
         	rdi5=rdi3*rdi2;
         	rindi2=x^2+y^2+(zscatindi)^2;	%indirect distances
         	rindi=sqrt(rindi2);
         	rindi3=rindi2*rindi;
         	rindi5=rindi3*rindi2; 
            
            Theta(I,J)=a_nf*((3*x*x/rdi2-1)/rdi3 -((3*x*x/rindi2-1)/rindi3)*e2pm1)*alfa_xx_yy; 		%xx
      		Theta(I,J+N)=a_nf*(3*x*y/rdi5 -(3*x*y/rindi5)*e2pm1)*alfa_xx_yy;    							%xy
      		Theta(I,J+2*N)=a_nf*(3*x*zscatindi/rindi5*e2pm1)*alfa_zz;  										%xz
            Theta(I+N,J)=Theta(I,J+N);%a_nf*(3*y*x/rdi5 -(3*y*x/rindi5)*e2pm1)*alfa_xx_yy;     		%yx %the same as Theta(I,J+N)
            Theta(I+N,J+N)=a_nf*((3*y*y/rdi2-1)/rdi3 -((3*y*y/rindi2-1)/rindi3)*e2pm1)*alfa_xx_yy;   	%yy
      		Theta(I+N,J+2*N)=a_nf*(3*y*zscatindi/rindi5*e2pm1)*alfa_zz;										%yz
      		Theta(I+2*N,J)=a_nf*(-3*x*zscatindi/rindi5*e2pm1)*alfa_xx_yy;									%zx    
      		Theta(I+2*N,J+N)=a_nf*(-3*y*zscatindi/rindi5*e2pm1)*alfa_xx_yy;								%zy
      		Theta(I+2*N,J+2*N)=a_nf*(-1/rdi3+((3*zscatindi*zscatindi/rindi2-1)/rindi3)*e2pm1)*alfa_zz;%zz
         else   
         	h=a_zz*exp_ik_1z2Rscat*besselh(0,kappa*rdi);
      		Theta(I,J)=h*k_1z_kappa2*x*x/rdi2*alfa_xx_yy;      			%xx
      		Theta(I,J+N)=h*k_1z_kappa2*x*y/rdi2*alfa_xx_yy;    			%xy
      		Theta(I,J+2*N)=h*k_1z_kappa*x/rdi*alfa_zz;      				%xz
            Theta(I+N,J)=Theta(I,J+N);%h*k_1z_kappa2*y*x/rdi/rdi*alfa_xx_yy;     %yx	%the same as Theta(I,J+N)
            Theta(I+N,J+N)=h*k_1z_kappa2*y*y/rdi2*alfa_xx_yy;     		%yy
      		Theta(I+N,J+2*N)=h*k_1z_kappa*y/rdi*alfa_zz;						%yz
      		Theta(I+2*N,J)=-h*k_1z_kappa*x/rdi*alfa_xx_yy;					%zx    
      		Theta(I+2*N,J+N)=-h*k_1z_kappa*y/rdi*alfa_xx_yy;				%zy
      		Theta(I+2*N,J+2*N)=-h*alfa_zz;										%zz
         end   
     	end
   end
   I  
end
ThetaINV=Theta^(-1);
clear Theta;
ThetaINVr=real(ThetaINV); save ThetaINVr.dat ThetaINVr -ascii; clear ThetaINVr;
ThetaINVi=imag(ThetaINV); save ThetaINVi.dat ThetaINVi -ascii; clear ThetaINVi;
clear ThetaINV;
end	%end of if file...
time_for_theta=toc

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 2==2,
for pp=1:N
 	x_vec=-coord_x+coord_x(pp);				%vectors with distances from all scatterers to the pp'te scatter (the observation points)
 	y_vec=-coord_y+coord_y(pp);
  	r2=(x_vec).^2+(y_vec).^2+0.00001;			%(as seen from far-field), add 0.0001 to awoid (0/0)-situation...  
   r=sqrt(r2);
         
   rdi=sqrt((x_vec).^2+(y_vec).^2+(zdi)^2);	%obtained at a height z=Rscat+Rprobe.
   rdi2=rdi.^2;
   rdi3=rdi2.*rdi;
   rdi5=rdi3.*rdi2;
   rindi2=x_vec.^2+y_vec.^2+(zindi)^2;			%indirect distances
   rindi=sqrt(rindi2);
   rindi3=rindi2.*rindi;
   rindi5=rindi3.*rindi2;
       
  	%%%%%the near-field green propagator%%%%%
  	gxxnf(:,pp)=-a_nf*((3*x_vec.*x_vec./rdi2-1)./rdi3 -((3*x_vec.*x_vec./rindi2-1)./rindi3)*e2pm1)*alfa_xx_yy;	%xx	%check signs for all these elements, oposit as below with theta and gxxff
 	gxynf(:,pp)=-a_nf*(3*x_vec.*y_vec./rdi5 -(3*x_vec.*y_vec./rindi5)*e2pm1)*alfa_xx_yy;    				%xy
 	gxznf(:,pp)=-a_nf*(3*x_vec.*zdi./rdi5 +(3*x_vec*zindi./rindi5)*e2pm1)*alfa_zz;  			  				%xz
  	gyxnf(:,pp)=gxynf(:,pp);%a_nf*(3*y_vec.*x_vec./rdi5 -(3*y_vec.*x_vec./rindi5)*e2pm1)*alfa_xx_yy;     	%yx	%the same as gxynf
  	gyynf(:,pp)=-a_nf*((3*y_vec.*y_vec./rdi2-1)./rdi3 -((3*y_vec.*y_vec./rindi2-1)./rindi3)*e2pm1)*alfa_xx_yy; 	%yy
 	gyznf(:,pp)=-a_nf*(3*y_vec.*zdi./rdi5 +(3*y_vec*zindi./rindi5)*e2pm1)*alfa_zz;								%yz
   gzxnf(:,pp)=-a_nf*(3*zdi.*x_vec./rdi5 -(3*x_vec*zindi./rindi5)*e2pm1)*alfa_xx_yy;							%zx    
   gzynf(:,pp)=-a_nf*(3*zdi.*y_vec./rdi5 -(3*y_vec*zindi./rindi5)*e2pm1)*alfa_xx_yy;							%zy
   gzznf(:,pp)=-a_nf*((3*zdi*zdi./rdi2-1)./rdi3+((3*zindi*zindi./rindi2-1)./rindi3)*e2pm1)*alfa_zz;			%zz
              
   %%%%%the far-field green propagator%%%%%  
   h=a_zz*exp_ik_1z3Rscat*besselh(0,kappa*r);
   gxxff(:,pp)=-h.*k_1z_kappa2.*x_vec.*x_vec./r2*alfa_xx_yy;				%xx
   gxyff(:,pp)=-h.*k_1z_kappa2.*x_vec.*y_vec./r2*alfa_xx_yy;				%xy
   gxzff(:,pp)=-h.*k_1z_kappa.*x_vec./r*alfa_zz;								%xz
   gyxff(:,pp)=gxyff(:,pp);%-h.*k_1z_kappa2.*x_vec.*y_vec./r./r*alfa_xx_yy;	%yx	%the same as gxyff
   gyyff(:,pp)=-h.*k_1z_kappa2.*y_vec.*y_vec./r2*alfa_xx_yy;				%yy
   gyzff(:,pp)=-h.*k_1z_kappa.*y_vec./r*alfa_zz;								%yz
   gzxff(:,pp)=h.*k_1z_kappa.*x_vec./r*alfa_xx_yy;							%zx
   gzyff(:,pp)=h.*k_1z_kappa.*y_vec./r*alfa_xx_yy;							%zy
   gzzff(:,pp)=h*alfa_zz;
            
  	for q=1:N,					%maybe select between r before calculation, but then calc as seperate values instead of vector..?
  		if r(q)<nflimit
  			gxx(q,pp)=gxxnf(q,pp);
  			gyx(q,pp)=gyxnf(q,pp);
   		gzx(q,pp)=gzxnf(q,pp);
   		gxy(q,pp)=gxynf(q,pp);
   		gyy(q,pp)=gyynf(q,pp);
   		gzy(q,pp)=gzynf(q,pp);
   		gxz(q,pp)=gxznf(q,pp);
   		gyz(q,pp)=gyznf(q,pp);
   		gzz(q,pp)=gzznf(q,pp);
   	else
   		gxx(q,pp)=gxxff(q,pp);
   		gyx(q,pp)=gyxff(q,pp);
   		gzx(q,pp)=gzxff(q,pp);
   		gxy(q,pp)=gxyff(q,pp);
   		gyy(q,pp)=gyyff(q,pp);
    		gzy(q,pp)=gzyff(q,pp);
     		gxz(q,pp)=gxzff(q,pp);
   		gyz(q,pp)=gyzff(q,pp);
    		gzz(q,pp)=gzzff(q,pp);
   	end
   end
end
gxxr=real(gxx); save gxxr.dat gxxr -ascii; clear gxxr;
gyxr=real(gyx); save gyxr.dat gyxr -ascii; clear gyxr;
gzxr=real(gzx); save gzxr.dat gzxr -ascii; clear gzxr;
gxyr=real(gxy); save gxyr.dat gxyr -ascii; clear gxyr;
gyyr=real(gyy); save gyyr.dat gyyr -ascii; clear gyyr;
gzyr=real(gzy); save gzyr.dat gzyr -ascii; clear gzyr;
gxzr=real(gxz); save gxzr.dat gxzr -ascii; clear gxzr;
gyzr=real(gyz); save gyzr.dat gyzr -ascii; clear gyzr;
gzzr=real(gzz); save gzzr.dat gzzr -ascii; clear gzzr;

gxxi=imag(gxx); save gxxi.dat gxxi -ascii; clear gxxi;
gyxi=imag(gyx); save gyxi.dat gyxi -ascii; clear gyxi;
gzxi=imag(gzx); save gzxi.dat gzxi -ascii; clear gzxi;
gxyi=imag(gxy); save gxyi.dat gxyi -ascii; clear gxyi;
gyyi=imag(gyy); save gyyi.dat gyyi -ascii; clear gyyi;
gzyi=imag(gzy); save gzyi.dat gzyi -ascii; clear gzyi;
gxzi=imag(gxz); save gxzi.dat gxzi -ascii; clear gxzi;
gyzi=imag(gyz); save gyzi.dat gyzi -ascii; clear gyzi;
gzzi=imag(gzz); save gzzi.dat gzzi -ascii; clear gzzi;

clear gxx; clear gyx; clear gzx; clear gxy; clear gyy; 
clear gzy; clear gxz; clear gyz; clear gzz;
end			%end of if file...
time_for_gxx=toc

%%%%%%%%%%%%%%%%Second-harmonic constants%%%%%%%%%%%%%%

lambda=lambda/2;					%Wavelength in nm
eps1=1;
filmmateriale=1;						%Gold 1, Silver 2
%EpsValue;								%find eps2
eps2=eps2;								%value from interpolation in program EpsValue;

%Rscat=40;								%Radius of the scatterer in nm
%Rprobe=0;
zscatindi=2*Rscat;					%The indirect distance in z from one scatter to another

zdi=Rscat+Rprobe;						%When thought of as spherical scatterers placed on a bulk substrate at z=0 (mirror plane)
														%which gives the indirect scattering
zindi=3*Rscat+Rprobe;				%The indirect distance in z (from probe to zms)

k0=2*pi/lambda;
k0_2=(2*pi/lambda)^2;
kappa=k0*sqrt(eps1*eps2/(eps1+eps2));					%SPP wavevector, kappa_p;
a_zz=kappa/2*i*(1/(sqrt(-eps1*eps2)*(1-(eps1^2/eps2^2))*((eps1+eps2)/(eps1*eps2)))); %seen from eq(25) PRB69-045422, eq(6) PRB67_165405 and a_zz= -3.2e-5+9.32e-4i at 750nm
k_1z=sqrt(k0_2-kappa^2);									%kappa_z
k_1z_kappa=k_1z/kappa;										%used (7*N^2+8*Nsize^2) times - gains 4 seconds with N=200
k_1z_kappa2=k_1z_kappa^2;									%used 3*N^2 + 4*Nsize^2 times
exp_ik_1z2Rscat=exp(i*k_1z*2*Rscat);					%used in Gspp - see below eq 6 in PRB67 165405 (p.2) (2003)
exp_ik_1z3Rscat=exp(i*k_1z*3*Rscat);								%height of source point plus height of observation point
alfa0=4*pi*(Rscat^3)*(eps2-1)/(eps2+2);				%Rayleigh polarizability;
renorm_xx_yy=1/(1-1/8*(eps2-1)/(eps2+2)*(eps2-1)/(eps2+1));
renorm_zz=1/(1-1/4*(eps2-1)/(eps2+2)*(eps2-1)/(eps2+1));
alfa_zz=alfa0*renorm_zz*k0_2;								%the k0^2 originates in from of the sum in eq(1) and (5) in PRB67 165405 
alfa_xx_yy=alfa0*renorm_xx_yy*k0_2;

nflimit=lambda/4;												%below this distance (nm) the nearfield propergator is used
a_nf= -1/(4*pi*k0_2);										%the front (-c^2/(4*pi*omega^2))-constant in D_nf (eq.3.8)from "wiley-paper"
e2pm1=(eps2-1)/(eps2+1);


%%%%%%%%%%%%%%%%%%The Second-harmonic matrices%%%%%%%%%%%%%%%
tic

if 2==2,		%select true if theta should be recalculated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Theta(1:3*N,1:3*N)=0;					%For calculation of elements of green between scatterers	 
for I=1:N
   for J=1:N
      if I==J
      	Theta(I,J)=1;					%xx											
      	Theta(I+N,J+N)=1;				%yy
      	Theta(I+2*N, J+2*N)=1;		%zz
      	Theta(I,J+N)=0;				%xy
      	Theta(I,J+2*N)=0;				%xz
      	Theta(I+N,J)=0;				%yx
      	Theta(I+N,J+2*N)=0;			%yz
      	Theta(I+2*N,J)=0;				%zx	
      	Theta(I+2*N,J+N)=0;			%zy
   	else
      
      	x=(coord_x(I)-coord_x(J));
         y=(coord_y(I)-coord_y(J));
         
         rdi2=x^2+y^2;						%direct distances between scatterers
         rdi=sqrt(rdi2);
         
         if rdi<nflimit
            
         	rdi3=rdi2*rdi;
         	rdi5=rdi3*rdi2;
         	rindi2=x^2+y^2+(zscatindi)^2;	%indirect distances
         	rindi=sqrt(rindi2);
         	rindi3=rindi2*rindi;
         	rindi5=rindi3*rindi2; 
            
            Theta(I,J)=a_nf*((3*x*x/rdi2-1)/rdi3 -((3*x*x/rindi2-1)/rindi3)*e2pm1)*alfa_xx_yy; 		%xx
      		Theta(I,J+N)=a_nf*(3*x*y/rdi5 -(3*x*y/rindi5)*e2pm1)*alfa_xx_yy;    							%xy
      		Theta(I,J+2*N)=a_nf*(3*x*zscatindi/rindi5*e2pm1)*alfa_zz;  										%xz
            Theta(I+N,J)=Theta(I,J+N);%a_nf*(3*y*x/rdi5 -(3*y*x/rindi5)*e2pm1)*alfa_xx_yy;     		%yx %the same as Theta(I,J+N)
            Theta(I+N,J+N)=a_nf*((3*y*y/rdi2-1)/rdi3 -((3*y*y/rindi2-1)/rindi3)*e2pm1)*alfa_xx_yy;   	%yy
      		Theta(I+N,J+2*N)=a_nf*(3*y*zscatindi/rindi5*e2pm1)*alfa_zz;										%yz
      		Theta(I+2*N,J)=a_nf*(-3*x*zscatindi/rindi5*e2pm1)*alfa_xx_yy;									%zx    
      		Theta(I+2*N,J+N)=a_nf*(-3*y*zscatindi/rindi5*e2pm1)*alfa_xx_yy;								%zy
      		Theta(I+2*N,J+2*N)=a_nf*(-1/rdi3+((3*zscatindi*zscatindi/rindi2-1)/rindi3)*e2pm1)*alfa_zz;%zz
         else   
         	h=a_zz*exp_ik_1z2Rscat*besselh(0,kappa*rdi);
      		Theta(I,J)=h*k_1z_kappa2*x*x/rdi2*alfa_xx_yy;      			%xx
      		Theta(I,J+N)=h*k_1z_kappa2*x*y/rdi2*alfa_xx_yy;    			%xy
      		Theta(I,J+2*N)=h*k_1z_kappa*x/rdi*alfa_zz;      				%xz
            Theta(I+N,J)=Theta(I,J+N);%h*k_1z_kappa2*y*x/rdi/rdi*alfa_xx_yy;     %yx	%the same as Theta(I,J+N)
            Theta(I+N,J+N)=h*k_1z_kappa2*y*y/rdi2*alfa_xx_yy;     		%yy
      		Theta(I+N,J+2*N)=h*k_1z_kappa*y/rdi*alfa_zz;						%yz
      		Theta(I+2*N,J)=-h*k_1z_kappa*x/rdi*alfa_xx_yy;					%zx    
      		Theta(I+2*N,J+N)=-h*k_1z_kappa*y/rdi*alfa_xx_yy;				%zy
      		Theta(I+2*N,J+2*N)=-h*alfa_zz;										%zz
         end   
     	end
   end
   I  
end
ThetaINV=Theta^(-1);
clear Theta;
ThetaINVr=real(ThetaINV); save ThetaINVrSH.dat ThetaINVr -ascii; clear ThetaINVr;
ThetaINVi=imag(ThetaINV); save ThetaINViSH.dat ThetaINVi -ascii; clear ThetaINVi;
clear ThetaINV;
end	%end of if file...
time_for_thetaSH=toc

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 2==2,
for pp=1:N
 	x_vec=-coord_x+coord_x(pp);				%vectors with distances from all scatterers to the pp'te scatter (the observation points)
 	y_vec=-coord_y+coord_y(pp);
  	r2=(x_vec).^2+(y_vec).^2+0.00001;			%(as seen from far-field), add 0.0001 to awoid (0/0)-situation...  
   r=sqrt(r2);
         
   rdi=sqrt((x_vec).^2+(y_vec).^2+(zdi)^2);	%obtained at a height z=Rscat+Rprobe.
   rdi2=rdi.^2;
   rdi3=rdi2.*rdi;
   rdi5=rdi3.*rdi2;
   rindi2=x_vec.^2+y_vec.^2+(zindi)^2;			%indirect distances
   rindi=sqrt(rindi2);
   rindi3=rindi2.*rindi;
   rindi5=rindi3.*rindi2;
       
  	%%%%%the near-field green propagator%%%%%
  	gxxnf(:,pp)=-a_nf*((3*x_vec.*x_vec./rdi2-1)./rdi3 -((3*x_vec.*x_vec./rindi2-1)./rindi3)*e2pm1)*alfa_xx_yy;	%xx	%check signs for all these elements, oposit as below with theta and gxxff
 	gxynf(:,pp)=-a_nf*(3*x_vec.*y_vec./rdi5 -(3*x_vec.*y_vec./rindi5)*e2pm1)*alfa_xx_yy;    				%xy
 	gxznf(:,pp)=-a_nf*(3*x_vec.*zdi./rdi5 +(3*x_vec*zindi./rindi5)*e2pm1)*alfa_zz;  			  				%xz
  	gyxnf(:,pp)=gxynf(:,pp);%a_nf*(3*y_vec.*x_vec./rdi5 -(3*y_vec.*x_vec./rindi5)*e2pm1)*alfa_xx_yy;     	%yx	%the same as gxynf
  	gyynf(:,pp)=-a_nf*((3*y_vec.*y_vec./rdi2-1)./rdi3 -((3*y_vec.*y_vec./rindi2-1)./rindi3)*e2pm1)*alfa_xx_yy; 	%yy
 	gyznf(:,pp)=-a_nf*(3*y_vec.*zdi./rdi5 +(3*y_vec*zindi./rindi5)*e2pm1)*alfa_zz;								%yz
   gzxnf(:,pp)=-a_nf*(3*zdi.*x_vec./rdi5 -(3*x_vec*zindi./rindi5)*e2pm1)*alfa_xx_yy;							%zx    
   gzynf(:,pp)=-a_nf*(3*zdi.*y_vec./rdi5 -(3*y_vec*zindi./rindi5)*e2pm1)*alfa_xx_yy;							%zy
   gzznf(:,pp)=-a_nf*((3*zdi*zdi./rdi2-1)./rdi3+((3*zindi*zindi./rindi2-1)./rindi3)*e2pm1)*alfa_zz;			%zz
              
   %%%%%the far-field green propagator%%%%%  
   h=a_zz*exp_ik_1z3Rscat*besselh(0,kappa*r);
   gxxff(:,pp)=-h.*k_1z_kappa2.*x_vec.*x_vec./r2*alfa_xx_yy;				%xx
   gxyff(:,pp)=-h.*k_1z_kappa2.*x_vec.*y_vec./r2*alfa_xx_yy;				%xy
   gxzff(:,pp)=-h.*k_1z_kappa.*x_vec./r*alfa_zz;								%xz
   gyxff(:,pp)=gxyff(:,pp);%-h.*k_1z_kappa2.*x_vec.*y_vec./r./r*alfa_xx_yy;	%yx	%the same as gxyff
   gyyff(:,pp)=-h.*k_1z_kappa2.*y_vec.*y_vec./r2*alfa_xx_yy;				%yy
   gyzff(:,pp)=-h.*k_1z_kappa.*y_vec./r*alfa_zz;								%yz
   gzxff(:,pp)=h.*k_1z_kappa.*x_vec./r*alfa_xx_yy;							%zx
   gzyff(:,pp)=h.*k_1z_kappa.*y_vec./r*alfa_xx_yy;							%zy
   gzzff(:,pp)=h*alfa_zz;
            
  	for q=1:N,					%maybe select between r before calculation, but then calc as seperate values instead of vector..?
  		if r(q)<nflimit
  			gxx(q,pp)=gxxnf(q,pp);
  			gyx(q,pp)=gyxnf(q,pp);
   		gzx(q,pp)=gzxnf(q,pp);
   		gxy(q,pp)=gxynf(q,pp);
   		gyy(q,pp)=gyynf(q,pp);
   		gzy(q,pp)=gzynf(q,pp);
   		gxz(q,pp)=gxznf(q,pp);
   		gyz(q,pp)=gyznf(q,pp);
   		gzz(q,pp)=gzznf(q,pp);
   	else
   		gxx(q,pp)=gxxff(q,pp);
   		gyx(q,pp)=gyxff(q,pp);
   		gzx(q,pp)=gzxff(q,pp);
   		gxy(q,pp)=gxyff(q,pp);
   		gyy(q,pp)=gyyff(q,pp);
    		gzy(q,pp)=gzyff(q,pp);
     		gxz(q,pp)=gxzff(q,pp);
   		gyz(q,pp)=gyzff(q,pp);
    		gzz(q,pp)=gzzff(q,pp);
   	end
   end
end
gxxr=real(gxx); save gxxrSH.dat gxxr -ascii; clear gxxr;
gyxr=real(gyx); save gyxrSH.dat gyxr -ascii; clear gyxr;
gzxr=real(gzx); save gzxrSH.dat gzxr -ascii; clear gzxr;
gxyr=real(gxy); save gxyrSH.dat gxyr -ascii; clear gxyr;
gyyr=real(gyy); save gyyrSH.dat gyyr -ascii; clear gyyr;
gzyr=real(gzy); save gzyrSH.dat gzyr -ascii; clear gzyr;
gxzr=real(gxz); save gxzrSH.dat gxzr -ascii; clear gxzr;
gyzr=real(gyz); save gyzrSH.dat gyzr -ascii; clear gyzr;
gzzr=real(gzz); save gzzrSH.dat gzzr -ascii; clear gzzr;

gxxi=imag(gxx); save gxxiSH.dat gxxi -ascii; clear gxxi;
gyxi=imag(gyx); save gyxiSH.dat gyxi -ascii; clear gyxi;
gzxi=imag(gzx); save gzxiSH.dat gzxi -ascii; clear gzxi;
gxyi=imag(gxy); save gxyiSH.dat gxyi -ascii; clear gxyi;
gyyi=imag(gyy); save gyyiSH.dat gyyi -ascii; clear gyyi;
gzyi=imag(gzy); save gzyiSH.dat gzyi -ascii; clear gzyi;
gxzi=imag(gxz); save gxziSH.dat gxzi -ascii; clear gxzi;
gyzi=imag(gyz); save gyziSH.dat gyzi -ascii; clear gyzi;
gzzi=imag(gzz); save gzziSH.dat gzzi -ascii; clear gzzi;

clear gxx; clear gyx; clear gzx; clear gxy; clear gyy; 
clear gzy; clear gxz; clear gyz; clear gzz;
end			%end of if file...

time_for_gxxSH=toc

%%%%%%%%%Load saved matrices%%%%%%%%

tic

load ThetaINVr.dat; load ThetaINVi.dat;                 
ThetaINV=ThetaINVr+i*ThetaINVi;
clear ThetaINVr; clear ThetaINVi;

toc

load ThetaINVrSH.dat; load ThetaINViSH.dat;                 
ThetaINVSH=ThetaINVrSH+i*ThetaINViSH;
clear ThetaINVrSH; clear ThetaINViSH;

toc

load gxxr.dat; load gxxi.dat; gxx=gxxr+i*gxxi; clear gxxr; clear gxxi;
load gyxr.dat; load gyxi.dat; gyx=gyxr+i*gyxi; clear gyxr; clear gyxi;  
load gzxr.dat; load gzxi.dat; gzx=gzxr+i*gzxi; clear gzxr; clear gzxi;  
load gxyr.dat; load gxyi.dat; gxy=gxyr+i*gxyi; clear gxyr; clear gxyi;
load gyyr.dat; load gyyi.dat; gyy=gyyr+i*gyyi; clear gyyr; clear gyyi;  
load gzyr.dat; load gzyi.dat; gzy=gzyr+i*gzyi; clear gzyr; clear gzyi;  
load gxzr.dat; load gxzi.dat; gxz=gxzr+i*gxzi; clear gxzr; clear gxzi;
load gyzr.dat; load gyzi.dat; gyz=gyzr+i*gyzi; clear gyzr; clear gyzi;  
load gzzr.dat; load gzzi.dat; gzz=gzzr+i*gzzi; clear gzzr; clear gzzi;   

toc

load gxxrSH.dat; load gxxiSH.dat; gxxSH=gxxrSH+i*gxxiSH; clear gxxrSH; clear gxxiSH;
load gyxrSH.dat; load gyxiSH.dat; gyxSH=gyxrSH+i*gyxiSH; clear gyxrSH; clear gyxiSH;  
load gzxrSH.dat; load gzxiSH.dat; gzxSH=gzxrSH+i*gzxiSH; clear gzxrSH; clear gzxiSH;  
load gxyrSH.dat; load gxyiSH.dat; gxySH=gxyrSH+i*gxyiSH; clear gxyrSH; clear gxyiSH;
load gyyrSH.dat; load gyyiSH.dat; gyySH=gyyrSH+i*gyyiSH; clear gyyrSH; clear gyyiSH;  
load gzyrSH.dat; load gzyiSH.dat; gzySH=gzyrSH+i*gzyiSH; clear gzyrSH; clear gzyiSH;  
load gxzrSH.dat; load gxziSH.dat; gxzSH=gxzrSH+i*gxziSH; clear gxzrSH; clear gxziSH;
load gyzrSH.dat; load gyziSH.dat; gyzSH=gyzrSH+i*gyziSH; clear gyzrSH; clear gyziSH;  
load gzzrSH.dat; load gzziSH.dat; gzzSH=gzzrSH+i*gzziSH; clear gzzrSH; clear gzziSH;   

toc

%%%%%%%%%%%%%%%The Scanning Process%%%%%%%%%%%%%%%%%%%

%INCIDENT FIELD - A Gaussian Beam Scanning the sample%%%%%%
ww=ww;				%beam waist ~ww in nm
ww2=ww^2;			%squared beam waist in nm. 
  
x_start=x_start;		%nm
x_end=x_end;
y_start=-y_start;
y_end=y_end;
   
Nscansteps=Nscansteps;
  
x=1:Nscansteps;
x=x/Nscansteps;
x=x-0.5;
x=x*(x_end-x_start);
y=x;

E0y(1:N)=0.0; E0y=transpose(E0y);
%E0x(1:N)=0.0; E0x=transpose(E0x);
E0z(1:N)=0.0; E0z=transpose(E0z);
alphax=100;
alphay=0;
%alphaz=0;

intensityX(1:Nscansteps,1:Nscansteps)=0;
intensityY(1:Nscansteps,1:Nscansteps)=0;
%intensityZ(1:Nscansteps,1:Nscansteps)=0;	%not detected for FH

E0xSH(1:N)=0.0; E0xSH=transpose(E0xSH);
E0ySH(1:N)=0.0; E0ySH=transpose(E0ySH);

intensityXSH(1:Nscansteps,1:Nscansteps)=0;
intensityYSH(1:Nscansteps,1:Nscansteps)=0;

tic
   
for p=1:Nscansteps
    x_s=x(p);
    for q=1:Nscansteps
        y_s=y(q);									%Position of center of the scanning beam
        E0x=exp(-((coord_x-x_s).^2+(coord_y-y_s).^2)/ww2)*reflec_1;	%E components at the position of scatterers
        E0xyz=[E0x;E0y;E0z];
        
        Exyz=ThetaINV*E0xyz;		
		  Ex=Exyz(1:N);							%x,y,z-components of the electric field at the site of scatterers
		  Ey=Exyz(N+1:2*N);
        Ez=Exyz(2*N+1:3*N);
         
        for pp=1:N
            Gtx=sum(gxx(:,pp).*Ex)+sum(gxy(:,pp).*Ey)+sum(gxz(:,pp).*Ez);	
            Gty=sum(gyx(:,pp).*Ex)+sum(gyy(:,pp).*Ey)+sum(gyz(:,pp).*Ez);
            Gtz=sum(gzx(:,pp).*Ex)+sum(gzy(:,pp).*Ey)+sum(gzz(:,pp).*Ez);	%used to estimate incident (generated main component) of SH-field
            
            Etx(pp)=Gtx;		%used to plot x-polarised FH-light
            Ety(pp)=Gty;		%used to plot y-polarised FH-light
            Etz(pp)=Gtz;		%used to estimate incident (generated main component) of SH-field at surface under each scatterer.
        end
         
        ETOTx=alphax+sum(Etx);
        ETOTy=alphay+sum(Ety);
        %ETOTz=alphaz+sum(Etz);	%not plotted/used (never detected) for FH-light
        intensityX(p,q)=sum(ETOTx*conj(ETOTx));	%with intensity(I,J) x and y axes are properbly interchanged
        intensityY(p,q)=sum(ETOTy*conj(ETOTy));
        %intensityZ(p,q)=sum(ETOTz*conj(ETOTz));	%not plotted/used (never detected) for FH-light
        
        %%%%%%%%%%%%%%Second-harmonic generation%%%%%%%%%%%%%%%%
        
        E0zSH=Etz.^2;		%relation should be fitted with an additional constant in front instead of only squaring
        E0xyz=[E0xSH;E0ySH;E0zSH'];
        
        Exyz=ThetaINVSH*E0xyz;		
		  Ex=Exyz(1:N);							%x,y,z-components of the SH-electric fields at the site of scatterers
		  Ey=Exyz(N+1:2*N);
        Ez=Exyz(2*N+1:3*N);

        for pp=1:N
            Gtx=sum(gxxSH(:,pp).*Ex)+sum(gxySH(:,pp).*Ey)+sum(gxzSH(:,pp).*Ez);	%used to calculate new field distribution of SH
            Gty=sum(gyxSH(:,pp).*Ex)+sum(gyySH(:,pp).*Ey)+sum(gyzSH(:,pp).*Ez);
            %Gtz=sum(gzxSH(:,pp).*Ex)+sum(gzySH(:,pp).*Ey)+sum(gzzSH(:,pp).*Ez);
            
            Etx(pp)=Gtx;		%used to plot x-polarised SH-light
            Ety(pp)=Gty;		%used to plot y-polarised SH-light
            %Etz(pp)=Gtz;		%used to plot z-polarised SH-light
        end
        
        ETOTx=sum(Etx);
        ETOTy=sum(Ety);
        %ETOTz=alphaz+sum(Etz);	%not plotted/used (never detected) for SH-light
        intensityXSH(p,q)=sum(ETOTx*conj(ETOTx));	%with intensity(I,J) x and y axes are properbly interchanged
        intensityYSH(p,q)=sum(ETOTy*conj(ETOTy));

    end
    p
    time_left=toc*(Nscansteps)-toc*p
    tic
end   
if 1==1
figure
h=pcolor(x,y,intensityX);shading interp; colorbar;
colormap;
hold on
plot (coord_y,coord_x,'w.')
axis image;

figure
h=pcolor(x,y,intensityY);shading interp; colorbar;
colormap;
hold on
plot (coord_y,coord_x,'w.') 
axis image;

figure
h=pcolor(x,y,intensityXSH);shading interp; colorbar;
colormap;
hold on
plot (coord_y,coord_x,'w.')
axis image;

figure
h=pcolor(x,y,intensityYSH);shading interp; colorbar;
colormap;
hold on
plot (coord_y,coord_x,'w.') 
axis image;
else
end
if 1==2
% Create the first subplot and plot the first figure
subplot(2, 2, 1)
h = pcolor(x, y, intensityX);
shading interp
colorbar
colormap
hold on
plot(coord_y, coord_x, 'w.')
axis image
title('IntensityX')

% Create the second subplot and plot the second figure
subplot(2, 2, 2)
h = pcolor(x, y, intensityY);
shading interp
colorbar
colormap
hold on
plot(coord_y, coord_x, 'w.')
axis image
title('IntensityY')

% Create the third subplot and plot the third figure
subplot(2, 2, 3)
h = pcolor(x, y, intensityXSH);
shading interp
colorbar
colormap
hold on
plot(coord_y, coord_x, 'w.')
axis image
title('IntensityXSH')

% Create the fourth subplot and plot the fourth figure
subplot(2, 2, 4)
h = pcolor(x, y, intensityYSH);
shading interp
colorbar
colormap
hold on
plot(coord_y, coord_x, 'w.')
axis image
title('IntensityYSH')
else
end
% Set the same y-axis limits for all subplots
%ylim([1, 10])








save AU_N100_Cross_y_x.dat intensityX -ascii;		%implement filnames with automatic numbering everywhere
save AU_N100_Cross_y_y.dat intensityY -ascii;
%save AU_N100_Cross_y_z.dat intensityZ -ascii;
   
save AU_N100_Cross_y_xSH.dat intensityXSH -ascii;		%implement filnames with automatic numbering everywhere
save AU_N100_Cross_y_ySH.dat intensityYSH -ascii;
   

