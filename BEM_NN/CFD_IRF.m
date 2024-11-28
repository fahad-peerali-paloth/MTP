function [] = CFD_IRF(fn,lamdabyL,Hbylamda,CFD)
% fn = 0;
% lamdabyL = 1;
% CFD =1;


clc

[TP,nmax,del_t,amp] = creat_inputs(fn,lamdabyL,Hbylamda);
%particulars=[1/B L/B,B/T,CB,CWL)
particulars=[0.06 6.61 3.1781 0.57 0.698];
hydrostaticP=1025*9.81*0.3153;
L=6.6067;
nirf= 100;
% del_t = 0.01;
% disp(strcat('delt == ',num2str(del_t)))%0.01; % time_step
n = nmax; % #time_steps
omegaF = 2*pi()/TP ;
M = 1000*1.1; 
Ry=1.652; %radius of gyration from centroid
Iyy = M*(Ry^2); 
c33 = 41897.461*1.2; %rho*g*Awp;
c35=5251.379; %rho*g*first_moment_of_area
c53=5251.379; %rho*g*first_moment_of_area
c55 =76080.547*1.2; %M*g**GM_L
kf=1;
fHz=1/TP;
b_dist=[0.28186 0.253765 0.224765 0.19486 0.163925 0.13158 0.103515 0.096465 0.05174 0.0010005]; %breadth distribution at the Pressure view points ##PUT THE CURRECT VALUES
h_betwn_b=0.125625606;


% c55
% load('inputdata.mat');
fid1=fopen('excitingF.txt');  %'excitingF'
fid2=fopen('addedmass.txt'); %'addedmass.txt'
fid3=fopen('damp.txt');
num1=fscanf(fid1,'%4d',1);
num2=fscanf(fid2,'%4d',1);
num3=fscanf(fid3,'%4d',1);
F33 = zeros(num1,1);
F55 = zeros(num1,1);
time= zeros(num1,1);
for j = 1:num1
    A = fscanf(fid1,['%f' ''],3);  %3
    %B = fscanf(fid1,['%f' ''],3);
    F33(j,1) = A(2,1);
    time(j,1) = A(1,1);
    F55(j,1) = A(3,1) ; 
end
a33 = zeros(num2,1);
b33 = zeros(num3,1);
a35 = zeros(num2,1);
b35 = zeros(num3,1);
a53 = zeros(num2,1);
b53 = zeros(num3,1);
a55 = zeros(num2,1);
b55 = zeros(num3,1);
omega = zeros(num2,1);
v3 = zeros(nmax+1,1);
v5 = zeros(nmax+1,1);
z3 = zeros(nmax+1,1);
z5 = zeros(nmax+1,1);
for j = 1:num2
    C = fscanf(fid2,['%f' ''],5);  %5
    omega(j,1) = C(1,1); 
    a33(j,1) = C(2,1);
    a35(j,1) = C(3,1);
    a53(j,1) = C(4,1);
    a55(j,1) = C(5,1);
end

for j = 1:num3
B = fscanf(fid3,['%f' ''],5);  %5
b33(j,1) = B(2,1); % B(2,1)
b35(j,1) = B(3,1);
b53(j,1)= B(4,1);
b55(j,1) = B(5,1);
end

ndata = num3;
time2  = time(1:nirf);
B33 = IRF(ndata, nirf,time2, b33,del_t, omega);
B35 = IRF(ndata, nirf,time2, b35,del_t, omega);
B53 = IRF(ndata, nirf,time2, b53,del_t, omega);
B55 = IRF(ndata, nirf,time2, b55,del_t, omega);
A33 = GetAddedMass(omega, del_t,b33,a33,ndata);
A35 = GetAddedMass(omega, del_t,b35,a35,ndata);
A53 = GetAddedMass(omega, del_t,b53,a53,ndata);
A55 = GetAddedMass(omega, del_t,b55,a55,ndata);

z3(1) = 0;
z5(1)=0;
v3(1) = 0;
v5(1)=0;
t(1)=0;
temp_trapz=0;
% % figure
% plot(omega,a35)
% title('a35 vs w')
% % figure
% plot(omega,b35)
% title('b35 vs w')
% figure
% plot(time,F33)
% title('F3 vs t')
% figure
% plot(time,F55)
% title('F5 vs t')
% needed for radiation restoring forces

if fn~=0
    % needed for radiation restoring forces
    a33w = interp1(omega,a33,omegaF);  % calculating a33 freq domain for omegaF by interpolation
    a35w = interp1(omega,a35,omegaF);
    a55w = interp1(omega,a55,omegaF);
    a53w = interp1(omega,a53,omegaF);
    Crad33  = radiation_restoration(a33w,B33,omegaF,time2,A33);
    Crad35  = radiation_restoration(a35w,B35,omegaF,time2,A35); % no need to write twice as it doesnot take velocity as input
    Crad55  = radiation_restoration(a55w,B55,omegaF,time2,A55);
    Crad53 = radiation_restoration(a53w,B53,omegaF,time2,A53);
else
    Crad33  = 0;
    Crad35  = 0; % no need to write twice as it doesnot take velocity as input
    Crad55  = 0;
    Crad53 = 0;
end

%%%%% Area calculation for slammimg force computaion
n_b=length(b_dist);
bow_A=zeros(1,n_b);
for i_b=1:n_b-1
    if (i_b < n_b-1)
        bow_A(i_b)=(h_betwn_b/12)*(5*b_dist(i_b)+8*b_dist(i_b+1)-b_dist(i_b+2));
    else
        bow_A(i_b)=(h_betwn_b/12)*(5*b_dist(i_b)+8*b_dist(i_b-1)-b_dist(i_b-2));
    end
end
%%
fidout = 'RESULTS_OUT.txt';
fopen(fidout ,'wt+');
dlmwrite(fidout,'t  z3  v3  z5  v5  P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13','newline', 'pc')


% xlswrite('results.xlsx',{'t','z3','v3','z5','v5','P'},strcat('fn=',num2str(fn),',lbyL=',num2str(lamdabyL)),'A1')

for i=1:n-1
    t(i) = (i-1)*del_t;
    Damp33 = DampingCoeff(i,nirf,B33,v3, del_t);
    Damp35 = DampingCoeff(i,nirf,B35,v5, del_t);
    Damp53 = DampingCoeff(i,nirf,B53,v3, del_t);
    Damp55 = DampingCoeff(i,nirf,B55,v5, del_t);
    
  

    

 %%
 if CFD ~= 1
      P1(i) = 0;
      P2(i) = 0;
      P3(i) = 0;
      P4(i) = 0;
      P5(i) = 0;
      P6(i) = 0;
      P7(i) = 0;
      P8(i) = 0;
      P9(i) = 0;
      P10(i) = 0;
      P11(i) = 0;
      P12(i) = 0;
      P13(i) = 0;
      cfdforce(i) = 0; % comment while running CFD coupling
 else
      pressure_dat=nnmodel(hydrostaticP,particulars,z3(i),v3(i),z5(i),v5(i),fn,lamdabyL,Hbylamda);
      P1(i) = pressure_dat(1);
      P2(i) = pressure_dat(2);
      P3(i) = pressure_dat(3);
      P4(i) = pressure_dat(4);
      P5(i) = pressure_dat(5);
      P6(i) = pressure_dat(6);
      P7(i) = pressure_dat(7);
      P8(i) = pressure_dat(8);
      P9(i) = pressure_dat(9);
      P10(i) = pressure_dat(10);
      P11(i) = pressure_dat(11);
      P12(i) = pressure_dat(12);
      P13(i) = pressure_dat(13);
      
 % %  fclose(fresult);
      %disp(P(i));
     % cfdforce(i) = (1/0.2*L)*(P1(i)*bow_A(1)+P2(i)*bow_A(2)+P3(i)*bow_A(3)+P4(i)*bow_A(4)+P5(i)*bow_A(5)+P6(i)*bow_A(6)+P7(i)*bow_A(7)+P8(i)*bow_A(8)+P9(i)*bow_A(9));
       cfdforce(i) = (1/(0.2*L))*2*(P1(i)*bow_A(1)+P2(i)*bow_A(2)+P3(i)*bow_A(3)+P4(i)*bow_A(4)+P5(i)*bow_A(5)+P6(i)*bow_A(6)+P7(i)*bow_A(7)+P8(i)*bow_A(8));
 end   
 %%    
     %%%%%%%%%%
     
     
    b3 = F33(i)-c33*z3(i)+ A33*omegaF^2*z3(i)- Damp33+ A35*omegaF^2*z5(i) - del_t*Damp35-c35*z5(i)+cfdforce(i)    - Crad33*z3(i) - Crad35*z5(i);
    b5 = F55(i)-c55*z5(i)+ A55*omegaF^2*z5(i)- Damp55+ A53*omegaF^2*z3(i) - del_t*Damp53-c53*z3(i)+(cfdforce(i)*(2.8945))    - Crad55*z5(i) - Crad53*z3(i);
    
    v3(i+1) = v3(i) + del_t*b3/(M);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z3(i+1) = z3(i) + del_t*v3(i+1);
    
    v5(i+1) = v5(i) + del_t*b5/(Iyy);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z5(i+1) = z5(i) + del_t*v5(i+1);
    
    
    %%%%%%%%%%%%%%%%%%%
    
    Damp33 = DampingCoeff(i+1,nirf,B33,v3, del_t);
    Damp35 = DampingCoeff(i+1,nirf,B35,v5, del_t);
    Damp53 = DampingCoeff(i+1,nirf,B53,v3, del_t);
    Damp55 = DampingCoeff(i+1,nirf,B55,v5, del_t);
    
    b3 = F33(i+1)-c33*z3(i+1)+ A33*omegaF^2*z3(i+1) - Damp33+ A35*omegaF^2*z5(i) - del_t*Damp35-c35*z5(i)+cfdforce(i)  - Crad33*z3(i+1) - Crad35*z5(i+1);
    b5 = F55(i+1)-c55*z5(i+1)+ A55*omegaF^2*z5(i+1) - Damp55+ A53*omegaF^2*z3(i) - del_t*Damp53-c53*z3(i)+(cfdforce(i)*(2.8945))  - Crad55*z5(i+1) - Crad53*z3(i+1);
    v3(i+1) = v3(i) + del_t*b3/(M);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z3(i+1) = z3(i) + del_t*v3(i+1);
    
    v5(i+1) = v5(i) + del_t*b5/(Iyy);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z5(i+1) = z5(i) + del_t*v5(i+1);
    
   %if(i==99)
      % dummy = 1.0;
   %end
    
      a=z3(i+1);
    %disp(z3(i+1));
    
      b=z5(i+1);
    %disp(z5(i+1));  
% xlswrite('results.xlsx',[t(i),z3(i),v3(i),z5(i),v5(i),cfdforce(i)],strcat('fn=',num2str(fn),',lbyL=',num2str(lamdabyL)),strcat('A',num2str(i+1)))
dlmwrite(fidout,[t(i),z3(i),v3(i),z5(i),v5(i),P1(i),P2(i),P3(i),P4(i),P5(i),P6(i),P7(i),P8(i),P9(i),P10(i),P11(i),P12(i),P13(i)],'delimiter','\t','newline', 'pc','precision',18,'-append')
disp(strcat('PERCENT COMPLETED == ',num2str((i/n)*100)))
end

nondz3=z3./amp;
nondz5=z5./amp;
figure(11)
plot(t,z3(1:length(t)))
title('z3')
figure(22)
plot(t,z5(1:length(t)))
title('z5')
figure(1101)
plot(t,nondz3(1:length(t)))
title('nondz3')
figure(2202)
plot(t,nondz5(1:length(t)))
title('nondz5')


end

