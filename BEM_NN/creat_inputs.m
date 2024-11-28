function [TP,nmax,del_t,amp] = creat_inputs(fn,lamdabyL,Hbylamda)
L  = 6.6067;
lamda=lamdabyL*L;
H=Hbylamda*lamda;
v_o=fn*sqrt(9.81*L);
amp = (0.5)*H;

filename = strcat('max.xlsx');
alldata = xlsread(filename,strcat('fn',num2str(fn)),'A2:M41');
omegaval = alldata(1:end,1);
A = alldata(1:end,2:5);
B = alldata(1:end,6:9);
fid = 'addedmass.txt';
fopen(fid ,'wt+');
dlmwrite(fid,length(omegaval),'delimiter','\n','newline', 'pc','precision',3)
dlmwrite(fid,[omegaval,A],'delimiter','\t','newline', 'pc','precision',6,'-append')
fid = 'damp.txt';
fopen(fid ,'wt+');
dlmwrite(fid,length(omegaval),'delimiter','\n','newline', 'pc','precision',3)
dlmwrite(fid,[omegaval,B],'delimiter','\t','newline', 'pc','precision',6,'-append')

F3val = abs(alldata(1:end,10))*1000;
P3val = alldata(1:end,11);

F5val = abs(alldata(1:end,12))*1000;
P5val = alldata(1:end,13);

TP = sqrt((lamdabyL*L)/1.56);
omega = 2*pi/TP;
% omega0 = omega;
% K = omega0^2/9.81;
% U = fn*sqrt(L*9.81);
% omega = omega0*(1 + U*K)
% TP = 2*pi/omega;
% omegaval = omegaval.*(1 + U*(omegaval.^2/9.81));

% if halftimestep==1
%     dt_factor = 60;
% else
%     dt_factor = 30;
% end
dt = TP/25;
del_t = dt;
% while dt <= 0.01
%     dt_factor = dt_factor+10;
%     dt = TP/dt_factor;
% end
T = 0:dt:4*TP;
nmax = length(T);

F3a = interp1(omegaval,F3val,omega);
P3a = interp1(omegaval,P3val,omega);
F33 = amp*F3a * sin(omega*T + P3a);



F5a = interp1(omegaval,F5val,omega);
P5a = interp1(omegaval,P5val,omega);
F55 = amp*F5a * sin(omega*T + P5a);
% ramp
for i = 1:length(T)
    if T(i) < TP
%fmod = 0.5*(1 - cos(w*t/2))
    F33(i) = F33(i) * 0.5*(1 - cos(omega*T(i)/2));
    F55(i) = F55(i) * 0.5*(1 - cos(omega*T(i)/2));
    
    else
        
    F33(i) = F33(i);
    F55(i) = F55(i);
        
        
    end
    
end
fid = 'excitingF.txt';
fopen(fid ,'wt+');

dlmwrite(fid,length(T),'delimiter','\n','newline', 'pc','precision',6)
dlmwrite(fid,[T',F33',F55'],'delimiter','\t','newline', 'pc','precision',6,'-append')

disp('CREATED INPUTS FOR THE FOLLOWING PARAMETERS')

disp(strcat('fn  =',num2str(fn)))
disp(strcat('lamdabyL  =',num2str(lamdabyL)))
disp(strcat('Hbylamda  =',num2str(Hbylamda)))
disp(strcat('TP =',num2str(TP)))

%editjou(H,v_o,lamda,dt);%edit the jounal file to set up the CFD wave steepness and flow velocity
fclose('all');
end

