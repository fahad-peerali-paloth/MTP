function [ p ] = nnmodel(hydrostaticP,particulars,x3,v3,x5,v5,fn,lamdabyL,Hbylamda)
url = 'http://127.0.0.1:5000/predict';
inputData = struct('input', [x3,v3,x5,v5,fn,lamdabyL,Hbylamda]);
options = weboptions('Timeout', 20, 'MediaType', 'application/json');
response = webwrite(url, inputData, options);
yk = response.predictions;
%y1=prediction1(1,:);
%y2=prediction1(2,:);
%y3=prediction1(3,:);
%y4=prediction1(4,:);
%x=[1/Bs L/B,B/T,CB,CWL)
hp(1,1)=1025*9.81*0.3345;
hp(2,1)=1025*9.81*0.37;
hp(3,1)=1025*9.81*0.32;
hp(4,1)=1025*9.81*0.4167;
anfis_output=zeros(1,13);
for i=1:13
    y1=yk(:,i).*hp;
    x1=[0.03106 7.14 2.99 0.649 0.845;0.0394 6.88 2.7 0.558 0.692;0.23 7.31 3.129 0.458 0.711;3.333 10 2.4 0.388 0.665];
    optg=genfisOptions('GridPartition');
    optg.NumMembershipFunctions = 2;
    input_mf_types = ["gaussmf","trimf", "trimf", "gaussmf", "gaussmf"];
    optg.InputMembershipFunctionType = input_mf_types;
    inFIS=genfis(x1,y1,optg); \()
    options=anfisOptions('InitialFIS',inFIS,'EpochNumber',30);
    FIS=anfis([x1,y1],options);
    anfis_output(i)=evalfis(FIS,particulars);
end

p = anfis_output;
end

