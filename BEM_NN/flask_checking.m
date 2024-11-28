% Set the URL of your Colab server
url = 'http://127.0.0.1:5000/predict';
particulars=[0.06 6.6067 3.1781 0.57 0.698];
% Define and send the input data
inputData = struct('input', [[[0.10447516789025962,-0.01331899873839977,0.053671428801788108,0.15888160862936418,0.2,1,1/35]]]);
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
    inFIS=genfis(x1,y1,optg);
    options=anfisOptions('InitialFIS',inFIS,'EpochNumber',30);
    FIS=anfis([x1,y1],options);
    anfis_output(i)=evalfis(FIS,particulars);
end

p = anfis_output.*(1025*9.81*0.1525*4/(1025*9.81*0.3345+1025*9.81*0.37+1025*9.81*0.32+1025*9.81*0.4167));
disp(p);
