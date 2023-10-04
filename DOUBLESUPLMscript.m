clear all
close all
clc


sz=20;


simIn = Simulink.SimulationInput("DOUBLESUPLMsim");
set_param("DOUBLESUPLMsim",'FastRestart','on');
simIn = setModelParameter(simIn,"TimeOut",1e+6);

vR1=linspace(0.1,200,3);%o ideal aqui é separar em casos
vR12=linspace(0.1,200,sz);
vR23=linspace(0.1,200,sz);
vR3=linspace(0.1,200,3);%o ideal aqui é separar em casos


%calores especificos

%massas

t1=zeros(3,sz,sz,3);
t2=zeros(3,sz,sz,3);
t3=zeros(3,sz,sz,3);



p1=zeros(3,sz,sz,3);
p2=zeros(3,sz,sz,3);

for i=1:3
    R1 = vR1(i);
    for j=1:length(vR1)
        R12 = vR12(j);
        for k=1:length(vR1)
            R23 = vR23(k);
            for l=1:3
                R3 = vR3(l);
                simOut=sim(simIn);
                t1(i,j,k,l)=simOut.t1(1,1,end);
                t2(i,j,k,l)=simOut.t2(1,1,end);
                t3(i,j,k,l)=simOut.t3(1,1,end);
                p1(i,j,k,l)=simOut.p1(1,1,end);
                p2(i,j,k,l)=simOut.p2(1,1,end);
            end
        end
    end
end




