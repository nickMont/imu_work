%clear;clc;
loadline=44;
initstate=[];

rcg2p=[0;0;0];

load dat999.mat
data=squeeze(trainDat(loadline,:,1,1));
ewxv0=data(1:12);
x0=ewxv0(7:9);
v0=ewxv0(10:12);
RBI0=euler2dcm(ewxv0(1:3));
gpsMeasStore={};
imuMeasStore={};

for ijk=1:400
    gpsMeas={data(12+(ijk-1)*4+1:12+4*ijk)'};
    gpsMeasStore{ijk}=gpsMeas;
end
imuData=data(12+4*400+1:end)';
times=imuData(1:7:end);

tlast=0;
for ijk=1:400
    gpsTime=gpsMeasStore{ijk}{1}(1);
    imuInd=find((times(:)>(tlast-1e-9)) & (times(:)<(gpsTime+1e-9)));
    for j=1:length(imuInd)
        tt=imuInd(j);
        imuMeas{j}=imuData((tt-1)*7+1:tt*7);
    end
    imuMeasStore{ijk}=imuMeas;
    
end
