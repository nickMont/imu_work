function state_out = runUKF(imuStruct,gpsStruct, state0, RBI0, P0)


%Qimu=diag([10^-3*ones(3,1); 10^-7*ones(3,1); 10^-3*ones(3,1); 10^-8*ones(3,1)]);

state_out=[];

persistent state Pk RBI tLast Limu
if isempty(state)
    state=state0;
    RBI=RBI0;
    Pk=P0;
    tLast = 0;
    Limu=zeros(3,1);
end

if ~isempty(gpsStruct)
    dtt = gpsStruct{1}(1)-tLast;
    tLast = gpsStruct{1}(1);
   [state,Pk,RBI,Limu] = ukfUnknownLever(x0,Pk,imuStruct,gpsStruct,state,Limu,systemParams,RBI);
end

state_out=state;

end

