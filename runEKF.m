function state_out = runEKF(imuStruct,gpsStruct, state0, RBI0, P0, systemParams)

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
   [state,Pk,RBI,Limu] = ekfUnknownLever(state,Pk,imuStruct,gpsStruct,Limu,systemParams,RBI);
end

state_out=state;

end

