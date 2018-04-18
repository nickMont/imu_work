function [state_out,P_out,RBI_out,Limu_out] = runEKF(imuStruct,gpsStruct, state0, RBI0, P0, systemParams)

state_out=[];
P_out=[];

persistent state Pk RBI tLast
if isempty(state)
    state=state0;
    RBI=RBI0;
    Pk=P0;
    tLast = 0;
end

if ~isempty(gpsStruct)
    dtt = gpsStruct{1}(1)-tLast;
    tLast = gpsStruct{1}(1);
   [state,Pk,RBI] = ekfUnknownLever(state,Pk,imuStruct,gpsStruct,systemParams,RBI);
end

state_out=state;
P_out=Pk;
RBI_out=RBI;
Limu_out=state(16:18);

end

