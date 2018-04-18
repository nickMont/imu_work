function [state_out,P_out,RBI_out] = runUKF15(imuStruct,gpsStruct, state0, RBI0, P0, systemParams,Limu)

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
   [state,Pk,RBI] = ukfFixedLever(state,Pk,imuStruct,gpsStruct,Limu,systemParams,RBI);
end

state_out=state;
P_out=Pk;
RBI_out=RBI;

end

