function [state_out,P_out,RBI_out,Sk_out,nuj_out] = runUKF15(imuStruct,gpsStruct, state, RBI, Pk, systemParams,Limu)

state_out=[];
P_out=[];

if ~isempty(gpsStruct)
   [state,Pk,RBI,Sk_out,nuj_out] = ukfFixedLever(state,Pk,imuStruct,gpsStruct,Limu,systemParams,RBI);
end

state_out=state;
P_out=Pk;
RBI_out=RBI;

end

