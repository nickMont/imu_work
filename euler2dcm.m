function RBI = euler2dcm(e)
%vB  = RBI * vI
  
cPhi = cos(e(1)); sPhi = sin(e(1));
cThe = cos(e(2)); sThe = sin(e(2));
cPsy = cos(e(3)); sPsy = sin(e(3));
%Psi -> Psy to avoid P_i confusion with Phi

RBI = [ cThe*cPsy, cThe*sPsy, -sThe;
        -cPhi*sPsy + sPhi*sThe*cPsy,   cPhi*cPsy + sPhi*sThe*sPsy,  sPhi*cThe;
         sPhi*sPsy + cPhi*sThe*cPsy,  -sPhi*cPsy + cPhi*sThe*sPsy,  cPhi*cThe];







