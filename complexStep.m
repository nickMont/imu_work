function df = complexStep(funHandle,x0,h)
for ii=1:length(x0)
    xd = x0;
    xd(ii) = x0(ii) + 1j*h;
    df(:,ii) = imag(funHandle(xd))./h;
end
end

