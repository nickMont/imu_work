function v = vectorSaturationF(v,vlim,vlim2)
if nargin==2
    j=length(v);
    for i=1:j
        if v(i)<-vlim(i)
            v(i)=-vlim(i);
        elseif v(i)>vlim(i)
            v(i)=vlim(i);
        end
    end
else
    j=length(v);
    for i=1:j
        if v(i)<vlim(i)
            v(i)=vlim(i);
        elseif v(i)>=vlim2(i)
            v(i)=vlim2(i);
        end
    end
end
end

