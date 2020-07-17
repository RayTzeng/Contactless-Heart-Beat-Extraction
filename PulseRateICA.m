function BVP_I=PulseRateICA(BVP,Px,Py) 
    M=zeros(length(BVP),3);

    M(:,1)=BVP;
    M(:,2)=Px;
    M(:,3)=Py;

    [W,Sx] = ica(M',3);   %JADE

    S=transpose(Sx);

    MinPx=zeros(1,3);
    for c=1:3
        cc = abs(corrcoef(S(:,c),M(:,2)))+abs(corrcoef(S(:,c),M(:,3)));
        MinPx(1,c)=cc(1,2);
    end

    BVP_I = S(:,argmin(MinPx,2));
end