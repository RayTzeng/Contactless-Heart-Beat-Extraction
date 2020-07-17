function [MAE,SR,CC]=Quantitate(GT,PR)
    AE=abs(GT-PR);
    MAE=mean(AE);
    success=zeros(length(AE),1);
    for i=1:length(AE)
        success(i)=AE(i)<=5;
    end
    SR=mean(success);
    cc=corrcoef(PR,GT);
    CC=cc(1,2);
end