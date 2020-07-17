function PR=PulseRateCal(BVP_I,PRres,WSZPR,FS)
    TPR=floor(length(BVP_I)/FS);
    NPR=floor((length(BVP_I)-WSZPR)/(FS*PRres));
    WSZPRres=(FS*PRres);
    PR = zeros(NPR);
    threshold=10;
 
    movAvg=0;
    mu=0.2;

    for i = 1 : NPR
        SPR=(i-1)*WSZPRres+1;
        EPR=SPR+WSZPR-1;
    
        %fprintf('\n%d %d %d',SPR,EPR,NPR);
    
        PR(i,1) = prpsd(BVP_I((SPR:EPR)),FS,40,240,false); 
    
        if i==1
            movAvg=PR(i,1);
            PR(i,2)=movAvg;
        else
            differ=PR(i,1)-PR(i-1,2);
            if abs(differ)>threshold
                beta=exp(-abs(differ-threshold));
            PR(i,2)=(1-beta)*movAvg+beta*PR(i,1);
            else
                PR(i,2)=PR(i,1);
            end
            movAvg=movAvg*(1-mu)+PR(i,2)*mu;
        end 
    end
end