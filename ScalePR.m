function PR=ScalePR(refPR,scale,offset)
    PR=zeros(length(refPR)/scale,2);
    for i=1:length(refPR)/scale
        PR(i,1)=refPR(scale*i-scale+1,1);
        PR(i,2)=refPR(scale*i-scale+1,2);
    end
end