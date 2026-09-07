function lVal = SampleLofF(forceLengthData,inpForce,delF)
    %Find the corresponding bin to given unifVal
%     if inpForce~=0.1 && inpForce~=5
    lVal = interp1(forceLengthData(:,1),forceLengthData(:,2),inpForce);
%     residue=(inpForce-forceLengthData(1,1))./delF
%     lbInd = floor(residue);
%     ubInd = lbInd+1;
%     xlb = forceLengthData(lbInd,1)
%     ylb = forceLengthData(lbInd,2);
%     xub = forceLengthData(ubInd,1);
%     yub = forceLengthData(ubInd,2);
%     lVal = (ylb+(inpForce-xlb).*(yub-ylb)./(xub-xlb));
%     elseif inpForce==0.1
%         lVal = forceLengthData(1,2);
%     elseif inpForce==5
%         lVal = forceLengthData(end,2);   
%     end
end