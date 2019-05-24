function ampSur = makeAmpRPhiSurr(R,ampSegCol,relativePhiCol,phiBin)
for n = 1:500
    for i = 1:4
        relativePhiColSur = randbetween(-pi,pi,1,length(relativePhiCol));
        [shiftPhiCol phipeak]= findAmpPhi(R,ampSegCol(i,:),relativePhiColSur,phiBin);
        [ampBinMuSur(n,i,:) ampBinSEMSur ampBinDatSur] = binstats(shiftPhiCol,ampSegCol(i,:),phiBin);
    end
end
for i=1:4
    X = squeeze(ampBinMuSur (:,i,:));
    ampSur(:,i,1) = prctile(X,95);
    ampSur(:,i,2) = prctile(X,5);
end