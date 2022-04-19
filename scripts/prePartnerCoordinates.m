function PreCoordinates = PrePartnerCoordinates(PrePSDSynapseSegID,df)
% *PrePartnerCoordinates* is used to obtain the x,y,z coordinates in NG
% framework for the listed PrePartnerSynapseID

PreX = [];
PreY = [];
PreZ = [];

    for i = 1: length(PrePSDSynapseSegID)
         PreX = [PreX; df.presyn_x(df.psd_segid==PrePSDSynapseSegID(i))];
         PreY = [PreY; df.presyn_y(df.psd_segid==PrePSDSynapseSegID(i))];
         PreZ = [PreZ; df.presyn_z(df.psd_segid==PrePSDSynapseSegID(i))];
    end
    PreCoordinates = [PreX,PreY,PreZ];
end
