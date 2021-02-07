function vHex = SpringRelax_func( vHex, eHex, hHex )


%%
% spring model

nV = size(vHex,1);
nE = size(eHex,1);

%%
% edge lengths
eHexLen = sum((vHex(eHex(:,1),:)-vHex(eHex(:,2),:)).^2,2).^0.5;

% figure(2)
% subplot(1,3,1)
% hist(eHexLen)
% title( sprintf('RMS dLen=%.3f',mean((eHexLen-hHex).^2).^0.5) )


% Forces
Cx = zeros(nV,nV);
for ii=1:nE
    Cx(eHex(ii,1),eHex(ii,2)) = (eHexLen(ii)-hHex(ii)) * (vHex(eHex(ii,2),1)-vHex(eHex(ii,1),1)) / hHex(ii);
    Cx(eHex(ii,2),eHex(ii,1)) = (eHexLen(ii)-hHex(ii)) * (vHex(eHex(ii,1),1)-vHex(eHex(ii,2),1)) / hHex(ii);
end
Fx = sum(Cx,2);

Cy = zeros(nV,nV);
for ii=1:nE
    Cy(eHex(ii,1),eHex(ii,2)) = (eHexLen(ii)-hHex(ii)) * (vHex(eHex(ii,2),2)-vHex(eHex(ii,1),2)) / hHex(ii);
    Cy(eHex(ii,2),eHex(ii,1)) = (eHexLen(ii)-hHex(ii)) * (vHex(eHex(ii,1),2)-vHex(eHex(ii,2),2)) / hHex(ii);
end
Fy = sum(Cy,2);

Cz = zeros(nV,nV);
for ii=1:nE
    Cz(eHex(ii,1),eHex(ii,2)) = (eHexLen(ii)-hHex(ii)) * (vHex(eHex(ii,2),3)-vHex(eHex(ii,1),3)) / hHex(ii);
    Cz(eHex(ii,2),eHex(ii,1)) = (eHexLen(ii)-hHex(ii)) * (vHex(eHex(ii,1),3)-vHex(eHex(ii,2),3)) / hHex(ii);
end
Fz = sum(Cz,2);

F = (Fx.^2+Fy.^2+Fz.^2).^0.5;



% update positions based on forces
Fmax = max(abs(F));
if Fmax>0
    scl = 0.1;
    vHex(:,1) = vHex(:,1) + scl * Fx/Fmax;
    vHex(:,2) = vHex(:,2) + scl * Fy/Fmax;
    vHex(:,3) = vHex(:,3) + scl * Fz/Fmax;
end

