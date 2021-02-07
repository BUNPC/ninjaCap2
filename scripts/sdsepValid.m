function [sdseps] = sdsepValid(grommets)
%SDSEPVALID checks the source detectors separation of all relevant optodes
%on the atlas and on the panels to enable a validation of the projection

%identify emitters
elist = arrayfun(@(x) x.optType==1, grommets);
%identify detectors
dlist = arrayfun(@(x) x.optType==2, grommets);
% init results struct
sdseps = struct('EmitterID', {}, 'DetectorID', {}, 'DistAtlas', {}, 'DistPanel', {}, 'DistError', {});
% for all detectors
for dd = find(dlist)
       % and corresponding sources 
       for ss = grommets(dd).MList'
           elist = find(elist);
           % are on the same panel?
           
           if ~isempty(grommets(dd).flags)
               str = grommets(dd).flags{1};
           else
               str = 'empty';
           end
           
           if strcmp(grommets(elist(ss)).panel, grommets(dd).panel) && ~strcmp(str, 'short-separation')
               % calculate euclidean distance on Atlas
               vA=grommets(elist(ss)).posHead-grommets(dd).posHead;
               distA = norm(vA);
               % calculate euclidean distance on panel
               vP=grommets(elist(ss)).posPanel-grommets(dd).posPanel;
               distP = norm(vP);
               % save
               sdsep = struct('EmitterID', ss, 'DetectorID', dd-numel(elist), 'DistAtlas', distA, 'DistPanel', distP, 'DistError', distA-distP);
               sdseps = [sdseps sdsep];
           end
       end
end

errList = arrayfun(@(x) x.DistError, sdseps);
figure
plot(errList, 'xr')
xlabel('SD pair')
ylabel('SDSep deviation / mm')
title('SD separations: deviation 2D panel from Atlas registered probe')
hold on
plot(1:numel(errList), mean(errList)*ones(numel(errList),1), '--k')
legend('single probe','avg probes')
dim = [.2 .3 .1 .1];
str = ['min: ' num2str(min(errList),'%2.2f') 'mm | max: ' num2str(max(errList), '%2.2f')...
    'mm | avg/std: ' num2str(mean(errList),'%2.2f') '/' num2str(std(errList), '%2.2f' ) 'mm'];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on')
end

