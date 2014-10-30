function computeSmrDualDynamic_pairwisePotts_clearGlobal()
%computeSmrDualDynamic_pairwisePotts_clearGlobal clears global variables created by computeSmrDualDynamic_pairwisePotts
%
% Anton Osokin (firstname.lastname@gmail.com),  23.05.2013

global computeSmrDualDynamic_pairwisePotts_graphHandle
global computeSmrDualDynamic_pairwisePotts_lastPoint
computeSmrDualDynamic_pairwisePotts_lastPoint = [];

if ~isempty(computeSmrDualDynamic_pairwisePotts_graphHandle)
        for iLabel = 1 : length(computeSmrDualDynamic_pairwisePotts_graphHandle);
            if ~isempty(computeSmrDualDynamic_pairwisePotts_graphHandle{iLabel})
                deleteGraphCutDynamicMex( computeSmrDualDynamic_pairwisePotts_graphHandle{iLabel} )
            end
        end
end

clear global computeSmrDualDynamic_pairwisePotts_lastPoint
clear global computeSmrDualDynamic_pairwisePotts_graphHandle

end
