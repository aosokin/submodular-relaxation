function computeSmrDualDynamic_highOrderPotts_clearGlobal()
%computeSmrDualDynamic_highOrderPotts_clearGlobal clears global variables created by computeSmrDualDynamic_highOrderPotts
%
% Anton Osokin (firstname.lastname@gmail.com),  27.10.2013

global computeSmrDualDynamic_highOrderPotts_graphHandle
global computeSmrDualDynamic_highOrderPotts_lastPoint
computeSmrDualDynamic_highOrderPotts_lastPoint = [];

if ~isempty(computeSmrDualDynamic_highOrderPotts_graphHandle)
        for iLabel = 1 : length(computeSmrDualDynamic_highOrderPotts_graphHandle);
            if ~isempty(computeSmrDualDynamic_highOrderPotts_graphHandle{iLabel})
                deleteGraphCutDynamicMex( computeSmrDualDynamic_highOrderPotts_graphHandle{iLabel} )
            end
        end
end

clear global computeSmrDualDynamic_highOrderPotts_lastPoint
clear global computeSmrDualDynamic_highOrderPotts_graphHandle

end
