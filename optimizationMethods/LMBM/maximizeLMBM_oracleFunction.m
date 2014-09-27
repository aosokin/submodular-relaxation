function [fValNegative, subgradientNegative] = maximizeLMBM_oracleFunction(x)

    global globalOraclePointer
    global globalPrimalFuncPointer
    global globalIterPrintFlag
    global globalDualVarsScale


    [fVal, subgradient, primalEstimate] = globalOraclePointer(x / globalDualVarsScale);
    fValNegative = -fVal;
    subgradientNegative = -subgradient * globalDualVarsScale;
    
    curStart = tic;
    
    global globalBestPoint
    global globalBestVal
    global globalBestPrimal
    global globalBestLabeling

    global globalTimePlot
    global globalFuncPlot
    global globalPrimalPlot

    global globalTimeStart
    global globalTimeGarbage
    global globalNumOracleCalls

    globalNumOracleCalls = globalNumOracleCalls + 1;
    globalFuncPlot( globalNumOracleCalls ) = fVal;
    if fVal > globalBestVal
        globalBestVal = fVal;
        globalBestPoint = x;
    end
    
    [curPrimal, newLabeling] = globalPrimalFuncPointer( primalEstimate );
    if curPrimal < globalBestPrimal
        globalBestPrimal = curPrimal;
        globalBestLabeling = newLabeling;
    end
    globalPrimalPlot( globalNumOracleCalls ) = curPrimal;
    
    globalTimePlot( globalNumOracleCalls ) = nan;
    
    curGarbage = toc( curStart );
    globalTimeGarbage = globalTimeGarbage + curGarbage;
    
    globalTimePlot( globalNumOracleCalls ) = toc( globalTimeStart ) - globalTimeGarbage;
    
    if globalIterPrintFlag 
        fprintf('Oracle call: %d, dual: %f, primal: %f\n', globalNumOracleCalls, fVal, curPrimal);
    end
end



