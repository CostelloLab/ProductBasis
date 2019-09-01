numCells = 1;

onecell = cell(1, 2);

onecell{1} = { [], [ 0 ], [] };
onecell{2} = { [ 13, 1, 15, 2 ], [ 0 0 0 0 0 1 0 0 0 1 0 0 1 1 0 0 ], [] };     % 13?
onecell{3} = { [ 2 ], [ 0 1 ], [] };
onecell{4} = { [ 3-15, 3+15, 1 ], [ 0 0 1 0 1 0 1 0 ], [] };
onecell{5} = { [ 4 ], [ 0 1 ], [] };
onecell{6} = { [ 5, 15 ], [ 0 0 1 0 ], [] };
onecell{7} = { [ 6 ], [ 0 1 ], [] };
onecell{8} = { [ 14, 5, 15 ], [ 0 0 0 0 1 0 0 0 ], [] };
onecell{9} = { [ 8, 9, 7-15, 7+15 ], [ 0 0 0 0 1 0 0 0 1 1 1 1 1 1 1 1 ], [] };
onecell{10} = { [ 9, 7-15, 7+15 ], [ 0 0 0 0 0 1 1 1 ], [] };
onecell{11} = { [ 9, 7-15, 7+15 ], [ 1 1 1 1 0 1 1 1 ], [] };
onecell{12} = { [ 5 ], [ 1 0 ], [] };
onecell{13} = { [ 12 ], [ 0 1 ], [] };
onecell{14} = { [ 13, 11, 6-15, 6+15 ], [ 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 ], [] };
onecell{15} = { [ 13, 11, 6-15, 6+15 ], [ 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ], [] };

logicTables = cell(1, 15*numCells);
for loopCell = 1:numCells
    moduleRange = 15*loopCell + (-14:0);
    logicTables(moduleRange) = onecell;
    if (mod(loopCell, 4) == 0) || (mod(loopCell, 4) == 3)
        logicTables{moduleRange(1)}{2}(1) = 1;
    end
    for loopModule = 1:15
        logicTables{moduleRange(loopModule)}{1} = logicTables{moduleRange(loopModule)}{1} + 15*(loopCell-1);
        logicTables{moduleRange(loopModule)}{1} = mod(logicTables{moduleRange(loopModule)}{1}-1, 15*numCells)+1;
    end
end

M = polynomialModelT(logicTables);          % using time-dependent test