
function bemio(filename)
    current_directory = pwd;
    cd('C:\Users\shusain\OneDrive - NREL\Documents\GitHub\WEC-Sim')
    run('addWecSimSource.m')
    cd(current_directory)

    hydro = struct();

    hydro = readCAPYTAINE(hydro, strcat(filename, '.nc')); % <-- uses the passed filename
    hydro = radiationIRF(hydro,30,3000,[],[],5);
    hydro = excitationIRF(hydro,30,3000,[],[],5);
    writeBEMIOH5(hydro)
    plotBEMIO(hydro)
end
