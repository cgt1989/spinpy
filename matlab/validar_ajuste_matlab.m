%% =====================================================================
%  validar_ajuste_matlab.m — ejecuta el AJUSTE REAL de la app, sin interfaz
%
%  Corre `fitSpinodoidToVOI('slow')` sobre VOIs reales usando el arnes de
%  Validacion_Anexo (optlib sombrea applib con stubs de interfaz) y vuelca el
%  estado final del optimizador para contrastarlo con el port a Python.
%
%  QUE SE PUEDE CONTRASTAR Y QUE NO
%  --------------------------------
%  El generador es estocastico y MATLAB y NumPy tienen flujos de numeros
%  aleatorios distintos, asi que cada plataforma evalua REALIZACIONES
%  DISTINTAS en los mismos puntos de la rejilla. Eso parte la validacion en
%  dos mitades de naturaleza opuesta:
%
%    DETERMINISTA (debe coincidir exactamente):
%      - las rejillas de candidatos de las etapas A, B y C
%      - los presets de angulos conicos y su orden
%      - la matriz de rotacion R_fit derivada del eje principal del VOI
%      - el numero de evaluaciones previstas
%    Ninguna de esas cosas depende del azar: son formulas sobre BV/TV y sobre
%    la direccion principal del VOI. Si difieren, el port esta mal.
%
%    ESTOCASTICO (solo comparable en distribucion):
%      - el error alcanzado y los parametros ganadores
%
%  Se fijan explicitamente los valores de partida de los "deslizadores"
%  (thetas y numero de ondas) para que las dos plataformas arranquen del mismo
%  sitio; si no, las rejillas ya no serian comparables.
%
%  Salida: Port_Python/resultados/ajuste_matlab.json
% =====================================================================

clear; clc;

AQUI  = fileparts(mfilename('fullpath'));
PROY  = fileparts(AQUI);
CODE  = fullfile(PROY, 'Validacion_Anexo', 'code');

addpath(CODE);
ctx = opt_setup('Quiet', false);   %#ok<NASGU>

SALIDA = fullfile(AQUI, 'resultados');
if ~exist(SALIDA,'dir'), mkdir(SALIDA); end

% --- Valores de partida, identicos en los dos lados ----------------------
THETAS0   = [15 15 15];
NUM_WAVES = 700;
MODO      = 'slow';

VOIDIR = fullfile(PROY, 'H4', 'Segmentadas');
casos  = {'VOI_proximal_cubico.vtk', 'VOI_medio_cubico.vtk', 'VOI_distal_cubico.vtk'};

global APPCTX VALRES %#ok<GVMIS>

res = struct('archivo',{},'dims',{},'spacing',{},'metricasVOI',{}, ...
             'rejillaDensidad',{},'rejillaOnda',{},'presetsTheta',{}, ...
             'R',{},'densidad',{},'waveNum',{},'numWaves',{},'thetas',{}, ...
             'errorFinal',{},'metricasGanador',{},'nEvaluaciones',{}, ...
             'nEvalPrevistas',{},'tiempo_s',{});

fprintf('=========================================================\n');
fprintf(' AJUSTE REAL DE LA APP (fitSpinodoidToVOI, modo %s)\n', MODO);
fprintf(' thetas de partida = [%d %d %d], numWaves = %d\n', THETAS0, NUM_WAVES);
fprintf('=========================================================\n\n');

for k = 1:numel(casos)
    f = fullfile(VOIDIR, casos{k});
    fprintf('--- %s ---\n', casos{k});

    [VOI, spacing] = readVTKVOI(f);
    if ~islogical(VOI), VOI = VOI > 0; end

    % Metricas del VOI por el mismo motor que usa la app
    mVOI = localMorphometryFromVoxels(VOI, spacing, true);

    opt_contexto(VOI, spacing, ...
                 'Thetas',   THETAS0, ...
                 'NumWaves', NUM_WAVES, ...
                 'WaveNum',  15, ...
                 'Metrics',  mVOI);

    VALRES = struct();
    t0 = tic;
    fitSpinodoidToVOI(MODO);
    t = toc(t0);

    if ~isfield(VALRES,'errorFinal')
        fprintf('   SIN RESULTADO (VALRES incompleto)\n\n');
        continue;
    end

    fprintf('   BV/TV VOI=%.4f  DA=%.4f\n', mVOI.BVTV, mVOI.DA);
    fprintf('   ganador: dens=%.4f  wave=%g  nw=%d  thetas=[%d %d %d]\n', ...
            VALRES.densidad, VALRES.waveNum, VALRES.numWaves, VALRES.thetas);
    fprintf('   error=%.6f   %d evaluaciones en %.1f s\n\n', ...
            VALRES.errorFinal, VALRES.nEvaluaciones, t);

    res(end+1) = struct( ...
        'archivo', casos{k}, 'dims', size(VOI), 'spacing', spacing(:)', ...
        'metricasVOI', VALRES.metricasVOI, ...
        'rejillaDensidad', VALRES.rejillaDensidad(:)', ...
        'rejillaOnda',     VALRES.rejillaOnda(:)', ...
        'presetsTheta',    {VALRES.presetsTheta}, ...
        'R', VALRES.R, ...
        'densidad', VALRES.densidad, 'waveNum', VALRES.waveNum, ...
        'numWaves', VALRES.numWaves, 'thetas', VALRES.thetas(:)', ...
        'errorFinal', VALRES.errorFinal, ...
        'metricasGanador', VALRES.metricasGanador, ...
        'nEvaluaciones', VALRES.nEvaluaciones, ...
        'nEvalPrevistas', VALRES.nEvalPrevistas, ...
        'tiempo_s', t); %#ok<AGROW>
end

meta = struct('descripcion', ...
    'fitSpinodoidToVOI ejecutado sin interfaz mediante el arnes optlib/applib.', ...
    'modo', MODO, 'thetas0', THETAS0, 'numWaves0', NUM_WAVES, ...
    'matlabVersion', version, ...
    'fecha', datestr(now,'yyyy-mm-dd HH:MM:SS'), 'resultados', res); %#ok<TNOW1,DATST>

fid = fopen(fullfile(SALIDA,'ajuste_matlab.json'),'w');
fwrite(fid, jsonencode(meta,'PrettyPrint',true), 'char');
fclose(fid);

fprintf('=========================================================\n');
fprintf(' Escrito: %s\n', fullfile(SALIDA,'ajuste_matlab.json'));
fprintf('=========================================================\n');
