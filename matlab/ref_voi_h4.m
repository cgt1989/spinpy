%% =====================================================================
%  ref_voi_h4.m — REFERENCIA MATLAB SOBRE LOS VOIs REALES DE H4
%
%  A diferencia de gen_referencia.m, aqui NO HAY ESTOCASTICIDAD: el VOI es un
%  dato fijo leido de disco. Python y MATLAB deben coincidir hasta la
%  precision de la aritmetica, salvo por el variante de marching cubes.
%
%  Es por tanto un contraste MUCHO MAS EXIGENTE que el del generador: no hay
%  varianza donde esconder una discrepancia. Cualquier diferencia en BV/TV,
%  Tb.Sp, DA o la direccion principal es un error real del port.
%
%  Se usa readVTKVOI + localMorphometryFromVoxels directamente en lugar de
%  computeVOIMetrics para no arrastrar el cacheo mecanico (J1), que depende
%  del struct `app` de la interfaz y no aporta nada a esta comparacion.
%
%  Salida: Port_Python/resultados/referencia_voi_h4.json
% =====================================================================

clear; clc;

AQUI   = fileparts(mfilename('fullpath'));
PROY   = fileparts(AQUI);
GIBBON = 'C:\Users\carlo\OneDrive\Escritorio\Adds-On Matlab\GIBBON-master\lib';
APPLIB = fullfile(PROY, 'Validacion_Anexo', 'applib');
VOIDIR = fullfile(PROY, 'H4', 'Segmentadas');

addpath(GIBBON); addpath(APPLIB);

SALIDA = fullfile(AQUI,'resultados');
if ~exist(SALIDA,'dir'), mkdir(SALIDA); end

global APPCTX APPFIG %#ok<GVMIS>
APPCTX = struct(); APPCTX.VOI = []; APPFIG = [];

archivos = dir(fullfile(VOIDIR,'*.vtk'));
campos = {'BVTV','DA','DA2','BSBV','BSTV','BS','BV','TV','TbTh','TbSp','TbN', ...
          'PoTot','FracPort'};

fprintf('=========================================================\n');
fprintf(' REFERENCIA MATLAB — VOIs REALES DE H4  (%d archivos)\n', numel(archivos));
fprintf('=========================================================\n\n');

res = struct('archivo',{},'dims',{},'spacing',{},'metricas',{}, ...
             'dir_principal',{},'eigenvalues',{},'MIL_valid',{}, ...
             'MIL_clamped',{},'tiempo_s',{});

for k = 1:numel(archivos)
    f = fullfile(VOIDIR, archivos(k).name);
    fprintf('--- %s ---\n', archivos(k).name);

    % Un VOI truncado no es motivo para abortar el lote: se avisa y se sigue.
    % (VOI_0060_y0065_x0449_z0193.vtk esta al 59%, exportacion interrumpida.)
    t0 = tic;
    try
        [VOI, spacing] = readVTKVOI(f);
    catch ME
        fprintf('   OMITIDO: %s\n\n', ME.message);
        continue;
    end
    if ~islogical(VOI), VOI = VOI > 0; end
    m = localMorphometryFromVoxels(VOI, spacing, true);
    t = toc(t0);

    vals = nan(1, numel(campos));
    for q = 1:numel(campos)
        if isfield(m, campos{q}) && isscalar(m.(campos{q}))
            vals(q) = double(m.(campos{q}));
        end
    end

    fprintf('   dims=%dx%dx%d  spacing=%.6f mm\n', size(VOI,1), size(VOI,2), size(VOI,3), spacing(1));
    fprintf('   BV/TV=%.5f  BS/BV=%.4f  Tb.Th=%.6f  Tb.Sp=%.6f  Tb.N=%.4f\n', ...
            vals(1), vals(4), vals(9), vals(10), vals(11));
    fprintf('   DA=%.5f  dir=[%.4f %.4f %.4f]  clamped=%d  (%.1f s)\n\n', ...
            vals(2), m.dir_principal(1), m.dir_principal(2), m.dir_principal(3), ...
            m.MIL_clamped, t);

    res(end+1) = struct( ...
        'archivo', archivos(k).name, ...
        'dims', size(VOI), 'spacing', spacing(:)', ...
        'metricas', vals, ...
        'dir_principal', m.dir_principal(:)', ...
        'eigenvalues', m.eigenvalues(:)', ...
        'MIL_valid', m.MIL_valid, 'MIL_clamped', m.MIL_clamped, ...
        'tiempo_s', t); %#ok<AGROW>
end

meta = struct('descripcion', ...
    'Morfometria MATLAB sobre los VOIs reales de H4 (deterministica).', ...
    'campos', {{campos}}, 'matlabVersion', version, ...
    'fecha', datestr(now,'yyyy-mm-dd HH:MM:SS'), 'resultados', res); %#ok<TNOW1,DATST>

fid = fopen(fullfile(SALIDA,'referencia_voi_h4.json'),'w');
fwrite(fid, jsonencode(meta,'PrettyPrint',true), 'char');
fclose(fid);

fprintf('=========================================================\n');
fprintf(' Escrito: %s\n', fullfile(SALIDA,'referencia_voi_h4.json'));
fprintf('=========================================================\n');
