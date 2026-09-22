%% =====================================================================
%  validar_elastico_matlab.m — lado MATLAB de la validacion elastica
%
%  Homogeneiza con `localHomogenizeVoxel` (copia verbatim de AppFinal_V2.m via
%  applib) las mismas estructuras que evalua el port, y vuelca los tensores.
%
%  La homogeneizacion es DETERMINISTA: no interviene el generador aleatorio en
%  ningun punto. A diferencia del ajuste, aqui si cabe exigir coincidencia
%  numerica entre las dos plataformas.
% =====================================================================

clear; clc;

AQUI   = fileparts(mfilename('fullpath'));
PROY   = fileparts(AQUI);
GIBBON = 'C:\Users\carlo\OneDrive\Escritorio\Adds-On Matlab\GIBBON-master\lib';
APPLIB = fullfile(PROY, 'Validacion_Anexo', 'applib');
addpath(GIBBON); addpath(APPLIB);

global APPCTX APPFIG %#ok<GVMIS>
APPCTX = struct(); APPCTX.VOI = []; APPFIG = [];

SALIDA = fullfile(AQUI, 'resultados');
S = load(fullfile(SALIDA, 'casos_elastico.mat'));

nombres = S.nombres;
if ~iscell(nombres), nombres = cellstr(nombres); end
E_s     = double(S.E_s);
nu_s    = double(S.nu_s);
voxSize = double(S.voxSize);

fprintf('=========================================================\n');
fprintf(' HOMOGENEIZACION MATLAB — %d casos\n', numel(nombres));
fprintf(' E_s=%g  nu_s=%g  voxSize=%g\n', E_s, nu_s, voxSize);
fprintf('=========================================================\n\n');

res = struct('nombre',{},'dims',{},'rho',{},'C',{},'ok',{},'solver',{},'tiempo_s',{});

for k = 1:numel(nombres)
    nm  = nombres{k};
    BW  = logical(S.(nm));
    fprintf('--- %-14s %dx%dx%d  rho=%.4f ---\n', nm, size(BW), nnz(BW)/numel(BW));

    t0 = tic;
    [Ch, info] = localHomogenizeVoxel(BW, E_s, nu_s, voxSize, []);
    t = toc(t0);

    if ~info.ok
        fprintf('   FALLO: %s\n\n', info.msg);
        continue;
    end
    ec = localEngineeringConstants(Ch);
    fprintf('   %s\n', info.solver);
    fprintf('   Ex=%.6g  Ey=%.6g  Ez=%.6g   (%.1f s)\n\n', ...
            ec.Ex, ec.Ey, ec.Ez, t);

    res(end+1) = struct('nombre', nm, 'dims', size(BW), ...
        'rho', nnz(BW)/numel(BW), 'C', Ch, 'ok', info.ok, ...
        'solver', info.solver, 'tiempo_s', t); %#ok<AGROW>
end

meta = struct('descripcion', ...
    'localHomogenizeVoxel sobre las estructuras de casos_elastico.mat.', ...
    'E_s', E_s, 'nu_s', nu_s, 'voxSize', voxSize, ...
    'matlabVersion', version, ...
    'fecha', datestr(now,'yyyy-mm-dd HH:MM:SS'), 'resultados', res); %#ok<TNOW1,DATST>

fid = fopen(fullfile(SALIDA,'elastico_matlab.json'),'w');
fwrite(fid, jsonencode(meta,'PrettyPrint',true), 'char');
fclose(fid);

fprintf('=========================================================\n');
fprintf(' Escrito: %s\n', fullfile(SALIDA,'elastico_matlab.json'));
fprintf('=========================================================\n');
