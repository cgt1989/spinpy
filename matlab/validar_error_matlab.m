%% =====================================================================
%  validar_error_matlab.m — lado MATLAB de la validacion de la funcion de error
%
%  Lee los pares generados por `python validar_error.py --exportar`, evalua
%  computeMorphometricError (la copia verbatim de applib) sobre cada uno y
%  escribe los resultados para que el lado Python los contraste.
%
%  Los campos con valor `null` en el JSON se OMITEN del struct, que es
%  exactamente lo que ocurre cuando una metrica no se puede calcular. Ese es
%  el caso que ejercita la invariante C6.
% =====================================================================

clear; clc;

AQUI   = fileparts(mfilename('fullpath'));
PROY   = fileparts(AQUI);
APPLIB = fullfile(PROY, 'Validacion_Anexo', 'applib');
addpath(APPLIB);

global APPCTX APPFIG %#ok<GVMIS>
APPCTX = struct(); APPCTX.VOI = []; APPFIG = [];

SALIDA = fullfile(AQUI,'resultados');
d = jsondecode(fileread(fullfile(SALIDA,'pares_error.json')));
campos = cellstr(d.campos);
pares  = d.pares;

n = numel(pares);
fprintf('Evaluando computeMorphometricError sobre %d pares...\n', n);

errores  = nan(n,1);
n_usados = nan(n,1);

for k = 1:n
    P = pares(k);
    a = P.a; b = P.b;
    if iscell(a), a = cell2mat(cellfun(@(x) toNum(x), a, 'UniformOutput', false)); end
    if iscell(b), b = cell2mat(cellfun(@(x) toNum(x), b, 'UniformOutput', false)); end

    mA = struct(); mB = struct();
    for q = 1:numel(campos)
        if q <= numel(a) && isfinite(a(q)), mA.(campos{q}) = a(q); end
        if q <= numel(b) && isfinite(b(q)), mB.(campos{q}) = b(q); end
    end

    [e, nu] = computeMorphometricError(mA, mB);
    errores(k)  = e;
    n_usados(k) = nu;

    if mod(k, 200) == 0, fprintf('  %d/%d\n', k, n); end
end

out = struct('errores', errores, 'n_usados', n_usados, ...
             'matlabVersion', version, ...
             'fecha', datestr(now,'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
fid = fopen(fullfile(SALIDA,'error_matlab.json'),'w');
fwrite(fid, jsonencode(out,'PrettyPrint',true), 'char');
fclose(fid);

fprintf('Escrito: %s\n', fullfile(SALIDA,'error_matlab.json'));


function v = toNum(x)
    if isempty(x) || (~isnumeric(x) && ~islogical(x))
        v = NaN;
    else
        v = double(x);
    end
end
