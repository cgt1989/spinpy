%% =====================================================================
%  ref_vois_todos.m — REFERENCIA MATLAB SOBRE TODOS LOS VOIs DISPONIBLES
%
%  Extiende ref_voi_h4.m a todo el material disponible:
%
%    - 24 VOIs equinos  (H1..H4 x {proximal, medio, distal} x {cubico, PCA})
%                       97^3 a 51.489 um/voxel, sesamoideos, SkyScan
%    - 10 VOIs porcinos (C1..C5, V1..V5) del conjunto publico de Koria,
%                       Mengoni & Brockett (2020), 188^3 a 16 um/voxel
%
%  Los porcinos son el control externo: otra especie, otro escaner, otro
%  tamano de voxel y otro tamano de VOI. Si el port coincidiera solo en los
%  equinos podria ser casualidad de un unico protocolo de adquisicion.
%
%  Como el VOI es un dato fijo, la morfometria es DETERMINISTA: MATLAB y
%  Python deben coincidir hasta la aritmetica. No hay varianza estocastica
%  donde esconder una discrepancia.
%
%  Salida: Port_Python/resultados/referencia_vois_todos.json
% =====================================================================

clear; clc;

AQUI   = fileparts(mfilename('fullpath'));
PROY   = fileparts(AQUI);
GIBBON = 'C:\Users\carlo\OneDrive\Escritorio\Adds-On Matlab\GIBBON-master\lib';
APPLIB = fullfile(PROY, 'Validacion_Anexo', 'applib');
BANCO  = 'C:\Users\carlo\OneDrive\Escritorio\VOIs Caballos';

addpath(GIBBON); addpath(APPLIB);

SALIDA = fullfile(AQUI,'resultados');
if ~exist(SALIDA,'dir'), mkdir(SALIDA); end

global APPCTX APPFIG %#ok<GVMIS>
APPCTX = struct(); APPCTX.VOI = []; APPFIG = [];

% --- Inventario -----------------------------------------------------------
lista = struct('ruta',{},'etiqueta',{},'grupo',{},'espec',{},'sitio',{},'tipo',{});

sitios = {'proximal','medio','distal'};
tipos  = {'cubico','PCAaligned'};
for h = 1:4
    for s = 1:numel(sitios)
        for t = 1:numel(tipos)
            f = fullfile(BANCO, sprintf('H%d',h), ...
                         sprintf('VOI_%s_%s.vtk', sitios{s}, tipos{t}));
            if exist(f,'file')
                lista(end+1) = struct('ruta',f, ...
                    'etiqueta',sprintf('H%d %s %s',h,sitios{s},tipos{t}), ...
                    'grupo','equino','espec',sprintf('H%d',h), ...
                    'sitio',sitios{s},'tipo',tipos{t}); %#ok<AGROW>
            end
        end
    end
end

dPorc = dir(fullfile(BANCO,'Publico_Porcino','VOI_*.mat'));
for k = 1:numel(dPorc)
    [~,nm] = fileparts(dPorc(k).name);
    tag = strrep(nm,'VOI_','');
    lista(end+1) = struct('ruta',fullfile(dPorc(k).folder,dPorc(k).name), ...
        'etiqueta',sprintf('Porcino %s',tag), 'grupo','porcino', ...
        'espec',tag, 'sitio','talar', 'tipo','cubico'); %#ok<AGROW>
end

campos = {'BVTV','PoTot','BSBV','BSTV','BS','BV','TV','TbTh','TbSp','TbN', ...
          'DA','DA2','FracPort'};

fprintf('=========================================================\n');
fprintf(' REFERENCIA MATLAB — %d VOIs\n', numel(lista));
fprintf('=========================================================\n\n');

res = struct('etiqueta',{},'grupo',{},'espec',{},'sitio',{},'tipo',{}, ...
             'archivo',{},'dims',{},'spacing',{},'metricas',{}, ...
             'dir_principal',{},'eigenvalues',{},'MIL_valid',{}, ...
             'MIL_clamped',{},'tiempo_s',{});

for k = 1:numel(lista)
    L = lista(k);
    fprintf('[%2d/%2d] %-28s ', k, numel(lista), L.etiqueta);

    t0 = tic;
    try
        [~,~,ext] = fileparts(L.ruta);
        if strcmpi(ext,'.mat')
            S = load(L.ruta);
            VOI = S.VOI;
            if isfield(S,'spacing'), spacing = double(S.spacing(:)'); ...
            else, spacing = [1 1 1]; end
        else
            [VOI, spacing] = readVTKVOI(L.ruta);
        end
    catch ME
        fprintf('OMITIDO (%s)\n', ME.message);
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

    fprintf('%3dx%3dx%3d  BV/TV=%.4f  DA=%.4f  (%.1f s)\n', ...
            size(VOI,1), size(VOI,2), size(VOI,3), vals(1), vals(11), t);

    res(end+1) = struct('etiqueta',L.etiqueta,'grupo',L.grupo, ...
        'espec',L.espec,'sitio',L.sitio,'tipo',L.tipo, ...
        'archivo',L.ruta,'dims',size(VOI),'spacing',spacing(:)', ...
        'metricas',vals,'dir_principal',m.dir_principal(:)', ...
        'eigenvalues',m.eigenvalues(:)', ...
        'MIL_valid',m.MIL_valid,'MIL_clamped',m.MIL_clamped, ...
        'tiempo_s',t); %#ok<AGROW>
end

meta = struct('descripcion', ...
    'Morfometria MATLAB (applib/localMorphometryFromVoxels) sobre todos los VOIs disponibles.', ...
    'campos',{{campos}}, 'nVOIs',numel(res), 'matlabVersion',version, ...
    'fecha',datestr(now,'yyyy-mm-dd HH:MM:SS'), 'resultados',res); %#ok<TNOW1,DATST>

fid = fopen(fullfile(SALIDA,'referencia_vois_todos.json'),'w');
fwrite(fid, jsonencode(meta,'PrettyPrint',true), 'char');
fclose(fid);

fprintf('\n=========================================================\n');
fprintf(' %d VOIs medidos. Escrito: %s\n', numel(res), ...
        fullfile(SALIDA,'referencia_vois_todos.json'));
fprintf('=========================================================\n');
