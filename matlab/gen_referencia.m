%% =====================================================================
%  gen_referencia.m — LADO MATLAB DE LA COMPARACION CON EL PORT A PYTHON
%
%  Genera la referencia contra la que se contrasta Port_Python/spinpy.
%
%  DISENO DEL EXPERIMENTO
%  ----------------------
%  Para que la comparacion mida el GENERADOR y no el camino de medida, ambos
%  lados recorren la MISMA ruta voxel-primero:
%
%      spinodoid(...)  ->  GRF, levelset  ->  BW = GRF <= levelset
%                      ->  localMorphometryFromVoxels(BW, spacing)
%
%  `spinodoid` de GIBBON devuelve el campo GRF (salida 4) y el umbral
%  (salida 8), de modo que se puede umbralizar directamente sin pasar por la
%  malla triangular. Asi cualquier diferencia con Python es atribuible al
%  muestreo de las ondas y a las fases, no al mallado ni a la voxelizacion.
%
%  Se registra ADEMAS la ruta de malla completa (evaluateSpinodoidCandidate),
%  que es la que usa realmente la app, para cuantificar cuanto cuesta el
%  rodeo malla -> voxelizacion frente a umbralizar el campo directamente.
%
%  Los puntos P1-P3 son los mismos de run_ruido_estocastico.m, de modo que la
%  banda de ruido ya medida (resultados/ruido_estocastico.csv) es directamente
%  aplicable. P4 y P5 tienen conos de ANGULOS DISTINTOS: son los casos donde
%  el esquema de muestreo de GIBBON (rechazo) y el del repositorio del
%  profesor (reparto equitativo) divergen.
%
%  Salida: Port_Python/resultados/referencia_matlab.json
% =====================================================================

clear; clc;
warning('off','MATLAB:MKDIR:DirectoryExists');

AQUI    = fileparts(mfilename('fullpath'));
PROY    = fileparts(AQUI);
GIBBON  = 'C:\Users\carlo\OneDrive\Escritorio\Adds-On Matlab\GIBBON-master\lib';
APPLIB  = fullfile(PROY, 'Validacion_Anexo', 'applib');

addpath(GIBBON);
addpath(APPLIB);

SALIDA = fullfile(AQUI, 'resultados');
if ~exist(SALIDA,'dir'), mkdir(SALIDA); end

% --- Contexto de la app: VOI de 64^3 a 0.02 mm/voxel (1.28 mm de lado) ----
global APPCTX APPFIG %#ok<GVMIS>
APPCTX = struct();
APPCTX.VOI              = false(64,64,64);
APPCTX.VOI_spacing      = [1 1 1]*0.02;
APPCTX.VOI_physicalSize = [1.28 1.28 1.28];
APPFIG = [];

% --- Configuracion --------------------------------------------------------
RESOLUCION = 64;                 % rejilla del GRF = rejilla de medida
LADO_MM    = 1.28;               % tamano fisico del VOI
SPACING    = (LADO_MM/RESOLUCION) * [1 1 1];
NREP       = getenv_num('SPIN_NREP', 6);      % realizaciones por punto
NREP_MALLA = getenv_num('SPIN_NREP_MALLA', 2);% cuantas repiten la ruta de malla
SEMILLA0   = 20260805;

puntos = { ...
  struct('nombre','P1_baja_dens',   'dens',0.25,'wave',12,'nw', 600,'th',[15 15 15]), ...
  struct('nombre','P2_media_dens',  'dens',0.35,'wave',15,'nw', 700,'th',[15 15 45]), ...
  struct('nombre','P3_alta_dens',   'dens',0.55,'wave',20,'nw',1000,'th',[15 15 60]), ...
  struct('nombre','P4_conos_desig', 'dens',0.35,'wave',15,'nw', 800,'th',[15 45  0]), ...
  struct('nombre','P5_cono_ancho',  'dens',0.30,'wave',14,'nw', 800,'th',[15 90  0])};

campos = {'BVTV','DA','DA2','BSBV','BS','TbTh','TbSp','TbN','PoTot','FracPort'};

fprintf('=========================================================\n');
fprintf(' REFERENCIA MATLAB PARA EL PORT A PYTHON\n');
fprintf(' resolucion=%d  spacing=%.6g mm  %d puntos x %d realizaciones\n', ...
        RESOLUCION, SPACING(1), numel(puntos), NREP);
fprintf('=========================================================\n\n');

ref = struct('punto',{},'parametros',{},'voxel',{},'malla',{});

for ip = 1:numel(puntos)
    P = puntos{ip};
    fprintf('--- %s : dens=%.2f wave=%gpi nw=%d thetas=[%d %d %d] ---\n', ...
            P.nombre, P.dens, P.wave, P.nw, P.th);

    Mvox = nan(NREP, numel(campos));
    Mmal = nan(NREP_MALLA, numel(campos));
    tvox = nan(NREP,1);

    iS = struct('isocap',true, 'domainSize',1, 'resolution',RESOLUCION, ...
                'waveNumber',P.wave*pi, 'numWaves',P.nw, ...
                'relativeDensity',P.dens, 'thetas',P.th, 'R',eye(3));

    for r = 1:NREP
        rng(SEMILLA0 + 1000*ip + r, 'twister');
        t0 = tic;

        % --- Ruta VOXEL-PRIMERO: se umbraliza el campo directamente -------
        [~,~,~,GRF,~,~,~,levelset] = spinodoid(iS);
        BW = GRF <= levelset;
        m  = localMorphometryFromVoxels(BW, SPACING, true);
        tvox(r) = toc(t0);

        for q = 1:numel(campos)
            if isfield(m, campos{q}) && isscalar(m.(campos{q}))
                Mvox(r,q) = double(m.(campos{q}));
            end
        end
        fprintf('   vox rep %d/%d  BVTV=%.4f  DA=%.3f  BSBV=%.2f  (%.1f s)\n', ...
                r, NREP, Mvox(r,1), Mvox(r,2), Mvox(r,4), tvox(r));
    end

    % --- Ruta de MALLA (la que usa realmente la app) -----------------------
    for r = 1:NREP_MALLA
        rng(SEMILLA0 + 1000*ip + r, 'twister');
        try
            cm = evaluateSpinodoidCandidate(iS);
            for q = 1:numel(campos)
                if isfield(cm, campos{q}) && isscalar(cm.(campos{q}))
                    Mmal(r,q) = double(cm.(campos{q}));
                end
            end
            fprintf('   malla rep %d/%d  BVTV=%.4f  DA=%.3f  BSBV=%.2f\n', ...
                    r, NREP_MALLA, Mmal(r,1), Mmal(r,2), Mmal(r,4));
        catch ME
            fprintf('   malla rep %d/%d  FALLO: %s\n', r, NREP_MALLA, ME.message);
        end
    end

    ref(end+1) = struct( ...
        'punto', P.nombre, ...
        'parametros', P, ...
        'voxel', struct('campos',{{campos}}, 'valores',Mvox, ...
                        'media',nanmean_(Mvox), 'std',nanstd_(Mvox), ...
                        'tiempo_s',nanmean_(tvox)), ...
        'malla', struct('campos',{{campos}}, 'valores',Mmal, ...
                        'media',nanmean_(Mmal), 'std',nanstd_(Mmal))); %#ok<AGROW>
    fprintf('\n');
end

meta = struct('descripcion', ...
    ['Referencia MATLAB voxel-primero (GRF de GIBBON umbralizado) mas la ' ...
     'ruta de malla de la app, para contrastar con Port_Python/spinpy.'], ...
    'resolucion', RESOLUCION, 'lado_mm', LADO_MM, 'spacing_mm', SPACING(1), ...
    'nrep', NREP, 'nrep_malla', NREP_MALLA, 'semilla0', SEMILLA0, ...
    'esquema_muestreo', 'rechazo (GIBBON spinodoid.m)', ...
    'campos', {{campos}}, 'matlabVersion', version, ...
    'fecha', datestr(now,'yyyy-mm-dd HH:MM:SS'), 'resultados', ref); %#ok<TNOW1,DATST>

fid = fopen(fullfile(SALIDA,'referencia_matlab.json'),'w');
fwrite(fid, jsonencode(meta, 'PrettyPrint', true), 'char');
fclose(fid);

fprintf('=========================================================\n');
fprintf(' Escrito: %s\n', fullfile(SALIDA,'referencia_matlab.json'));
fprintf('=========================================================\n');


% -------------------------------------------------------------------------
function v = getenv_num(nombre, porDefecto)
    s = getenv(nombre);
    if isempty(s), v = porDefecto; else, v = str2double(s); end
    if ~isfinite(v) || v < 1, v = porDefecto; end
end

function m = nanmean_(A)
    m = nan(1, size(A,2));
    for j = 1:size(A,2)
        v = A(:,j); v = v(isfinite(v));
        if ~isempty(v), m(j) = mean(v); end
    end
end

function s = nanstd_(A)
    s = nan(1, size(A,2));
    for j = 1:size(A,2)
        v = A(:,j); v = v(isfinite(v));
        if numel(v) > 1, s(j) = std(v); end
    end
end
