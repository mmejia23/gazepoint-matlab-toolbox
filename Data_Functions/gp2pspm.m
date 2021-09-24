function gp2pspm(num_sub, new_sr, experimento)
%
% Funcion para convertir a formato PsPM archivos de Gazepoint Biometrics.
% A partir de archivos .txt que el codigo de Matlab brinda.
%
% ej. gp2pspm(num_sub, hz)
%
% Donde num_sub = numero de participante.
% Y hz = la frecuencia de muestreo del equipo (usualmente 60 o 63 hz).
% La variable 'experimento' permite aplicar la division de trials
% dependiendo del codigo que se uso en el USER data para cada trial.
% - experimento = 1, es '01IAPS' USER = STIMULI_01_POS_XXXX.BMP
%   i.e. STIMULI_#trial_tipodeimagen_codigodeIAPS.BMP
%   y para los inicios de trials: TRIAL_START_01
% - experimento = 4.1, es '04rep' USER = STIMULI_01_POSITIVE_XXXX.BMP
%   i.e. 
%   y para los inicios de trials es: TRIALSTRAT_01_POSITIVE_XXXX.BMP
% - experimento = 4.1, es '04caras'
% 
% Requiere un archivo .txt por participante: Ejemplo:
% - subject_4_SCR_01_IAPS_1_SCR_results.txt: 
% Requiere, ademas, dos archivos template en una carpeta inspeccionada:
% - pspm_SCR_template.mat: 
% - onsets_pspm_template.mat: 
%
% Exporta dos archivos por participante:
% - pspm_exp1_subject_#.mat: 
% - onsets_exp1_subject_#.mat: 
% num_sub = 3;
% hz = 63;
% experimento = 1;
% new_sr = 60;

% Dependiendo del experimento, el nombre del archivo:
if experimento == 1
	txt_filename = 'subject_%i_SCR_01_IAPS_1_SCR_results.txt';
	% export filenames:
	pspmfilename = 'pspm_exp1_subject_%i.mat';
	onsetfilename = 'onsets_exp1_subject_%i.mat';
elseif experimento == 4.1
    txt_filename = 'subject_%i_SCR_04_rep_1_SCR_results.txt';
	% export filenames:
	pspmfilename = 'pspm_exp4rep_subject_%i.mat';
	onsetfilename = 'onsets_exp4rep_subject_%i.mat';
elseif experimento == 4.2
    txt_filename = 'subject_%i_SCR_04_caras_1_SCR_results.txt';
	% export filenames:
	pspmfilename = 'pspm_exp4caras_subject_%i.mat';
	onsetfilename = 'onsets_exp4caras_subject_%i.mat';
end


%% Definicion de variables generales:
sourcefile = sprintf(txt_filename, num_sub);
% Archivos
% Nombres de archivos a exportar: -----------------------------
PsPMSCRfile = sprintf(pspmfilename, num_sub);
PsPMonsetsfile = sprintf(onsetfilename, num_sub);

% Abrir archivos: --------------------------------------------
% Abrir archivos templates
load('pspm_SCR_template.mat')
load('onsets_pspm_template.mat')



%% Seccion de leer archivo txt de Gazepoint a matlab.

% Abrir archivo y leerlo. Se guarda en var gpdata:
fid = fopen(sourcefile);
headers = fgetl(fid);
c = 0;

% Get header names depending on the first row
ff = strfind(headers, '	');
headernames{1} = headers(1:ff(1)-1);
headernames{2} = headers(ff(1)+1:ff(2)-1);
for i = 3:length(ff)
    headernames{i} = headers(ff(i-1)+2:ff(i)-1);
end

while ~feof(fid)
    templine = fgetl(fid);
    if ~strcmp(templine(1:3), 'MSG')
        c = c + 1;
        ff = strfind(templine, '	');
        index = str2num(templine( ff(1)+1 : ff(2)-1));
        if c == 1
            start_CNT = index;
        end
        index = index - start_CNT + 1;
        vals = [templine(ff(2):ff(3)),' ', templine(ff(3):ff(4)), ' ', templine(ff(7):end)];
        
        for i = 2:length(headernames)
            if ~strcmp(headernames{1,i}, 'USER')
                gpdata(index).(headernames{1, i}) = str2num(templine( ff(i-1) : ff(i) ));
            else
                gpdata(index).(headernames{1, i}) = templine( ff(i-1)+1 : ff(i) );
            end
        end
        
        
%         gpdata(c).CNT = str2num(templine(ff(1):ff(2)));
%         gpdata(c).TIME = str2num(templine(ff(2):ff(3)));
%         gpdata(c).GSR = str2num(templine(ff(3):ff(4)));
%         gpdata(c).GSRV = str2num(templine(ff(4):ff(5)));
%         gpdata(c).HR = str2num(templine(ff(5):ff(6)));
%         gpdata(c).HRV = str2num(templine(ff(6):ff(7)));
%         gpdata(c).USER = templine(ff(7)+1:end);
    end
end
fclose(fid);
clearvars c ff fid headers kk templine vals


% La var gpdata.USER, tiene informacion de condiciones, segun la
% programacion del experimento en particular. Aqui se descompone la
% variable USER, y se integran nuevas columnas segun la informacion del
% experimento en particular.
% El dato de GSR en Ohms se convierte a microsiemens:
if experimento == 1 % Exp01 de imagenes de IAPS, afectivas vs neutrales
    % Determine length of session:
    cut_end = sum((~isinf(1./[gpdata.GSR])));
    for i = 1:cut_end
        temprow = gpdata(i).USER;
        if ~isempty(temprow)
            ff = strfind(temprow, '_');
            firstword = temprow(1:ff(1)-1);
            if strcmp(firstword, 'START')
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'CLIENT2')
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'TRIAL')
                gpdata(i).TRIAL = str2num(temprow(ff(2)+1:end));
                gpdata(i).TRIALTYPE = 'BASE';
                gpdata(i).STIMULI = 'BASE';
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'STIMULI')
                gpdata(i).TRIAL = str2num(temprow(ff(1)+1:ff(2)-1));
                gpdata(i).TRIALTYPE = temprow(ff(2)+1:ff(3)-1);
                gpdata(i).STIMULI = temprow(ff(3)+1:end);
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            else
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            end
        else
            % Aqui agregar codigo si los valores vacios sera mejor tenerlos
            % como NaN!!
        end
        if isinf(gpdata(i).CONDUCT) %%%%%%%%%%%%%%%%%%%% Corregir:
            gpdata(i).CONDUCT = NaN;
        end
    end
    
    clearvars ff firstword i temprow
end


if experimento == 4.1 % Exp04 de imagenes de IAPS, replica de Bach et al. 2013
    % Determine length of session:
%     cut_end = sum((~isinf(1./[gpdata.GSR])));
    cut_end = length(gpdata);
    for i = 1:cut_end
        temprow = gpdata(i).USER;
        if ~isempty(temprow)
            ff = strfind(temprow, '_');
            firstword = temprow(1:ff(1)-1);
            if strcmp(firstword, 'START')
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'CLIENT2')
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'TRIALSTART')
                gpdata(i).TRIAL = str2num(temprow(ff(1)+1:ff(2)-1));
                gpdata(i).TRIALTYPE = 'BASE';
                gpdata(i).STIMULI = 'BASE';
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'STIMULI')
                gpdata(i).TRIAL = str2num(temprow(ff(1)+1:ff(2)-1));
                gpdata(i).TRIALTYPE = temprow(ff(2)+1:ff(3)-1);
                gpdata(i).STIMULI = temprow(ff(3)+1:end);
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            else
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            end
        else
%            % Aqui agregar codigo si los valores vacios sera mejor tenerlos
%            % como NaN!! No meter NaNs.
            gpdata(i).USER = 'empty';
           gpdata(i).CONDUCT = NaN;
           gpdata(i).TIME = NaN;
        end
        if isinf(gpdata(i).CONDUCT)
            gpdata(i).CONDUCT = NaN;
        end
    end
    
    clearvars ff firstword i temprow
end


if experimento == 4.2 % Exp04 de imagenes de IAPS, replica de Bach et al. 2013
    % Determine length of session:
%     cut_end = sum((~isinf(1./[gpdata.GSR])));
    cut_end = length(gpdata);
    for i = 1:cut_end
        temprow = gpdata(i).USER;
        if ~isempty(temprow)
            ff = strfind(temprow, '_');
            firstword = temprow(1:ff(1)-1);
            if strcmp(firstword, 'START')
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'CLIENT2')
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'TRIALSTART')
                gpdata(i).TRIAL = str2num(temprow(ff(1)+1:ff(2)-1));
                gpdata(i).TRIALTYPE = 'BASE';
                gpdata(i).STIMULI = 'BASE';
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            elseif strcmp(firstword, 'STIMULI')
                gpdata(i).TRIAL = str2num(temprow(ff(1)+1:ff(2)-1));
                gpdata(i).TRIALTYPE = temprow(ff(2)+1:ff(3)-1);
                gpdata(i).STIMULI = temprow(ff(3)+1:end);
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            else
                gpdata(i).TRIAL = [];
                gpdata(i).TRIALTYPE = [];
                gpdata(i).STIMULI = [];
                gpdata(i).CONDUCT = 1/gpdata(i).GSR * 1000000;
                gpdata(i).TIMEREL = gpdata(i).TIME - gpdata(1).TIME;
            end
        else
%            % Aqui agregar codigo si los valores vacios sera mejor tenerlos
%            % como NaN!! No meter NaNs.
            gpdata(i).USER = 'empty';
           gpdata(i).CONDUCT = NaN;
           gpdata(i).TIME = NaN;
        end
        if isinf(gpdata(i).CONDUCT)
            gpdata(i).CONDUCT = NaN;
        end
    end
    
    clearvars ff firstword i temprow
end


fprintf('\nCantidad de datos cargados: %.0f\n', length(gpdata));



%% Plot data to debug Borrar eventualmente
% subplot(2, 1, 1); plot(([gpdata.CNT]-25054), [gpdata.GSR])
% subplot(2, 1, 2); plot([gpdata.TIME]-393.9983, [gpdata.GSR])
% 
% subplot(2, 1, 1); plot(([gpdata.CNT]-25054)/63, [gpdata.GSR])
% subplot(2, 1, 2); plot([gpdata.TIME]-393.9983, [gpdata.GSR])
% xlim([0 60])
% 
% % data has 63.5891 Hz????
% subplot(2, 1, 1); plot([gpdata.CNT]/63.5891, [gpdata.GSR])
% subplot(2, 1, 2); plot([gpdata.TIME], [gpdata.GSR])
% 
% plot([gpdata.CNT] ./ [gpdata.TIME])
% mean([gpdata.CNT] ./ [gpdata.TIME])


%% Estimar sampling rate original con base en el reloj del equipo:
distances = [];
datapoints = [];
for i = 1:length(gpdata)-1
    if ~isempty(gpdata(i+1).TIME) & ~isempty(gpdata(i).TIME)
        distances(i) = [gpdata(i+1).TIME] - [gpdata(i).TIME];
    else
        distances(i) = NaN;
    end
    if ~isempty(gpdata(i+1).CNT) & ~isempty(gpdata(i).CNT)
    datapoints(i) = [gpdata(i+1).CNT] - [gpdata(i).CNT];
    else
        datapoints(i) = NaN;
    end
end
% plot(distances)
% 0.0159 equal to 62.8931 Hz
fprintf('\n----------------------------------------------\n');
fprintf('\nSampling rate original: %.05f Hz\n', 1/mean(distances, 'omitnan'));

% Revisar esta parte... da numeros demasiado grandes, superiores a la
% cantidad de datapoints que hay....
% Estimar datos o paquetes perdidos segun el conteo del equipo: CNT
% plot(datapoints) % straight line on y=1 means no lost datapoints
fprintf('\nDatos perdidos en el registro: %.0f\n', sum(isnan(datapoints)));

% Segundos segun el reloj del equipo:
fprintf('\nDuracion de la sesion segun reloj: %.03f seg\n', max([gpdata.TIME]) - min([gpdata.TIME]));
fprintf('\n----------------------------------------------\n');


%% Ajustar sample rate (que era aprox de 62.8931 Hz, sampling interval 0.0159)
% Si aplicar esta seccion!!! Downsample (interpolando) a 60 Hz exactos, para poder pasar
% esta data a pspm y otros downsamples!!
% Upsample?????
vartime = transpose([gpdata(1:cut_end).TIME] *100000); % Para tener numeros enteros en TIME
vartime = vartime - min(vartime); % Se ajusta inicio de TIME a cero.
varconduct = transpose([gpdata(1:cut_end).CONDUCT]); % Para plotear onda, y verificar cosas.
fs = new_sr/100000; % 60 muestras cada 100000 clockmarks (cada segundo)=60Hz.
[ds_data, ty] = resample(varconduct, vartime, fs, 1, 1, 'linear'); % quitar 'spline', y dejar linear???
ty = ty / 100000; % Reconvertir tiempo a segundos
% subplot(2,1,1); plot(ty, y)
% En su caso, cambiar valor de hz por el nuevo
hz = new_sr;
fprintf('\nNuevo sampling rate despues de downsample: %.05f Hz\n', new_sr);
fprintf('\nDuracion de la sesion segun downsample: %.03f seg\n', length(ds_data)/new_sr);
fprintf('\n----------------------------------------------\n');


%% Take care of NaN values on the original data after resampling:
% In this section, we find the previously NaN values in the data,
% then, we look for groups of them. And then, we look for the beginning and
% end of each group, to index the time of those groups, and convert the
% downsampled data to NaN again. Therefore, the downsample doesnt imput
% those values interpolating them, since we have lots of NaN values.
% %{
nanvalues = find(isnan([gpdata.CONDUCT]));
if ~isempty(nanvalues)
    % This is used to select if the extended NaN conversion of missing values
    % is applied. Right now, it is applied always.
    nan_extended_trials = 1;
    group = 1;
    groups(group).values = nanvalues(1);
    for i = 2:length(nanvalues)
        if nanvalues(i) == nanvalues(i-1)+1
            groups(group).values = [groups(group).values nanvalues(i)];
        else
            group = group + 1;
            groups(group).values = nanvalues(i);
        end
    end

    for i = 1:length(groups)
        if groups(i).values(1) ~= 1
            groups(i).startofgroup = find(round(gpdata(groups(i).values(1)-1).TIMEREL, 2) == round(ty, 2));
            groups(i).endofgroup = find(round(gpdata(groups(i).values(end)+1).TIMEREL, 2) == round(ty, 2));
            ds_data(groups(i).startofgroup:groups(i).endofgroup) = NaN;

            if nan_extended_trials == 1
                % Find trial to delete previous to NaN block and posterior:
                trial_to_delete_start = gpdata(groups(i).values(1)-1).TRIAL;
                trial_to_delete_end = gpdata(groups(i).values(end)+1).TRIAL;
                % Find the datapoint where the trial to delete at the start, starts:
                trial_to_delete_start_position = min(find(strncmp({gpdata.USER}, sprintf('TRIALSTART_%02d', trial_to_delete_start), 13)));
                trial_to_delete_end_position = max(find(strncmp({gpdata.USER}, sprintf('STIMULI_%02d', trial_to_delete_end), 10)));
                % Find the corresponding position of the cutoffs in the downsampled data:
                groups(i).startofgroup_extended = find(round(gpdata(trial_to_delete_start_position).TIMEREL, 2) == round(ty, 2));
                groups(i).endofgroup_extended = find(round(gpdata(trial_to_delete_end_position).TIMEREL, 2) == round(ty, 2));
                % Set the group to NaN:
                ds_data(groups(i).startofgroup_extended:groups(i).endofgroup_extended) = NaN;
            end
        end
    end
end
% %}


%% Codigo para armar variables de datos de onsets:

% Armar tabla con datos por trial: num, onset time, type, stimuli data.
uniquenametrials = unique({gpdata.USER});
uniquenametrials = uniquenametrials(strncmp('STIMULI', uniquenametrials, 7));
for i = 1:length(uniquenametrials)
    indexoftrial = find(strcmp({gpdata.USER}, uniquenametrials(i)));
    startoftrial = min(indexoftrial);
    trial(i).CNT = gpdata(startoftrial).CNT;
    trial(i).TIME = gpdata(startoftrial).TIME;
    trial(i).TIMEREL = gpdata(startoftrial).TIMEREL;
    trial(i).TRIAL = gpdata(startoftrial).TRIAL;
    trial(i).TRIALTYPE = gpdata(startoftrial).TRIALTYPE;
    trial(i).STIMULI = gpdata(startoftrial).STIMULI;
end

conditions = unique({trial.TRIALTYPE});
howmanyconditions = length(conditions);
for i = 1:howmanyconditions
%    onsets{i} = [trial(find(strcmp(conditions(i), {trial.TRIALTYPE}))).TRIAL];
    % Lets try to run this line instead of the one above, for data where
    % there are missing trials. Since the trial number wont correspond to
    % the trial "index".
    onsets{i} = find(strcmp(conditions(i), {trial.TRIALTYPE}));
    names{i} = [conditions{i}];
end


%% 
% Guardar en template, los datos SCR --------------------------
data{1, 1}.data = ds_data;

% Guardar en template, los datos de onsets --------------------
data{2, 1}.data = transpose([trial.TIMEREL]);

% Guardar variables generales del registro: frecuencia de muestreo (hz),
% duracion de la sesion (datapoints / hz), unidad de medida (microsiemens
% por conductancia), tipo de equipo, nombres de archivos originales y
% exportados: -------------------------------------------------
infos.duration = length([data{1, 1}.data]) / hz;
data{1, 1}.header.sr = hz;
data{1, 1}.header.units = 'microsiemens';
infos.source.type = 'Gazepoint Biometrics';
infos.source.file = sourcefile;
infos.importfile = PsPMSCRfile;
infos.importdate = date;


%% 

% Guardar nuevos archivo .mat con los datos del participante. ------
eval(['save ' PsPMSCRfile ' data infos gpdata trial'])
eval(['save ' PsPMonsetsfile ' names onsets'])

clearvars experimento gpdata data infos names onsets num_sub sourcefile
