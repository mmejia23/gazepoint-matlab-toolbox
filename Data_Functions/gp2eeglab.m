function [gpdata, events, EEG] = gp2eeglab(sourcefile, new_sr)
%
% Funcion para convertir a formato PsPM archivos de Gazepoint Biometrics.
% A partir de archivos .txt que el codigo de Matlab brinda.
%
% ej. gp2eeglab('CF6_4_1.txt', hz)
% Donde:
% sourcefile = el nombre del archivo a convertir.
% new_sr = la frecuencia de muestreo del equipo (usualmente 60 o 63 hz).
%
% Requiere tener instalado EEGlab, particularmente la funcion eeg_emptyset.
%
%

% Nombres de archivos a exportar: -----------------------------
eeglabfile = [sourcefile(1:end-4) '.m'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Seccion de leer archivo txt de Gazepoint a matlab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Abrir archivo y leerlo. Se guarda en var gpdata:
fid = fopen(sourcefile);
headers = fgetl(fid);
c = 0;
event_counter = 0;

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
	else
		event_counter = event_counter + 1;
		ff = strfind(templine, '	');
		index = str2num(templine( ff(1)+1 : ff(2)-1));
		events(event_counter).type = templine( ff(3)+1 : end);
		events(event_counter).latencySecs = str2num(templine( ff(2)+1 : ff(3)-1));
    end
end
fclose(fid);
clearvars c ff fid headers kk templine vals



fprintf('\nCantidad de datos cargados: %.0f\n', length(gpdata));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load to EEGlab structure:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create empty EEG structure.
EEG = eeg_emptyset;
 
% Define basic items of EEG structure.
EEG        = eeg_emptyset();
EEG.data   = [gpdata.GSR; gpdata.LPMM; gpdata.RPMM];
EEG.times  = [gpdata.TIME];
EEG.xmin   = EEG.times(1);
EEG.xmax   = EEG.times(end); 
EEG.srate  = round(1/((EEG.xmax-EEG.xmin)/length(EEG.times))); % Rounded actual sampling rate. Note that the unit of the time must be in second.
EEG.nbchan = size(EEG.data,1);
EEG.pnts   = size(EEG.data,2);
 
% Define event information.
eventStructure = struct('type', [], 'latency', []);
latencyValues  = num2cell(([events.latencySecs]-EEG.xmin)*EEG.srate);
[eventStructure(1:length([events.latencySecs])).latency] = latencyValues{:};
eventTypes     = {events.type};
[eventStructure(1:length([events.latencySecs])).type] = eventTypes{:};
EEG.event = eventStructure;
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');


end