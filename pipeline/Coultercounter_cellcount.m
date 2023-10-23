%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This code processes Coulter counter .#m4 files to plot single-cell
%  volume data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all, clc
currentFolder = pwd;
% Input UI to grab path to FBM_metadata_assembly_instruction sheet
fprintf('\nGetting coulter counter cell count instruction...\n')
[input_info.instruction_filename, input_info.instruction_dir, exist_pmt] = uigetfile('../*.*','Select FBM metadata assembly instruction File',' ');

instruction_path = [input_info.instruction_dir,'\',input_info.instruction_filename];
opts = detectImportOptions(instruction_path,'ReadVariableNames',true,'VariableNamingRule','preserve','Delimiter',',');
opts = setvartype(opts,'string');
instruction = readtable(instruction_path,opts);

report=instruction;
%%
for i = 1:height(instruction)
    close all
% open the files
%     [fname,pname] = uigetfile('*.#m4',['Select file to process: ']); % selects file
    fname = instruction.filename(i);
    pname = instruction.path(i);

    cd(pname); % changes the current directory to the one with the file
    fid = fopen(fname); %opens file for reading
    if fid == -1 % if file can't be opened
        disp('File could not be opened');
        return % end program
    end

% declare constant
    countspervolt = 1/(4*298.02e-9);
    
% get the following data from file: nPulses, Kd, current,
% resistance, MaxHtCorr, pulse data    
    
    % get the number of pulses
    linenum = 135;
    aPulses = textscan(fid, '%*9s %n', 1, 'delimiter',...
        '\n', 'headerlines', linenum-1); 
    nPulses = aPulses{1};
    fseek(fid, 0, 'bof'); % resets pointer to bof
        
    % get the Kd value
    linenum = 127;
    aKd = textscan(fid, '%*4s %n', 1, 'delimiter',...
        '\n', 'headerlines', linenum-1);
    Kd = aKd{1}; % converts the cell array into an integer
    fseek(fid, 0, 'bof'); % resets pointer to the beginning of file

    % get the aperture current (mA)
    linenum = 178;
    aCurrent = textscan(fid, '%*9s %n', 1, 'delimiter',...
        '\n', 'headerlines', linenum-1); 
    current = aCurrent{1}/1000; % gets current as an integer in mA
    fseek(fid, 0, 'bof'); % resets pointer to the beginning of file

    % get the gain and convert it to resistance
    linenum = 184;
    gain = textscan(fid, '%*6s %n', 1, 'delimiter',...
        '\n', 'headerlines', linenum-1); 
    resistance = 25*gain{1}; % convert gain to equiv resistance (kohms)
    fseek(fid, 0, 'bof'); % resets pointer to the beginning of file

    % get the MaxHeight Correction
    linenum = 187;
    MaxHtCorr = textscan(fid, '%*11s %n', 1, 'delimiter',...
        '\n', 'headerlines', linenum-1); 
    fseek(fid, 0, 'bof'); % resets pointer to the beginning of file

    % get the pulse data
    linenum = 1533;
    pulsearray = textscan(fid,'%s %*s %*s %*s %*s', nPulses, ...
        'delimiter', ',', 'headerlines', linenum-1); 

    st = fclose(fid); % closes file after getting all useful data

    % convert the pulse data to volume
    heightt = (hex2dec(pulsearray{1}) + MaxHtCorr{1})';
    diameter = Kd*((heightt./(countspervolt*resistance*current)).^(1/3));
    volume = 4/3*pi*(diameter/2).^3;
    
    % get cells
    volume_cells = volume(volume>double(instruction.cutoff_low_fL(i))&volume<double(instruction.cutoff_high_fL(i)));
    report.cell_count(i) = length(volume_cells);
    
% Plot the pulse data
    histogram(volume,'BinWidth',20);
    hold on
    histogram(volume_cells,'BinWidth',20);
    xlim([50,10000])
    set(gca,'xscale','log')
    title(strrep(instruction.filename(i),"_"," "))
%     input_pause = input('pass');
end

%%
cd(input_info.instruction_dir)
out_file_name = ['Report_coulter_cellcount.txt'];
writetable(report,out_file_name, 'delimiter', '\t');
cd(currentFolder)