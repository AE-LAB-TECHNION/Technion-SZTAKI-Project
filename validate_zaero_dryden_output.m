function comparison = validate_zaero_dryden_output(pltFile,outFile)
%VALIDATE_ZAERO_DRYDEN_OUTPUT Check DRY.PLT against RMS values in ZAERO .out.
% This is a ZAERO-internal consistency check; no MATLAB plant is involved.

if nargin < 1, pltFile = 'DRY.PLT'; end
if nargin < 2, outFile = 'ASE_ANALYSIS_DRYDEN.out'; end

data = load_zaero_dryden_response(pltFile);
f = data.frequency;
signals = [data.modalX,data.modalV];
types = [repmat("MODALX",1,3),repmat("MODALV",1,3)];
modes = [1:3,1:3];

outText = fileread(outFile);
pattern = ['L\s*A\s*B\s*E\s*L\s*=\s*(RE|IM)\(\s*(MODAL[XV])\s*\)' ...
    '.*?I\s*K\s*E\s*Y\s*=\s*(\d+).*?' ...
    'R\s*M\s*S\s+V\s*A\s*L\s*U\s*E\s*=\s*([+\-0-9.Ee]+)'];
tokens = regexp(outText,pattern,'tokens');
if isempty(tokens)
    error('No modal RMS entries were found in %s.',outFile);
end

n = numel(signals(1,:));
response = strings(2*n,1);
component = strings(2*n,1);
rmsFromPlt = zeros(2*n,1);
rmsFromOut = nan(2*n,1);
row = 0;
for j = 1:n
    for part = 1:2
        row = row+1;
        if part == 1
            component(row) = "RE";
            spectralComponent = real(signals(:,j));
        else
            component(row) = "IM";
            spectralComponent = imag(signals(:,j));
        end
        response(row) = sprintf('%s %d',types(j),modes(j));
        rmsFromPlt(row) = sqrt(trapz(f,spectralComponent.^2));
        for k = 1:numel(tokens)
            if strcmp(tokens{k}{1},component(row)) && ...
                    strcmp(tokens{k}{2},types(j)) && ...
                    str2double(tokens{k}{3}) == modes(j)
                rmsFromOut(row) = str2double(tokens{k}{4});
                break
            end
        end
    end
end

comparison = table(response,component,rmsFromPlt,rmsFromOut, ...
    rmsFromPlt./rmsFromOut,100*(rmsFromPlt-rmsFromOut)./rmsFromOut, ...
    'VariableNames',{'Response','Component','RmsFromDryPlt','RmsPrintedByZaero', ...
    'Ratio','DifferencePercent'});

fprintf('\nZAERO internal Dryden RMS consistency check\n');
fprintf('DRY.PLT: %s\nZAERO output: %s\n\n',pltFile,outFile);
disp(comparison);
end
