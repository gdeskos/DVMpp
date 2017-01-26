% fParseForces(outdir, stamp) return the force time history for a run in
% output directory outdir with timestamp stamp.
%
% MAB

function out = fParseForces(outdir, stamp)

% General prefix
run = fullfile(outdir, fCheckStampEnd(stamp));

% Forces
force_file = [run, 'loads.dat'];
res = dlmread(force_file, '', 1, 0);

out.Fx = res(:,2);
out.Fz = res(:,3);

% Time
xml_file = [run, 'xml_in.xml'];
xml = xmlread(xml_file);
dt = str2double(xml.getElementsByTagName('dt').item(0).getAttribute('val'));

out.t = 0:dt:dt*(length(out.Fx)-1);

end

% fCheckStampEnd ensures that the timestamp has a trailing _
function stamp = fCheckStampEnd(stamp)

if stamp(end) ~= '_'
    stamp = [stamp, '_'];
end

end
