% fBlobVisualse(outdir, stamp) visualises the vortex blobs for runs in
% output directory outdir with timestamp stamp.
%
% Optional Arguments
% plot_vorticity - plot the vorticity field (default true)
% plot_contours  - do  contour plot (default false)
% save_avi       - save a video of the visualisation (default false)
% avi_name       - file to write to (default stamp + flowfield.avi)
%
% GD, MAB

function fBlobVisualise(outdir, stamp, varargin)

% Deal with the optional arguments
plot_vorticity = true;
plot_contours = false;
save_avi = false;
avi_name = '';
if ~isempty(varargin)
    if mod(length(varargin), 2) ~= 0
        error('Optional arguments must be in name-value pairs')
    end
    
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'plot_vorticity'
                plot_vorticity = varargin{i+1};
            case 'plot_contours'
                plot_contours = varargin{i+1};
            case 'save_avi'
                save_avi = varargin{i+1};
            case 'avi_name'
                avi_name = varargin{i+1};
            otherwise
                error('Argument not recognised: %s', varargin{i})
        end
    end
end

% Start the processing
run = fullfile(outdir, fCheckStampEnd(stamp));
coord_file = fGetCoordFile(run);

vortex_file = [run, 'vortex.mat'];
if exist(vortex_file, 'file')
    fprintf('Loading the mat file\n');
    load(vortex_file);
else
    fid = fOpenFile([run, 'vortex.dat']);
    
    nodes       = cell2mat(textscan(fid,'%f',1,'CommentStyle','#'));
    step        = cell2mat(textscan(fid,'%f,',1,'CommentStyle','#'));
    steps       = cell2mat(textscan(fid,'%f,',1,'CommentStyle','#'));
    while ~feof(fid)
        A = textscan(fid,'%f %f %f %f','headerlines',1);
        t = A{1};
        x = A{2};
        z = A{3};
        circ = A{4};
    end
    fclose(fid);
    save(vortex_file,'t','x','z','circ','nodes','steps','step')
end

fidNum = fOpenFile([run,'vortex_num.dat']);
fidBody = fOpenFile(coord_file);
B = textscan(fidNum,'%f %f','headerlines',1);
C = textscan(fidBody,'%f %f','headerlines',1);
fclose(fidNum);
fclose(fidBody);

numvb=B{2};
bodyx=C{1};
bodyz=C{2};

if(plot_vorticity)
    scrsz = get(0,'ScreenSize');
    fig=figure('Position',[100 scrsz(4)/2-100 scrsz(3)*2/3 scrsz(4)/2]);
    colormap(fig,'jet')
    caxis([-0.001 0.001])
    colorbar
    if (save_avi)
        if isempty(avi_name)
            avi_name = [run, flowstream.avi];
        end
        aviobj = VideoWriter(avi_name);
        open(aviobj);
    end
    
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    
    index_start = 1;
    
    for n=1:steps-1
        index_end = index_start + numvb(n);
        plot(bodyx,bodyz,'k-');hold on;
        scatter(x(index_start:index_end),z(index_start:index_end),4,circ(index_start:index_end),'filled');
        title(['t = ',num2str(t(index_start),'%5.3f')]);
        ylabel('z[m]');
        xlim([-2*max(bodyx) 20*max(bodyx)]);
        ylim([-7*max(bodyz) 7*max(bodyz)]);
        
        index_start = index_end;
        
        
        caxis([-0.001 0.001])
        colorbar
        
        frame = getframe;
        hold off;
        if save_avi
            writeVideo(aviobj,frame);
        end
    end
    if save_avi
        close(aviobj);
    end
    
end

R=0.5;
Center(1)=0;
Center(2)=0;
dx=0.1;
dy=0.1;
[X,Y]=meshgrid(-2*R+Center(1):dx:15*R+Center(1),-3*R+Center(2):dy:3*R+Center(2));
Z=X+1i*Y;
psi_inf=imag(2.5*(Z+R^2./Z));

if(plot_contours)
    scrsz = get(0,'ScreenSize');
    
    fig=figure('Position',[100 scrsz(4)/2-100 scrsz(3)*2/3 scrsz(4)/2]);
    colormap(fig,'jet')
    colorbar
    if (save_avi)
        if isempty(avi_name)
            avi_name = [run, flowstream_contours.avi];
        end
        aviobj = VideoWriter(avi_name);
        open(aviobj);
    end
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    
    index_start = 1;
    
    for n=1:steps-1
        index_end = index_start + numvb(n);
        % Compute the stream function on the Grid
        psi_vort=streamfunction2D(X,Y,x(index_start:index_end),z(index_start:index_end),circ(index_start:index_end),0,numvb(n));
        psi=psi_inf+psi_vort;
        a=sqrt((X-Center(1)).^2+(Y-Center(2)).^2);
        psi(a<R)=NaN;
        plot(bodyx,bodyz,'k-');hold on;
        contour(X,Y,psi,60);
        title(['t = ',num2str(t(index_start),'%5.3f')]);
        ylabel('z[m]');
        ylim([-3*max(bodyz) 3*max(bodyz)]);
        index_start = index_end;
        xlim([-2*max(bodyx) 15*max(bodyx)]);
        %caxis([-0.001 0.001])
        axis equal
        colorbar
        frame = getframe;
        hold off;
        if save_avi
            writeVideo(aviobj,frame);
        end
    end
    if save_avi
        close(aviobj);
    end
    
end

end

% fOpenFile opens a file and checks the ID is valid
function fid = fOpenFile(file)
fid = fopen(file);

if fid < 0
    error('cannot open file: %s', file)
end

end

% fCheckStampEnd ensures that the timestamp has a trailing _
function stamp = fCheckStampEnd(stamp)

if stamp(end) ~= '_'
    stamp = [stamp, '_'];
end

end

% fGetCoordFile get the full path of the input coordinate file from the
% xml file written in the output directory
function file = fGetCoordFile(run)
xmlfile = [run, 'xml_in.xml'];
xml = xmlread(xmlfile);
fGetXmlString = @(tag) char(xml.getElementsByTagName(tag).item(0).getAttribute('string'));

file = fullfile(fGetXmlString('input_dir'), fGetXmlString('domain_file'));

end

