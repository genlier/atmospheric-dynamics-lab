function run_dense_bve(labels)
if nargin==0, labels=1:5; end
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'outputs','dense');
if ~exist(out,'dir'), mkdir(out); end
source=fileread(fullfile(root,'models','matlab','bve.m'));
source=strrep(source,'clc; clear all; close all;','');
first=strfind(source,'figure(1);'); last=strfind(source,'tic;');
source=[source(1:first(1)-1) source(last(1):end)];
source=strrep(source,'if (mod(t, p.output_fre) == 0 || t == p.Step)','if false');
split=strfind(source,'function [dealias]');
main=source(1:split(1)-1); helpers=source(split(1):end);
setup=sprintf(['animation_steps=unique(round(linspace(0,p.Step,601)));\n' ...
 'animation_times=animation_steps*p.dt;\n' ...
 'animation_snapshots=zeros(p.N1,p.N2,numel(animation_steps));\n' ...
 'animation_snapshots(:,:,1)=Iz0; animation_index=2;\n' ...
 'animation_spectra=zeros(floor(p.N1/3),numel(animation_steps));\n' ...
 '[animation_k,animation_spectra(:,1)]=compute_spectrum(Iz1_hat,inv_ksq,kx_array,ky_array,k_array,p.N1);\n']);
main=strrep(main,'tic;',[setup newline 'tic;']);
capture=sprintf(['\n    if animation_index<=numel(animation_steps) && t==animation_steps(animation_index)\n' ...
 '        animation_snapshots(:,:,animation_index)=real(ifft2(Iz1_hat))*(p.N1*p.N2);\n' ...
 '        [animation_k,animation_spectra(:,animation_index)]=compute_spectrum(Iz1_hat,inv_ksq,kx_array,ky_array,k_array,p.N1);\n' ...
 '        animation_index=animation_index+1;\n    end\n']);
main=strrep(main,'Iz1_hat = Iz2_hat;',['Iz1_hat = Iz2_hat;' capture]);
for label=labels
    name=sprintf('dense_case%d',label);
    body=strrep(main,'p.model_case = 1;',sprintf('p.model_case = %d;',label));
    target=strrep(fullfile(out,sprintf('case%d.mat',label)),'''','''''');
    finish=sprintf(['\nassert(animation_index==numel(animation_steps)+1);\n' ...
      'assert(all(isfinite(animation_snapshots(:))));\n' ...
      'save(''%s'',''animation_snapshots'',''animation_times'',''animation_steps'',''animation_spectra'',''animation_k'',''p'',''xpos'',''ypos'',''-v7'');\n' ...
      'fprintf(''Completed case %d: %%d actual frames, %%d integration steps.\\n'',numel(animation_steps),p.Step);\n'],target,label);
    file=fullfile(out,[name '.m']);fid=fopen(file,'w','n','UTF-8');
    fprintf(fid,'function %s\n%s\n%s\nend\n%s',name,body,finish,helpers);fclose(fid);
    addpath(out);feval(name);
end
end
