function SMD=loadH5(DataDir,FileName)

if ~strcmp('h5',FileName(end-1:end)) && ~strcmp('f5',FileName(end-1:end))
    error('File type has to be given at the end of file name.'); 
end
SMD=h5read(fullfile(DataDir,FileName),'/locs');

end