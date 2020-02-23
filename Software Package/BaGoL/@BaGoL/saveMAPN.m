function saveMAPN(Directory,FileType,MAPN)
%Saving the MAPN-results in either mat-file or h5-file. 
switch FileType
    case 'mat'
        save(fullfile(Directory,'MAPN'),'MAPN');
    case 'h5'
        SZ = size(MAPN.X);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.X',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN/Y',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.X_SE',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.Y_SE',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.Nmean',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.AlphaX',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.AlphaY',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.AlphaX_SE',SZ);
        h5create(fullfile(Directory,'MAPN.h5'),'/MAPN.AlphaY_SE',SZ);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/X',MAPN.X);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/Y',MAPN.Y);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/X_SE',MAPN.X_SE);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/Y_SE',MAPN.Y_SE);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/Nmean',MAPN.Nmean);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaX',MAPN.AlphaX);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaY',MAPN.AlphaY);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaX_SE',MAPN.AlphaX_SE);
        h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaY_SE',MAPN.AlphaY_SE);
        if ~isempty(MAPN.Z) && isfield(MAPN,'Z')
            h5create(fullfile(Directory,'MAPN.h5'),'/MAPN/Z',SZ);
            h5create(fullfile(Directory,'MAPN.h5'),'/MAPN/Z_SE',SZ);
            h5create(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaZ',SZ);
            h5create(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaZ_SE',SZ);
            h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/Z',MAPN.Z);
            h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/Z_SE',MAPN.Z_SE);
            h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaZ',MAPN.AlphaZ);
            h5write(fullfile(Directory,'MAPN.h5'),'/MAPN/AlphaZ_SE',MAPN.AlphaZ_SE);
        end
    case 'hdf5'
        %Getting dtructure Info
        FileName = 'MAPN.hdf5';
        Dataset = 'locs';
        MAPN=rmfield(MAPN,'N');
        Fields = fieldnames(MAPN);
        for ii = length(Fields):-1:1
             MAPN=setfield(MAPN,Fields{ii},single(getfield(MAPN,Fields{ii})));
        end
        NFields = length(Fields);
        Dims = length(getfield(MAPN,Fields{1}));
        %Creating an empty hdf5-file in the given directory
        File = H5F.create(fullfile(Directory,FileName), 'H5F_ACC_TRUNC', ...
                'H5P_DEFAULT', 'H5P_DEFAULT');
        %Creating single data types
        SingleType = H5T.copy('H5T_NATIVE_FLOAT');
        SZt = H5T.get_size(SingleType);
        SZ = SZt*ones(1,NFields);
        %Computer offsets to each field. The first offset is always zero.
        offset(1) = 0;
        offset(2:NFields) = cumsum(SZ(1:NFields-1));
        %Creating compound datatype for memory
        Memtype = H5T.create('H5T_COMPOUND',sum(SZ));
        for ii = 1:NFields
            H5T.insert(Memtype,Fields{ii},offset(ii),SingleType); 
        end
        %Create the compound datatype for the file.
        FileType = H5T.create('H5T_COMPOUND', sum(SZ));
        for ii = 1:NFields
            H5T.insert(FileType,Fields{ii},offset(ii),SingleType); 
        end
        % Create dataspace. Setting the maximum size to [] sets 
        %  the maximum size to be current size
        Space = H5S.create_simple(1,fliplr(Dims),[]);
        %Creat the dataset and write the compound data to it
        Dset = H5D.create(File, Dataset, FileType, Space, 'H5P_DEFAULT');
        H5D.write(Dset,Memtype,'H5S_ALL','H5S_ALL','H5P_DEFAULT',MAPN);
        %Close and release resources
        H5D.close(Dset);
        H5S.close(Space);
        H5T.close(FileType);
        H5F.close(File);
end
end
