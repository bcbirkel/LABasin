
% Translated from the SEISAN bindings to ObsPy core module.

function type=seisan_getVersion(data)

%    Extracts SEISAN version from given data chunk.
%
%    Parameters
%    ----------
%    data : string
%        Data chunk.
%
%    Returns
%    -------
%    type = [ '<' | '>' ,  '32' | '64', ' 6' | '7']
%%
%        Byte order (little endian '<' or big endian '>'), architecture (32 or
%        64) and SEISAN version (6 or 7).%
%
%    From the SEISAN documentation::
%
%        When Fortran writes a files opened with "form=unformatted", additional
%        data is added to the file to serve as record separators which have to
%        be taken into account if the file is read from a C-program or if read
%        binary from a Fortran program. Unfortunately, the number of and meaning
%        of these additional characters are compiler dependent. On Sun, Linux,
%        MaxOSX and PC from version 7.0 (using Digital Fortran), every write is
%        preceded and terminated with 4 additional bytes giving the number of
%        bytes in the write. On the PC, Seisan version 6.0 and earlier using
%        Microsoft Fortran, the first 2 bytes in the file are the ASCII
%        character "KP". Every write is preceded and terminated with one byte
%        giving the number of bytes in the write. If the write contains more
%        than 128 bytes, it is blocked in records of 128 bytes, each with the
%        start and end byte which in this case is the number 128. Each record is
%        thus 130 bytes long. All of these additional bytes are transparent to
%        the user if the file is read as an unformatted file. However, since the
%        structure is different on Sun, Linux, MacOSX and PC, a file written as
%        unformatted on Sun, Linux or MacOSX cannot be read as unformatted on PC
%        or vice versa.
%
%        The files are very easy to write and read on the same computer but
%        difficult to read if written on a different computer. To further
%        complicate matters, the byte order is different on Sun and PC. With 64
%        bit systems, 8 bytes is used to define number of bytes written. This
%        type of file can also be read with SEISAN, but so far only data written
%        on Linux have been tested for reading on all systems.
%
%        From version 7.0,the Linux and PC file structures are exactly the same.
%        On Sun the structure is the same except that the bytes are swapped.
%        This is used by SEISAN to find out where the file was written. Since
%        there is always 80 characters in the first write, character one in the
%        Linux and PC file will be the character P (which is represented by 80)
%        while on Sun character 4 is P.
%    
%% ************************************************************************
%    # check size of data chunk
    load bsval.mat % blank space value
    if (length(data) < 12 * 80);
        type='none';
        return 
    end
    
    if(strcmp(data(1:2),'KP') & strcmp(data(83),'P'))
        type='<326';
    end
        
    if(strcmp(data(1:8),[bs bs bs bs bs bs bs 'P']) & strcmp(data(89:96),[ bs bs bs bs bs bs bs 'P']))
        type='>647';
        return
    end
                
    if(strcmp(data(1:8)',['P' bs  bs  bs bs bs bs bs']) &  strcmp(data(89:96)',[bs bs bs bs bs bs bs 'P']))
        type='<647';
        return 
    end
    
    if(strcmp(data(1:4)',[bs bs bs 'P']) & strcmp(data(85:88)',[bs bs bs 'P']))
        type='>327';
        return
    end
    
    if (strcmp(data(1:4)',['P' bs bs bs]) & strcmp(data(85:88)',['P' bs bs bs]))
        type='<327';
        return
    end
    
    type='none';
%    return None

