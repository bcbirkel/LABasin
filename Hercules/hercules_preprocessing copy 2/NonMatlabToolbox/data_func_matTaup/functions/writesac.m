% writesac.m
% function to write a sac file. Returns 0 on failure.
% syntax:
% ierr = writesac('output_file_name',data,header);
% Because I typically work with two sac readers, this checks whether header
% is a structured array or a vector and makes assumptions about the format
% readsac() only returns a subset of the sac header, but load_sac returns most headers.
% Therefore the result is a more complete header with load_sac.
%
% The return value is the number of points written to the output file. Clearly, it fails if 0 is returned
% 
%
function ierr = writesac(output_name,data,header)
    ierr = 0;
    %initialize with null values
    
    fheader=zeros(1,70);
    for i=1:70
        fheader(i) = -12345.;
    end
    iheader=zeros(1,35);
    for i=1:35
        iheader(i)=-12345;
    end
    lheader=zeros(1,5);
    for i=1:5
        lheader(i) = -12345;
    end
    
    %kheader = ['abcdefgh'; 'abcdefgh'; 'abcdefgh'; 'abcdefgh'; 'abcdefgh'; 'abcdefgh'; 'abcdefgh'; 'abcdefgh'; 'abcdefgh'; 'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh';'abcdefgh'];
    %kheader =  ['        '; '        '; '        '; '        '; '        '; '        '; '        '; '        '; '        '; '        ';'        ';'        ';'        ';'        ';'        ';'        ';'        ';'        ';'        ';'        ';'        ';'        ';'        ';'        '];
    kheader =  ['-12345  '; '-12345  '; '-12345  '; '-12345  '; '-12345  '; '-12345  '; '-12345  '; '-12345  '; '-12345  '; '-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  ';'-12345  '];
    
    a=size(header);
    if (a(2) == 1)
    %    form header matrices from structured array
    
        %floats
        fheader(1) = header.delta;
        fheader(2) = min(data);
        fheader(3) = max(data);
        fheader(4) = header.scale;
        fheader(5) = header.odelta;
        fheader(6) = header.b;
        fheader(7) = header.e;
        fheader(8) = header.o;
        fheader(9) = header.a;
        fheader(11) = header.t0;
        fheader(12) = header.t1;
        fheader(13) = header.t2;
        fheader(14) = header.t3;
        fheader(15) = header.t4;
        fheader(16) = header.t5;
        fheader(17) = header.t6;
        fheader(18) = header.t7;
        fheader(19) = header.t8;
        fheader(20) = header.t9;
        fheader(22) = header.resp0;
        fheader(23) = header.resp1;
        fheader(24) = header.resp2;
        fheader(25) = header.resp3;
        fheader(26) = header.resp4;
        fheader(27) = header.resp5;
        fheader(28) = header.resp6;
        fheader(29) = header.resp7;
        fheader(30) = header.resp8;
        fheader(31) = header.resp9;
        fheader(32) = header.stla;
        fheader(33) = header.stlo;
        fheader(34) = header.stel;
        fheader(35) = header.stdp;
        fheader(36) = header.evla;
        fheader(37) = header.evlo;
        fheader(38) = header.evel;
        fheader(39) = header.evdp;
        fheader(41) = header.user0;
        fheader(42) = header.user1;
        fheader(43) = header.user2;
        fheader(44) = header.user3;
        fheader(45) = header.user4;
        fheader(46) = header.user5;
        fheader(47) = header.user6;
        fheader(48) = header.user7;
        fheader(49) = header.user8;
        fheader(50) = header.user9;
        fheader(51) = header.dist;
        fheader(52) = header.az;
        fheader(53) = header.baz;
        fheader(54) = header.gcarc;
        fheader(57) = header.depmen;
        fheader(58) = header.cmpaz;
        fheader(59) = header.cmpinc;
  
        %integers
        iheader(1) = header.nzyear;
        iheader(2) = header.nzjday;
        iheader(3) = header.nzhour;
        iheader(4) = header.nzmin;
        iheader(5) = header.nzsec;
        iheader(6) = header.nzmsec;
        iheader(7) = header.nvhdr;
        iheader(10) = header.npts;
        iheader(16) = header.iftype;
        iheader(17) = header.idep;
        iheader(18) = header.iztype;
        iheader(19) = header.iqual;
        iheader(20) = header.isynth;

        %longs
        lheader(1) = header.leven;
        lheader(2) = header.lpspol;
        lheader(3) = 1;
        lheader(4) = header.lcalda;

        %strings    
        for i=1:8
            header.kevnm1(i) = header.kevnm(i);
            header.kevnm2(i) = header.kevnm(i+8);
        end
        kheader = [header.kstnm; header.kevnm1; header.kevnm2; header.khole; header.ko; header.ka; header.kt0; header.kt1; header.kt2; header.kt3; header.kt4; header.kt5; header.kt6; header.kt7; header.kt8; header.kt9; header.kf; header.kuser0; header.kuser1; header.kuser2; header.kcmpnm; header.knetwk; header.kdatrd; header.kinst];
   
        
    else
        %floats
        fheader(1) = header(1);
        fheader(2) = min(data);
        fheader(3) = max(data);
        fheader(6) = header(4);
        fheader(7) = header(5);
        fheader(8) = header(6);
        fheader(11:20) = header(7:16);
        fheader(32:39) = header(17:24);
        fheader(51:54) = header(25:28);
        fheader(41) = header(41);
  
        
        %integers
        iheader(1:6) = header(29:34);
        iheader(10) = header(35);
        iheader(16:18) = header(36:38);
        iheader(23) = header(39);
        iheader(25) = header(40);
        

        %longs
        lheader(3) = 1;
        
    end
    
    %set the version number which I use to check endian
    iheader(7) = 6;
    
    fid = fopen(output_name,'w');
    if (fid < 0) 
        disp('could not write file')
        ierr = 0;
        return
    end
    
    
    c1 = fwrite(fid,fheader,'single');
    c2 = fwrite(fid,iheader,'int');
    c3 = fwrite(fid,lheader,'long');
    c4 = fwrite(fid,kheader','char');
    c5 = fwrite(fid,data,'single');
    ierr = c1 + c2 + c3 + c4 + c5;
    
    fclose(fid);
    
    
    return
    
    