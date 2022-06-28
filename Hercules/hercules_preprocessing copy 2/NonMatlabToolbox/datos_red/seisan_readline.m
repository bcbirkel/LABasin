function data=seisan_readline(fp,length)
        
    dataT=fread(fp,length+8,'char=>char');
    
    endL   = length + 4;
    startL = 5;
    data = dataT(startL:endL);

return
