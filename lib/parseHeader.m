function [headers, nHeaderLines] = parseHeader(fpath)
fidP = fopen(fpath);
nextLine = fgets(fidP);
nHeaderLines = 0;
while ischar(nextLine) && contains(nextLine, "=")
    if startsWith(nextLine, "VOX_DIM=")
        VOX_DIM=sscanf(nextLine, 'VOX_DIM=%d,%d,%d'); 
        headers.VOX_DIM=VOX_DIM;
    end
    if startsWith(nextLine, "VOX_SIDELEN=")
        headers.VOX_SIDELEN=sscanf(nextLine, 'VOX_SIDELEN=%f'); 
    end
    if startsWith(nextLine, "VOX_ZLEN=")
        headers.VOX_ZLEN=sscanf(nextLine, 'VOX_ZLEN=%f'); 
    end
    if startsWith(nextLine, "NUM_ABS=")
        headers.NUM_ABS = sscanf(nextLine, 'NUM_ABS=%d');
        arrFmt = "%d";
        for i = 1:headers.NUM_ABS-1
            arrFmt = append(arrFmt, ",%d");
        end
    end
    if startsWith(nextLine, "ABSVOX_ROW")
        absRowFmt = append("ABSVOX_ROW=", arrFmt);
        headers.ABSVOX_ROW = sscanf(nextLine, absRowFmt);
    end
    if startsWith(nextLine, "ABSVOX_COL")
        absColFmt = append("ABSVOX_COL=", arrFmt);
        headers.ABSVOX_COL = sscanf(nextLine, absColFmt);
    end
    if startsWith(nextLine, "ABSVOX_Z")
        absZFmt = append("ABSVOX_Z=", arrFmt);
        headers.ABSVOX_Z = sscanf(nextLine, absZFmt);
    end
    
    if startsWith(nextLine, "SRC_DIM=")
        SRC_DIM=sscanf(nextLine, 'SRC_DIM=%d,%d'); 
        headers.SRC_DIM=SRC_DIM;
    end
    if startsWith(nextLine, "SRC_ORIG=")
        SRC_ORIG=sscanf(nextLine, 'SRC_ORIG=%f,%f,%f'); 
        headers.SRC_ORIG=SRC_ORIG;
    end
    if startsWith(nextLine, "SRC_SEP=")
        headers.SRC_SEP=sscanf(nextLine, 'SRC_SEP=%f'); 
    end
    if startsWith(nextLine, "SENS_DIM=")
        SENS_DIM=sscanf(nextLine, 'SENS_DIM=%d,%d'); 
        headers.SENS_DIM=SENS_DIM;
    end
    if startsWith(nextLine, "SENS_ORIG=")
        SENS_ORIG=sscanf(nextLine, 'SENS_ORIG=%f,%f,%f'); 
        headers.SENS_ORIG=SENS_ORIG;
    end
    if startsWith(nextLine, "SENS_PIXW=")
        headers.SENS_PIXW=sscanf(nextLine, 'SENS_PIXW=%f'); 
    end
    if startsWith(nextLine, "SENS_SEP=")
        headers.SENS_SEP=sscanf(nextLine, 'SENS_SEP=%f'); 
    end
    if startsWith(nextLine, "VOX_DIM=")
        VOX_DIM=sscanf(nextLine, 'VOX_DIM=%d,%d,%d');
        headers.VOX_DIM = VOX_DIM;
    end
    if startsWith(nextLine, "VOX_ORIG=")
        VOX_ORIG=sscanf(nextLine, 'VOX_ORIG=%f,%f,%f'); 
        headers.VOX_ORIG=VOX_ORIG;
    end
    if startsWith(nextLine, "NUM_BINS=")
        headers.NBINS=sscanf(nextLine, 'NUM_BINS=%d');
    end
    if startsWith(nextLine, "TIME_MIN")
        headers.TIME_MIN = sscanf(nextLine, "TIME_MIN=%f");
    end
    if startsWith(nextLine, "TIME_MAX")
        headers.TIME_MAX = sscanf(nextLine, "TIME_MAX=%f");
    end
    if startsWith(nextLine, "ABSVOX_UA")
        headers.ABSVOX_UA = sscanf(nextLine, "ABSVOX_UA=%f");
    end
    if startsWith(nextLine, "NUM_SAMPS")
        headers.NUM_SAMPS = sscanf(nextLine, "NUM_SAMPS=%d");
    end
    if startsWith(nextLine, "NUM_ITERS")
        headers.NUM_ITERS = sscanf(nextLine, "NUM_ITERS=%d");
    end
    if (startsWith(nextLine, "NORMALIZETPSF"))
        headers.NORMALIZETPSF = sscanf(nextLine, "NORMALIZETPSF=%d");
    end
    
    if startsWith(nextLine, "NUMLAYERS")
        configs.NUMLAYERS = sscanf(nextLine, "#NUMLAYERS=%d");
        layerFmt = "{%f";
        for i = 1:configs.NUMLAYERS-1
            layerFmt = append(layerFmt, ",%f");
        end
        layerFmt = append(layerFmt, "}");
    end
    if startsWith(nextLine, "SCAT_VEC")
        scatVecFmt = append("SCAT_VEC=", layerFmt);
        headers.SCAT_VEC = sscanf(nextLine, scatVecFmt);
    end
    
    nextLine = fgets(fidP);
    nHeaderLines = nHeaderLines + 1;
end

headers.VOX_L = VOX_DIM(1); headers.VOX_W = VOX_DIM(2); headers.VOX_H = VOX_DIM(3);
headers.SRC_L = SRC_DIM(1); headers.SRC_W = SRC_DIM(2);
headers.SRC_ORIGX = SRC_ORIG(1); headers.SRC_ORIGY = SRC_ORIG(2); headers.SRC_ORIGZ = SRC_ORIG(3);
headers.SENS_L = SENS_DIM(1); headers.SENS_W = SENS_DIM(2);
headers.SENS_ORIGX = SENS_ORIG(1); headers.SENS_ORIGY = SENS_ORIG(2); headers.SENS_ORIGZ = SENS_ORIG(3);
headers.VOX_ORIGX = VOX_ORIG(1); headers.VOX_ORIGY = VOX_ORIG(2); headers.VOX_ORIGZ = VOX_ORIG(3);
fclose all;

end