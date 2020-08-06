function [M,S] = DistGSDF(Src)
% Estimate GSDF Distribution

nm = 1;
for ns = 1:size(Src,2) %Source number
    for nr = 1:size(Src(ns).Rec,2)
        for np = 1:size(Src(ns).Rec(nr).GSDF,1)
            for nf = 1:size(Src(ns).Rec(nr).GSDF,2)
                Meas(nm) = Src(ns).Rec(nr).GSDF(np,nf,2);
            end
            nm = nm+1;
        end
    end
end

M = mean(Meas);
S = std(Meas);



