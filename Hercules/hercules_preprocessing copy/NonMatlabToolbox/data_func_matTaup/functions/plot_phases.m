% plot_phases.m
% function to plot expected arrival times of various phases.
% uses the taup toolkit to calculated expected times.
% usage:
% plot_phases(
%     time_array, 
%     data_array,
%     origin_time_relative_to_start_of_time_array,
%     model - 'iasp91', 'prem', or 'qdt'
%     depth
%     phase (array of strings - eg 'P,S,PmP'
%     'deg' or 'km'
%     distance in degrees or kilometers as specified above
%       )
%
% if unsure of the distance, use the distance([lat lon], [lat lon])
% function or the degreedistance([lat1,long1;lat2,long2]) function
%
% requires the installation of the matTaup library. To see if it is working
% type: help taupTimes
% if you get help information on taupTimes, this will work
%
% addition of phase 'all'
% allows getting all possible phases
function plot_phases(time, data, origin, model, depth, phase, distance_flag, distance)
    % check if the user wants all phases
    if strcmp(phase,'all')
        phase = 'p,s,P,S,pP,sS,Pn,Sn,PcP,ScS,Pdiff,Sdiff,PKP,SKS,PKiKP,SKiKS,PKIKP,SKIKS';
    end

    % call tauP
    tt = taupTime(model, depth, phase, distance_flag, distance);
    
    % get the number of phases returned
    tmp = size(tt);
    nphases = tmp(2);
    
    % get the phase arrivals in an array
    % Because this can have 0 to 18 entries and it returns in an annoying
    % way, we have to use a trick to get it into an array
    if (nphases == 0)
        disp('No phases returned')
        return
    elseif (nphases == 1)
        [pp] = tt.time;
        [aaa] = tt.phaseName;
        pn=char(aaa);
    elseif (nphases == 2)
        [pp(1),pp(2)] = tt.time;
        [aaa, bbb] = tt.phaseName;
        pn=char(aaa,bbb);
    elseif (nphases == 3)
        [pp(1),pp(2),pp(3)] = tt.time;
        [aaa,bbb,ccc] = tt.phaseName;
        pn=char(aaa,bbb,ccc);
    elseif (nphases == 4)
        [pp(1), pp(2), pp(3), pp(4)] = tt.time;
        [aaa, bbb, ccc, ddd] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd);
    elseif (nphases == 5)
        [pp(1), pp(2), pp(3), pp(4), pp(5)] = tt.time;
        [aaa, bbb, ccc, ddd, eee] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee);
    elseif (nphases == 6)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff);
    elseif (nphases == 7)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg);
    elseif (nphases == 8)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh);
    elseif (nphases == 9)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii);
    elseif (nphases == 10)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj);
    elseif (nphases == 11)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk);
    elseif (nphases == 12)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11), pp(12)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll);
    elseif (nphases == 13)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11), pp(12),pp(13)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm);
    elseif (nphases == 14)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11), pp(12),pp(13),pp(14)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn] = tt.phaseName;
         pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn);
    elseif (nphases == 15)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11), pp(12),pp(13),pp(14),pp(15)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo] = tt.phaseName;
         pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo);
    elseif (nphases == 16)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11), pp(12),pp(13),pp(14),pp(15),pp(16)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo,ppp] = tt.phaseName;
         pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo,ppp);
    elseif (nphases == 17)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11), pp(12),pp(13),pp(14),pp(15),pp(16),pp(17)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo,ppp,qqq] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo,ppp,qqq);
    elseif (nphases == 18)
        [pp(1),pp(2),pp(3),pp(4),pp(5),pp(6), pp(7), pp(8),pp(9),pp(10),pp(11), pp(12),pp(13),pp(14),pp(15),pp(16),pp(17),pp(18)] = tt.time;
        [aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo,ppp,qqq,rrr] = tt.phaseName;
        pn=char(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,iii,jjj,kkk,lll,mmm,nnn,ooo,ppp,qqq,rrr);
    end
    
    
    % plot the data
    plot(time,data);
    
    % keep hold on to overlay figures
    hold on
    
    % format the axis using information about the origin and maximum phases
    a = min(origin-100,min(time));
    b = max(pp(nphases)+100,max(time));
    time_axis = [a,b];
    xlim(time_axis);
    
    % add the origin
    plot(origin,0,'bd');
    text(origin,min(data)/5,'O');
    
    % add each of the phases
    for (j=1:nphases)
        if mod(j,2) == 1
            flip_flag = -1;
        else
            flip_flag = 1;
        end
        offset = mod(j,4)+1;
        plot(pp(j)+origin,0,'rd');
        text(pp(j)+origin,flip_flag * offset * min(data)/5,pn(j,:))
    end

    hold off
    
    return
