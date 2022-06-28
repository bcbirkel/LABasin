function [str, rake, dip]=fpsol(nu,u)

dip=acos(-1*nu(3));
  if(nu(1)==0. & nu(2)==0.)
   str=0.;
  else
   str=atan2(-1*nu(1),nu(2));
  end
  

sstr=sin(str);
cstr=cos(str);
sdip=sin(dip);
cdip=cos(dip);


  if(abs(sdip) > 0.)
    rake=asin(-1*u(3)/sin(dip));
  else
    arg1=1.;
    arg2=u(3);
    arg=sign(arg2)
        if(arg < 0.)
          rake=pi;
        else
          rake=0.;
      end
  end
  
  slambda=sin(rake);
  cdsl=cdip*slambda;
  
  if(abs(sstr) > abs(cstr))
        clambda=(u(2)+cdsl*cstr)/sstr;
  else
        clambda=(u(1)-cdsl*sstr)/cstr;
    end
      if(slambda == 0. & clambda == 0.)
        slip=0.;
      else
        slip=atan2(slambda,clambda)+pi;%pi added
    end
      
  if(dip > pi/2)
        dip=pi-dip;
        str=str+pi;
        slip=2*pi-slip;
    end
  
    if(str < 0.) 
        str=str+2*pi;
    end
    if(slip>= pi) 
        slip=slip-2*pi;
    end
    str=str*180/pi;
    rake=slip*180/pi;
    dip=dip*180/pi;