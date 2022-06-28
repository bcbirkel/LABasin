%    
%  [maxdisp, maxvel, maxacc]= max_osc_response
%               
%                    Computes the maximum response of an oscilator
%                    given the accelration, period and damping.
%
%   INPUT:
%
%         acceleration -- input signal
%                   dt -- time step
%                  csi -- damping ratio
%               period -- oscilator's period
%          initialdisp -- initial displacement 
%           initialvel -- initial velocity
%
%
%   OUTPUT:maximum values of displacment, velocity, acceleration
%                    
%                     maxdisp 
%                     maxvel
%                     maxacc
%
%
%   Copyright (c) 2006 Leonardo Ramirez-Guzman
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 2
%   of the License, or any later version. 
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 
%   51 Franklin Street,
%   Fifth Floor, Boston, MA 
%   02110-1301, USA.
%
%   Contact:
%   Leonardo Ramirez-Guzman 
%   Civil and Environmental Engineering
%   Carnegie Mellon University
%   5000 Forbes Avenue
%   Pittsburgh, PA 15213
%   lramirez@andrew.cmu.edu

function [maxdisp, maxvel, maxacc] = max_osc_response_2011( acceleration, dt, csi, ...
                                     period,initialdisp,initialvel)     

    signalsize=size(acceleration);
    np=length(acceleration);
    
    d(1)=initialdisp;
    v(1)=initialvel;
    
    w=2*pi/period;
	ww=w^2;
    csicsi=csi^2;
	dcsiw=2.*csi*w;	
    
    rcsi=sqrt(1-csicsi);
    csircs=csi/rcsi; 
    wd=w*rcsi;
    ueskdt=-1./(ww*dt);
    dcsiew=2.*csi/w;
    um2csi=(1.-2*csicsi)/wd;
    e=exp(-w*dt*csi);
    s=sin(wd*dt);
    c0=cos(wd*dt);
    aa(1)=-ww*d(1)-dcsiw*v(1);
    
    ca=e*(csircs*s+c0);
    cb=e*s/wd;
    cc=(e*((um2csi-csircs*dt)*s-(dcsiew+dt)*c0)+dcsiew)*ueskdt;
    cd=(e*(-um2csi*s+dcsiew*c0)+dt-dcsiew)*ueskdt;
    cap=-cb*ww;
    cbp=e*(c0-csircs*s);
    ccp=(e*((w*dt/rcsi+csircs)*s+c0)-1.)*ueskdt;
    cdp=(1.-ca)*ueskdt;

    for  i=2:np
        d(i) =ca*d(i-1)+cb*v(i-1)+cc*acceleration(i-1)+cd*acceleration(i);
        v(i) =cap*d(i-1)+cbp*v(i-1)+ccp*acceleration(i-1)+cdp*acceleration(i);
       
    end
    aa=-ww*d-dcsiw*v;
    maxdisp =max(abs(d));
    maxvel  =max(abs(v));
    maxacc  =max(abs(aa));           

    return; 
