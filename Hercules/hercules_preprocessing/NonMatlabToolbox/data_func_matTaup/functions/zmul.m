% function zmul.m
% replicates zmul function in evalresp source code file calc_functions.c
% essentially a simple routine to multiply complex values
%%  ORIGINAL C CODE
%/*==================================================================
%  *    Complex multiplication:  complex version of val1 *= val2;
%  *=================================================================*/
% void zmul(struct complex *val1, struct complex *val2) {
%     double r, i;
%     r = val1->real*val2->real - val1->imag*val2->imag;
%     i = val1->imag*val2->real + val1->real*val2->imag;
%     val1->real = r;
%     val1->imag = i;
% }
%% USAGE
% [complex_out] = zmul(complex_1, complex_2)
%% STATMENT OF REDUDANCY
% This function is redundant as matlab happily handles complex values
%%

function [complex_out] = zmul(complex_1, complex_2)
    complex_out = (real(complex_1) .* real(complex_2) - imag(complex_1) .* imag(complex_2)) + (imag(complex_1) .* real(complex_2) + real(complex_1) .* imag(complex_2))*i;
%%

 return

