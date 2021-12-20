% This function consists of objective function f and constraint
% functions gj, not necessarily convex;
% There are 25 constraints in the Minimization problem,
% P: min f ,s.t 0 <= gj <= 1;
% The number of variabls is 20 and the highest degree is 2;

function [F,G,I,J,d,k] = test20_2

d=1;
k=1;

F=[ 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
    0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1;
    1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; 
    0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1;        
    0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1;
    0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 -1];  %f=x1^2-x2^2+x3^2-x4^2+x5^2-x6^2+x7^2-x8^2+x9^2-x10^2+x11^2-x12^2+x13^2-x14^2+x15^2-x16^2+x17^2-x18^2+x19^2-x20^2+x1-x2;
 
 g1=[2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 2]; %g1=2*x1^2+3*x2^2+2*x1*x2+2*x3^2+3*x4^2+2*x3*x4+2*x5^2+3*x6^2+2*x5*x6+2*x7^2+3*x8^2+2*x7*x8+2*x9^2+3*x10^2+2*x9*x10+2*x11^2+3*x12^2+2*x11*x12+2*x13^2+3*x14^2+2*x13*x14+2*x15^2+3*x16^2+2*x15*x16+2*x17^2+3*x18^2+2*x17*x18+2*x19^2+3*x20^2+2*x19*x20;
 
 g2=[2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 -4]; %g2=3*x1^2+2*x2^2-4*x1*x2+3*x3^2+2*x4^2-4*x3*x4+3*x5^2+2*x6^2-4*x5*x6+3*x7^2+2*x8^2-4*x7*x8+3*x9^2+2*x10^2-4*x9*x10+3*x11^2+2*x12^2-4*x11*x12+3*x13^2+2*x14^2-4*x13*x14+3*x15^2+2*x16^2-4*x15*x16+3*x17^2+2*x18^2-4*x17*x18+3*x19^2+2*x20^2-4*x19*x20;
 
 g3=[2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6;
     1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6;
     0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6;
     0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 6;
     0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 6;
     0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 6;
     0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 6;
     0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 6;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 6;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 -4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 6;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 -4];  %g3=x1^2+6*x2^2-4*x1*x2+x3^2+6*x4^2-4*x3*x4+x5^2+6*x6^2-4*x5*x6+x7^2+6*x8^2-4*x7*x8+x9^2+6*x10^2-4*x9*x10+x11^2+6*x12^2-4*x11*x12+x13^2+6*x14^2-4*x13*x14+x15^2+6*x16^2-4*x15*x16+x17^2+6*x18^2-4*x17*x18+x19^2+6*x20^2-4*x19*x20;

 g4=[2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4;
     1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3;
     0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4;
     0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3;
     0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4;
     0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3;
     0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 4;
     0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 -3;
     0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 4;
     0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 -3;
     0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 4;
     0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 -3;
     0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 4;
     0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 -3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 -3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 -3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 4;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 -3];  %g4=x1^2+4*x2^2-3*x1*x2+x3^2+4*x4^2-3*x3*x4+x5^2+4*x6^2-3*x5*x6+x7^2+4*x8^2-3*x7*x8+x9^2+4*x10^2-3*x9*x10+x11^2+4*x12^2-3*x11*x12+x13^2+4*x14^2-3*x13*x14+x15^2+4*x16^2-3*x15*x16+x17^2+4*x18^2-3*x17*x18+x19^2+4*x20^2-3*x19*x20;
 
 g5=[2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5;
     1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5;
     0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5;
     0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 5;
     0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 5;
     0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 5;
     0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 5;
     0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 5;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 5;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 3;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 2;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 5;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 3]; %g5=2*x1^2+5*x2^2+3*x1*x2+2*x3^2+5*x4^2+3*x3*x4+2*x5^2+5*x6^2+3*x5*x6+2*x7^2+5*x8^2+3*x7*x8+2*x9^2+5*x10^2+3*x9*x10+2*x11^2+5*x12^2+3*x11*x12+2*x13^2+5*x14^2+3*x13*x14+2*x15^2+5*x16^2+3*x15*x16+2*x17^2+5*x18^2+3*x17*x18+2*x19^2+5*x20^2+3*x19*x20;
 
 g6 =[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g7 =[0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g8 =[0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g9 =[0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g10=[0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g11=[0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g12=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g13=[0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1];
 g14=[0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1];
 g15=[0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1];
 g16=[0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1];
 g17=[0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1];
 g18=[0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1];
 g19=[0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1];
 g20=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1];
 g21=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1];
 g22=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1];
 g23=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1];
 g24=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];
 g25=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1];

 G={g1;g2;g3;g4;g5;g6;g7;g8;g9;g10;g11;g12;g13;g14;g15;g16;g17;g18;g19;g20;g21;g22;g23;g24;g25};  
 I={[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]};
 J={[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]};
