clc; clear all; close all

%Tapered Beam Properties 
b1 = 0.2;      % Width of the beam at node 1 (m)
h1 = 0.2;      % Height of the beam at node 1 (m)
b2 = 0.1;      % Width of the beam at node 2 (m)
h2 = 0.1;      % Height of the beam at node 2 (m)
E = 200e9;    % Young's Modulus of the steel (Pa)
L = 0.8;        % Lenght of the beam (m)
I1 = b1*h1^3/12; % Area moment of inertia at node 1 (m^4)
I2 = b2*h2^3/12; % Area moment of inertia at node 2 (m^4)
A1 = b1*h1;      % Area of the beam at node 1 (m^2)
A2 = b2*h2;      % Area of the beam at node 2 (m^2)
G = 75e9;     % Shear Modulus of the steel (Pa)
J1 = (b1*h1/12)*(b1^2+h1^2);  % Polar area moment of inertia (m^4)
J2 = (b2*h2/12)*(b2^2+h2^2);  % Polar area moment of inertia (m^4)


% The Stifness matrix of the rod

Krod = (E*(A1+A2)/2*L) *[1    -1;
                        -1     1];    

% The Stifness matrix of the beam in xy plane

Kbeam_xy =[(6*E*(I1 + I2))/L^3     (2*E*(2*I1 + I2))/L^2    -(6*E*(I1 + I2))/L^3   (2*E*(I1 + 2*I2))/L^2;
           (2*E*(2*I1 + I2))/L^2       (E*(3*I1 + I2))/L  -(2*E*(2*I1 + I2))/L^2         (E*(I1 + I2))/L;
          -(6*E*(I1 + I2))/L^3    -(2*E*(2*I1 + I2))/L^2     (6*E*(I1 + I2))/L^3  -(2*E*(I1 + 2*I2))/L^2;
          (2*E*(I1 + 2*I2))/L^2         (E*(I1 + I2))/L  -(2*E*(I1 + 2*I2))/L^2       (E*(I1 + 3*I2))/L];
   
% The Stifness matrix of the beam in xz plane   
Kbeam_xz =[  (6*E*(I1 + I2))/L^3  -(2*E*(2*I1 + I2))/L^2   -(6*E*(I1 + I2))/L^3  -(2*E*(I1 - 4*I2))/L^2;
          -(2*E*(2*I1 + I2))/L^2       (E*(3*I1 + I2))/L  (2*E*(2*I1 + I2))/L^2       (E*(3*I1 - I2))/L;
            -(6*E*(I1 + I2))/L^3   (2*E*(2*I1 + I2))/L^2    (6*E*(I1 + I2))/L^3   (2*E*(I1 - 4*I2))/L^2;
          -(2*E*(I1 - 4*I2))/L^2       (E*(3*I1 - I2))/L  (2*E*(I1 - 4*I2))/L^2    (E*(9*I1 + 19*I2))/L]
                   
% The Torsional stifness matrix for tapered beam 
Kt = (G * (J1+J2)/2*L) * [1  -1;
                         -1   1];
    
% 3D Stifness matrix for tapered beam             
Ke= zeros(12,12); 
Ke([1,7],[1,7]) = Krod;
Ke([4,10],[4,10]) = Kt;
Ke([2,6,8,12],[2,6,8,12]) = Kbeam_xy;
Ke([3,5,9,11],[3,5,9,11]) = Kbeam_xz               
 
          


               
  