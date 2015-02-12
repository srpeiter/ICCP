%% Argon gas simulation
% This file is for trying out ways to model an argon gas
% Feed back for testing methods is faster with matlab than with C++ For
% now


function main
%% Parameters
NR = 3;                 % Define the number of particles per line (Total number = NR^3)
U = 10;                 % Dimensions of the two dimensional universe (U by U)
C1 = 10;                % repelling force constant
C2 = 1;                 % attractive force constant
t = 0.02;                % discrete time unit
E0 = 0.001;               % energy of the initial system
R = 1;              % Number of iterations

%% Pre-calculations
N = NR^3;               % Calculate the number of particles in the system

%% File
P = InitiatePosition(N,NR,U);           % Initiate [N] particles position [P] that are in lanes of [NR] in universe size [U]
V = InitiateVelocity(N,E0);             % Initiate velocity [V] of [N] particles with energy [E0]
A = zeros(N,3);
for i = 1:R
    P1 = UpdatePosition(P,V,t,U);             % Updated positions [P1] according to the velocities [V] after time period [t]
    for j = 1:N
        [d,p] = Distance(j,P,N);            % calculate the distance [d] between particle [j] and all others. [p] is the vector between two points [j] and others
        u = UnitVector(p,d,N);              % Unit vecotr [u] is nessesary to give directionality to the acceleration. it is determined using the distant veotors [p],the distance [d] between the particles and the number of particles [N] 
        A(j,:) = Acceleration(u,d,N,C1,C2); % Calculate the accelartion [A] for for an element, [u] for directionality, [d],[C1] and [C2] are used to deterimine the magnitute of the force 
    end
    V = V + A * t;                          % Update the velocities [V] according to the accelaration [A] applied during time [t]
    P = P1;                                 % The new postitions [P] are the updated positions [P1]
    Ek = KineticEnergy(V);                  % Calculate the total kinetic energy [Ek] 
    
    %% Plotting the inforation
    subplot(1,2,1)
    scatter3(P(:,1),P(:,2),P(:,3))
    axis([0 U 0 U 0 U])
    pause(0.001)
    bla(i) = Ek;
    subplot(1,2,2)
    plot(bla)
    
%     %% saving the information
%     save('temp','V','P')
    
end

end

function P = InitiatePosition(N,NR,U)
%% Initial particle position
% FCC structure is het best for high density packing
% This function is meant to create such a structure
% Inputs:   N (Number of particles)
%           NR: Number of particle row
%           U (Universe size)
% Outputs:  P (N by 2 matrix that uncludes the x and y position, size of U by U)

% Number of particles per row

% Prealocating P
P = zeros(N,3);                     %
k = 0;
for i = 0:NR-1
    for j = 0:NR-1
        for h = 0:NR-1
            k = k + 1;
            x = 3 * i +  j +  h;        % x position of the particle
            y = 3 * j + h;                    % y position of the particle
            z = 3 * h ;                    % z position of the particle
            
            l1 = 3*NR * floor(x/(3*NR));       % if the x position of the particle is outside the box, it will start at the other side again
            l2 = 3*NR * floor(y/(3*NR));
            
            P(k,:) = [(x-l1)/3; (y-l2)/3;z/3];
        end
    end
end
% convert the total space units to a square of size one.
P = P / NR * U;                   % divided by particles per row, multiply with the size of the universe


figure(1)
% plot(P(:,1),P(:,2),'.')
scatter3(P(:,1),P(:,2),P(:,3))
end

function V = InitiateVelocity(N,E0)
%% Initial particle velocity
% In this file I give the particles a initial velocity
% Gaussian distribution used
% Note: the implementation of energy is shit atm
% Inputs:   N (Number of particles)
%           E0 Energy of the begin system
% Outputs:  V (N by 2 matrix giving x and y velocities, gaussian
%              distributed)

V = randn(N,3)*E0;

% figure(1)
% plot(V(:,1),V(:,2),'.')
end

function [d,p] = Distance(i,P,N)
%% Distance between element i and the rest
% This is a very simple script that calculates the distance between element
% i and all the other elements
% Note: the way to get rid of interaction with the same point is probably
% not optimal.
% Inputs:   i
%           P
%           N
% Outputs:  d (distance between particles) (N vector)
%           p (distance vector between particle i and every other) (N by 3 matrix)

p = zeros(N,3);
p(:,1) = P(:,1) - P(i,1);
p(:,2) = P(:,2) - P(i,2);
p(:,3) = P(:,3) - P(i,3);
p(i,:) = 100;               % cheap way to get rid of effect of the same element
d = (p(:,1).*p(:,1)+p(:,2).*p(:,2)+p(:,3).*p(:,3)).^(1/2);
end

function u = UnitVector(p,d,N)
%% Unit vector
% The directionality of the force is of course also crusial
% Inputs:   p (distance vector between element i and every other element) (N by 3 matrix)
%           d (distance between elements)
%           N (number of particles)
% Outputs:  u (unitary vectors) (N by 3)
u = zeros(N,3);
u(:,1) = p(:,1)./d;
u(:,2) = p(:,2)./d;
u(:,3) = p(:,3)./d;
end

function A = Acceleration(u,d,N,C1,C2)
%% Acceleration
% The acceleration can now be determined
% Note: I did not pay much attention in getting the force right
% Inputs:   u: unitary vectors) (N by 3)
%           d: distance between elements)
%           N: number of particles)
%           C1: repelling force constant
%           C2: atractive force constant
% Outputs:  a: accelerations in all 3 dirrections) (N by 3 matrix)
%

r6 = d.^6;
f = (C1./r6.^2 - C2./r6).*d;       % force equation, probably wrong

a = zeros(N,3);
a(:,1) = u(:,1).*f;
a(:,2) = u(:,2).*f;
a(:,3) = u(:,3).*f;

% Total accelation
A = -sum(a);
end

function Ek = KineticEnergy(V)
%% In this file I calculate the kinetic energy because why not right?
% Inputs:   V: Speeds of all the particles in the three directions of space
% Outputs:  Ek: Total kinetic energy
Ek = sum(0.5 * (sqrt(V(:,1).^2+V(:,2).^2+V(:,3).^2)));

end

function P1 = UpdatePosition(P,V,t,U)
%% The position gets updated according to the speed V and the time distance t
% Inputs:   P: Old positions
%           V: speeds of the particles
%           t: time distance of the particles
%           U: Size universe
% Outputs:  P1: New positions
          
P1 = P + V * t;
% P1 = P1 - U * floor(P1./U);         % respawn a point at the other side
% after leaving the universe kinda works, but mainly not. Maybe due to the
% lack of repeated boundary conditions

end

% function Potential energy
% % 
% 
% 
% end

