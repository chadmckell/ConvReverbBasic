%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% PBMMI Assignment 2 - Image Source Method (Equal Dimensions)
% 
% Author: Chad McKell
% Date: 06.03.17
%
% Description: This script computes a simulated impulse response from a
% three-dimensional room with flat walls positioned perpendicular to each
% other.  The dimensions of the room (i.e. Lx, Ly, and Lz) are equal.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
tic; close; clc; clear;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

Fs = 44100; % sample rate [samples/sec]
L = 10; % length of room along each dimension [meter]
alpha = 0.9; % absorption coefficient 
beta = sqrt(1-alpha); % reflection coefficient 
T60 = 0.161*L/(6*alpha); % reverberation time [sec]
v = 343; % speed of sound [meter/sec]
N = ceil(v*T60/L); % maximum integer image source number 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define source location S = (p, q, r)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

p = 2; % x displacement [meter]
q = 7; % y displacement [meter]
r = 1; % z displacement [meter]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define receiver location R = (a, b, c)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

a = 2; % x displacement [meter]
b = 4; % y displacement [meter]
c = 5; % z displacement [meter]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Initialize the impulse response vector 'y'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

Q = (N+1)*L; % maximum dimensional length between R and virtual source
Lmax = sqrt(3*Q*Q); % maximum total distance between R and virtual source
tmax = ceil(Fs*Lmax/v); % maximum arrival time [samples] 
y = zeros(tmax,1); % initialized impulse response vector

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Error handling: Terminate code for the following errors
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Fs or alpha is not a real number
if ~isfloat(Fs) || ~isfloat(alpha) || ~isreal(Fs) || ~isreal(alpha)  
   error('Fs and alpha must each be a real number.')
end

% L is not a real number
if ~isfloat(L) || ~isreal(L)
   error('L must be a real number.')
end

% Fs has a non-zero decimal
if rem(Fs,1) ~= 0
   error('Fs must be zero-decimal.')
end

% Fs, L, or alpha is negative
if Fs < 0 || L < 0 || alpha < 0
   error('Fs, L, and alpha must each be positive.')
end

% alpha is out of bounds
if alpha > 1 
   error('alpha cannot be greater than 1.')
end

% p or a is out of bounds
if p > L || q > L || r > L || a > L || b > L || c > L
   error('p, q, r, a, b, and c must each be less than L.')
end

% p, q, r, a, b, or c is negative
if p < 0 || q < 0 || r < 0 || a < 0 || b < 0 || c < 0 
   error('p, q, r, a, b, and c must each be positive.')
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Compute impulse response vector 'y'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

for d = -N:N
    for e = -N:N
        for f = -N:N
            
            % Compute distance of virtual source from R along x-dimension 
            if rem(d,2) == 0
                A = abs(d)*L + p - a; 
            else
                A = (abs(d)+1)*L - p - a;
            end
            
            % Compute distance of virtual source from R along y-dimension    
            if rem(e,2) == 0
                B = abs(e)*L + q - b;
            else
                B = (abs(e)+1)*L - q - b;
            end
            
            % Compute distance of virtual source from R along z-dimension     
            if rem(f,2) == 0
                C = abs(f)*L + r - c;
            else
                C = (abs(f)+1)*L - r - c;
            end
            
            % Find the total distance between receiver and virtual source
            l = sqrt(A*A + B*B + C*C);
            
            % Find the number of collisions between walls and sound wave
            w = abs(d+e+f);
            
            % Calculate the magnitude of the impulse
            g = beta^w/l;
            
            % Compute the time of arrival 
            t = l/v;
            
            % Find the time step 
            n = floor(t*Fs);
            
            % Calculate the total impulse reponse at time step 'n'
            y(n) = g + y(n);
        end
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Normalize output signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
y = y/max(abs(y));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Listen to output signal 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soundsc(y,Fs);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Save simulation as .wav file
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define filename. Include parameter values in filename.
filename = ['IR_' num2str(L) 'x' num2str(L) 'x' num2str(L)...
    '_s1669455_McKell.wav'];

% Write .wav file to MATLAB pwd at 16 bits
audiowrite(filename, y, Fs, 'BitsPerSample', 16);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define plotting parameters 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

n = 0:tmax-1; % length bins [samples]
tbin = (n/Fs)'; % time bins [sec]
font = 14; % font size for plots

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Generate plot of impulse response waveform 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Input signal  
figure(1)
plot(tbin, y); 
title(['Impulse Response Waveform (\alpha='  num2str(alpha) ')']);
xlabel('Time (sec)'); 
ylabel('Normalized Magnitude');
set(gca,'fontsize',font)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Check code efficiency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
toc % print elapsed time

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% References
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% 1. Allen, J. B., Berkley, D.A. "Image method for efficiently simulating
%    small?room acoustics." The Journal of the Acoustical Society of
%    America, 65.4 (1979): 943-950.
% 2. http://www.umiacs.umd.edu/~ramani/cmsc828d_audio/828d_l20.pdf

