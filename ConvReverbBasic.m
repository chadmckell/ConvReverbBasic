%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Convolution Reverberation (Equal Dimensions)
% 
% Author: Chad McKell
% Date: 06.03.17
%
% Description: This script convolves an input audio signal with a simulated 
% impulse response from a three-dimensional room with equal wall dimensions. 
% The convolution is performed using a vector-based implementation.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
tic; close; clc; clear;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Read in the input .WAV file 'x' and sample rate 'Fx'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
[filename, pathname]=uigetfile('*.wav', 'Select a wave file');
[x,Fx] = audioread(filename);

% If user inserts stereo audio, only read the left channel 
if size(x,2) > 1
    x = x(:,1); 
    warning('File x has two channels. Only left channel was read.');
end

% Determine length of 'x'
Lx = length(x);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Read in .wav file of simulated impulse response 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

[h,Fh] = audioread('IR_10x10x10_s1669455_McKell.wav'); 

% If user inserts stereo audio, only read the left channel 
if size(x,2) > 1
    h = h(:,1); 
    warning('File h has two channels. Only left channel was read.');
end

% Determine length of 'h'
Lh = length(h);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Convolve input audio signal with simulated impulse response 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Compute length of the output signal
N = Lh+Lx;

% Initialize the output signal
y = zeros(N,1);

% Zero-pad the beginning and end of the input signal 
x = [zeros(Lh,1); x]; 
x(end+Lh) = 0; 

% Flip the simulated impulse response vector
h = fliplr(h');

% Make 'x' a row vector and 'h' a column vector
x = x';
h = h';

% Use vector-based implementation to perform convolution
for n = 1:N-1
    y(n) = x(n+1:n+Lh)*h(1:Lh);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Normalize output signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
y = y/max(abs(y));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Listen to output signal 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soundsc(y,Fx);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Save simulation as .wav file
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define filename. Include parameter values in filename.
filename = 'ConvReverb_s1669455_McKell.wav';

% Write .wav file to MATLAB pwd at 16 bits
audiowrite(filename, y, Fx, 'BitsPerSample', 16);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Check code efficiency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
toc % print elapsed time

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% References
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% 1. Smith, S. W. "The scientist and engineer's guide to digital signal 
% processing," Chapter 6: Convolution (1997).

