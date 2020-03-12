% This function MSFM3D calculates the shortest distance from a list of
% points to all other pixels in an image volume, using the  
% Multistencil Fast Marching Method (MSFM). This method gives more accurate 
% distances by using second order derivatives and cross neighbours.
% 
% T=msfm3d(F, SourcePoints, UseSecond, UseCross)
%
% inputs,
%   F: The 3D speed image. The speed function must always be larger
%			than zero (min value 1e-8), otherwise some regions will
%			never be reached because the time will go to infinity. 
%   SourcePoints : A list of starting points [3 x N] (distance zero)
%   UseSecond : Boolean Set to true if not only first but also second 
%                order derivatives are used (default)
%   UseCross : Boolean Set to true if also cross neighbours 
%                are used (default)
% outputs,
%   T : Image with distance from SourcePoints to all pixels
%
% Note:
%   First you need to compile the c file "mex msfm3d.c"
%
% Literature : M. Sabry Hassouna et Al. Multistencils Fast Marching 
%   Methods: A Highly Accurate Solution to the Eikonal Equation on
%   Cartesian Domains
%
% Example,
%   % compile the needed mex file
%   mex msfm3d.c
%
%   SourcePoint = [21; 21; 21];
%   SpeedImage = ones([41 41 41]);
%   [X,Y,Z] = ndgrid(1:41, 1:41, 1:41);
%   T1 = sqrt((X-SourcePoint(1)).^2 + (Y-SourcePoint(2)).^2 + (Z-SourcePoint(3)).^2);
%
%   % Run fast marching 1th order, 1th order multi stencil 
%   % and 2th orde and 2th orde multi stencil
%
%   tic; T1_FMM1 = msfm3d(SpeedImage, SourcePoint, false, false); toc;
%   tic; T1_MSFM1 = msfm3d(SpeedImage, SourcePoint, false, true); toc;
%   tic; T1_FMM2 = msfm3d(SpeedImage, SourcePoint, true, false); toc;
%   tic; T1_MSFM2 = msfm3d(SpeedImage, SourcePoint, true, true); toc;
%
%   % Show results
%   fprintf('\nResults with T1 (Matlab)\n');
%   fprintf('Method   L1        L2        Linf\n');
%   Results = cellfun(@(x)([mean(abs(T1(:)-x(:))) mean((T1(:)-x(:)).^2) max(abs(T1(:)-x(:)))]), {T1_FMM1(:) T1_MSFM1(:) T1_FMM2(:) T1_MSFM2(:)}, 'UniformOutput',false);
%   fprintf('FMM1:   %9.5f %9.5f %9.5f\n', Results{1}(1), Results{1}(2), Results{1}(3));
%   fprintf('MSFM1:  %9.5f %9.5f %9.5f\n', Results{2}(1), Results{2}(2), Results{2}(3));
%   fprintf('FMM2:   %9.5f %9.5f %9.5f\n', Results{3}(1), Results{3}(2), Results{3}(3));
%   fprintf('MSFM2:  %9.5f %9.5f %9.5f\n', Results{4}(1), Results{4}(2), Results{4}(3));

