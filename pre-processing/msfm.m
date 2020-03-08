function [T,Y]=msfm(F, SourcePoints, UseSecond, UseCross)
% This function MSFM calculates the shortest distance from a list of
% points to all other pixels in an image volume, using the  
% Multistencil Fast Marching Method (MSFM). This method gives more accurate 
% distances by using second order derivatives and cross neighbours.
% 
%   [T,Y]=msfm(F, SourcePoints, UseSecond, UseCross)
%
% inputs,
%   F: The 2D or 3D speed image. The speed function must always be larger
%			than zero (min value 1e-8), otherwise some regions will
%			never be reached because the time will go to infinity. 
%   SourcePoints : A list of starting points [2 x N ] or  [3 x N] (distance zero)
%   UseSecond : Boolean Set to true if not only first but also second 
%                order derivatives are used (default)
%   UseCross : Boolean Set to true if also cross neighbours 
%                are used (default)
% outputs,
%   T : Image with distance from SourcePoints to all pixels
%   Y : Image for augmented fastmarching with, euclidian distance from 
%       SourcePoints to all pixels. (Used by skeletonize method)
%
% Note:
%   Run compile_c_files.m to allow 3D fast marching and for cpu-effective 
%	registration of 2D fast marching.
%
% Note(2):
%   Accuracy of this method is enhanced by just summing the coefficients
% 	of the cross and normal terms as suggested by Olivier Roy.
%
% Literature : M. Sabry Hassouna et Al. Multistencils Fast Marching 
%   Methods: A Highly Accurate Solution to the Eikonal Equation on
%   Cartesian Domains
%
%
% Example 2D,
%   SourcePoint = [51; 51];
%   SpeedImage = ones([101 101]);
%   [X Y] = ndgrid(1:101, 1:101);
%   T1 = sqrt((X-SourcePoint(1)).^2 + (Y-SourcePoint(2)).^2);
%
%   % Run fast marching 1th order, 1th order multi stencil 
%   % and 2th orde and 2th orde multi stencil
%
%   tic; T1_FMM1 = msfm(SpeedImage, SourcePoint, false, false); toc;
%   tic; T1_MSFM1 = msfm(SpeedImage, SourcePoint, false, true); toc;
%   tic; T1_FMM2 = msfm(SpeedImage, SourcePoint, true, false); toc;
%   tic; T1_MSFM2 = msfm(SpeedImage, SourcePoint, true, true); toc;
%
%   % Show results
%   fprintf('\nResults with T1 (Matlab)\n');
%   fprintf('Method   L1        L2        Linf\n');
%   Results = cellfun(@(x)([mean(abs(T1(:)-x(:))) mean((T1(:)-x(:)).^2) max(abs(T1(:)-x(:)))]), {T1_FMM1(:) T1_MSFM1(:) T1_FMM2(:) T1_MSFM2(:)}, 'UniformOutput',false);
%   fprintf('FMM1:   %9.5f %9.5f %9.5f\n', Results{1}(1), Results{1}(2), Results{1}(3));
%   fprintf('MSFM1:  %9.5f %9.5f %9.5f\n', Results{2}(1), Results{2}(2), Results{2}(3));
%   fprintf('FMM2:   %9.5f %9.5f %9.5f\n', Results{3}(1), Results{3}(2), Results{3}(3));
%   fprintf('MSFM2:  %9.5f %9.5f %9.5f\n', Results{4}(1), Results{4}(2), Results{4}(3));
%
%
% Example 2D, multiple starting points,
%
%   SourcePoint=rand(2,100)*255+1;
%   SpeedImage = ones([256 256]);
%   tic; T1_MSFM2 = msfm(SpeedImage, SourcePoint, true, true); toc;
%   figure, imshow(T1_MSFM2,[]); colormap(hot(256));
%
%
% Example 3D,
%   SourcePoint = [21; 21; 21];
%   SpeedImage = ones([41 41 41]);
%   [X,Y,Z] = ndgrid(1:41, 1:41, 1:41);
%   T1 = sqrt((X-SourcePoint(1)).^2 + (Y-SourcePoint(2)).^2 + (Z-SourcePoint(3)).^2);
%
%   % Run fast marching 1th order, 1th order multi stencil 
%   % and 2th orde and 2th orde multi stencil
%
%   tic; T1_FMM1 = msfm(SpeedImage, SourcePoint, false, false); toc;
%   tic; T1_MSFM1 = msfm(SpeedImage, SourcePoint, false, true); toc;
%   tic; T1_FMM2 = msfm(SpeedImage, SourcePoint, true, false); toc;
%   tic; T1_MSFM2 = msfm(SpeedImage, SourcePoint, true, true); toc;
%
%   % Show results
%   fprintf('\nResults with T1 (Matlab)\n');
%   fprintf('Method   L1        L2        Linf\n');
%   Results = cellfun(@(x)([mean(abs(T1(:)-x(:))) mean((T1(:)-x(:)).^2) max(abs(T1(:)-x(:)))]), {T1_FMM1(:) T1_MSFM1(:) T1_FMM2(:) T1_MSFM2(:)}, 'UniformOutput',false);
%   fprintf('FMM1:   %9.5f %9.5f %9.5f\n', Results{1}(1), Results{1}(2), Results{1}(3));
%   fprintf('MSFM1:  %9.5f %9.5f %9.5f\n', Results{2}(1), Results{2}(2), Results{2}(3));
%   fprintf('FMM2:   %9.5f %9.5f %9.5f\n', Results{3}(1), Results{3}(2), Results{3}(3));
%   fprintf('MSFM2:  %9.5f %9.5f %9.5f\n', Results{4}(1), Results{4}(2), Results{4}(3));
%
% Function is written by D.Kroon University of Twente (Oct 2010)
add_function_paths();

if(nargin<3), UseSecond=false; end
if(nargin<4), UseCross=false; end

if(nargout>1)
    if(size(F,3)>1)
        [T,Y]=msfm3d(F, SourcePoints, UseSecond, UseCross);        
    else
        [T,Y]=msfm2d(F, SourcePoints, UseSecond, UseCross);
    end
else
    if(size(F,3)>1)
        T=msfm3d(F, SourcePoints, UseSecond, UseCross);
    else
        T=msfm2d(F, SourcePoints, UseSecond, UseCross);
    end
end

function add_function_paths()
try
    functionname='msfm.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions'])
    addpath([functiondir '/shortestpath'])
catch me
    disp(me.message);
end
