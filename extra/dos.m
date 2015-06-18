%
% Extra function to compare DOS, with get from tight-binding.
% Please be sure eigenvalues data (eigen.dat file ) exists.
%
% Create eigen.dat example, and run this script:
%   $ cat input/cadenaPt.xyz | ./tight-biding -e > eigen.dat
%   $ octave --silent dos.m 
%
% Gonzalo Aguirre <graguirre@gmail.com>
%

e=load('eigen.dat');

e=e(:,2); % get 2nd column
eta=5e-3; % delta

xe = [min(e):0.001:max(e)]; % X-axis

DOS=zeros(1,length(xe));    % Y-axis init

i=1;
for ii=[min(e):0.001:max(e)]
	d=0;
	for jj=[1:length(e)]
		d = d + 1./(ii - e(jj) + j*eta) ;
	end
	DOS(i) = d;
	i += 1;
end

for i=[1:length(xe)]
	printf("%.3f %f\n",xe(i),-imag(DOS(i))/pi)
end
%plot(xe,-imag(DOS),'-');

