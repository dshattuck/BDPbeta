% 
% BDP BrainSuite Diffusion Pipeline
% 
% Copyright (C) 2023 The Regents of the University of California and
% the University of Southern California
% 
% Created by Chitresh Bhushan, Divya Varadarajan, Justin P. Haldar, Anand A. Joshi,
%            David W. Shattuck, and Richard M. Leahy
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; version 2.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
% USA.
% 


function rand_str = randstr(n)
% generates a random string of lower case letters and numbers of length n


if usejava('jvm')
   tmp_name = strrep(char(java.util.UUID.randomUUID),'-','');
   while length(tmp_name)<n
      tmp_name = [tmp_name strrep(char(java.util.UUID.randomUUID),'-','')];
   end
else
   tmp_name = num2str(feature('timing','cpucount'));
   while length(tmp_name)<n
      tmp_name = [tmp_name num2str(feature('timing','cpucount'))];
   end
end

rand_str = tmp_name(1:n);

end
