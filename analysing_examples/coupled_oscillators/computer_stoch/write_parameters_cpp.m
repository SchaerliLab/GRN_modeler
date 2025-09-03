function write_parameters_cpp(Mobj)
%WRITE_PARAMETERS_CPP Generate parameters.cpp implementation file
% Mobj: simple SimBiology model object

fid = fopen('parameters.cpp','w');
if fid == -1
    error('Could not open parameters.cpp for writing.');
end

fprintf(fid, '#include "parameters.h"\n\n');
fprintf(fid, '// define shared parameters\n');
fprintf(fid, 'parameters p;\n');
fprintf(fid, 'double t;\n');
fprintf(fid, 'double dt;\n\n');

fprintf(fid, '// load parameters\n');
fprintf(fid, 'void load_parameters(const mxArray *prhs[]) {\n');

% read every parameters
fprintf(fid, '\n\t//Reaction kinetics parameters\n');
for i = 1:numel(Mobj.Parameters)
    fprintf(fid, '\tp.%s = mat2eig_double(prhs[0], "%s");\n',...
        Mobj.Parameters(i).Name,Mobj.Parameters(i).Name);
end

fprintf(fid, '\n\t//Simulation parameters\n');
% other simulation parameters
fprintf(fid, '\tp.t_sim = mat2eig_double(prhs[0], "t_sim");\n');
% fprintf(fid, '\tp.omega   = mat2eig_double(prhs[0], "omega");\n');
fprintf(fid, '\tp.sigma_inf = mat2eig_double(prhs[0], "sigma_inf");\n');
fprintf(fid, '\tp.tau = mat2eig_double(prhs[0], "tau");\n\n');

fprintf(fid, '\tp.n_step = mat2eig_size_t(prhs[0], "n_step");\n');
fprintf(fid, '\tp.n_data = mat2eig_size_t(prhs[0], "n_data");\n');
fprintf(fid, '}\n');

fclose(fid);


end
