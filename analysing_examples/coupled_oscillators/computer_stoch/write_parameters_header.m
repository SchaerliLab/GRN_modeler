function write_parameters_header(Mobj)
%WRITE_PARAMETERS_HEADER Generate parameters.h header file
% Mobj: simple SimBiology model object

fid = fopen('parameters.h','w');
if fid == -1
    error('Could not open parameters.h for writing.');
end

fprintf(fid, '#ifndef PARAMETERS_H\n');
fprintf(fid, '#define PARAMETERS_H\n\n');

fprintf(fid, '#include "mex.h"\n');
fprintf(fid, '#include "mex_eig_dense.h"\n\n');

fprintf(fid, '// number of species\n');
fprintf(fid, 'inline constexpr size_t NSPEC = %d;\n',numel(Mobj.Species));

fprintf(fid, '// number of reactions\n');
fprintf(fid, 'inline constexpr size_t NREACT = %d;\n\n',numel(Mobj.Reactions));

fprintf(fid, 'struct parameters {\n\n');
fprintf(fid, '\t// reaction kinetics paramertes\n');

% write reaction kinetics parameter names
fprintf(fid, '\tdouble ');
for i = 1:numel(Mobj.Parameters)
    if i == 1
        fprintf(fid, '%s',Mobj.Parameters(i).Name);
    else
        fprintf(fid, ', %s',Mobj.Parameters(i).Name);
    end
end
fprintf(fid, ';\n');

% extra simulation parameters
fprintf(fid, '\t// simulation parameters\n');
fprintf(fid, '\tdouble t_sim, sigma_inf, tau;\n');
fprintf(fid, '\tsize_t n_step, n_data;\n');
fprintf(fid, '};\n\n');

fprintf(fid, '// parameters\n');
fprintf(fid, 'extern parameters p;\n\n');

fprintf(fid, '// time step\n');
fprintf(fid, 'extern double t,dt;\n\n');

fprintf(fid, '// load the parameters\n');
fprintf(fid, 'void load_parameters(const mxArray *prhs[]);\n\n');

fprintf(fid, '#endif // PARAMETERS_H\n');

fclose(fid);

end
