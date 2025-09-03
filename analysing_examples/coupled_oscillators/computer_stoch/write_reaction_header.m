function [] = write_reaction_header(Mobj,text)
%WRITE_REACTION_HEADER creates the header file with reactions for
%stochastic simulations
% Mobj: simple SimBiology model object
% text: a cell, column vector shape with strings, its content will be
% written at the beginning of the "reactions" function

% check that the model not contains extra elements
if ~isempty(Mobj.Events)
    error('Events are not implemented!')
end
if ~isempty(Mobj.Rules)
    error('Rules are not implemented!')
end

fid = fopen('reactions.h','w');
if fid == -1
    error('Could not open reactions.h for writing.');
end

fprintf(fid, '#ifndef REACTIONS_H\n');
fprintf(fid, '#define REACTIONS_H\n\n');

fprintf(fid, '#include <cmath>\n');
fprintf(fid, '#include <random>\n\n');

fprintf(fid, '#include <iostream>\n');
fprintf(fid, 'using std::pow;\n');
fprintf(fid, 'using std::sqrt;\n\n');

fprintf(fid, '#include "mex_eig_dense.h"\n');
fprintf(fid, '#include "parameters.h"\n\n\n');

fprintf(fid, '// ********  REACTIONS ********\n');
fprintf(fid, 'inline void reactions(double c[],std::normal_distribution<double> &distribution,std::mt19937 &generator) {\n\n');

fprintf(fid, '\tdouble r[NREACT]; // reatction rates\n');
fprintf(fid, '\tstatic double W[NREACT] = {0}; // Wiener process\n');
fprintf(fid, '\tdouble decay = std::exp(-dt/p.tau); // decay in the noise\n');

if nargin > 1
    fprintf(fid, '\n');
    for i = 1:size(text,1)
        fprintf(fid, '\t%s\n',text{i});
    end
end

% ---------- write the reactions -------------

% print reactions
% rename the parameters for structure format
Mobj = rename_parameters(Mobj);
% rename species
for i = 1:numel(Mobj.Species)
    rename(Mobj.Species(i),['c' int2str(i-1)]);
end
fprintf(fid, '\n\t// Reactions\n');
for i = 1:numel(Mobj.Reactions)
    reaction = Mobj.Reactions(i).ReactionRate;
    % remove [] brackets around the p.ParameterName
    reaction = strrep(reaction,'[','');
    reaction = strrep(reaction,']','');
    % Rename c# to c[#]
    reaction = regexprep(reaction, 'c(\d+)', 'c[$1]');
    % print the reaction rates 
    fprintf(fid, '\tr[%d] = %s;\n',i-1,reaction);
end
fprintf(fid, '\n');

M = getstoichmatrix(Mobj);

% generate wiener noise
fprintf(fid, '\n\t// Wiener process\n');
fprintf(fid, '\tfor (size_t nreact=0;nreact<NREACT;++nreact) {\n\t\tW[nreact] *= decay;\n\t\tW[nreact] += distribution(generator);\n\t}\n');

% change reactions according to the wiener process
fprintf(fid, '\n\t// Change reactions with Wiener process\n');
fprintf(fid, '\tfor (size_t nreact=0;nreact<NREACT;++nreact) {\n\t\tr[nreact] = r[nreact]*dt + sqrt(r[nreact])*W[nreact]*sqrt(dt);\n\t}\n');

fprintf(fid, '\n\tswitch (SIMTYPE) {\n');
fprintf(fid, '\t\tcase Euler:\n');
fprintf(fid, '\t\t\t// Euler–Maruyama method\n');

% ---- Euler–Maruyama ----
% print concentration changes
for nspec = 1:size(M,1)
    % if there is no reaction for the species
    if all((M(nspec,:)==0))
        continue;
    end
    fprintf(fid, '\t\t\tc[%d] += (',nspec-1);
    % reaction terms
    for nreact = 1:size(M,2)
        if M(nspec,nreact) ~= 0
            if full(M(nspec,nreact))==1
                fprintf(fid, '+');
            elseif full(M(nspec,nreact))==-1
                fprintf(fid, '-');
            else
                fprintf(fid, '%+d.*',full(M(nspec,nreact)));
            end
            fprintf(fid, 'r[%d]',nreact-1);
        end
    end
    fprintf(fid, ');\n');
end

fprintf(fid, '\n\t\t\tbreak;\n');
fprintf(fid, '\n\t\tcase Heun:\n');
fprintf(fid, '\t\t\tdouble c_pred[NSPEC];\n');

% --- Heun method ---

% Heun prediction
% print concentration changes
fprintf(fid, '\n\t\t\t// Stochastic Heun method, prediction\n');
for nspec = 1:size(M,1)
    % if there is no reaction for the species
    if all((M(nspec,:)==0))
        continue;
    end
    fprintf(fid, '\t\t\tc_pred[%d] = c[%d] + (',nspec-1,nspec-1);
    % reaction terms
    for nreact = 1:size(M,2)
        if M(nspec,nreact) ~= 0
            if full(M(nspec,nreact))==1
                fprintf(fid, '+');
            elseif full(M(nspec,nreact))==-1
                fprintf(fid, '-');
            else
                fprintf(fid, '%+d.*',full(M(nspec,nreact)));
            end
            fprintf(fid, 'r[%d]',nreact-1);
        end
    end
    fprintf(fid, ');\n');
end

fprintf(fid, '\n\t\t\t// Preserve positivity\n');
fprintf(fid, '\t\t\tfor (size_t nspec=0;nspec<NSPEC;++nspec) {\n\t\t\t\tif(c_pred[nspec]<0.0) {\n\t\t\t\t\tc_pred[nspec] *= -1.0;\n\t\t\t\t}\n\t\t\t}\n');

fprintf(fid, '\n\t\t\t// Stochastic Heun correction\n');

fprintf(fid, '\n\t\t\t// Reactions\n');
for i = 1:numel(Mobj.Reactions)
    reaction = Mobj.Reactions(i).ReactionRate;
    % remove [] brackets around the p.ParameterName
    reaction = strrep(reaction,'[','');
    reaction = strrep(reaction,']','');
    % Rename c# to c[#]
    reaction = regexprep(reaction, 'c(\d+)', 'c_pred[$1]');
    % print the reaction rates 
    fprintf(fid, '\t\t\tr[%d] = %s;\n',i-1,reaction);
end
fprintf(fid, '\n');
% change reactions according to the wiener process
fprintf(fid, '\n\t\t\t// Change reactions with Wiener process\n');
fprintf(fid, '\t\t\tfor (size_t nreact=0;nreact<NREACT;++nreact) {\n\t\t\t\tr[nreact] = r[nreact]*dt + sqrt(r[nreact])*W[nreact]*sqrt(dt);\n\t\t\t}\n');

% Heun correction
% print concentration changes
fprintf(fid, '\n\t\t\t// Stochastic Heun method, correction\n');
for nspec = 1:size(M,1)
    % if there is no reaction for the species
    if all((M(nspec,:)==0))
        continue;
    end
    fprintf(fid, '\t\t\tc[%d] = 0.5*(c_pred[%d] + c[%d] + (',nspec-1,nspec-1,nspec-1);
    % reaction terms
    for nreact = 1:size(M,2)
        if M(nspec,nreact) ~= 0
            if full(M(nspec,nreact))==1
                fprintf(fid, '+');
            elseif full(M(nspec,nreact))==-1
                fprintf(fid, '-');
            else
                fprintf(fid, '%+d.*',full(M(nspec,nreact)));
            end
            fprintf(fid, 'r[%d]',nreact-1);
        end
    end
    fprintf(fid, '));\n');
end

fprintf(fid, '\t\t\tbreak;\n');
fprintf(fid, '\t}\n');

% positivity
fprintf(fid, '\n\t// Preserve positivity\n');
fprintf(fid, '\tfor (size_t nspec=0;nspec<NSPEC;++nspec) {\n\t\t\tif(c[nspec]<0.0) {\n\t\t\tc[nspec] *= -1.0;\n\t\t}\n\t}\n');

% close the function
fprintf(fid, '}\n');
fprintf(fid, '#endif //REACTIONS_H\n');

fclose(fid);

% change back parameter names from the struct from "p.XY" -> "XY"
% rename the parameters to store them in "p." struct later
for i = 1:numel(Mobj.Parameters)
    parname = Mobj.Parameter(i).Name;
    rename(Mobj.Parameter(i),parname(3:end));
end

end

