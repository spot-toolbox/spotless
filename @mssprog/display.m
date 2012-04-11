function display(p)
% function display(p)
%
% display mss program as a structure

fprintf('\n mss program: ')
fprintf('%d free, ',length(p.o(p.t==1)))
fprintf('%d pos, ',length(p.o(p.t==2)))
fprintf('%d lor, ',length(p.o(p.t==3)))
fprintf('%d rlor, ',length(p.o(p.t==4)))
fprintf('%d psd, ',length(p.o(p.t==5)))
fprintf('%d eqs, ',length(p.e))
if ~isempty(p.x),
    fprintf('solution available');
else
    fprintf('no solution');
end
fprintf('\n\n')