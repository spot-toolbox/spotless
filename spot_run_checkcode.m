function map = spot_run_checkcode(path)
    if nargin < 1,
        path = fileparts(mfilename('fullpath'));
    end
    paths = vertcat({path},list_subdirectories(path));

    result = containers.Map;
    for j = 1:length(paths)
        D = dir([paths{j} filesep '*.m']);
        for f = 1:length(D)
            fn = [paths{j} filesep D(f).name];
            check = checkcode(fn);
            if size(check, 1) > 0
                result(fn) = check;
            end
        end
    end
    
    keys = result.keys;
    for k = 1:length(keys)
        print_report(keys{k}, result(keys{k}));
    end
    if nargout > 0,
        map = result;
    end
end

function [] = print_report(fn, r)
    disp(repmat('<', 1, 60))
    disp(fn)
    for i = 1:length(r),
        disp(sprintf('%5d: %s', r(i).line, r(i).message));
    end
    disp(repmat('>', 1, 60))
end

function paths = list_subdirectories(path)
    D = dir(path);
    names = {D.name};
    isdir = {D.isdir};
    paths = {};
    for i = 1:length(names)
        if names{i}(1) == '.',
            continue;
        elseif isdir{i},
            sub_path = [path filesep names{i}];
            paths{end+1,1} = sub_path;
            paths = vertcat(paths, list_subdirectories(sub_path));
        end                
    end
end