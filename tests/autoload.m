%"""
%autoloader.m
%Autoloader to load the dependency packages. This file must be called by the startup (or main) 
%file to load biotracs libraries (the core application, and other applications).
%* Date: 2019
%* Author:  D. A. Ouattara
%* License: BIOASTER License
%
%Omics Hub, Bioinformatics team
%BIOASTER Technology Research Institute (http://www.bioaster.org)
%"""

function autoload( varargin )
    p = inputParser();
    p.addParameter('Paths', struct(), @isstruct);
    p.addParameter('Variables', struct(), @isstruct);
    p.addParameter('PkgPaths', {}, @(x)iscellstr(x)); %#ok<ISCLSTR>
    p.addParameter('Dependencies', {}, @(x)iscellstr(x)); %#ok<ISCLSTR>
    p.addParameter('Production', false, @islogical);  
    p.addParameter('WorkingDirectory', '', @ischar);
    p.addParameter('Verbose', true, @islogical);
    p.KeepUnmatched = true;
    p.parse(varargin{:});

    if p.Results.Verbose
        fprintf('BIOTRACS\n\n');
    end
    
    % load dependencies
    if ~isdeployed
        [ depPaths, depVars ] = createDepPaths(p.Results.PkgPaths, p.Results.Dependencies);
        
        if isempty(depPaths)
            if p.Results.Verbose
                if isempty(p.Results.Dependencies)
                    fprintf('No dependencies provided.\n');
                else
                    error('The dependencies could not be loaded. Please check PkgPaths.');
                end
            end
        else
            if p.Results.Verbose
                fprintf('Loading dependencies ...\n');
            end
            
            loadDep(depPaths, p.Results.Verbose);
            
            if p.Results.Verbose
                fprintf('The dependencies have been successfully loaded.\n\n');
            end
        end
    end
    
    % Set environment variables
    biotracs.core.env.Env.depPaths(depPaths);
    
    if ~isempty(p.Results.Variables)
        f = fieldnames(p.Results.Variables);
        for i=1:length(f)
            depVars.(f{i}) = p.Results.Variables.(f{i});
        end
    end
    biotracs.core.env.Env.vars( depVars );

    depNames = cell(size(depPaths));
    tokens = cell(size(depPaths));
    for i=1:length(depPaths)
        tab = strsplit(fullfile(depPaths{i}), filesep);
        depNames{i} = tab{end};
        tokens{i} = strrep(upper(tab{end}),'-','_');
        tokens{i} = strcat(tokens{i},'_DIR');
        biotracs.core.env.Env.vars(struct(tokens{i}, depPaths{i}));
    end
    biotracs.core.env.Env.depNames(depNames);
    biotracs.core.env.Env.depPathTokens(tokens);
    
    try
        if ischar(p.Results.WorkingDirectory) && ~strcmp( p.Results.WorkingDirectory, '' )
            workingDir = p.Results.WorkingDirectory;
        else
            workingDir = biotracs.core.env.Env.userDir();
        end
        workingDir = fullfile(workingDir, biotracs.core.env.Env.company(), biotracs.core.env.Env.name());
        biotracs.core.env.Env.workingDir( workingDir );
        biotracs.core.env.Env.logDir( workingDir );
    catch err
        error('Biocode has not been loaded properly. Please check module directories.\n%s', err.message);
    end    
end

function [ oDepPaths, oDepVars ] = createDepPaths( iRootPaths, iDeps, iExceptions )
    oDepPaths = {};
    foundDeps = {};
    oDepVars = struct();
    if isempty(iDeps)
        return;
    end

    if nargin <= 2
        iExceptions = {};
    end
    
%     if nargin <= 3
%         oDepVars = struct();
%     else
%         oDepVars = iDepVars;
%     end
    
    isRelativeDir = ~cellfun(@isempty,  regexp(iRootPaths, '^\.', 'once'));
    if isRelativeDir
        error('The RootPaths must be absolutes paths');
    end
        
    for i=1:length(iRootPaths)
        for j=1:length(iDeps)
            if any(ismember(foundDeps, iDeps{j}))
                continue;
            end
            root = iRootPaths{i};
            dep = iDeps{j};
            depPath = fullfile(root, dep);
            if exist(depPath, 'dir') == 7
                oDepPaths{end+1} = depPath; %#ok<AGROW>
                foundDeps{end+1} = dep; %#ok<AGROW>
                
                %search 'package.json' file
                appFile = fullfile(depPath, 'package.json');
                if exist(appFile, 'file') == 2
                    
                    try
                        data = jsondecode(fileread(appFile));
                    catch exception
                        error('BIOTRACS:Autoload:InvalidPkgFile', 'An error occured while loading the pkg file ''package.json''. Please check.\n %s', exception.message)
                    end
                    
                    if isfield(data, 'variables')
                        oDepVars = data.variables;
                    end
                    
                    if isfield(data, 'dependencies')
                        subDeps = {};
                        for k=1:length(data.dependencies)
                            if ~any(ismember(foundDeps, data.dependencies{k})) && ...
                                ~any(ismember(iExceptions, data.dependencies{k}))
                                subDeps{end+1} = data.dependencies{k}; %#ok<AGROW>
                            end
                        end
                        
                        [ dPaths, dVars ] = createDepPaths( iRootPaths, subDeps, foundDeps );
                        
                        %concatenation of dependency paths
                        oDepPaths = [ oDepPaths, dPaths ]; %#ok<AGROW>
                        
                        %concatenation of dependency variables                        
                        f = fieldnames(dVars);
                        for kk=1:length(f)
                            oDepVars.(f{kk}) = dVars.(f{kk});
                        end
                    end

                end
            end
        end
    end
           
    oDepPaths = unique(oDepPaths);
end

%
% Load dependencies
function loadDep( iDep, iVerbose ) 
    for i=1:length( iDep )
        if iVerbose
            tab = strsplit(fullfile(iDep{i}), filesep);
            fprintf(' %s -> %s\n', tab{end}, iDep{i});
        end
        
        moduleDir = iDep{i};

        if exist(moduleDir, 'dir') == 7
            addpath(moduleDir);

            %check if the module contains an 'externs' subdirectory
            externDir = fullfile(moduleDir, 'externs/matlab');
            if exist( externDir, 'dir' ) == 7
                loadRecursive( externDir );
            end

            %check if the module contains an 'backcomp' subdirectory
            backcompDir = fullfile(moduleDir, 'backcomp');
            if exist( backcompDir, 'dir' ) == 7
                loadBackcomp( backcompDir );
            end
        else
            error('BIOTRACS:InvalidModule', 'The module ''%s'' is not found', moduleDir);
        end
    end

end

% Load paths recursively
function loadRecursive( moduleDir )
    list = subdir(moduleDir);
    p = arrayfun( @(x)(fullfile(x.folder, x.name)), list, 'UniformOutput', false);
    addpath( p{:} );
end

% Load backward compatibility paths
function loadBackcomp( backcompDir )
    list = subdir(backcompDir);
    p = arrayfun( @(x)(fullfile(x.folder, x.name)), list, 'UniformOutput', false);

    currentVer = strcat('R',version('-release'));
    for j=1:length(p)
        folderVer = regexprep(p{j}, '.*(R\d+\w)([/\\\.]+)?$', '$1');
        list = sort({currentVer, folderVer});
        isFolderVersionNewerThanCurrentVersion = strcmp(list{2}, folderVer);     
        if isFolderVersionNewerThanCurrentVersion
            addpath( p{j} );
        end
    end
end

% List all sub directories
function [ list ] = subdir( path )
    p = dir(path);
    list = [];
    for i=1:length(p)
       if strcmp(p(i).name, '.')
           list = p(i);
       elseif regexp(p(i).name, '(^\.\.)|(^\..+)', 'once')
           continue;
       elseif p(i).isdir && isempty(regexp(p(i).name, '(^\+.*)|(^@.*)', 'once'))
           d = fullfile(p(i).folder, p(i).name);
           list = [list, subdir(d)]; %#ok<AGROW>
       end
    end
end