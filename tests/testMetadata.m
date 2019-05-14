%"""
%Unit tests for biotracs.metadata.*
%* License: BIOASTER License
%* Create: 2016
%Bioinformatics team, Omics Hub, BIOASTER Technology Research Institute (http://www.bioaster.org)
%"""

function testMetadata( cleanAll )
    if nargin == 0 || cleanAll
        clc; close all force;
        restoredefaultpath();
    end

    autoload( ...
        'PkgPaths', {fullfile(pwd, '../../')}, ...
        'Dependencies', {...
            'biotracs-m-metadata', ...
        }, ...
        'Variables',  struct(...
        ) ...
    );

    %% Tests
    import matlab.unittest.TestSuite;
    Tests = TestSuite.fromFolder('./', 'IncludingSubfolders', true);
    Tests.run;
end