%"""
%biotracs.metadata.model.TableAnnotator
%Process to annotate a DataTable using experimental metadata.
%* License: BIOASTER License
%* Created by: Bioinformatics team, Omics Hub, BIOASTER Technology Research Institute (http://www.bioaster.org), 2017
%* See: biotracs.metadata.model.TableAnnotatorConfig
%"""

classdef TableAnnotator < biotracs.core.mvc.model.Process
    
    properties(Constant)
    end
    
    properties(SetAccess = protected)
        abbreviationTable = struct();
    end
    
    % -------------------------------------------------------
    % Public methods
    % -------------------------------------------------------
    
    methods
        
        % Constructor
        function this = TableAnnotator()
            this@biotracs.core.mvc.model.Process();            
            % enhance inputs specs
            this.addInputSpecs({...
                struct(...
                'name', 'DataTable',...
                'class', 'biotracs.data.model.DataTable' ...
                ),...
                struct(...
                'name', 'MetaTable',...
                'class', {{'biotracs.data.model.DataTable','biotracs.core.mvc.model.ResourceSet'}} ...
                )...
                });
            
            % enhance outputs specs
            this.addOutputSpecs({...
                struct(...
                'name', 'DataTable',...
                'class', 'biotracs.data.model.DataTable' ...
                )...
                });
        end
        
        %-- R --

        %-- S --
        
    end
    
    % -------------------------------------------------------
    % Protected methods
    % -------------------------------------------------------
    
    methods(Access = protected)
        
        % Prettify sample names in the FeatureSet
        function doRun( this, varargin )
            dataTable = this.getInputPortData('DataTable'); 
			dataTable = dataTable.copy();
            sampleTable = this.getInputPortData('MetaTable');

            if isa(sampleTable, 'biotracs.core.mvc.model.ResourceSet')
                if this.config.getParamValue('Verbose')
                    fprintf('\t Merging set of meta tables ...\n');
                end
                sampleTable = vertmerge( sampleTable.elements{:}, 'Force', true );  %Merge tables
            end
            
            %dataTable
            %sampleTable.summary
            
            sampleSelector = this.getConfig().getParamValue('SampleSelectorQuery');
            if ~isempty(sampleSelector)
                sampleTable = sampleTable.select(sampleSelector{:});
            end
            
            [ prettySampleNames ] = this.doBuildPrettyNamesFromSampleTable(sampleTable, varargin{:});
            if isempty(prettySampleNames), return; end
            
            %replace the old sample names by the pretty names
            if this.getConfig().getParamValue('Verbose')
                fprintf('Prettify row names\n');
            end
            
            rowNames = dataTable.getRowNames();
            namesOfRowsToRemove = {};
            for i=1:length(rowNames)
                originalSampleName = rowNames{i};
                originalSampleName = regexprep(originalSampleName, '.csv', '');
                if isKey(prettySampleNames, originalSampleName)
                    prettyName = prettySampleNames(originalSampleName);
                    if this.getConfig().getParamValue('Verbose')
                        fprintf('%s is now %s\n', originalSampleName, prettyName);
                    end
                else
                    prettyName = strrep(originalSampleName, ' ', '');
                    if this.getConfig().getParamValue('Verbose')
                        fprintf('%s is now %s\n', originalSampleName, prettyName);
                    end
                    fprintf('%s is not found in the sample table\n', originalSampleName);
                    if this.config.getParamValue('RemoveNonAnnotated')
                        namesOfRowsToRemove{end+1} = originalSampleName; %#ok<AGROW>
                    end
                end
                dataTable.setRowName( i, prettyName );
                dataTable.setRowTag( i, 'OriginalSampleName', originalSampleName);
            end
            
            if ~isempty(originalSampleName)
                fprintf('\nwarning: the following rows are removed for the data\n: %s\n', strjoin(namesOfRowsToRemove,', '));
                dataTable = dataTable.removeByRowName(namesOfRowsToRemove);
            end
            this.setOutputPortData('DataTable', dataTable);
        end
        
        % Give the pretty sample names using information in the sample
        % table
        function [ oPrettyNames ] = doBuildPrettyNamesFromSampleTable( this, iSampleTable, varargin )
            userIdHeader = this.config.getParamValue('UserIdHeader');
            headersToUse = this.config.getParamValue('HeadersToUse');
            tokenDelimiter = this.config.getParamValue('TokenDelimiter');
            keyValDelimiter = this.config.getParamValue('KeyValueDelimiter');
            useRegExp = this.config.getParamValue('UseRegExp');
            abbrvTable = this.config.getParamValue('AbbreviationTable');
   
            if isempty(headersToUse)
                error('Please give headers to use');
            end
            
            if ~isempty(userIdHeader)
                headersToUse = unique([{userIdHeader}, headersToUse], 'stable');
            else
                error('Parameter UserIdHeader is required');
            end
                
            %cellect indexes of headers to use
            nbHeadersToUse = length(headersToUse);
            headerIndexes = cell(nbHeadersToUse,1);
            for i=1:nbHeadersToUse
                header = headersToUse{i};
                if useRegExp
                    %custom regexp can be used
                    k =  iSampleTable.getColumnIndexesByName(header);
                else
                    %match against the whole word by default
                    k =  iSampleTable.getColumnIndexesByName(['^',header,'$']);
                end
                if ~isempty(k)
                    headerIndexes{i} = k(:);
                end
            end
            headerIndexes = cell2mat(headerIndexes);

            if isempty(headerIndexes)
                error('The headers are found in the sample table.');
            end
                
            %create names names
            [nbSamples,~] = getSize(iSampleTable);
            oPrettyNames = containers.Map();
            userIds = iSampleTable.getDataByColumnName(userIdHeader);
            if isempty(userIds)
                error('No column name matches with the UserFileIdHeader ''%s'' in the meta data file', userIdHeader)
            end

            for i=1:nbSamples
                currentSampleName = userIds{i};
                if isempty( currentSampleName )
                    fprintf('Sample name %d is an empty string. Skip it!\n', i);
                    continue;
                end
                
                if isKey(oPrettyNames, currentSampleName)
                    error('UniqueFileName ''%s'' is not unique (position %d). Please check the sample table', currentSampleName, i);
                end
                
                %Update the abbrv. table
                keys = fieldnames(this.abbreviationTable);
                for j = 1:length(keys)
                    key = keys{j};
                    if ~isfield(abbrvTable, key)
                        abbrvTable.(key) = this.abbreviationTable.(key);
                    end
                end
                
                nbHeaders = length(headerIndexes);
                tokens = cell(1,nbHeaders);
                %tokenKeys = cell(1,nbHeaders);
                %tokenVals = cell(1,nbHeaders);
                for j=1:nbHeaders
                    value = iSampleTable.getDataAt( i, headerIndexes(j) );
                    if isnumeric(value)
                        if isnan(value)
                            value = '';
                        else
                            value = num2str(value);
                        end
                    end
                    header = iSampleTable.getColumnName( headerIndexes(j) );
                    %clean headers and values
                    tokenKey = regexprep( header, {'(^\s+)|(\s+$)', '\s+', [keyValDelimiter,'|',tokenDelimiter]}, {'', ' ', '-'} );
                    tokenVal = regexprep( value, {'(^\s+)|(\s+$)', '\s+', [keyValDelimiter,'|',tokenDelimiter]}, {'', ' ', '-'} );

                    %abbreviate tokenKeys
                    keys = fieldnames(abbrvTable);
                    for k=1:length(keys)
                        key = keys{k};
                        tokenKey = regexprep( tokenKey, ['^',key,'$'], abbrvTable.(key) );
                    end
                
                    tokens{j} = [ tokenKey, ':', tokenVal ];
                end

                %remove empty tokens
                idx = cellfun(@isempty, tokens);
                tokens(idx) = [];
                
                uniqueIdVal = regexprep( iSampleTable.rowNames{i}, {'(^\s+)|(\s+$)', '\s+', [keyValDelimiter,'|',tokenDelimiter]}, {'', ' ', '-'} );
                tokens = [['UniqueID:',uniqueIdVal], tokens];
                prettyName = strjoin(tokens,tokenDelimiter);
                if isempty(prettyName)
                    frpintf('Warning: cannot build the pretty name of %s (position %d)\n', currentSampleName, i);
                    continue;
                end

                oPrettyNames(currentSampleName) = prettyName;
            end
            
        end

    end
end
