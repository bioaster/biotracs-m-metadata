classdef TableAnnotatorTests < matlab.unittest.TestCase
    
    properties (TestParameter)
    end
    
    properties
        workingDir = fullfile(biotracs.core.env.Env.workingDir(), '/biotracs/metadata/TableAnnotatorTests');
    end

    methods (Test)
        
        function testTableAnnotation(testCase)
            metaTable = biotracs.data.model.DataTable.import( '../../../testdata/sample_metadata.xlsx' );
            dataSet = biotracs.data.model.DataSet.import( '../../../testdata/small_amix_data.csv'  );
            
            process = biotracs.metadata.model.TableAnnotator();
            c = process.getConfig();
            c.updateParamValue('HeadersToUse', {'Acronym','BiologicalReplicate'});
            c.updateParamValue('UserIdHeader', 'UniqueFileName');
            c.updateParamValue('WorkingDirectory', testCase.workingDir);
            process.setInputPortData('MetaTable', metaTable);
            process.setInputPortData('DataTable', dataSet);
            process.run();
            
            annotatedDataTable = process.getOutputPortData('DataTable');
            
            testCase.verifyEqual( metaTable.rowNames{1}, 'ID1' );
            testCase.verifyEqual( metaTable.rowNames{3}, 'ID3' );
            testCase.verifyEqual( metaTable.rowNames{4}, 'ID4' );
            testCase.verifyEqual( metaTable.rowNames{10}, 'ID10' );
            testCase.verifyEqual( dataSet.rowNames{1}, 'S1' );
            testCase.verifyEqual( dataSet.rowNames{3}, 'S3' );
            testCase.verifyEqual( dataSet.rowNames{4}, 'S4' );
            testCase.verifyEqual( dataSet.rowNames{10}, 'S10' );
            testCase.verifyEqual( annotatedDataTable.rowNames{1}, 'UniqueID:ID1_UniqueFileName:S1_Acronym:CEX_BiologicalReplicate:10' );
            testCase.verifyEqual( annotatedDataTable.rowNames{3}, 'UniqueID:ID3_UniqueFileName:S3_Acronym:CEX_BiologicalReplicate:12' );
            testCase.verifyEqual( annotatedDataTable.rowNames{4}, 'UniqueID:ID4_UniqueFileName:S4_Acronym:_BiologicalReplicate:13' );
            testCase.verifyEqual( annotatedDataTable.rowNames{10}, 'UniqueID:ID10_UniqueFileName:S10_Acronym:DOR_BiologicalReplicate:7' );
        end
        
    end


end