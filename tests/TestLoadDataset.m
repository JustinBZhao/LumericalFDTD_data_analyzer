classdef TestLoadDataset < matlab.unittest.TestCase
    % Test suite for loadLumDataset method.

    properties
        ds_matrix % normal matrix dataset
        ds_recti  % normal rectilinear dataset
    end

    properties (TestParameter)
        dataset = loadTestParameter();
        % Define a set of non-struct inputs
        non_struct_scalar = struct(...
            "char_array", 'a char array', ...
            "string", "a string", ...
            "string_array", ["string array 1", "string array 2"], ...
            "cell_array", {{1, 2, 'cell3'}}, ...
            "logical", true, ...
            "numeric_array", [1, 2, 3], ...
            "numeric_scalar", 15, ...
            "numeric_muldim", rand(3, 5, 4), ...
            "complex", 4 + 5i, ...
            "anonymous_function", @(x) x, ...
            "struct_vector", struct('field1', {1, 'char', "string"}, 'field2', {"string", 2, 'char'}) ...
            );
    end

    methods (Test, TestTags = {'1'})
        function testDatasetType(testCase, non_struct_scalar)
            % Test that non-struct-scalar dataset type throws an error.
            testCase.verifyErrorMessage(@() loadLumDataset(non_struct_scalar), ...
                'Input dataset must be a struct scalar!');
        end

        function testDatasetFields1(testCase, dataset)
            % Test exist field 'Lumerical_dataset'
            ds_rm = rmfield(dataset, "Lumerical_dataset");
            testCase.verifyErrorMessage(@() loadLumDataset(ds_rm), ...
                'Input dataset does not have the field ''Lumerical_dataset''!');
        end

        function testDatasetFields2(testCase, dataset, non_struct_scalar)
            % Test field 'Lumerical_dataset' is struct scalar
            ds_wrong = dataset;
            ds_wrong.Lumerical_dataset = non_struct_scalar;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_wrong), ...
                'Field ''Lumerical_dataset'' is not a struct scalar!');
        end

        function testDatasetFields3(testCase, dataset)
            % Test field 'Lumerical_dataset' subfield 'attributes'
            ds_attributes_rm = dataset;
            ds_attributes_rm.Lumerical_dataset = rmfield(ds_attributes_rm.Lumerical_dataset, "attributes");
            testCase.verifyErrorMessage(@() loadLumDataset(ds_attributes_rm), ...
                'Field ''Lumerical_dataset'' does not have the ''attributes'' subfield!');
        end

        function testDatasetFields4(testCase, dataset)
            % Test field 'Lumerical_dataset' subfield 'parameters'
            ds_parameters_rm = dataset;
            ds_parameters_rm.Lumerical_dataset = rmfield(ds_parameters_rm.Lumerical_dataset, "parameters");
            testCase.verifyErrorMessage(@() loadLumDataset(ds_parameters_rm), ...
                'Field ''Lumerical_dataset'' does not have the ''parameters'' subfield!');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Category results hard to test, but can test throwing errors
    methods (Test, TestTags = {'recti_special'})
        function testDatasetClass(testCase)
            % Test wrong 'geometry' label
            ds_recti_ = testCase.ds_recti;
            ds_recti_.Lumerical_dataset.geometry = 2; % only case we test
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_), ...
                'Wrong label in ''lum_dataset.geometry'' for the rectilinear dataset!');
        end

        function testMissingXYZ(testCase)
            % Test missing x,y,z field
            ds_recti_rm_x = rmfield(testCase.ds_recti, 'x');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_rm_x), ...
                'No x field data in the rectilinear dataset!');
            ds_recti_rm_y = rmfield(testCase.ds_recti, 'y');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_rm_y), ...
                'No y field data in the rectilinear dataset!');
            ds_recti_rm_z = rmfield(testCase.ds_recti, 'z');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_rm_z), ...
                'No z field data in the rectilinear dataset!');
        end
    end

    methods (Test, TestTags = {'3'})
        % function 
    end

    methods (TestClassSetup)
        function loadDatasetOnce(testCase)
            data = load("normal_test_dataset.mat");
            testCase.ds_matrix = data.ds_matrix;
            testCase.ds_recti = data.ds_recti;
        end
    end

    methods (Access = private)
        function verifyErrorMessage(testCase, actual, errorMsg)
            % Test to verify correct error message thrown
            % Interface similar to verifyError, errorMsg must be char array
            try
                actual();
                testCase.verifyFail('An expected error was not thrown.');
            catch ME
                testCase.verifyEqual(ME.message, errorMsg, ...
                    'Incorrect error message was thrown! See below.');
            end
        end
    end
end

function dataset = loadTestParameter()
% Load dataset for the parameterized test
data = load("normal_test_dataset.mat");
dataset.ds_matrix = data.ds_matrix;
dataset.ds_recti = data.ds_recti;
end
