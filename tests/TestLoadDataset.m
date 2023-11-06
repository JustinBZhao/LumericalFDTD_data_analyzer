classdef TestLoadDataset < matlab.unittest.TestCase
    % Test suite for loadLumDataset method.

    properties
        ds_matrix % normal matrix dataset
        ds_recti  % normal rectilinear dataset
    end

    properties (TestParameter)
        dataset = loadTestParameter();

        non_struct_scalar = struct(...
            "char_array", 'a char array', ...
            "string", "a string", ...
            "string_array", ["string array 1", "string array 2"], ...
            "cell_array", {{1, 2, 'cell3'}}, ...
            "numeric_array", [1, 2, 3], ...
            "struct_array", struct('field1', {1, 'char', "string"}, 'field2', {"string", 2, 'char'}) ...
            );

        non_geometry = struct(...
            "wrong_word", "matrix", ...
            "string_array", ["string array 1", "string array 2"], ...
            "numeric_array", [1, 2, 3] ...
            );

        non_numeric = struct(...
            "char_array", 'a char array', ...
            "string", "a string", ...
            "string_array", ["string array 1", "string array 2"], ...
            "cell_array", {{1, 2, 'cell3'}}, ...
            "logical", true, ...
            "struct_scalar", struct('field1', 2, 'field2', "haha"), ...
            "anonymous_function", @(x) x, ...
            "struct_array", struct('field1', {1, 'char', "string"}, 'field2', {"string", 2, 'char'}) ...
            );

        non_cell_column = struct(...
            "cell_row", {{1, "a", 'xyz'}}, ...
            "cell_matrix", {{1, "a"; [1, 2, 3], 'xyz'}}, ...
            "struct_array", struct('field1', {1, 'char', "string"}, 'field2', {"string", 2, 'char'}) ...
            );

        non_struct = struct(...
            "char_array", 'a char array', ...
            "string", "a string", ...
            "string_array", ["string array 1", "string array 2"], ...
            "cell_array", {{1, 2, 'cell3'}}, ...
            "numeric_array", [1, 2, 3] ...
            );

        non_text_scalar = struct(...
            "string_array", ["string1", "string2"], ...
            "char_array", {{'string1', 'string2'}} ...
            );

        non_var_name = struct(...
            "leading_digit", "1name", ...
            "space", "name 1", ...
            "too_long", repmat('a', 1, 100) ...
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
    % Category results (non-throwing) hard to test, but can test throwing errors
    methods (Test, TestTags = {'recti_special'})
        function testDatasetClass(testCase, non_geometry)
            % Test wrong 'geometry' label
            ds_recti_ = testCase.ds_recti;
            ds_recti_.Lumerical_dataset.geometry = non_geometry;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_), ...
                'Wrong label in ''lum_dataset.geometry'' for the rectilinear dataset!');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % May parameterize x,y,z
        function testMissingXYZ(testCase)
            % Test missing x,y,z field
            ds_recti_rm_x = rmfield(testCase.ds_recti, 'x');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_rm_x), ...
                'No x data in the rectilinear dataset!');
            ds_recti_rm_y = rmfield(testCase.ds_recti, 'y');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_rm_y), ...
                'No y data in the rectilinear dataset!');
            ds_recti_rm_z = rmfield(testCase.ds_recti, 'z');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_rm_z), ...
                'No z data in the rectilinear dataset!');
        end

        function testBadXYZNonNumeric(testCase, non_numeric)
            % Test bad x,y,z values: non-numeric
            ds_recti_ = testCase.ds_recti;
            ds_recti_.x = non_numeric;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_), ...
                'x data must be a numeric vector!');
            ds_recti_ = testCase.ds_recti;
            ds_recti_.y = non_numeric;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_), ...
                'y data must be a numeric vector!');
            ds_recti_ = testCase.ds_recti;
            ds_recti_.z = non_numeric;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_), ...
                'z data must be a numeric vector!');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % You cannot test a warning unless you can terminate the execution
        % at the point of warning generation
        function testBadXYZMuldim(testCase)
            % Test bad x,y,z values: multi-dimensional matrix
            % ds_recti_ = testCase.ds_recti;
            % ds_recti_.x = rand(3, 5, 2);
            % testCase.verifyWarning(@() loadLumDataset(ds_recti_), ...
            %     'TestID');
            testCase.verifyFail;
        end

        function testBadXYZImag(testCase)
            % ds_recti_ = testCase.ds_recti;
            % ds_recti_.x = complex(ds_recti_.x); %%%%%%%%%%%%% test??
            testCase.verifyFail;
        end
    end

    methods (Test, TestTags = {'parameters'})
        function testParametersType(testCase, dataset, non_cell_column)
            % 'parameters' type cell, column vector
            ds = dataset;
            ds.Lumerical_dataset.parameters = non_cell_column;
            testCase.verifyErrorMessage(@() loadLumDataset(ds), ...
                'Field ''Lumerical_dataset.parameters'' should be a cell column vector!');
        end

        %%%%%%% Only test one parameter in MATRIX dataset
        function testParameterContent1(testCase, non_struct)
            % Test the cell in 'parameters' with correct format
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1} = non_struct;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
        end

        function testParameterContent2(testCase)
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1} = rmfield(ds_matrix_.Lumerical_dataset.parameters{1}, 'variable');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
            ds_matrix_.Lumerical_dataset.parameters{1} = rmfield(ds_matrix_.Lumerical_dataset.parameters{1}, 'name');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
        end

        function testSingleParameterContent1(testCase, non_text_scalar)
            % Test one interdependent parameter
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).variable = non_text_scalar;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
        end

        function testSingleParameterContent2(testCase, non_var_name)
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).variable = non_var_name;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
        end

        function testIllegalCharacters(testCase)
            %%%%%%%%%%%%%%%%% Same warning problem
            testCase.verifyFail;
        end

        function testParameterInDataset(testCase)
            % Test parameter data exists
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_rm = rmfield(ds_matrix_, parameter_name);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_rm), ...
                ['Parameter field ''', parameter_name, ''' data not found!']);
        end
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
