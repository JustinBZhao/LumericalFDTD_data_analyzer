classdef TestLoadDataset < matlab.unittest.TestCase
    % Test suite for loadLumDataset method.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % How do you test if something does NOT generate the error?
    % How do you test warnings?

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
            "empty", [], ...
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

        non_numeric_vector = struct(...
            "empty", [], ...
            "matrix", [1, 2; 3, 4], ...
            "cell_array", {{1, 2, 3}} ...
            );
    end

    methods (Test, TestTags = {'Golden master'})
        function testCorrectDataset(testCase)
            % Test correct dataset output
            golden = load("golden_master_normal_dataset.mat");
            testCase.verifyEqual(LumericalDataset.createObject(testCase.ds_matrix), golden.dataset_matrix);
            testCase.verifyEqual(LumericalDataset.createObject(testCase.ds_recti), golden.dataset_recti);
        end
    end

    methods (Test, TestTags = {'Special dataset'})
        function testEmptyDataset(testCase)
            ds_test = load("special_test_dataset.mat", 'ds_matrix_empty');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_test.ds_matrix_empty), ...
                'Empty dataset is not supported!');
        end

        function testDatasetNoAttribute(testCase)
            % Test dataset with no attribute triggers this error
            ds_test = load("special_test_dataset.mat", 'ds_matrix_no_attribute', 'ds_recti_no_attribute', 'ds_recti_xyzonly');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_test.ds_matrix_no_attribute), ...
                ['Field ''Lumerical_dataset'' does not have the ''attributes'' subfield! ' ...
                'Maybe the dataset does not have any attribute. This type of dataset is not supported.']);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_test.ds_recti_no_attribute), ...
                ['Field ''Lumerical_dataset'' does not have the ''attributes'' subfield! ' ...
                'Maybe the dataset does not have any attribute. This type of dataset is not supported.']);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_test.ds_recti_xyzonly), ...
                ['Field ''Lumerical_dataset'' does not have the ''attributes'' subfield! ' ...
                'Maybe the dataset does not have any attribute. This type of dataset is not supported.']);
        end

        function testDatasetComplexMuldimParameter(testCase)
            % Test generated real-valued parameters (matrix dataset) and
            % warning generation
            % Test multi-dimensional parameter (x,y,z) get correctly
            % converted, size correctly calculated
            testCase.verifyFail;
        end
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
                ['Field ''Lumerical_dataset'' does not have the ''attributes'' subfield! ' ...
                'Maybe the dataset does not have any attribute. This type of dataset is not supported.']);
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

        function testParameterNames1(testCase, non_text_scalar)
            % Test one interdependent parameter names
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).variable = non_text_scalar;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 names cannot be resolved!');
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).name = non_text_scalar;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 names cannot be resolved!');
        end

        function testParameterNames2(testCase, non_var_name)
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).variable = non_var_name;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The interdependent parameter set 1 names are not valid variable names!');
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

        function testParameterNumericVector(testCase, non_numeric_vector)
            % Test parameter data format (non-empty numeric vector)
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.(parameter_name) = non_numeric_vector;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                ['Parameter field ''', parameter_name, ''' data is not a numeric vector!']);
        end

        function testInterdepParameterSameLength(testCase)
            % Test interdependent parameter data same length
            ds_matrix_ = testCase.ds_matrix;
            if length(ds_matrix_.Lumerical_dataset.parameters{1}) == 1 % make sure it has interdependent parameter
                error("Test error: The testing dataset do not have an interdependent parameter set at the 1st index! Cannot proceed with this test.");
            end
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.(parameter_name)(end+1) = 5; % append a value at the end of the data
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'Interdependent parameters data do not have the same length!');
        end
    end

    methods (Test, TestTags = {'attributes'})
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Only test MATRIX dataset
        function testAttributesFieldType1(testCase, non_struct)
            % Test 'attributes' field is struct and has 'variable' and
            % 'name'fields
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes = non_struct;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'The attributes are not properly defined!');
        end

        function testAttributesFieldType2(testCase)
            ds_matrix_rm = testCase.ds_matrix;
            ds_matrix_rm.Lumerical_dataset.attributes = rmfield(ds_matrix_rm.Lumerical_dataset.attributes, 'variable');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_rm), ...
                'The attributes are not properly defined!');
            ds_matrix_rm = testCase.ds_matrix;
            ds_matrix_rm.Lumerical_dataset.attributes = rmfield(ds_matrix_rm.Lumerical_dataset.attributes, 'name');
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_rm), ...
                'The attributes are not properly defined!');
        end

        function testAttributeNames1(testCase, non_text_scalar)
            % Test attribute names
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes(1).variable = non_text_scalar;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'One or more attribute names cannot be resolved!');
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes(1).name = non_text_scalar;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'One or more attribute names cannot be resolved!');
        end

        function testAttributeNames2(testCase, non_var_name)
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes(1).variable = non_var_name;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                'One or more attribute names are not valid variable names!');
        end

        function testAttributesIllegalCharacters(testCase)
            %%%%%%%%%%%%%%%%% Same warning problem
            testCase.verifyFail;
        end

        function testAttributeInDataset(testCase)
            % Test attribute data exists
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_rm = rmfield(ds_matrix_, attribute_name);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_rm), ...
                ['Attribute field ''', attribute_name, ''' data not found!']);
        end

        function testAttributeDataNumeric(testCase, non_numeric)
            % Test attribute data format (non-empty numeric)
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_.(attribute_name) = non_numeric;
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                ['Attribute field ''', attribute_name, ''' data must be numeric!']);
        end

        function testAttributeDataSize1(testCase)
            % Test attribute data size
            % 1st dimension (xyz)
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz(1) = sz(1) + 1;
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                ['Unexpected size for attribute field ''', attribute_name, ''' data at 1st dimension!']);
            ds_recti_ = testCase.ds_recti;
            attribute_name = ds_recti_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_recti_.(attribute_name));
            sz(1) = sz(1) + 1;
            ds_recti_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_recti_), ...
                ['Unexpected size for attribute field ''', attribute_name, ''' data at 1st dimension!']);
        end

        function testAttributeDataSize2(testCase)
            % neither scalar nor vector
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz(2) = 2;
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                ['Unexpected size for attribute field ''', attribute_name, ''' data at 2nd dimension!']);
        end

        function testAttributeDataSize3(testCase)
            % Test first parameter dimension (3rd)
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz(3) = sz(3) + 1;
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                ['Unexpected size for attribute field ''', attribute_name, ''' data at dimension 3 !']);
        end

        function testAttributeDataSize4(testCase)
            % Test two many dimensions for attribute data
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz = [sz, 5];
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() loadLumDataset(ds_matrix_), ...
                ['Too many dimensions for attribute field ''', attribute_name, ''' data!']);
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
