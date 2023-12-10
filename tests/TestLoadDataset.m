classdef TestLoadDataset < matlab.unittest.TestCase
    % Test suite for LumericalDataset constructor.

    properties
        ds_matrix % normal matrix dataset
        ds_recti  % normal rectilinear dataset
        ds_special % collection of special test datasets
    end

    properties (TestParameter)
        dataset = loadTestParameter(); % both matrix and rectilinear as test parameter

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

    methods (Test, TestTags = {'NormalDataset'})
        function testNormalDatasetGeneration(testCase)
            golden = load("golden_master_normal_dataset.mat");
            testCase.verifyEqual(LumericalDataset.createObject(testCase.ds_matrix), golden.dataset_matrix);
            testCase.verifyEqual(LumericalDataset.createObject(testCase.ds_recti), golden.dataset_recti);
        end
    end

    methods (Test, TestTags = {'SpecialDataset'})
        function testEmptyDataset(testCase)
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(testCase.ds_special.ds_matrix_empty), ...
                'Empty dataset is not supported!');
        end

        function testDatasetNoAttribute(testCase)
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(testCase.ds_special.ds_matrix_no_attribute), ...
                ['Field ''Lumerical_dataset'' is missing the ''attributes'' subfield! ' ...
                'Dataset without any attribute is not supported.']);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(testCase.ds_special.ds_recti_no_attribute), ...
                ['Field ''Lumerical_dataset'' is missing the ''attributes'' subfield! ' ...
                'Dataset without any attribute is not supported.']);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(testCase.ds_special.ds_recti_xyzonly), ...
                ['Field ''Lumerical_dataset'' is missing the ''attributes'' subfield! ' ...
                'Dataset without any attribute is not supported.']);
        end
    end

    methods (Test, TestTags = {'DatasetStructure'})
        function testDatasetType(testCase, non_struct_scalar)
            % Test that non-struct-scalar dataset type throws an error.
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(non_struct_scalar), ...
                'Input dataset must be a struct scalar!');
        end

        function testDatasetFields1(testCase, dataset)
            % Test field 'Lumerical_dataset' exists
            ds_rm = rmfield(dataset, "Lumerical_dataset");
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_rm), ...
                'Input dataset does not have the field ''Lumerical_dataset''!');
        end

        function testDatasetFields2(testCase, dataset, non_struct_scalar)
            % Test field 'Lumerical_dataset' is struct scalar
            ds_wrong = dataset;
            ds_wrong.Lumerical_dataset = non_struct_scalar;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_wrong), ...
                'Field ''Lumerical_dataset'' is not a struct scalar!');
        end

        function testDatasetFields3(testCase, dataset)
            % Test field 'Lumerical_dataset' subfield 'attributes'
            ds_attributes_rm = dataset;
            ds_attributes_rm.Lumerical_dataset = rmfield(ds_attributes_rm.Lumerical_dataset, "attributes");
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_attributes_rm), ...
                ['Field ''Lumerical_dataset'' is missing the ''attributes'' subfield! ' ...
                'Dataset without any attribute is not supported.']);
        end

        function testDatasetFields4(testCase, dataset)
            % Test field 'Lumerical_dataset' subfield 'parameters'
            ds_parameters_rm = dataset;
            ds_parameters_rm.Lumerical_dataset = rmfield(ds_parameters_rm.Lumerical_dataset, "parameters");
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_parameters_rm), ...
                'Field ''Lumerical_dataset'' is missing the ''parameters'' subfield!');
        end

        function testWrongDatasetType1(testCase)
            % Test directly constructing wrong dataset type gives an error
            testCase.verifyErrorMessage(@() MatrixDataset(testCase.ds_recti), ...
                'This is not a matrix dataset!');
        end

        function testWrongDatasetType2(testCase)
            testCase.verifyErrorMessage(@() RectilinearDataset(testCase.ds_matrix), ...
                'This is not a rectilinear dataset!');
        end
    end

    methods (Test, TestTags = {'LoadRectilinearError'})
        function testRectiLabel(testCase, non_geometry)
            % Test wrong 'geometry' label
            ds_recti_ = testCase.ds_recti;
            ds_recti_.Lumerical_dataset.geometry = non_geometry;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_), ...
                'Wrong label in ''lum_dataset.geometry'' for the rectilinear dataset!');
        end

        function testMissingXYZField(testCase)
            ds_recti_rm_x = rmfield(testCase.ds_recti, 'x');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_rm_x), ...
                'No x data in the rectilinear dataset!');
            ds_recti_rm_y = rmfield(testCase.ds_recti, 'y');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_rm_y), ...
                'No y data in the rectilinear dataset!');
            ds_recti_rm_z = rmfield(testCase.ds_recti, 'z');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_rm_z), ...
                'No z data in the rectilinear dataset!');
        end

        function testXYZDataNonNumeric(testCase, non_numeric)
            ds_recti_ = testCase.ds_recti;
            ds_recti_.x = non_numeric;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_), ...
                'x data must be a numeric vector!');
            ds_recti_ = testCase.ds_recti;
            ds_recti_.y = non_numeric;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_), ...
                'y data must be a numeric vector!');
            ds_recti_ = testCase.ds_recti;
            ds_recti_.z = non_numeric;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_), ...
                'z data must be a numeric vector!');
        end

        function testXYZDataMuldim(testCase)
            testCase.verifyWarning(@() LumericalDataset.createObject(testCase.ds_special.ds_recti_muldimxyz), ...
                "PositionalVector:DataIsMuldim");
        end

        function testXYZMuldimConvertedResult(testCase)
            warning('off', "PositionalVector:DataIsMuldim");
            dataset_recti_muldimxyz = LumericalDataset.createObject(testCase.ds_special.ds_recti_muldimxyz);
            warning('on', "PositionalVector:DataIsMuldim");
            testCase.verifyEqual(dataset_recti_muldimxyz.y, testCase.ds_special.ds_recti_muldimxyz.y(:));
        end

        function testXYZDataComplex(testCase)
            testCase.verifyWarning(@() LumericalDataset.createObject(testCase.ds_special.ds_recti_complexxyz), ...
                "PositionalVector:DataIsComplex");
        end

        function testXYZComplexConvertedResult(testCase)
            warning('off', "PositionalVector:DataIsComplex");
            dataset_recti_complexxyz = LumericalDataset.createObject(testCase.ds_special.ds_recti_complexxyz);
            warning('on', "PositionalVector:DataIsComplex");
            testCase.verifyEqual(dataset_recti_complexxyz.y, real(testCase.ds_special.ds_recti_complexxyz.y));
        end

        function testXYZDataNaNInf(testCase)
            ds_recti_ = testCase.ds_recti;
            ds_recti_.y(1) = NaN;
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_recti_), ...
                "PositionalVector:DataHasInvalidElement");
            ds_recti_ = testCase.ds_recti;
            ds_recti_.y(1) = -Inf;
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_recti_), ...
                "PositionalVector:DataHasInvalidElement");
        end

        function testRectiParameterNameClash(testCase)
            % Test rectilinear dataset parameters having name x, y or z
            % will result in an error
            ds_recti_ = testCase.ds_recti;
            ds_recti_.Lumerical_dataset.attributes(1).variable = 'x';
            ds_recti_.Lumerical_dataset.attributes(1).name = 'x';
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_), ...
                'Parameters, positional vectors (x,y,z) and attributes have duplicate names! Rejected!');
        end
    end

    methods (Test, TestTags = {'LoadParametersError'})
        function testParametersType(testCase, dataset, non_cell_column)
            % 'parameters' field should be a cell column vector
            ds = dataset;
            ds.Lumerical_dataset.parameters = non_cell_column;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds), ...
                'Field ''Lumerical_dataset.parameters'' should be a cell column vector!');
        end

        %%%%%%% Only test one parameter in MATRIX dataset
        function testParameterContent1(testCase, non_struct)
            % Test the cell in 'parameters' should be a struct
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1} = non_struct;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
        end

        function testParameterContent2(testCase)
            % Test the 'variable' and 'name' subfields
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1} = rmfield(ds_matrix_.Lumerical_dataset.parameters{1}, 'variable');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
            ds_matrix_.Lumerical_dataset.parameters{1} = rmfield(ds_matrix_.Lumerical_dataset.parameters{1}, 'name');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'The interdependent parameter set 1 is not properly defined!');
        end

        function testParameterNames1(testCase, non_text_scalar)
            % Test interdependent parameter names should be text
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).variable = non_text_scalar;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'The interdependent parameter set 1 names cannot be resolved!');
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).name = non_text_scalar;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'The interdependent parameter set 1 names cannot be resolved!');
        end

        function testParameterNames2(testCase, non_var_name)
            % Test interdependent parameter names should be valid variable
            % names
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).variable = non_var_name;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'The interdependent parameter set 1 names are not valid variable names!');
        end

        function testDuplicateParameterName(testCase)
            % Test no two parameters should have the same name, otherwise
            % throw parameter data not found error
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.Lumerical_dataset.parameters{2}(1).variable = parameter_name;
            ds_matrix_.Lumerical_dataset.parameters{2}(1).name = parameter_name;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'Parameters and attributes have duplicate names! Rejected!');
        end

        function testParameterNameIllegalCharacters(testCase)
            % When original parameter names from Lumerical contains illegal
            % characters, display a warning
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.Lumerical_dataset.parameters{1}(1).name = char(parameter_name + " ");
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_matrix_), ...
                "Parameter:NameHasIllegalCharacters");
        end

        function testParameterInDataset(testCase)
            % Test the specific parameter data exists
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_rm = rmfield(ds_matrix_, parameter_name);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_rm), ...
                ['Parameter field ''', parameter_name, ''' data not found!']);
        end

        function testParameterDataNumericVector(testCase, non_numeric_vector)
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.(parameter_name) = non_numeric_vector;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                ['Parameter field ''', parameter_name, ''' data is not a numeric vector!']);
        end

        function testParameterDataComplex(testCase)
            % Test warning from complex-valued parameter data in both
            % matrix and rectilinear dataset. xyz not tested here. No way
            % to test multidimensional parameter data either, because it is
            % already converted during MATLAB export.
            testCase.verifyWarning(@() LumericalDataset.createObject(testCase.ds_special.ds_matrix_complex), ...
                "Parameter:DataIsComplex");
            testCase.verifyWarning(@() LumericalDataset.createObject(testCase.ds_special.ds_recti_complex), ...
                "Parameter:DataIsComplex");
        end

        function testParameterComplexConvertedResult(testCase)
            warning('off', "Parameter:DataIsComplex"); % temporarily disable warning
            dataset_matrix_complex = LumericalDataset.createObject(testCase.ds_special.ds_matrix_complex);
            warning('on', "Parameter:DataIsComplex");
            testCase.verifyEqual(dataset_matrix_complex.parameters{1, 2}(:, 2), real(testCase.ds_special.ds_matrix_complex.linear_complex));
        end

        function testParameterDataNaNInf(testCase)
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.(parameter_name)(1) = NaN;
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_matrix_), ...
                "Parameter:DataHasInvalidElement");
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.(parameter_name)(1) = -Inf;
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_matrix_), ...
                "Parameter:DataHasInvalidElement");
        end

        function testInterdepParameterSameLength(testCase)
            % Test each interdependent parameter data has same length
            ds_matrix_ = testCase.ds_matrix;
            parameter_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.(parameter_name)(end+1) = 5; % append a value at the end of the data
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'Interdependent parameters data do not have the same length!');
        end
    end

    methods (Test, TestTags = {'LoadAttributesError'})
        %%%%%%% Only test one attribute in MATRIX dataset
        function testAttributesFieldType1(testCase, non_struct)
            % Test 'attributes' field is struct
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes = non_struct;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'The attributes are not properly defined!');
        end

        function testAttributesFieldType2(testCase)
            % Test 'attributes' field has 'variable' and 'name' fields
            ds_matrix_rm = testCase.ds_matrix;
            ds_matrix_rm.Lumerical_dataset.attributes = rmfield(ds_matrix_rm.Lumerical_dataset.attributes, 'variable');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_rm), ...
                'The attributes are not properly defined!');
            ds_matrix_rm = testCase.ds_matrix;
            ds_matrix_rm.Lumerical_dataset.attributes = rmfield(ds_matrix_rm.Lumerical_dataset.attributes, 'name');
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_rm), ...
                'The attributes are not properly defined!');
        end

        function testAttributeNames1(testCase, non_text_scalar)
            % Test attribute names should be text
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes(1).variable = non_text_scalar;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'One or more attribute names cannot be resolved!');
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes(1).name = non_text_scalar;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'One or more attribute names cannot be resolved!');
        end

        function testAttributeNames2(testCase, non_var_name)
            % Test attribute names should be valid variable names
            ds_matrix_ = testCase.ds_matrix;
            ds_matrix_.Lumerical_dataset.attributes(1).variable = non_var_name;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'One or more attribute names are not valid variable names!');
        end

        function testDuplicateAttributeName(testCase)
            % Test no two attributes should have the same name, otherwise
            % throw attribute data not found error
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_.Lumerical_dataset.attributes(2).variable = attribute_name;
            ds_matrix_.Lumerical_dataset.attributes(2).name = attribute_name;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'Parameters and attributes have duplicate names! Rejected!');
        end

        function testAttributeNameClashWithParameter(testCase)
            % Test an attribute with same name as a parameter, should throw
            % attribute data not found error
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.parameters{1}(1).variable;
            ds_matrix_.Lumerical_dataset.attributes(1).variable = attribute_name;
            ds_matrix_.Lumerical_dataset.attributes(1).name = attribute_name;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                'Parameters and attributes have duplicate names! Rejected!');
        end

        function testAttributesIllegalCharacters(testCase)
            % When original attribute names from Lumerical contains illegal
            % characters, display a warning
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_.Lumerical_dataset.attributes(1).name = char(attribute_name + " ");
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_matrix_), ...
                "Attribute:NameHasIllegalCharacters");
        end

        function testAttributeInDataset(testCase)
            % Test attribute data exists
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_rm = rmfield(ds_matrix_, attribute_name);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_rm), ...
                ['Attribute field ''', attribute_name, ''' data not found!']);
        end

        function testAttributeDataNumeric(testCase, non_numeric)
            % Attribute data can be multi-dimensional or complex
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_.(attribute_name) = non_numeric;
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                ['Attribute field ''', attribute_name, ''' data must be numeric!']);
        end

        function testAttributeDataNaNInf(testCase)
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_.(attribute_name)(1) = NaN;
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_matrix_), ...
                "Attribute:DataHasInvalidElement");
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            ds_matrix_.(attribute_name)(1) = -Inf;
            testCase.verifyWarning(@() LumericalDataset.createObject(ds_matrix_), ...
                "Attribute:DataHasInvalidElement");
        end

        function testAttributeDataSize1(testCase)
            % Test attribute data size at 1st dimension matches dataset type
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz(1) = sz(1) + 1; % increase size at 1st dimension by 1
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                ['In matrix dataset, the size for attribute field ''', ...
                attribute_name, ''' data at 1st dimension should be 1!']);
            ds_recti_ = testCase.ds_recti;
            attribute_name = ds_recti_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_recti_.(attribute_name));
            sz(1) = sz(1) + 1;
            ds_recti_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_recti_), ...
                ['In rectilinear dataset, the size for attribute field ''', ...
                attribute_name, ''' data at 1st dimension should match the product of x, y and z dimensions!']);
        end

        function testAttributeDataSize2(testCase)
            % Test attribute data size at 2nd dimension should be either 1
            % (scalar) or 3 (vector)
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz(2) = 2;
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                ['Unexpected size for attribute field ''', attribute_name, ''' data at 2nd dimension!']);
        end

        function testAttributeDataSize3(testCase)
            % Test attribute data size at 3rd dimension and beyond matches
            % corresponding parameter data length
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz(3) = sz(3) + 1; % test 3rd dimension
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                ['Unexpected size for attribute field ''', attribute_name, ''' data at dimension 3 !']);
        end

        function testAttributeDataSize4(testCase)
            % Test too many dimensions for attribute data
            ds_matrix_ = testCase.ds_matrix;
            attribute_name = ds_matrix_.Lumerical_dataset.attributes(1).variable;
            sz = size(ds_matrix_.(attribute_name));
            sz = [sz, 5]; % append a dimension of length 5 at the end
            ds_matrix_.(attribute_name) = rand(sz);
            testCase.verifyErrorMessage(@() LumericalDataset.createObject(ds_matrix_), ...
                ['Too many dimensions for attribute field ''', attribute_name, ''' data!']);
        end
    end

    methods (TestClassSetup)
        function loadDatasetOnce(testCase)
            % Load dataset into test class properties
            data = load("normal_test_dataset.mat");
            testCase.ds_matrix = data.ds_matrix;
            testCase.ds_recti = data.ds_recti;
            testCase.ds_special = load("special_test_dataset.mat");
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
                    'Incorrect error message was thrown! See below for details.');
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
