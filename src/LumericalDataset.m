classdef (Abstract) LumericalDataset < matlab.mixin.Copyable
    % Lumerical dataset base class

    % Parameter and attribute data array cannot be empty.
    %
    % Missing required data will result in an error. However, more
    % (unnecessary) data will not trigger any error, and will simply be
    % ignored.
    %
    % No two names in parameters, attributes and positional vectors (x,y,z)
    % can be the same.
    %
    % If a parameter or positional vector (x,y,z) data contains duplicate
    % elements or is not strictly monotonic, no error/warning will be
    % generated. However, 2D plot cannot be made. If it contains invalid
    % element (NaN or Inf), a warning will be issued, and no plot can be
    % made.
    %
    % If a parameter or positional vector (x,y,z) has complex data, the
    % real part will be taken during conversion.
    %
    % A dataset without any parameter is not accepted. A dataset without
    % any attribute is not accepted.

    properties (SetAccess = protected)
        parameters
        num_parameters
        parameters_indexes % selected the current parameter slice indexes

        attributes % common data but depends on the type of subclass
        attributes_component % select the current attribute component (x, y, z)
        % NaN-scalar 0-magnitude 1-x 2-y 3-z
        num_attributes
    end

    methods (Static)
        function obj = createObject(lum_dataset)
            % Factory method
            % Initialize the dataset object and choose the correct subclass
            % to invoke based on the type of the dataset

            % No parameter or attribute (empty matrix dataset)
            if isequal(lum_dataset, 'data type not supported')
                error("Empty dataset is not supported!");
            end

            % Decide which class to instantiate
            dataset_type = LumericalDataset.parseDatasetStructure(lum_dataset);
            if dataset_type == "rectilinear"
                obj = RectilinearDataset(lum_dataset);
            elseif dataset_type == "matrix"
                obj = MatrixDataset(lum_dataset);
            end
        end

        function converted_obj = createObjectFromMat(mat_name, varargin)
            % Create dataset from the .mat file
            % By default, attempt to convert all structs
            % Can specify specific files to convert

            % Check name of mat_file, check varargin
            try
                info = whos(matfile(mat_name));
            catch ME
                error("File name is invalid!");
            end
            % Can only have one dataset in one file
            if isempty(info)
                error("File is not found is empty!");
            end
            % Need to find how many structs are in there
            for i = 1:numel(info)
                variable_name = info(i).name;
                if isempty(varargin)
                    data = load(mat_name, variable_name);
                    converted_obj.(variable_name) = ...
                        LumericalDataset.createObject(data.(variable_name));
                end
                % varargin should be a cell array of char vector
                if strcmp(variable_name, varargin)
                    data = load(mat_name, variable_name);
                    converted_obj.(variable_name) = ...
                        LumericalDataset.createObject(data.(variable_name));
                end
            end
        end
        % Get the dataset in the current file
    end

    methods
        function obj = LumericalDataset()
            % Load Lumerical exported MATLAB dataset into custom class object
            % Empty. Definitions in derived classes.
        end

        function showInformation(obj)
            % Display the information for an overview of this dataset
            % Print attributes information
            fprintf("%d attribute(s):\n", obj.num_attributes);
            attributes_fields = fieldnames(obj.attributes);
            max_name_length = max([cellfun(@length, attributes_fields); 4]); % max name length
            fprintf("| %-*s |  Type  | Plot component |\n", max_name_length, 'Name');
            fprintf("+%s+--------+----------------+\n", repmat('-', 1, max_name_length + 2));
            for i = 1:obj.num_attributes
                field = attributes_fields{i};
                if isnan(obj.attributes_component.(field)) % scalar or vector
                    type = "scalar";
                else
                    type = "vector";
                end
                % Print attribute line
                if isnan(obj.attributes_component.(field))
                    component = "              ";
                else
                    switch obj.attributes_component.(field)
                        case 0
                            component = "  magnitude   ";
                        case 1
                            component = "      x       ";
                        case 2
                            component = "      y       ";
                        case 3
                            component = "      z       ";
                    end
                end
                fprintf("| %-*s | %s | %s |\n", max_name_length, field, type, component);
            end

            % Print parameters information
            fprintf("%d parameter(s):\n", obj.num_parameters);
            params = obj.parameters(:, 1);
            max_digits_n_points = max(cellfun(@(x) length(num2str(x)), obj.parameters(:, 3))); % max number of points length
            max_digits_slice_index = max([arrayfun(@(x) length(num2str(x)), obj.parameters_indexes); 11]);
            fprintf("| %-*s | Slice index | Name\n", max_digits_n_points + 9, 'N_points');
            fprintf("+%s+%s+------\n", repmat('-', 1, max_digits_n_points + 11), ...
                repmat('-', 1, max_digits_slice_index + 2));
            for i = 1:obj.num_parameters
                fprintf("| %*d point(s) | %*d | ", max_digits_n_points, obj.parameters{i, 3}, ...
                    max_digits_slice_index, obj.parameters_indexes(i)); % print parameter line
                for param = params{i} % independent parameter names
                    fprintf("%s ", param);
                end
                fprintf('\n');
            end
        end

        function result = getParameterData(obj, parameter_name) % non-virtual
            % Get the data of a parameter
            para_loc = obj.iCheckAndFindParameter(parameter_name);
            result = obj.parameters{para_loc(1), 2}(:, para_loc(2));
        end

        function result = getAttributeData(obj, attribute_name) % non-virtual
            % Get the data for an attribute
            obj.iCheckAttributeExist(attribute_name);
            result = obj.attributes.(attribute_name);
        end

        function setAttributeComponent(obj, attribute_name, component) % non-virtual
            % set the component (x,y,z,magnitude) for an attribute
            obj.iCheckAttributeExist(attribute_name);
            if isnan(obj.attributes_component.(attribute_name))
                error("Scalar attribute cannot be modified!");
            end

            LumericalDataset.validateTextScalar(component, "Component name must be text!");
            switch component
                case "x"
                    obj.attributes_component.(attribute_name) = 1;
                case "y"
                    obj.attributes_component.(attribute_name) = 2;
                case "z"
                    obj.attributes_component.(attribute_name) = 3;
                case "magnitude"
                    obj.attributes_component.(attribute_name) = 0;
                otherwise
                    error("Invalid component name!");
            end
        end

        function hPlot = plotData1D(obj, parameter_name, attribute_name, ax)  % non-virtual
            % This function makes a 1D plot based on the parameter name and
            % the attribute name.

            % If 'ax' is a valid axes handle, plot on that axes object
            if nargin < 4
                ax = gca(); % plot on current axis
            elseif  ~(isscalar(ax) && isgraphics(ax, 'axes')) % valid axes handle
                error("'ax' argument must be an valid axes handle!");
            end

            % Calls the respective polymorphic function for each object
            [xdata, ydata] = obj.getPlot1DData(parameter_name, attribute_name);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Problem (also for 2D):
            % complex numbers? real?
            % Unit conversion and scaling?
            hPlot = plot(ax, xdata, real(ydata));
            xlabel(ax, parameter_name, 'Interpreter', 'none');
            ylabel(ax, attribute_name, 'Interpreter', 'none');
            set(ax, 'FontSize', 16);
        end

        function [hSurf, hClb] = plotData2D(obj, parameter1_name, parameter2_name, attribute_name, ax)  % non-virtual
            % This function makes a 2D plot based on the parameter names
            % and the attribute name.

            % If 'ax' is a valid axes handle, plot on that axes object
            if nargin < 5
                ax = gca(); % plot on current axis
            elseif  ~(isscalar(ax) && isgraphics(ax, 'axes')) % valid axes handle
                error("'ax' argument must be an valid axes handle!");
            end

            % Calls the respective polymorphic function for each object
            [xdata, ydata, zdata] = obj.getPlot2DData(parameter1_name, parameter2_name, attribute_name);

            % Throw an error is xdata or ydata (parameters) is not monotonic
            if ~LumericalDataset.isRealVectorMonotonic(xdata)
                error("x data is not monotonic! Cannot make 2D plot.");
            end
            if ~LumericalDataset.isRealVectorMonotonic(ydata)
                error("y data is not monotonic! Cannot make 2D plot.");
            end

            % Adjust data to make true "2D plot"
            % If length(xdata) == 1 or length(ydata) == 1, this program
            % will not report an error but will have nothing plotted
            xdata = ([xdata(1); xdata(:)] + [xdata(:); xdata(end)])/2;
            ydata = ([ydata(1); ydata(:)] + [ydata(:); ydata(end)])/2;
            zdata_new = zeros(length(ydata), length(xdata));
            zdata_new(1:end-1, 1:end-1) = zdata;
            zdata = zdata_new;

            hSurf = surface(ax, xdata*1e9, ydata*1e9, real(zdata), 'EdgeColor', 'none');
            xlim(ax, [min(xdata*1e9), max(xdata*1e9)]);
            ylim(ax, [min(ydata*1e9), max(ydata*1e9)]);
            xlabel(ax, parameter1_name + " (nm)", 'Interpreter', 'none');
            ylabel(ax, parameter2_name + " (nm)", 'Interpreter', 'none');
            colormap(ax, 'jet');
            hClb = colorbar;
            set(ax, 'FontSize', 16);
            set(ax, 'Layer', 'top');
            box(ax, 'on');
        end
    end

    methods (Abstract)
        obj = setParameterSliceIndex(obj, varargin);
        [xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name);
        [xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name);
        [x, y, z, data] = getPlot3DData(obj, parameter1_name, parameter2_name, parameter3_name, attribute_name);
        new_obj = removeDimensions(obj, varargin);
        new_obj = mergeDataset(obj, other_obj, varargin);
    end

    methods (Access = protected)
        function p = iAddParametersToParser(obj, p)
            % Add all parameters and their range validations to the parser
            for i = 1:size(obj.parameters, 1)
                for para = obj.parameters{i, 1}
                    p.addParameter(para, NaN, @(x) LumericalDataset.validateIndex(x, obj.parameters{i, 3}));
                end
            end
        end

        function iAnalyzeAndSetParsedParameter(obj, p)
            % You can supply the same parameter multiple times and it will
            % keep the last occurance of it
            for i = 1:size(obj.parameters, 1)
                interdep_indexes = [];
                for para = obj.parameters{i, 1}
                    interdep_indexes(end+1) = p.Results.(para);  % retrieve
                end
                unique_index = unique(interdep_indexes(~isnan(interdep_indexes)));
                if numel(unique_index) > 1
                    throwAsCaller(MException('', "Interdependent parameters for '" + para + "' should select the same index!"));
                elseif numel(unique_index) == 1
                    obj.parameters_indexes(i) = unique_index;
                end
            end
        end

        function [para_slice_indexes, para_value_list, para_remove_indexes] = ...
                iGenerateParametersSliceIndexAndData(obj, para_value_list, varargin)
            % Generate index list that only keep passed in parameters
            % (para_slice_indexes), or index list that only removes passed in
            % parameters (para_remove_indexes), and update the information of
            % the parameters passed in (para_value_list).
            % Certain elements in varargin can be empty, meaning that it
            % has already been parsed in the last step. In this case, the
            % corresponding element in para_value_list should already have
            % a value.
            % ---para_value_list---
            % First column: actual data for the parameter
            % Second column: the number of the parameter (1:num_parameters)
            try
                para_slice_indexes = num2cell(obj.parameters_indexes);
                for i = 1:length(varargin)
                    if isempty(varargin{i}) % already parsed
                        continue;
                    end
                    para_loc = obj.iCheckAndFindParameter(varargin{i});
                    para_slice_indexes{para_loc(1)} = ':';
                    % Extract found parameter data
                    para_value_list{i, 1} = obj.parameters{para_loc(1), 2}(:, para_loc(2));
                    para_value_list{i, 2} = para_loc(1);
                end

                para_remove_indexes = num2cell(obj.parameters_indexes);
                for i = 1:length(para_slice_indexes)
                    if isnumeric(para_slice_indexes{i})
                        para_remove_indexes{i} = ':';
                    end
                end
            catch ME
                throwAsCaller(ME)
            end
        end

        function location = iCheckAndFindParameter(obj, parameter_name)
            % Check if a parameter name exists and return location
            % Location: [row, column (within interdep)]
            try
                LumericalDataset.validateTextScalar(parameter_name, "Parameter name must be text scalar!");
            catch ME
                throwAsCaller(ME);
            end

            tf_found = false;
            location = zeros(1, 2);
            parameter_names = obj.parameters(:, 1);
            for i = 1:length(parameter_names)
                if ismember(parameter_name, parameter_names{i})
                    location(1) = i;
                    location(2) = find(parameter_names{i} == parameter_name);
                    tf_found = true;
                    break;
                end
            end
            if ~tf_found
                throwAsCaller(MException('', "Parameter name '" + parameter_name + "' not found!"));
            end
        end

        function iCheckParameterExist(obj, parameter_name)
            % Check if a parameter name exists
            try
                LumericalDataset.validateTextScalar(parameter_name, "Parameter name must be text scalar!");
            catch ME
                throwAsCaller(ME);
            end

            tf_found = false;
            parameter_names = obj.parameters(:, 1);
            for i = 1:length(parameter_names)
                if ismember(parameter_name, parameter_names{i})
                    tf_found = true;
                    break;
                end
            end
            if ~tf_found
                throwAsCaller(MException('', "Parameter name '" + parameter_name + "' not found!"));
            end
        end

        function iCheckAttributeExist(obj, attribute_name)
            % Check if an attribute name exists
            try
                LumericalDataset.validateTextScalar(attribute_name, "Attribute name must be text scalar!");
                LumericalDataset.validateFieldInStruct(obj.attributes, attribute_name, "This attribute '" + attribute_name + "' is not found!");
            catch ME
                throwAsCaller(ME);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deprecated, but this one takes in the value instead of the
        % indexes
        function parameters_value_list = parseParameterList(obj, parameters_value_list, arglist)
            % Check varargin (should be parameter name-value pair)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % In this version, only one interdependent parameter is allowed
            % to be declared. Further, pass in values rather than indexes.
            n = length(arglist);
            if mod(n, 2)
                error('Parameter name-value not in pairs!');
            end
            for i = 1:2:n
                parameter = arglist{i};
                value = arglist{i+1};
                para_loc = obj.iCheckAndFindParameter(parameter);
                if ~isnumeric(value) || ~isscalar(value)
                    error("The value for this parameter '" + parameter + "' is not a number!");
                end

                % Check if the parameter (or other interdep) has already been defined
                if ~isnan(parameters_value_list(para_loc(1), 1))
                    disp(parameter);
                    error('Repeat definition of this parameter!');
                end

                % Find the index of the closest value
                [~, value_index] = min(abs(obj.parameters{para_loc(1), 2}(:, para_loc(2)) - value));
                parameters_value_list(para_loc(1), 1) = para_loc(2);
                parameters_value_list(para_loc(1), 2) = value_index;
            end

            % Any not specified parameters are default to the first value
            parameters_value_list(isnan(parameters_value_list(:, 1)), :) = 1;
        end
    end

    methods (Static, Access = protected)
        % Helper functions to be shared with base and derived classes
        function validateTextScalar(input, errmsg)
            % Validate input as text scalar and throw if not
            if ~(ischar(input) && isrow(input)) && ~isStringScalar(input)
                throwAsCaller(MException('', errmsg));
            end
        end

        function validateStructScalar(input, errmsg)
            % Validate input as struct scalar and throw if not
            if ~isstruct(input) || ~isscalar(input)
                throwAsCaller(MException('', errmsg));
            end
        end

        function validateFieldInStruct(struct_in, field_in, errmsg)
            % Validate field exists in a struct
            if ~isfield(struct_in, field_in)
                throwAsCaller(MException('', errmsg));
            end
        end

        function validateNonEmptyNumericVector(input, errmsg)
            % Validate input as non-empty vector (N-by-1 or 1-by-N)
            if ~(isnumeric(input) && isvector(input) && ~isempty(input))
                throwAsCaller(MException('', errmsg));
            end
        end

        function validateIndex(idx, max_val)
            if ~isnumeric(idx) || floor(idx) ~= idx % test integer
                throwAsCaller(MException('', "Must be an integer!"));
            end
            if idx < 1 || idx > max_val
                throwAsCaller(MException('', "Must be between 1 and " + max_val + "!"));
            end
        end

        function result = sliceThroughAttributeVectorDim(attribute_data, component_index)
            % Slice off attribute data matrix based on the selected vector
            % component (NaN-scalar, 0-magnitude, 1-x, 2-y, 3-z)
            % on the second dimension
            if isnan(component_index) % NaN-scalar
                result = attribute_data;
            elseif component_index == 0 % 0-magnitude
                result = vecnorm(attribute_data, 2, 2);
            else % 1, 2, 3-x, y, z
                idx = repmat({':'}, 1, ndims(attribute_data)); % ugly
                idx{2} = component_index;
                result = attribute_data(idx{:});
            end
        end

        function tf = isequalWithinTol(first, second, absTol, relTol)
            % Compare two numeric arrays equal within a tolerance limit
            % Absolute OR relative tolerance satisfied
            if nargin <= 3
                relTol = 1e-10;
            end
            if nargin == 2
                absTol = 1e-12;
            end
            if ~isequal(size(first), size(second))
                tf = false;
                return;
            end
            absolute_error = abs(first - second);
            relative_error = abs((first - second) ./ first);
            tf = all((absolute_error <= absTol) | (relative_error <= relTol), 'all');
        end

        function tf = isRealVectorMonotonic(vec)
            % Returns true if a real-valued non-empty vector is strictly
            % monotonic (increasing or decreasing)
            tf =  all(diff(vec) > 0) || all(diff(vec) < 0);
        end
    end

    methods (Static, Access = protected)
        % Helper functions only shared with base class
        function dataset_type = parseDatasetStructure(lum_dataset)
            % Initially pase the dataset structure and determine dataset
            % type (matrix or rectilinear)

            % Input dataset has to be a struct scalar
            LumericalDataset.validateStructScalar(lum_dataset, "Input dataset must be a struct scalar!");

            % Check field: Lumerical_dataset
            LumericalDataset.validateFieldInStruct(lum_dataset, 'Lumerical_dataset', "Input dataset does not have the field 'Lumerical_dataset'!");
            LumericalDataset.validateStructScalar(lum_dataset.Lumerical_dataset, "Field 'Lumerical_dataset' is not a struct scalar!");
            % Check 'attribute' field (check the contents later)
            LumericalDataset.validateFieldInStruct(lum_dataset.Lumerical_dataset, 'attributes', "Field 'Lumerical_dataset' does not have the 'attributes' subfield! " + ...
                "Maybe the dataset does not have any attribute. This type of dataset is not supported.");
            % Check 'parameters' field (check the contents later)
            LumericalDataset.validateFieldInStruct(lum_dataset.Lumerical_dataset, 'parameters', "Field 'Lumerical_dataset' does not have the 'parameters' subfield!")

            % Determine dataset type (matrix or rectilinear) using
            % 'geometry' field
            if isfield(lum_dataset.Lumerical_dataset, 'geometry')
                dataset_type = 'rectilinear';
                if ~isequal(lum_dataset.Lumerical_dataset.geometry, "rectilinear")
                    error("Wrong label in 'lum_dataset.geometry' for the rectilinear dataset!");
                end
            else
                dataset_type = 'matrix';
            end
        end

        function [x, y, z, xyz_prod_size] = parseXYZ(lum_dataset)
            % Parse positional vectors (x,y,z) for a rectilinear dataset

            % Check if x, y and z data exist in the dataset
            LumericalDataset.validateFieldInStruct(lum_dataset, 'x', "No x data in the rectilinear dataset!");
            LumericalDataset.validateFieldInStruct(lum_dataset, 'y', "No y data in the rectilinear dataset!");
            LumericalDataset.validateFieldInStruct(lum_dataset, 'z', "No z data in the rectilinear dataset!");

            % Load x,y,z to xyz and remove those field from the dataset
            % Check x,y,z data.
            for axis = ['x', 'y', 'z']
                data.(axis) = lum_dataset.(axis);
                if ~isnumeric(data.(axis)) || isempty(data.(axis))
                    error(axis + " data must be a numeric vector!");
                end
                if ~isvector(data.(axis)) % 2+ dimensional matrix?
                    warning("PositionalVector:DataIsMuldim", ...
                        "Parameter " + axis + " is multi-dimensional! Stretched to one dimension!");
                end
                if any(imag(data.(axis))) % has imaginary part?
                    warning("PositionalVector:DataIsComplex", ...
                        "Parameter " + axis + " is complex! Takes the real part and proceed.");
                    data.(axis) = real(data.(axis));
                end
                if any(isnan(data.(axis)), 'all') || any(isinf(data.(axis)), 'all') % real part has NaN or Inf?
                    warning("PositionalVector:DataHasInvalidElement", "Parameter '" + axis + ...
                        "' data contains invalid (NaN or Inf) elements! Something to keep in mind.");
                end
            end
            x = data.x(:); % (vectorize) convert to column vector
            y = data.y(:);
            z = data.z(:);
            xyz_prod_size = length(x) * length(y) * length(z);
        end

        function parameters_info = parseParametersName(lum_dataset)
            % Parse parameter field names of the dataset

            % Load all parameters names and organize them
            parameters = lum_dataset.Lumerical_dataset.parameters;
            if ~(iscell(parameters) && iscolumn(parameters))
                error("Field 'Lumerical_dataset.parameters' should be a cell column vector!");
            end
            parameters_info = cell(length(parameters), 3);
            for i = 1:length(parameters)
                parameter = parameters{i};
                % Check struct
                if ~(isstruct(parameter) && isfield(parameter, 'variable') && isfield(parameter, 'name'))
                    error("The interdependent parameter set " + i + " is not properly defined!");
                end

                % Initialize variables
                parameter_names = strings(1, length(parameter));

                % Retrieve the interdependent parameter set
                for j = 1:length(parameter)
                    LumericalDataset.validateTextScalar(parameter(j).variable, "The interdependent parameter set " + i + " names cannot be resolved!");
                    interdep_parameter_name = parameter(j).variable;
                    parameter_names(j) = interdep_parameter_name; % char array also works
                    LumericalDataset.validateTextScalar(parameter(j).name, "The interdependent parameter set " + i + " names cannot be resolved!");
                    if ~isvarname(interdep_parameter_name) % must be legal variable name
                        error("The interdependent parameter set " + i + " names are not valid variable names!");
                    end
                    % Illegal characters in the names are converted to '_'
                    if ~isequal(parameter(j).variable, parameter(j).name)
                        warning("Parameter:NameHasIllegalCharacters", "Parameter name '" + parameter(j).name + ...
                            "' contains illegal characters. Converted to '" + parameter(j).variable + "'.");
                    end
                    LumericalDataset.validateFieldInStruct(lum_dataset, interdep_parameter_name, ...
                        "Parameter field '" + interdep_parameter_name + "' data not found!");
                end
                parameters_info{i, 1} = parameter_names;
            end
        end

        function parameters_info = parseParametersData(lum_dataset, parameters_info)
            % Parse parameter field data of the dataset

            for i = 1:size(parameters_info, 1)
                parameter_values = cell(1, length(parameters_info{i, 1}));
                parameter_length = nan(1, length(parameters_info{i, 1}));
                for j = 1:length(parameters_info{i, 1})
                    interdep_parameter_name = parameters_info{i, 1}(j);
                    value = lum_dataset.(interdep_parameter_name);
                    % parameter will always be non-empty N-by-1 column vector. Check that.
                    % In the check, can relax the condition to "non-empty vectors".
                    LumericalDataset.validateNonEmptyNumericVector(value, ... %% instead check iscolumn?
                        "Parameter field '" + interdep_parameter_name + "' data is not a numeric vector!");
                    value = value(:); % convert to column vector, if applicable
                    % Remove complex portion
                    if any(imag(value)) % has imaginary part?
                        warning("Parameter:DataIsComplex", "Parameter '" + interdep_parameter_name + ...
                            "' data is complex! Takes the real part and proceed.");
                        value = real(value);
                    end
                    if any(isnan(value)) || any(isinf(value)) % real part has NaN or Inf?
                        warning("Parameter:DataHasInvalidElement", "Parameter '" + interdep_parameter_name + ...
                            "' data contains invalid (NaN or Inf) elements! Something to keep in mind.");
                    end
                    parameter_length(j) = length(value);
                    parameter_values{j} = value; % each value could be different length
                end
                if ~isscalar(unique(parameter_length)) % lengths all the same?
                    error("Interdependent parameters data do not have the same length!");
                end

                parameters_info{i, 2} = cell2mat(parameter_values); % values same length, combine
                parameters_info{i, 3} = parameter_length(1);
            end
        end

        function attributes_names = parseAttributesName(lum_dataset)
            % Parse attribute field names of the dataset

            % Load attribute names
            attributes = lum_dataset.Lumerical_dataset.attributes;
            if ~(isstruct(attributes) && isfield(attributes, 'variable') && isfield(attributes, 'name'))
                error("The attributes are not properly defined!");
            end

            attributes_names = cell(length(attributes), 1);
            for i = 1:length(attributes)
                attribute = attributes(i);
                % Same here, illegal characters in the names are converted to '_'
                % Check names are text and legal variable name
                LumericalDataset.validateTextScalar(attribute.variable, "One or more attribute names cannot be resolved!");
                LumericalDataset.validateTextScalar(attribute.name, "One or more attribute names cannot be resolved!");
                if ~isvarname(attribute.variable)
                    error("One or more attribute names are not valid variable names!");
                end
                if ~isequal(attribute.variable, attribute.name)
                    warning("Attribute:NameHasIllegalCharacters", "Attribute name '" + attribute.name + ...
                        "' contains illegal characters. Converted to '" + attribute.variable + "'.");
                end
                % Check attribute data exists
                LumericalDataset.validateFieldInStruct(lum_dataset, attribute.variable, "Attribute field '" + attribute.variable + "' data not found!");
                attributes_names{i} = attribute.variable;
            end
        end


        function [attributes_info, attributes_component] = parseAttributesData(lum_dataset, attributes_names, parameters_info, dataset_type, total_xyz_size)
            % Parse attribute field data of the dataset

            for i = 1:length(attributes_names)
                attribute_name = attributes_names{i};
                attribute_value = lum_dataset.(attribute_name);
                % Verify attribute value non-empty numeric
                if ~isnumeric(attribute_value) || isempty(attribute_value)
                    error("Attribute field '" + attribute_name + "' data must be numeric!");
                end
                % Give warning if attribute data contains NaN or Inf
                if any(isnan(attribute_value), 'all') || any(isinf(attribute_value), 'all')
                    warning("Attribute:DataHasInvalidElement", "Attribute field '" + attribute_name + ...
                        "' data contains invalid (NaN or Inf) elements! Something to keep in mind.");
                end
                % Check first dimension: should equal to multiplied x,y,z lengths
                if isequal(dataset_type, 'rectilinear')
                    if size(attribute_value, 1) ~= total_xyz_size
                        error("Unexpected size for attribute field '" + attribute_name + "' data at 1st dimension!");
                    end
                else % should assume dataset_type == 'matrix'
                    if size(attribute_value, 1) ~= 1
                        error("Unexpected size for attribute field '" + attribute_name + "' data at 1st dimension!");
                    end
                end
                % Check second dimension: scalar or vector
                if size(attribute_value, 2) == 1
                    attributes_component.(attribute_name) = NaN;
                elseif size(attribute_value, 2) == 3
                    attributes_component.(attribute_name) = 0; % default-magnitude
                else
                    error("Unexpected size for attribute field '" + attribute_name + "' data at 2nd dimension!");
                end
                % Check remaining dimensions, should agree with each parameter length
                for k = 1:size(parameters_info, 3)
                    if size(attribute_value, k + 2) ~= parameters_info{k, 3}
                        error("Unexpected size for attribute field '" + attribute_name + "' data at dimension " + (k+2) + " !");
                    end
                end
                % If ndims-2 < number of parameters, that means there is no extra
                % dimension(s) in the attribute data
                if ndims(attribute_value) > size(parameters_info, 1) + 2
                    error("Too many dimensions for attribute field '" + attribute_name + "' data!");
                end

                attributes_info.(attribute_name) = attribute_value;
            end
        end
    end
end