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
            %
            % Input arguments are the same as 'load' function. First
            % argument specifies the MAT file name, and the remaining
            % arguments optionally specify the variables to load. If
            % variables not specified, all variables are loaded.
            %
            % For variables that are not convertible to dataset objects,
            % leave them as it is.

            % Legal arguments supplied to load the MAT file (may throw
            % exception from the 'load' function)
            data = load(mat_name, "-mat", varargin{:});
            % Specified file must be non-empty
            if isempty(fieldnames(data))
                error("Cannot process an empty MAT file!");
            end
            variable_names = fieldnames(data);
            % Attempt to convert each variable
            for i = 1:length(variable_names)
                variable_name = variable_names{i};
                try % if not convertable, leave it as it is and issue warning
                    converted_obj.(variable_name) = ...
                        LumericalDataset.createObject(data.(variable_name));
                catch ME
                    converted_obj.(variable_name) = data.(variable_name);
                    warning_msg = ['Variable ''', variable_name, ...
                        ''' in the .mat file was not a dataset and was left unconverted! Reason: ', ME.message];
                    warning(warning_msg);
                end
            end
        end

        function new_obj = convertDatasetArray(obj_array, new_parameter_name, new_parameter_values)
            % Convert an dataset array (from a manual parameter sweep) to
            % one single dataset array
            % Input dataset array can be 1D cell array of converted or
            % unconverted datasets, or a regular 1D homogeneous array of
            % converted datasets.

            % Check obj_array, size, content
            % Non-empty vector
            if isempty(obj_array) || ~isvector(obj_array)
                error("Must be a non-empty dataset vector!");
            end
            array_size = length(obj_array);
            % Cell array or 1D array
            if iscell(obj_array)
                % Check every element is the same dataset type
                type_name = unique(cellfun(@class, obj_array, 'UniformOutput', false));
                if ~isscalar(type_name) % all elements should be one type
                    error("Inhomogeneous type!");
                end
                if type_name == "RectilinearDataset"
                    isrecti = true;
                elseif type_name == "MatrixDataset"
                    isrecti = false;
                elseif type_name == "struct"
                    try
                        for ii = 1:array_size
                            obj_array{ii} = LumericalDataset.createObject(obj_array{ii});
                        end
                    catch ME
                        fprintf("Cannot process the datasets. At least one of the datasets cannot be interpreted as either the unconverted or converted dataset. See the error message for the reason.");
                        rethrow(ME);
                    end
                else
                    error("Wrong type!");
                end
            else % a homogeneous array **probably**
                type_name = class(obj_array);
                if type_name == "RectilinearDataset"
                    isrecti = true;
                elseif type_name == "MatrixDataset"
                    isrecti = false;
                else
                    error("Wrong type!");
                end
                obj_array = num2cell(obj_array);
            end
            %%%%%% Check everything the same
            % Check all existing parameter names the same
            % Check all parameters the same value within an error limit
            for ii = 1:array_size
                if ~isequal(obj_array{ii}.parameters(:, 1), obj_array{1}.parameters(:, 1))
                    error("Parameter names not the same!");
                end
                if isrecti
                    if ~LumericalDataset.isequalWithinTol(obj_array{ii}.x, obj_array{1}.x)
                        error("x not equal!");
                    end
                    if ~LumericalDataset.isequalWithinTol(obj_array{ii}.y, obj_array{1}.y)
                        error("y not equal!");
                    end
                    if ~LumericalDataset.isequalWithinTol(obj_array{ii}.z, obj_array{1}.z)
                        error("z not equal!");
                    end
                end
                for iP = 1:obj_array{1}.num_parameters
                    if ~LumericalDataset.isequalWithinTol(obj_array{ii}.parameters{iP, 2}, obj_array{1}.parameters{iP, 2})
                        error("parameter not equal!");
                    end
                end
            end
            % Check parameter slice indexes, only give off warning if not
            % matching
            for ii = 1:array_size
                if ~isequal(obj_array{ii}.parameters_indexes, obj_array{1}.parameters_indexes)
                    warning("Parameter slice index different!");
                end
                if isrecti
                    if ~isequal(obj_array{ii}.axes_indexes, obj_array{1}.axes_indexes)
                        warning("Axes slice index different!");
                    end
                end
            end
            % Check attribute fields (names) the same
            for ii = 1:array_size
                if ~isequal(fieldnames(obj_array{ii}.attributes), fieldnames(obj_array{1}.attributes))
                    error("Attribute sets not the same!");
                end
            end
            % Check attribute plot components the same
            for ii = 1:array_size
                % Attributes_component could be NaN
                if ~isequaln(obj_array{ii}.attributes_component, obj_array{1}.attributes_component)
                    error("Attribute component not the same!");
                end
            end
            %%%%%%%%%%%
            % Check parameter name validity (must not collide with existing
            % ones)
            if ~isvarname(new_parameter_name) % also checked text scalar
                error("New parameter name is not valid!");
            end
            % Name cannot collide with any existing parameter or attribute
            % name
            pnames = obj_array{1}.parameters(:, 1);
            if ismember(new_parameter_name, [pnames{:}])
                error("Name clash with existing parameter names!");
            end
            if ismember(new_parameter_name, fieldnames(obj_array{1}.attributes))
                error("Name clash with existing attribute names!");
            end
            if isrecti && ismember(new_parameter_name, ["x", "y", "z"])
                error("Name cannot be 'x', 'y' or 'z' for rectilinear datasets!");
            end

            % Check new parameter values (real, ...), size match
            LumericalDataset.validateNonEmptyNumericVector(new_parameter_values, ...
                "New parameter values is not a numeric vector!");
            new_parameter_values = new_parameter_values(:); % convert to column vector, if applicable
            if ~length(new_parameter_values) == array_size
                error("New parameter values length does not match the dataset count!");
            end
            % Remove complex portion
            if any(imag(new_parameter_values)) % has imaginary part?
                warning("Parameter:DataIsComplex", ...
                    "New parameter values data is complex! Takes the real part and proceed.");
                new_parameter_values = real(new_parameter_values);
            end
            if any(isnan(new_parameter_values)) || any(isinf(new_parameter_values)) % real part has NaN or Inf?
                warning("Parameter:DataHasInvalidElement", ...
                    "New parameter value data contains invalid (NaN or Inf) elements! Something to keep in mind.");
            end

            % Start conversion
            new_obj = obj_array{1}.copy();

            % Can be sure that same type (scalar, vector) for each
            % attribute because attributes_component has been checked to be
            % equal for all datasets
            % Concatenate attribute data
            attributes_names = fieldnames(new_obj.attributes);
            for iA = 1:length(attributes_names)
                name = attributes_names{iA};
                new_attr = cell(array_size, 1);
                for ii = 1:array_size
                    new_attr{ii} = obj_array{ii}.attributes.(name);
                end
                new_obj.attributes.(name) = cat(new_obj.num_parameters + 3, new_attr{:});
            end
            % Update parameters, num_parameters, parameters_indexes
            new_obj.num_parameters = new_obj.num_parameters + 1;
            new_obj.parameters{end+1, 1} = string(new_parameter_name);
            new_obj.parameters{end, 2} = new_parameter_values(:);
            new_obj.parameters{end, 3} = array_size;
            new_obj.parameters_indexes(end+1, 1) = 1;
        end
    end

    methods
        function obj = LumericalDataset()
            % Load Lumerical exported MATLAB dataset into custom class object
            % Empty. Definitions in derived classes.
        end

        function showInformation(obj)
            % Display the information for an overview of this dataset
            % Print attributes information (common to both dataset types)
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
        end

        function result = getAttributeData(obj, attribute_name) % non-virtual
            % Get the data for an attribute
            obj.iCheckAttributeExist(attribute_name);
            result = obj.attributes.(attribute_name);
        end

        function setParameterSliceIndex(obj, varargin) % non-virtual
            % Set parameters (including x, y, z for rectilinear datasets)
            % slice position based on input indexes. See function
            % 'setParameterSlice'.
            try
                obj.setParameterSlice("index", varargin{:});
            catch ME
                ME.throw();
            end
        end

        function setParameterSliceValue(obj, varargin) % non-virtual
            % Set parameters (including x, y, z for rectilinear dataset)
            % slice position based on input values. See function
            % 'setParameterSlice'.
            try
                obj.setParameterSlice("value", varargin{:});
            catch ME
                ME.throw();
            end
        end

        function setAttributeComponent(obj, attribute_name, component) % non-virtual
            % set the component (x,y,z,magnitude) for an attribute

            arguments
                obj
                attribute_name
                component {mustBeMember(component, ["x", "y", "z", "magnitude"])}
            end

            obj.iCheckAttributeExist(attribute_name);
            if isnan(obj.attributes_component.(attribute_name))
                error("Scalar attribute cannot be modified!");
            end

            switch component
                case "x"
                    obj.attributes_component.(attribute_name) = 1;
                case "y"
                    obj.attributes_component.(attribute_name) = 2;
                case "z"
                    obj.attributes_component.(attribute_name) = 3;
                case "magnitude"
                    obj.attributes_component.(attribute_name) = 0;
            end
        end

        function hPlot = plotData1D(obj, parameter_name, attribute_name, optargs, plotprops)  % non-virtual
            % This function makes a 1D plot based on the parameter name and
            % the attribute name.
            arguments
                obj
                parameter_name
                attribute_name
                optargs.ScalarOperation {mustBeMember(optargs.ScalarOperation, ["real", "imag", "abs", "angle"])} = "real"
                optargs.Ax = gca()
                optargs.XFactor (1, 1) {mustBeReal, mustBeNonNan} = 1
                optargs.YFactor (1, 1) {mustBeReal, mustBeNonNan} = 1
                optargs.LineSpec (1, 1) string = ""
                plotprops.?matlab.graphics.chart.primitive.Line % plot name-value pairs
            end

            % Convert received plot name-value pairs to a cell array
            plotpropscell = namedargs2cell(plotprops);

            % If 'ax' is a valid axes handle, plot on that axes object
            ax = optargs.Ax;
            if  ~(isscalar(ax) && isgraphics(ax, 'axes')) % valid axes handle
                error("'Ax' optional argument must be a valid axes handle!");
            end

            % Calls the respective polymorphic function for each object
            [xdata, ydata] = obj.getPlot1DData(parameter_name, attribute_name);

            ydata = applyScalarOperation(ydata, optargs.ScalarOperation);

            if any(isnan(xdata)) || any(isinf(xdata)) % xdata has NaN or Inf?
                error("Parameter:DataHasInvalidElement", ...
                    "Parameter data contains invalid (NaN or Inf) elements! Unable to make the plot.");
            end

            if any(isnan(ydata)) || any(isinf(ydata)) % ydata has NaN or Inf?
                error("Attribute:DataHasInvalidElement", ...
                    "Attribute data contains invalid (NaN or Inf) elements! Unable to make the plot.");
            end

            hPlot = plot(ax, xdata * optargs.XFactor, ydata * optargs.YFactor, ...
                optargs.LineSpec, plotpropscell{:});
            xlabel(ax, parameter_name, 'Interpreter', 'none');
            ylabel(ax, attribute_name, 'Interpreter', 'none');
        end

        function [hSurf, hClb] = plotData2D(obj, parameter1_name, parameter2_name, attribute_name, ...
                optargs, surfprops)  % non-virtual
            % This function makes a 2D plot based on the parameter names
            % and the attribute name.

            arguments
                obj
                parameter1_name
                parameter2_name
                attribute_name
                optargs.ScalarOperation {mustBeMember(optargs.ScalarOperation, ["real", "imag", "abs", "angle"])} = "real"
                optargs.Ax = gca()
                optargs.XFactor (1, 1) {mustBeReal, mustBeNonNan} = 1
                optargs.YFactor (1, 1) {mustBeReal, mustBeNonNan} = 1
                optargs.Mode {mustBeMember(optargs.Mode, ["none", "spatial"])} = "none"
                surfprops.?matlab.graphics.primitive.Surface % surface name-value pairs
            end

            % Convert received surface name-value pairs to a cell array
            surfpropscell = namedargs2cell(surfprops);

            % If 'ax' is a valid axes handle, plot on that axes object
            ax = optargs.Ax;
            if ~(isscalar(ax) && isgraphics(ax, 'axes')) % valid axes handle
                error("'Ax' optional argument must be a valid axes handle!");
            end

            % Calls the respective polymorphic function for each object
            [xdata, ydata, zdata] = obj.getPlot2DData(parameter1_name, parameter2_name, attribute_name);

            zdata = applyScalarOperation(zdata, optargs.ScalarOperation);

            % Throw an error is xdata or ydata (parameters) is singleton
            % dimension
            if isscalar(xdata) || isscalar(ydata)
                error("Cannot make 2D plot with singleton dimension in xdata or ydata!");
            end

            % Throw an error is xdata or ydata (parameters) is not monotonic
            if ~LumericalDataset.isRealVectorMonotonic(xdata)
                error("x data is not monotonic! Cannot make 2D plot.");
            end
            if ~LumericalDataset.isRealVectorMonotonic(ydata)
                error("y data is not monotonic! Cannot make 2D plot.");
            end

            if any(isnan(xdata)) || any(isinf(xdata)) || ...
                    any(isnan(ydata)) || any(isinf(ydata)) % xydata has NaN or Inf?
                error("Parameter:DataHasInvalidElement", ...
                    "Parameter data contains invalid (NaN or Inf) elements! Unable to make the plot.");
            end

            if any(isnan(zdata), 'all') || any(isinf(zdata), 'all') % zdata has NaN or Inf?
                error("Attribute:DataHasInvalidElement", ...
                    "Attribute data contains invalid (NaN or Inf) elements! Unable to make the plot.");
            end

            % Scale data
            if strcmp(optargs.Mode, "spatial")
                xdata = xdata * 1e9;
                ydata = ydata * 1e9;
            end
            xdata = xdata * optargs.XFactor;
            ydata = ydata * optargs.YFactor;

            % Adjust data to make true "2D plot"
            xdata = ([xdata(1); xdata(:)] + [xdata(:); xdata(end)])/2;
            ydata = ([ydata(1); ydata(:)] + [ydata(:); ydata(end)])/2;
            zdata_new = zeros(length(ydata), length(xdata));
            zdata_new(1:end-1, 1:end-1) = zdata;
            zdata = zdata_new;

            hSurf = surface(ax, xdata, ydata, zdata, 'EdgeColor', 'none', surfpropscell{:});
            xlim(ax, [min(xdata), max(xdata)]);
            ylim(ax, [min(ydata), max(ydata)]);
            if strcmp(optargs.Mode, "spatial")
                xlabel(ax, parameter1_name + " (nm)", 'Interpreter', 'none');
                ylabel(ax, parameter2_name + " (nm)", 'Interpreter', 'none');
            elseif strcmp(optargs.Mode, "none")
                xlabel(ax, parameter1_name + " (nm)", 'Interpreter', 'none');
                ylabel(ax, parameter2_name + " (nm)", 'Interpreter', 'none');
            end
            colormap(ax, 'jet');
            hClb = colorbar;
            set(ax, 'Layer', 'top');
            box(ax, 'on');
        end

        function interp_data = getInterpolatedPlot2DData(obj, parameter1_name, parameter2_name, attribute_name, ...
                Xq, Yq, method, extrapval, options)  % non-virtual
            % This function acquires 2D plot data and then interpolate
            % Z(X,Y) based on the query Xq and Yq points.
            % Can specify 'method' and 'extrapval' as optional positional
            % arguments as defined in 'interp2' function. Additionally, can
            % specify 'ScalarOperation' as an optional name-value pair
            % argument.
            % The requirements for Xq and Yq is the same as 'interp2'
            % function. They can take on different forms and do not need to
            % be monotonic. See documentation.

            arguments
                obj
                parameter1_name
                parameter2_name
                attribute_name
                Xq
                Yq
                method = 'linear' % make sure this is the default for 'interp2'
                extrapval = 'not specified'
                options.ScalarOperation (1, 1) string {mustBeMember(options.ScalarOperation, {'real', 'imag', 'abs', 'angle'})} = 'real'
            end

            [xdata, ydata, zdata] = obj.getPlot2DData(parameter1_name, parameter2_name, attribute_name);
            zdata = applyScalarOperation(zdata, options.ScalarOperation);

            % Use EAFP to let the 'interp2' function check these arguments
            % Do not check NaN, Inf or monotonicity of xdata or ydata.
            % zdata can have NaN or Inf, interp2 still works
            % Calls the respective polymorphic function for each object
            if strcmp(extrapval, 'not specified') % user does not give value
                interp_data = interp2(xdata, ydata, zdata, Xq, Yq, method);
            else % user specified extrapval
                interp_data = interp2(xdata, ydata, zdata, Xq, Yq, method, extrapval);
            end
        end
    end

    methods (Abstract)
        result = getParameterData(obj, parameter_name);
        [xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name);
        [xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name);
        [x, y, z, data] = getPlot3DData(obj, parameter1_name, parameter2_name, parameter3_name, attribute_name);
        new_obj = removeDimensions(obj, varargin);
        new_obj = mergeDataset(obj, other_obj, varargin);
    end

    methods (Abstract, Access = protected)
        setParameterSlice(obj, mode_flag, varargin);
    end

    methods (Access = protected)
        function p = iAddAllParametersToParser(obj, p, check_mode)
            % Add all parameters and their range validations to the parser
            for i = 1:size(obj.parameters, 1)
                for para = obj.parameters{i, 1}
                    if strcmp(check_mode, "index") % validate index
                        p.addParameter(para, NaN, ...
                            @(x) LumericalDataset.validateIndex(x, obj.parameters{i, 3}));
                    elseif strcmp(check_mode, "value") % validate value
                        p.addParameter(para, NaN, @LumericalDataset.mustBeRealNumericScalar);
                    end
                end
            end
        end

        function parsed_index_list = iAnalyzeParsedParameters(obj, p, parse_mode)
            % Parse mode can be either "index" or "value"
            parsed_index_list = NaN(obj.num_parameters, 1);
            for i = 1:obj.num_parameters % each interdep params set
                interdep_set_data = NaN(length(obj.parameters{i, 1}), 1); % either index or value
                interdep_indexes = NaN(length(obj.parameters{i, 1}), 1);
                for j = 1:length(obj.parameters{i, 1}) % loop through each interdep param
                    interdep_set_data(j) = p.Results.(obj.parameters{i, 1}(j));
                    if ~isnan(interdep_set_data(j)) && strcmp(parse_mode, "value")
                        try % find the index corresponding to the value
                            interdep_indexes(j) = LumericalDataset.findIndexFromValueWithinTol( ...
                                interdep_set_data(j), obj.parameters{i, 2}(:, j), ...
                                "Cannot find the value specified for '" + obj.parameters{i, 1}(j) + "'!");
                        catch ME
                            ME.throwAsCaller();
                        end
                    end
                end
                if strcmp(parse_mode, "index")
                    interdep_indexes = interdep_set_data;
                end
                % If multiple interdependent parameters are present, they
                % must resolve to the same index
                unique_index = unique(interdep_indexes(~isnan(interdep_indexes)));
                if numel(unique_index) > 1 % different indexes?
                    ME = MException('', "Multiple interdependent parameters in the same set were selected, " + ...
                        "but they do not resolve to the same index!");
                    ME.throwAsCaller();
                elseif numel(unique_index) == 1
                    parsed_index_list(i) = unique_index;
                end % no parameter might be specified. Simply skip
            end
        end

        function iUpdateParametersSliceIndex(obj, parsed_index_list)
            obj.parameters_indexes(~isnan(parsed_index_list)) = ...
                parsed_index_list(~isnan(parsed_index_list));
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
                ME.throwAsCaller();
            end
        end

        function location = iCheckAndFindParameter(obj, parameter_name)
            % Check if a parameter name exists and return location
            % Location: [row, column (within interdep)]
            try
                LumericalDataset.validateTextScalar(parameter_name, "Parameter name must be text scalar!");
            catch ME
                ME.throwAsCaller();
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
                ME = MException('', "Parameter name '" + parameter_name + "' not found!");
                ME.throwAsCaller();
            end
        end

        function iCheckParameterExist(obj, parameter_name)
            % Check if a parameter name exists
            try
                LumericalDataset.validateTextScalar(parameter_name, "Parameter name must be text scalar!");
            catch ME
                ME.throwAsCaller();
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
                ME = MException('', "Parameter name '" + parameter_name + "' not found!");
                ME.throwAsCaller();
            end
        end

        function iCheckAttributeExist(obj, attribute_name)
            % Check if an attribute name exists
            try
                LumericalDataset.validateTextScalar(attribute_name, "Attribute name must be text scalar!");
                LumericalDataset.validateFieldInStruct(obj.attributes, attribute_name, "This attribute '" + attribute_name + "' is not found!");
            catch ME
                ME.throwAsCaller();
            end
        end
    end

    methods (Static, Access = protected)
        % Helper functions to be shared with base and derived classes
        function mustBeRealNumericScalar(input)
            if ~(isnumeric(input) && isscalar(input) && isreal(input))
                error("Value must be real numeric scalar!");
            end
        end

        function validateTextScalar(input, errmsg)
            % Validate input as text scalar and throw if not
            if ~(ischar(input) && isrow(input)) && ~isStringScalar(input)
                ME = MException('MATLAB:NotTextScalar', errmsg);
                ME.throwAsCaller();
            end
        end

        function validateStructScalar(input, errmsg)
            % Validate input as struct scalar and throw if not
            if ~isstruct(input) || ~isscalar(input)
                ME = MException('MATLAB:NotStructScalar', errmsg);
                ME.throwAsCaller();
            end
        end

        function validateFieldInStruct(struct_in, field_in, errmsg)
            % Validate field exists in a struct
            if ~isfield(struct_in, field_in)
                ME = MException('MATLAB:FieldNotInStruct', errmsg);
                ME.throwAsCaller();
            end
        end

        function validateNonEmptyNumericVector(input, errmsg)
            % Validate input as non-empty vector (N-by-1 or 1-by-N)
            if ~(isnumeric(input) && isvector(input) && ~isempty(input))
                ME = MException('MATLAB:NotNonEmptyNumericVector', errmsg);
                ME.throwAsCaller();
            end
        end

        function validateIndex(idx, max_val)
            % Validate input index makes sense in the context
            if ~isnumeric(idx) || floor(idx) ~= idx % test integer
                ME = MException('MATLAB:InvalidIndex', "Must be an integer!");
                ME.throwAsCaller();
            end
            if idx < 1 || idx > max_val
                ME = MException('MATLAB:InvalidIndex', "Must be between 1 and " + max_val + "!");
                ME.throwAsCaller();
            end
        end

        function result = sliceThroughAttributeVectorDim(attribute_data, component_index)
            % Slice off attribute data matrix based on the selected vector
            % component (NaN-scalar, 0-magnitude, 1-x, 2-y, 3-z)
            % on the second dimension
            %
            % NOT checked
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
            %
            % NOT checked
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

        function index = findIndexFromValueWithinTol(value, array, errmsg)
            [~, index] = min(abs(array - value));
            % Is it within tol?
            if ~LumericalDataset.isequalWithinTol(value, array(index))
                if nargin < 3
                    ME = MException('', "Cannot find the value in the array within the tolerance level!");
                else
                    ME = MException('', errmsg);
                end
                ME.throwAsCaller();
            end
        end

        function tf = isRealVectorMonotonic(vec)
            % Returns true if a real-valued non-empty vector is strictly
            % monotonic (increasing or decreasing)
            %
            % NOT checked
            tf =  all(diff(vec) > 0) || all(diff(vec) < 0);
        end

        function printParametersInfo(param_names, lengths, indexes)
            % Print title
            num_parameters = length(param_names);
            fprintf("%d parameter set(s):\n", num_parameters);
            % Max parameter name length (multiple interdependent parameters names
            % separated by '|')
            max_strlength_names = max(cellfun(@(x) strlength(join(x, " | ")), param_names));
            max_strlength_names = max(max_strlength_names, strlength('Name(s)'));
            % Max number of digits for parameter data length
            max_digits_length = max(arrayfun(@(x) length(num2str(x)), lengths));
            max_digits_length = max(max_digits_length, strlength('Length'));
            % Max number of digits for parameter slice index
            max_digits_sliceindex = max(arrayfun(@(x) length(num2str(x)), indexes));
            max_digits_sliceindex = max(max_digits_sliceindex, strlength('Slice index'));
            % Print header
            fprintf("| %-*s | %-*s | %-*s |\n", ...
                max_strlength_names, 'Name(s)', ...
                max_digits_length, 'Length', ...
                max_digits_sliceindex, 'Slice index');
            % Print separator
            fprintf("+%s+%s+%s+\n", ...
                repmat('-', 1, max_strlength_names + 2), ...
                repmat('-', 1, max_digits_length + 2), ...
                repmat('-', 1, max_digits_sliceindex + 2));
            % Print contents
            for i = 1:num_parameters
                fprintf("| %-*s | %*d | %*d |\n", ...
                    max_strlength_names, join(param_names{i}, " | "), ...
                    max_digits_length, lengths(i), ...
                    max_digits_sliceindex, indexes(i));
            end
        end
    end

    methods (Static, Access = protected)
        % Helper functions for dataset loading
        function dataset_type = parseDatasetStructure(lum_dataset)
            % Initially pase the dataset structure and determine dataset
            % type (matrix or rectilinear)

            % Input dataset has to be a struct scalar
            LumericalDataset.validateStructScalar(lum_dataset, "Input dataset must be a struct scalar!");

            % Check field: Lumerical_dataset
            LumericalDataset.validateFieldInStruct(lum_dataset, 'Lumerical_dataset', "Input dataset does not have the field 'Lumerical_dataset'!");
            LumericalDataset.validateStructScalar(lum_dataset.Lumerical_dataset, "Field 'Lumerical_dataset' is not a struct scalar!");
            % Check 'attribute' field (check the contents later)
            LumericalDataset.validateFieldInStruct(lum_dataset.Lumerical_dataset, 'attributes', "Field 'Lumerical_dataset' is missing the 'attributes' subfield! " + ...
                "Dataset without any attribute is not supported.");
            % Check 'parameters' field (check the contents later)
            LumericalDataset.validateFieldInStruct(lum_dataset.Lumerical_dataset, 'parameters', "Field 'Lumerical_dataset' is missing the 'parameters' subfield!")

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
                    % Parameter data will never be multi-dimensional. Even
                    % if it is originally, it will be converted to 1D
                    % vector during MATLAB export.
                    % Check data to be non-empty vector (should always be
                    % column vector, actually)
                    LumericalDataset.validateNonEmptyNumericVector(value, ...
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
                if isequal(dataset_type, 'matrix')
                    if size(attribute_value, 1) ~= 1
                        error("In matrix dataset, the size for attribute field '" ...
                            + attribute_name + "' data at 1st dimension should be 1!");
                    end
                else % should assume dataset_type == 'rectilinear'
                    if size(attribute_value, 1) ~= total_xyz_size
                        error("In rectilinear dataset, the size for attribute field '" ...
                            + attribute_name + "' data at 1st dimension should match the product of x, y and z dimensions!");
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
                % Check remaining dimensions, should agree with each parameter data length
                for k = 1:size(parameters_info, 3)
                    if size(attribute_value, k + 2) ~= parameters_info{k, 3}
                        error("Unexpected size for attribute field '" + attribute_name + "' data at dimension " + (k+2) + " !");
                    end
                end
                % If ndims-2 <= number of parameters, that means there is
                % no extra dimension(s) in the attribute data
                if ndims(attribute_value) > size(parameters_info, 1) + 2
                    error("Too many dimensions for attribute field '" + attribute_name + "' data!");
                end
                % Finally write data to result
                attributes_info.(attribute_name) = attribute_value;
            end
        end
    end
end

%% Helper functions
function ydata = applyScalarOperation(ydata, scalar_operation)
% input already checked
switch scalar_operation
    case 'real'
        ydata = real(ydata);
    case 'imag'
        ydata = imag(ydata);
    case 'abs'
        ydata = abs(ydata);
    case 'angle'
        ydata = angle(ydata);
end
end