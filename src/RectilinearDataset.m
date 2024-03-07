classdef RectilinearDataset < LumericalDataset
    % Subclass for rectilinear dataset

    properties (SetAccess = protected)
        x
        y
        z
        axes_indexes % axes parameter indexes (x, y, z)
    end

    methods
        function obj = RectilinearDataset(lum_dataset)
            % Parse dataset structure, determine dataset type
            dataset_type = LumericalDataset.parseDatasetStructure(lum_dataset);
            if dataset_type ~= "rectilinear"
                error("This is not a rectilinear dataset!");
            end

            % Parse parameter and attribute names
            parameters_info = LumericalDataset.parseParametersName(lum_dataset);
            attributes_names = LumericalDataset.parseAttributesName(lum_dataset);
            % No duplicate names allowed
            allnames = [parameters_info{:, 1}, string(attributes_names).', "x", "y", "z"];
            if length(unique(allnames)) ~= length(allnames)
                error("Parameters, positional vectors (x,y,z) and attributes have duplicate names! Rejected!");
            end

            % Parse parameter data
            obj.parameters = LumericalDataset.parseParametersData(lum_dataset, parameters_info);
            obj.num_parameters = size(obj.parameters, 1);
            obj.parameters_indexes = ones(obj.num_parameters, 1); % initialize to all 1

            % Parse positional vectors
            [obj.x, obj.y, obj.z, xyz_prod_size] = LumericalDataset.parseXYZ(lum_dataset);
            obj.axes_indexes = ones(3, 1);

            % Parse attribute data
            [obj.attributes, obj.attributes_component] = LumericalDataset.parseAttributesData(lum_dataset, attributes_names, obj.parameters, dataset_type, xyz_prod_size);
            obj.num_attributes = length(fieldnames(obj.attributes));
        end

        function showInformation(obj)
            % Print rectilinear dataset header
            fprintf("This dataset is a rectilinear dataset.\n");
            fprintf('\n');
            showInformation@LumericalDataset(obj);
            fprintf('\n');
            % Print parameters information (including x, y, z)
            LumericalDataset.printParametersInfo(vertcat({"x"; "y"; "z"}, obj.parameters(:, 1)), ...
                [length(obj.x); length(obj.y); length(obj.z); cell2mat(obj.parameters(:, 3))], ...
                [obj.axes_indexes; obj.parameters_indexes]);
        end

        function [xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name)
            % Get x and y data for 1D plot
            para_value_list = cell(1, 2); % 1D, xdata
            % First check if they are x, y or z
            [arglist, axes_keep_indexes, para_value_list] = ...
                obj.iGenerateAxesSliceIndexAndData(para_value_list, parameter_name);
            [para_slice_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, arglist{:});
            xdata = para_value_list{1, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            attribute_data = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            % Expand first index into x,y,z
            sz = size(attribute_data);
            attribute_data = reshape(attribute_data, [length(obj.x), length(obj.y), length(obj.z), sz(2:end)]);

            ydata = squeeze(attribute_data(axes_keep_indexes{:}, :, para_slice_indexes{:}));
        end

        function [xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name)
            % Get x,y and z data for 2D plot
            para_value_list = cell(2, 2); % 2D, xdata & ydata
            % First check if they are x, y or z
            [arglist, axes_keep_indexes, para_value_list] = ...
                obj.iGenerateAxesSliceIndexAndData(para_value_list, parameter1_name, parameter2_name);
            [para_slice_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, arglist{:});
            xdata = para_value_list{1, 1};
            ydata = para_value_list{2, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            attribute_data = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            % Expand first index into x,y,z
            sz = size(attribute_data);
            attribute_data = reshape(attribute_data, [length(obj.x), length(obj.y), length(obj.z), sz(2:end)]);

            zdata = squeeze(attribute_data(axes_keep_indexes{:}, :, para_slice_indexes{:}));

            % Check for interdependent parameters and rearrange zdata based
            % on the order of parameters
            if (para_value_list{1, 2} == para_value_list{2, 2})
                error("Can't plot against two interdependent parameters!");
            elseif (para_value_list{1, 2} > para_value_list{2 ,2}) % flipped
                zdata = zdata.'; % flip
            end

            zdata = zdata.'; % swap x,y dimensions
        end

        function [x, y, z, data] = getPlot3DData(obj, parameter1_name, parameter2_name, parameter3_name, attribute_name)
            % Get x,y,z and data for 3D plot
            para_value_list = cell(3, 2); % 3D, x,y,z
            % First check if they are x, y or z
            [arglist, axes_keep_indexes, para_value_list] = ...
                obj.iGenerateAxesSliceIndexAndData(para_value_list, parameter1_name, parameter2_name, parameter3_name);
            [para_slice_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, arglist{:});
            x = para_value_list{1, 1};
            y = para_value_list{2, 1};
            z = para_value_list{3, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            attribute_data = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            % Expand first index into x,y,z
            sz = size(attribute_data);
            attribute_data = reshape(attribute_data, [length(obj.x), length(obj.y), length(obj.z), sz(2:end)]);

            data = squeeze(attribute_data(axes_keep_indexes{:}, :, para_slice_indexes{:}));

            % Check for interdependent parameters and rearrange data based
            % on the order of parameters
            if (para_value_list{1, 2} == para_value_list{2, 2} || ...
                    para_value_list{2, 2} == para_value_list{3, 2} || ...
                    para_value_list{1, 2} == para_value_list{3, 2})
                error("Can't plot against two interdependent parameters!");
            end

            [~, order] = sort(cell2mat(para_value_list(:, 2)));
            [~, rank] = sort(order); % sort twice to get params ranking
            data = permute(data, rank); % rearrange based on params ranking
            data = permute(data, [2, 1, 3]); % swap x,y dimensions
        end

        function new_obj = removeDimensions(obj, varargin)
            % Remove some dimensions (parameters or x,y,z) of the dataset
            % When removing x,y,z, reset them to be length of 1
            para_value_list = cell(nargin - 1, 2); % nargin includes obj
            % First check if they are x, y or z
            [arglist, ~, para_value_list, axes_remove_indexes] = ...
                obj.iGenerateAxesSliceIndexAndData(para_value_list, varargin{:});
            [~, para_value_list, para_remove_indexes] = obj.iGenerateParametersSliceIndexAndData(para_value_list, arglist{:});

            params_axes = cell2mat(para_value_list(:, 2));
            param_idx = params_axes(params_axes > 0);
            axes_idx = params_axes(params_axes < 0);

            new_obj = obj.copy();

            % Adjust attributes
            attributes_fields = fieldnames(obj.attributes);

            for i = 1:obj.num_attributes
                field = attributes_fields{i};
                attribute_data = new_obj.attributes.(field);

                % Expand first index into x,y,z
                sz = size(attribute_data); % this is a different sz only used on the next line
                attribute_data = reshape(attribute_data, [length(obj.x), length(obj.y), length(obj.z), sz(2:end)]);
                sz = size(attribute_data); % get the size after xyz reshape, before trimming
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%% This is dumb implementation of removing dimension
                attribute_data = attribute_data(axes_remove_indexes{:}, :, para_remove_indexes{:});
                % Should not squeeze (would remove other)
                sz(param_idx + 4) = [];
                sz(axes_idx + 4) = 1;
                sz2 = sz(3:end);
                sz2(1) = sz(1)*sz(2)*sz(3);
                attribute_data = reshape(attribute_data, sz2);
                new_obj.attributes.(field) = attribute_data;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Stupid implementation

            % Remove x,y or z
            for i = 1:length(axes_idx)
                switch axes_idx(i)
                    case -3
                        new_obj.x = obj.x(obj.axes_indexes(1));
                    case -2
                        new_obj.y = obj.y(obj.axes_indexes(2));
                    case -1
                        new_obj.z = obj.z(obj.axes_indexes(3));
                end
            end
            % Reset axes_indexes
            new_obj.axes_indexes(axes_idx + 4) = 1; % set sliced axis index to 1

            % Remove parameters
            new_obj.parameters(param_idx, :) = [];
            % Update parameters_indexes
            new_obj.parameters_indexes(param_idx) = [];
            % Update num_parameters
            new_obj.num_parameters = new_obj.num_parameters - length(param_idx);
        end

        function new_obj = mergeDataset(obj, other_obj, varargin)
            error("This method is not implemented.");
        end
    end

    methods (Access = private)
        function p = iAddAxesToParser(obj, p, check_mode)
            if strcmp(check_mode, "index") % validate index
                p.addParameter('x', NaN, @(x) LumericalDataset.validateIndex(x, length(obj.x)));
                p.addParameter('y', NaN, @(x) LumericalDataset.validateIndex(x, length(obj.y)));
                p.addParameter('z', NaN, @(x) LumericalDataset.validateIndex(x, length(obj.z)));
            elseif strcmp(check_mode, "value") % validate value
                p.addParameter('x', NaN, @LumericalDataset.mustBeRealNumericScalar);
                p.addParameter('y', NaN, @LumericalDataset.mustBeRealNumericScalar);
                p.addParameter('z', NaN, @LumericalDataset.mustBeRealNumericScalar);
            end
        end

        function iAnalyzeAndSetParsedAxes(obj, p, parse_mode)
            axisnames = ["x", "y", "z"];
            for iA = 1:3
                result = p.Results.(axisnames(iA));
                if ~isnan(result)
                    if strcmp(parse_mode, "value")
                        result = LumericalDataset.findIndexFromValueWithinTol( ...
                            result, obj.(axisnames(iA)), ...
                            "Cannot find the value specified for '" + axisnames(iA) + "'!");
                    end
                    obj.axes_indexes(iA) = result;
                end
            end
        end

        function [arglist, axes_keep_indexes, para_value_list, axes_remove_indexes] = ...
                iGenerateAxesSliceIndexAndData(obj, para_value_list, varargin)
            % Parse axes parameters input and generate results
            % axes_keep_indexes: mentioned axes to keep (':')
            % axes_remove_indexes: mentioned axes to be removed (slice off)
            axes_keep_indexes = num2cell(obj.axes_indexes);
            for i = 1:length(varargin)
                parameter_name = varargin{i};
                if isequal(parameter_name, "x") % parameter_name might not be text
                    axes_keep_indexes{1} = ':';
                    para_value_list{i, 1} = obj.x;
                    para_value_list{i, 2} = -3;
                    varargin{i} = [];
                elseif isequal(parameter_name, "y")
                    axes_keep_indexes{2} = ':';
                    para_value_list{i} = obj.y;
                    para_value_list{i, 2} = -2;
                    varargin{i} = [];
                elseif isequal(parameter_name, "z")
                    axes_keep_indexes{3} = ':';
                    para_value_list{i} = obj.z;
                    para_value_list{i, 2} = -1;
                    varargin{i} = [];
                end
            end

            % Get axes_remove_indexes from axes_keep_indexes
            axes_remove_indexes = num2cell(obj.axes_indexes);
            for i = 1:3
                if isnumeric(axes_keep_indexes{i})
                    axes_remove_indexes{i} = ':';
                end
            end

            arglist = varargin;
        end
    end

    methods (Access = protected)
        function setParameterSlice(obj, mode_flag, varargin)
            % Rectilinear version of set slice index that also accepts 'x',
            % 'y' and 'z' as parameters

            % Initialize inputParser, add all regular and axes parameters
            p = inputParser();
            p.PartialMatching = false;
            p = obj.iAddAllParametersToParser(p, mode_flag);
            p = obj.iAddAxesToParser(p, mode_flag); % add x, y, z

            % Parse input name-value pairs
            try
                p.parse(varargin{:});
            catch ME
                ME.throw();
            end

            % Analyze and set values for three axes input (x,y,z)
            obj.iAnalyzeAndSetParsedAxes(p, mode_flag);
            % Analyze indexes or values provided for the parameters If
            % multiple interdependent parameters are declared, they must
            % resolve to the same index. Otherwise, report an error.
            parsed_index_list = obj.iAnalyzeParsedParameters(p, mode_flag);
            obj.iUpdateParametersSliceIndex(parsed_index_list);
        end
    end
end