classdef MatrixDataset < LumericalDataset
    % Subclass for matrix dataset

    properties
    end

    methods
        function obj = MatrixDataset(lum_dataset)
            % Call superclass constructor
            obj = obj@LumericalDataset(lum_dataset);
        end

        function obj = setParameterSliceIndex(obj, varargin)
            % Set slice index of one or more parameters (specified as
            % name-value pair)

            % Initialize inputParser, add regular parameters
            p = inputParser();
            p.PartialMatching = false;
            p = obj.iAddParametersToParser(p);

            % Parse input argument
            try
                p.parse(varargin{:});
            catch ME
                throwAsCaller(ME);
            end

            % Analyze and set values for regular parameters
            % If multiple interdependent parameters are declared, their
            % indexes should be the same. Otherwise, report an error.
            obj.iAnalyzeAndSetParsedParameter(p);
        end

        function [xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name)
            % Get x and y data for 1D plot
            para_value_list = cell(1, 2); % 1D, xdata
            [para_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, parameter_name);
            xdata = para_value_list{1, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            ydata = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            ydata = squeeze(ydata(:, :, para_indexes{:}));
        end

        function [xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name)
            % Get x,y and z data for 2D plot
            para_value_list = cell(2, 2); % 2D, xdata & ydata
            [para_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, parameter1_name, parameter2_name);
            xdata = para_value_list{1, 1};
            ydata = para_value_list{2, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            zdata = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            zdata = squeeze(zdata(:, :, para_indexes{:}));

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
            [para_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, parameter1_name, parameter2_name, parameter3_name);
            x = para_value_list{1, 1};
            y = para_value_list{2, 1};
            z = para_value_list{3, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            data = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            data = squeeze(data(:, :, para_indexes{:}));

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
    end
end