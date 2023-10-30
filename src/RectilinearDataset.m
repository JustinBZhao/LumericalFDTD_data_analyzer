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
            % Call superclass constructor
            obj = obj@LumericalDataset(lum_dataset);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is very stupid
            [~, ~, ~, xyz] = loadLumDataset(lum_dataset);
            obj.x = xyz.x;
            obj.y = xyz.y;
            obj.z = xyz.z;
            obj.axes_indexes = ones(3, 1);
        end

        function showInformation(obj)
            % Also print x,y,z information
            fprintf("This dataset is a rectilinear dataset.\n");
            showInformation@LumericalDataset(obj);
            fprintf('x: [%.2e, %.2e], %d points\n', min(obj.x), max(obj.x), length(obj.x));
            fprintf('y: [%.2e, %.2e], %d points\n', min(obj.y), max(obj.y), length(obj.y));
            fprintf('z: [%.2e, %.2e], %d points\n', min(obj.z), max(obj.z), length(obj.z));
        end

        function obj = setParameterSliceIndex(obj, varargin)
            % Rectilinear version of set slice index that also accepts 'x',
            % 'y' and 'z' as parameters

            % Initialize inputParser, add regular and axes parameters
            p = inputParser();
            p.PartialMatching = false;
            p = obj.iAddParametersToParser(p);
            p = obj.iAddAxesToParser(p); % add x, y, z

            % Parse input argument
            try
                p.parse(varargin{:});
            catch ME
                throwAsCaller(ME);
            end

            % Analyze and set values for three axes input (x,y,z)
            obj.iAnalyzeAndSetParsedAxes(p);
            % Analyze and set values for regular parameters
            % If multiple interdependent parameters are declared, their
            % indexes should be the same. Otherwise, report an error.
            obj.iAnalyzeAndSetParsedParameter(p);
        end

        function [xdata, ydata] = getPlot1DData(obj, parameter_name, attribute_name)
            % Get x and y data for 1D plot
            para_value_list = cell(1, 2); % 1D, xdata
            % First check if they are x, y or z
            [arglist, axes_indexes, para_value_list] = ...
                obj.iGenerateAxesSliceIndexAndData(para_value_list, parameter_name);
            [para_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, arglist{:});
            xdata = para_value_list{1, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            attribute_data = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            % Expand first index into x,y,z
            sz = size(attribute_data);
            attribute_data = reshape(attribute_data, [length(obj.x), length(obj.y), length(obj.z), sz(2:end)]);
            
            ydata = squeeze(attribute_data(axes_indexes{:}, 1, para_indexes{:}));
        end

        function [xdata, ydata, zdata] = getPlot2DData(obj, parameter1_name, parameter2_name, attribute_name)
            % Get x,y and z data for 2D plot
            para_value_list = cell(2, 2); % 2D, xdata & ydata
            % First check if they are x, y or z
            [arglist, axes_indexes, para_value_list] = ...
                obj.iGenerateAxesSliceIndexAndData(para_value_list, parameter1_name, parameter2_name);
            [para_indexes, para_value_list] = ...
                obj.iGenerateParametersSliceIndexAndData(para_value_list, arglist{:});
            xdata = para_value_list{1, 1};
            ydata = para_value_list{2, 1};

            obj.iCheckAttributeExist(attribute_name);
            attribute_data = obj.attributes.(attribute_name);
            attribute_data = LumericalDataset.sliceThroughAttributeVectorDim(attribute_data, obj.attributes_component.(attribute_name));
            % Expand first index into x,y,z
            sz = size(attribute_data);
            attribute_data = reshape(attribute_data, [length(obj.x), length(obj.y), length(obj.z), sz(2:end)]);

            zdata = squeeze(attribute_data(axes_indexes{:}, 1, para_indexes{:}));

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
            [arglist, axes_indexes, para_value_list] = ...
                obj.iGenerateAxesSliceIndexAndData(para_value_list, parameter1_name, parameter2_name, parameter3_name);
            [para_indexes, para_value_list] = ...
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

            data = squeeze(attribute_data(axes_indexes{:}, 1, para_indexes{:}));

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

    methods (Access = private)
        function p = iAddAxesToParser(obj, p)
            p.addParameter('x', NaN, @(x) LumericalDataset.validateIndex(x, length(obj.x)));
            p.addParameter('y', NaN, @(x) LumericalDataset.validateIndex(x, length(obj.y)));
            p.addParameter('z', NaN, @(x) LumericalDataset.validateIndex(x, length(obj.z)));
        end

        function iAnalyzeAndSetParsedAxes(obj, p)
            if (~isnan(p.Results.x))
                obj.axes_indexes(1) = p.Results.x;
            end
            if (~isnan(p.Results.y))
                obj.axes_indexes(2) = p.Results.y;
            end
            if (~isnan(p.Results.z))
                obj.axes_indexes(3) = p.Results.z;
            end
        end

        function [arglist, axes_indexes, para_value_list] = ...
            iGenerateAxesSliceIndexAndData(obj, para_value_list, varargin)
            axes_indexes = num2cell(obj.axes_indexes);
            for i = 1:length(varargin)
                parameter_name = varargin{i};
                if parameter_name == "x"
                    axes_indexes{1} = ':';
                    para_value_list{i, 1} = obj.x;
                    para_value_list{i, 2} = -3;
                    varargin{i} = [];
                elseif parameter_name == "y"
                    axes_indexes{2} = ':';
                    para_value_list{i} = obj.y;
                    para_value_list{i, 2} = -2;
                    varargin{i} = [];
                elseif parameter_name == "z"
                    axes_indexes{3} = ':';
                    para_value_list{i} = obj.z;
                    para_value_list{i, 2} = -1;
                    varargin{i} = [];
                end
            end

            arglist = varargin;
    end
end
end