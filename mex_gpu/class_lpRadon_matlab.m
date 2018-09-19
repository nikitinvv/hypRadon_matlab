%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef class_lpRadon_matlab < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods (Access = private)
        %% Constructor - Create a new C++ class instance 
        function this = class_lpRadon_matlab(varargin)
            this.objectHandle = class_lpRadon('new', varargin{:});
        end
    end    
    methods (Static)
      function singleObj = getInstance(varargin)
         persistent localObj
         if isempty(localObj) || ~isvalid(localObj)
            localObj = class_lpRadon_matlab(varargin{:});
         end
         singleObj = localObj;
      end
    end  
    methods
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            class_lpRadon('delete', this.objectHandle);
        end
        function varargout = fftlp(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('fftlp', this.objectHandle, varargin{:});
        end
        function varargout = fftlpadj(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('fftlpadj', this.objectHandle, varargin{:});
        end        
        function varargout = convtauqgpu(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('convtauq', this.objectHandle, varargin{:});
        end        
        function varargout = convtxgpu(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('convtx', this.objectHandle, varargin{:});
        end                
        function varargout = convtauqadjgpu(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('convtauqadj', this.objectHandle, varargin{:});
        end                        
        function varargout = convtxadjgpu(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('convtxadj', this.objectHandle, varargin{:});
        end      
        function varargout = fwd(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('fwd', this.objectHandle, varargin{:});
        end            
        function varargout = adj(this, varargin)
            [varargout{1:nargout}] = class_lpRadon('adj', this.objectHandle, varargin{:});
        end                
    end
end
