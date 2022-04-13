function test_code_varargin(varargin)

if nargin ~=0
disp('Hello world')
batch = varargin{1};
disp(batch.a)

end

end