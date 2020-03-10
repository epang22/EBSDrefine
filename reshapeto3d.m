function dp_3d = reshapeto3d(dp, N)
%RESHAPETO3D
% Reshape vector to 3d cube array
% 12/13/19 (Edward Pang, MIT)
%
%%% Inputs:
% -dp: N^3x1 array of dot products
%%% Outputs:
% -dp_3d: NxNxN reshaped array
% Note: this functions reshapes in different order than built-in matlab
% 'reshape' function



% Reshape data to 3D
dp_3d = zeros(N,N,N);
for oo=1:N
    for pp=1:N
        for qq=1:N
            dp_3d(oo,pp,qq) = dp((oo-1)*N^2+(pp-1)*N+qq);
        end
    end
end



end




