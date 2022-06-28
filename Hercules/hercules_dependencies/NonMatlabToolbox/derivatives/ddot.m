function C = ddot(A, B)
    assert(~isvector(A) && ~isvector(B))
    idx = max(0, ndims(A) - 1);
    B_t = permute(B, circshift(1:ndims(A) + ndims(B), [0, idx - 1]));
    C = sum(squeeze(sum(squeeze(sum(squeeze(sum(bsxfun(@times, A, B_t), idx)), idx)), idx)));
% C = sum(dot(A,B));

% P = A.*B;
% C = sum(P(:));