function C = doubledot(A,B)
% C = doubledot(A,B) tensorial double dot product of tensors A and B,
% where ndims(A) >= 2 & ndims(B) >= 2
assert(ndims(A) >= 2 && ndims(B) >= 2)

AT = repmat(permute(A,[1:ndims(A)-2 ndims(A) ndims(A)-1]),[1 1 1:ndims(B)-2]);
C = squeeze(sum(sum(AT.*B,2),1));
