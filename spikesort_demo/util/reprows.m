% Return the matrix v but with each row repeated n times. 

function B = reprows(A,n)

[T N] = size(A);
B = reshape(repmat(A,1,n)',N,n*T)';