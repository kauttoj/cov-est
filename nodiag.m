function mat = nodiag(mat)
    mat = mat - diag(diag(mat));
end