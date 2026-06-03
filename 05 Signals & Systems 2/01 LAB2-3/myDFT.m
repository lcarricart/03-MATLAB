function X = myDFT(x)
    input_was_row = isrow(x);

    x = x(:);
    N = length(x);

    n = 0:N-1;
    k = (-(N-1)/2 : (N-1)/2).';

    nk_matrix = k * n;

    exp_matrix = exp((-1j) * 2 * pi / N * nk_matrix);

    X = exp_matrix * x;

    if input_was_row
        X = X.';
    end
end