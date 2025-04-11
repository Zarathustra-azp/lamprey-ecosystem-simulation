function stability = lyapunov_second_method()
    syms x y z M F s real
    Kx = 15000000;
    Ky = 0.1 * x;
    Kz = 0.15 * (M + F);
    Kl = 0.1 * s;
    Ks = 500000;

    r = 0.15;
    r1 = 0.2;

    qz = 0.24 * (1 - z / Kz);
    kpr = 1;

    mx = 0.0000005;
    mmf = 0.0001;
    ms = 0.00005;
    alpha = 0.56 + (0.78 - 0.56) / (1 + exp(0.6 * y / Ky));

    dy = 0.61;
    dg = 0.25;
    ds = 0.54;
    dz = 0.1;
    dm = 0.5;
    df = 0.5;

    x_dot = r * x * (1 - x / Kx) - y * (mx * x);
    y_dot = -dy * y - dg * y + 30 * F * exp(-y / Ky);
    z_dot = qz * z - dz * z;
    M_dot = alpha * dg * y * exp(-(M + F) / Kl) - (mmf * (M + F)) * z - dm * M;
    F_dot = y * (1 - alpha) * dg * exp(-(M + F) / Kl) - (mmf * (M + F)) * z - df * F;
    s_dot = r1 * exp(-s / Ks) - kpr * ds * (M + F) * (ms * s);
    
    V = x^2 + y^2 + z^2 + M^2 + F^2 + s^2;
    dV_matrix = [partial_derivative(V, x, x_dot); ...
                 partial_derivative(V, y, y_dot); ...
                 partial_derivative(V, z, z_dot); ...
                 partial_derivative(V, M, M_dot); ...
                 partial_derivative(V, F, F_dot); ...
                 partial_derivative(V, s, s_dot)];
    if all(isAlways(eig(dV_matrix) < 0))
        stability = 'Stable';
    else
        stability = 'Unstable';
    end
end

function derivative = partial_derivative(f, variable, dot_variable)
    derivative = diff(f, variable) * dot_variable;
end

