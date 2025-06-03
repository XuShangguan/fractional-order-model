% ========================================================================
% PSC fractional-order equivalent circuit model (PSC-FoECM)
% This MATLAB script is generated  from a Simulink model.
% It directly implements the following equation:
%   J(t) = J_ph - J_s * (exp((q (V_b + (J_1 + J_ca) R_s))/ n k T) - 1)            
%             - (V_b + (J_1 + J_ca) R_s)/R_p
%             - C_alpha * J_ca
% ========================================================================

function [J, alpha] = fomodel(Vb, I_forward, I_reverse, I_ref, dVb_dt, use_dynamic_alpha)
    dt = 1e-6;
    N = length(Vb);
    t = (0:N-1)' * dt;
    C_alpha = 100e-6 * ones(N,1);

    q = 1.6e-19;
    k = 1.38e-23;
    T = 25 + 273.15;
    n = 1.5;
    Rs = 21.5;
    Rp = 2860.2;
    Isc = 0.05528;
    Voc = 4.61;
    Jph = Isc;
    Js = 1e-6;
    J1 = 0.1;

    %% ----------- Fractional order a -----------
    if use_dynamic_alpha
        % dVb_dt is directly input from outside (measured or calculated)
        Ica = zeros(N,1);
        for i = 1:N
            if dVb_dt(i) >= 0
                Ica(i) = I_forward(i) - I_ref(i); % Forward scan
            else
                Ica(i) = I_reverse(i) - I_ref(i); % Reverse scan
            end
        end

        alpha = zeros(N,1);
        for i = 1:N
            if Ica(i) > 0 && dVb_dt(i) > 0
                alpha(i) = log(Ica(i)/C_alpha(i)) / log(dVb_dt(i));
            else
                alpha(i) = 1; % fallback: ordinary first-order derivative
            end
        end
    else
        alpha = 0.5 * ones(N,1); % Theoretical mode: set a fixed fractional order or use the provided vector directly
    end

    %% ----------- Current calculation -----------
    Jca = fractional_derivative(Vb, t, alpha);
    expo = q * (Vb + (J1 + Jca).*Rs) ./ (k * T * n);
    diode = Js .* (exp(expo) - 1);
    Rsh_term = (Vb + (J1 + Jca).*Rs) ./ Rp;
    J = Jph - diode - Rsh_term - C_alpha .* Jca;
end

function D = fractional_derivative(f, t, alpha)
    N = length(f);
    D = zeros(N,1);
    dt = t(2) - t(1);
    for k = 1:N
        a = alpha(k);
        gamma_arr = ones(k,1);
        for i = 2:k
            gamma_arr(i) = gamma_arr(i-1)*(a-(i-2))/(i-1);
        end
        sum_val = 0;
        for j = 0:(k-1)
            idx = k-j;
            sum_val = sum_val + gamma_arr(j+1)*f(idx);
        end
        D(k) = sum_val / dt^a;
    end
end
