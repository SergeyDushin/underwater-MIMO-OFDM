function  [a, bk_tild, bk_hat, ek, bit_error, y_k_m, x_k_m] =...
    Stocastic_grad(a, y_k_m, x_k_m, y_k_m_old, x_k_m_old, LMS, norm_grad, PIL_symbol)

Bk_hat = x_k_m./x_k_m_old;      %symbol estimates before MRC (we have M estimates)
bk_hat = x_k_m_old' * x_k_m / (x_k_m_old' * x_k_m_old);  %symbol estimate after MRC


if (LMS.is_pilot == 1)     %if pilots are availble (e.g. lower carriers of the first block)
    mu = LMS.mu_PIL;       %use the larger step size for when pilots are avaible
    bk_tild = PIL_symbol;  %use the pilot symbol
    
else                %if in decision directed mode
    mu = LMS.mu_decdir;     %use the smaller step size for decision directed mode
    bit_r   = sign(real (bk_hat * (1+1j)));     %make decision on real axis
    bit_i   = sign(imag (bk_hat * (1+1j)));     %make decision on imag axis
    bk_tild = (1-1j) * (bit_r + 1j*bit_i) / 2;  %regenerate the symbol for finding the ek
end

Ek   = bk_tild - Bk_hat;    %error before MRC (vector of length M)

%if thresholds are met, apply one iteration of the stochastic gradient method:
for m = 1:size(y_k_m,2)
    if (abs(Ek(m)) < LMS.Thrs_Ek)
        g = ((1/(x_k_m_old(m)^2))) * (y_k_m(:,m) * x_k_m_old(m) - x_k_m(m) * y_k_m_old(:,m))* conj(Ek(m));
        if (g < LMS.Thrs_g)
            a(:,m)=a(:,m)+abs(x_k_m_old(m)^norm_grad)*(g*mu);
        end
    end
end

%if the amplitude of the estimated symbol is too large, scale it down
%(neede to practical estimation of MSE)
if(abs(bk_hat)>2)
    bk_hat = 2*exp(1j*angle(bk_hat));
end

ek  = bk_hat-PIL_symbol;                    %error after MRC (used for performance analysis)
bit_error = 0.5*abs(bk_tild-PIL_symbol)^2;  %count the number of error bits (this relationship holds true for QPSK symbols)
