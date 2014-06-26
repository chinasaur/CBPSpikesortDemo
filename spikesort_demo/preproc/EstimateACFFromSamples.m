function acf = EstimateACFFromSamples(samples, num_lags)

acf = xcorr(samples, num_lags);
acf = acf(num_lags + 1 : end);
acf = acf - acf(end);
acf = acf ./ max(abs(acf)) .* var(samples);
