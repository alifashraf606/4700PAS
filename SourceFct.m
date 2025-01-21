function E = SourceFct(t,InputParas) % Calculation of electric field at a given time with specified Input Parameters

if isfield(InputParas,'rep') % Checks for repetition fields (rep), adjusting the time t such that it falls within one repetition period
    n = floor(t/InputParas.rep);
    t = t-n*InputParas.rep;
end

if ~isstruct(InputParas) % Directly use electric field value if input parameter is a struct, otherwise calculate it by using its parameters
    E = InputParas;
else
    E = InputParas.E0*exp(-(t-InputParas.t0)^2/InputParas.wg^2)*exp(1i*(InputParas.we*t + InputParas.phi));
end

end
