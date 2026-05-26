function eigvals = read_zaero_eigs(filename, nmodes)
%READ_ZAERO_EIGS Read FEM eigenvalues from a ZAERO output file.
%
%   eigvals = READ_ZAERO_EIGS(filename, nmodes)
%   eigvals = READ_ZAERO_EIGS(filename)
%
%   Compatibility wrapper around READ_ZAERO_MODAL_INFO. Without nmodes, it
%   returns the full eigenvalue lookup vector indexed by FEM mode number.

info = read_zaero_modal_info(filename);

if nargin < 2 || isempty(nmodes)
    eigvals = info.eigenvaluesByMode;
else
    eigvals = info.eigenvalues(1:nmodes);
end

end
