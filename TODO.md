# TO-DOs

## MUST-DOs`

- [ ] Inefficient use in splitting real part and imaginary part in `opf.acpf_consfcn_jac()`
- [ ] (Optional) Use an older version of `numpy` (Latest version for ARM debian is `1.13.3`)
- [x] Use an older version of `scipy` (Latest version for ARM debian is `0.19.1`)

## Code Simplification

- [ ] Use `:` instead of `ndarray.take()`
- [ ] Move `opf.makeYbus()` to `powerflow` module

## Performance Improvement

- [ ] Use sparse matrices when necessary
- [ ] Added `dg`, `d2g`, `dh`, `d2h` for better convergence

## Issues

### Selection of `method` for `scipy.optimize.minimize()`

- [ ] `'trust-ncg'` cannot handle `bounds` nor `constraints`
- [ ] `'COBYLA'` nor `'TNC'` cannot handle equality `constraints`