```@meta
DocTestSetup = quote
    using GaPSE
end
```

# Fast Fourier Transforms and Power Spectra


## The FFTLog package

```@docs
GaPSE.FFTLog.AbstractPlan
GaPSE.FFTLog._c_window
GaPSE.FFTLog._logextrap
GaPSE.FFTLog._zeropad
GaPSE.FFTLog._eval_cm!
GaPSE.FFTLog._eval_ηm!

GaPSE.FFTLog.SingleBesselPlan
GaPSE.FFTLog.HankelPlan
GaPSE.FFTLog._eval_gl
GaPSE.FFTLog._eval_y!
GaPSE.FFTLog._eval_gl_hm!
GaPSE.FFTLog.prepare_FFTLog!
GaPSE.FFTLog.prepare_Hankel!
GaPSE.FFTLog.get_y
GaPSE.FFTLog.evaluate_FFTLog
GaPSE.FFTLog.evaluate_Hankel
GaPSE.FFTLog.evaluate_Hankel!
```

## The Power Spectrum with FFTLog 

```@docs
GaPSE.FFTLog_PS_multipole
GaPSE.FFTLog_all_PS_multipole
```

## The Power Spectrum with TwoFAST

```@docs
GaPSE.TwoFAST_PS_multipole
GaPSE.TwoFAST_all_PS_multipole
```

## The Power Spectrum multipole computation

```@docs
GaPSE.PS_multipole
GaPSE.print_PS_multipole
GaPSE.all_PS_multipole
GaPSE.print_all_PS_multipole
```
