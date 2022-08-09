AWG_SHORT_PULSES = [
    [0.6, 1, 0.5, 0.1], # best low ripple pulse, low ripple, 2.1ns FWHM (no overshoot at peak)
    [0.5, 1, 0.7, 0.2], # pretty good, low ripple pulse : 2.3ns FWHM (no overshoot at peak)
    [0.75, 1, 0.3, 0.1], # this is pretty good, moderate ripple, 1.9ns FWHM (moderate ~5% amplitude overshoot)
    [0.7, 1, 0.4, 0.0] # also pretty good, moderate ripple, 2.1ns FWHM (roughly 5% amplitude overshoot at peak)
]

# ~11ns FWHM, ~5ns t_on, no ripple (similar to Keysight 33600A)
AWG_LONG_PULSE = [0.1, 0.3, 0.5, 0.7, 0.9, 1, 1, 1, 1, 1, 1, 0.9, 0.7, 0.5, 0.3, 0.1]