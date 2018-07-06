# adc-nonlinearity
waveform processing and fitter for the ADC nonlinearity calibration in protoDUNE-SP

1a) ADC bench test rawdata --> (save2ROOT.py) --> ROOT format
1b) (AdcSimulation.cc) --> ROOT format

2) commonFormat.cc (if necesary)

3) Saturation identifier.cc & tagger.cc

4) Pedestal subtractor.cc --> simplified waveform in ROOT format

5) eigenFitter.cc --> effective NL correction & response

6) NLReverser.C: convert the best-fit NL correction to an ideal-to-measurement map

7) predictor.C --> given an input charge, predict signal waveform

8) comparison.C: only for comparison
