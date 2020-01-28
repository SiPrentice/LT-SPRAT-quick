# LT-SPRAT-quick
Quick extractions of LT-SPRAT spectra 

Requires Pyraf, see Astroconda

Download the required spectrum from the LT repository and save in the same folder as the files.

In an IRAF-ready conda environment, type "python QuickSPRATPipeline.py" to run.

Applies a correction to compensate for the flux calibration of the official SPRAT pipeline

3/10/19 QuickSPRATPipeline.py updated to fix issue that calls the correction file in as the spectrum file

28/1/20 Accounts for mirror recorating on the LT by using new corrections file. Alsoaccoutns for object names beginning with a letter
