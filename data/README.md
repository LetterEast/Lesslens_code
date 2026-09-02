# Reconstruction input

Run `create_input_data` from the project root to generate:

```text
data/reconstruction_input.mat
```

The generated file contains registered intensity images, optical parameters,
distance steps and calibration geometry. Reconstruction code reads only this
MAT file and does not access the source image folder.
