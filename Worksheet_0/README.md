# Purpose

This is a simple variant of imagesc, a Matlab built-in function, to have a few desirable behaviours as default.

These behaviours are:

1. Always produce a new figure when called.

2. Use the Cartesian axis mode, instead of the "matrix" axis mode.

3. Use an aspect ratio that respects dimension.

4. Use the grayscale color map.

5. Include the color bar.

# Usage

To view a 2-D array, C, as an image, just call

```bash
  imagescxy(C)
```
